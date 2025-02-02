
# NB modules
import os
import re
import json
import sys
import collections
import hashlib
import shutil
import string
import glob
import copy # creating dict copies
import csv
import subprocess
import numpy as np
from collections import defaultdict
import argparse
from Bio import AlignIO
from collections import Counter
from Bio import SeqIO
#import pandas as pd


#16384 is current hard limit on the number of MSA records in the A3M AF input; we want to read in the paired sequences, perform our own pairing, then estimate 
# for now use interleaved as simple to implement and more balanced

# get the HPI mapping 
# basically, using a DB map between the host and the pathogen and will use this for our sequence pairing,
# similar to AF3, take the unpaired uniprot for both species, if there is info in the DB about the species pairing, pair 
# will make a dictionary of the HPI dt key will be the unique identifier (names joined together) and will reorder based on display name

def getMSAsFromJson(outDir, type='unpaired'):

    jsonPath = [f for f in glob.glob(outDir+'/*_data.json')]
    
    if len(jsonPath) == 0 or len(jsonPath) > 1:
        raise RuntimeError(f'Expected 1 json file in {output_dir} but found {len(jsonPath)}:\n{jsonPath}')
    
    with open(jsonPath[0], mode="r", encoding="utf-8") as json_file:
        json_data = json.load(json_file)

    
    # pull out the relevant info for each of the biomoleucles
    # focus on MSA and templates as longest to generate..
    for biomolecule in json_data.get("sequences", []):

      protein = biomolecule.get("protein", {})
      for bm in [protein]:
        with open(outDir+'/'+prot_name+'.unpaired.a3m', mode="w", encoding="utf-8") as unpairedMSAout:
            unpairedMSAout.write(bm['unpairedMsa'])
        with open(outDir+'/'+prot_name+'.paired.a3m', mode="w", encoding="utf-8") as pairedMSAout:
            pairedMSAout.write(bm['pairedMsa'])

# return dict of a3m and the output file
def readA3M(path):
  allSeqs = collections.defaultdict(str)

  with open(path) as fp:
    currentName = ""
    for line in fp.readlines():
      line = line.strip()
      if line.startswith("#"): #skip header if it contains '#'
        continue
      # start processing
      if line.startswith(">"):
        currentName = line
      else:
        assert currentName != "", "Unexpected format, empty sequence name or sequence data before first >"
        allSeqs[currentName] = allSeqs[currentName] + line
  return(allSeqs)

# get the total number of sequences, look for each sequence to see if it passes or fails the seq filter
# only keep those that pass the sequence filter
def filterGappedMSA(msaDict, maxSeqGap=0.9):

  print(f"Starting number of records in msa: {len(msaDict)}")
  print(f"Filtering out msa records with >= {maxSeqGap} % gapped characters...")

  filteredDict = collections.defaultdict(str)
  for k,v in msaDict.items():
     if int(v.count('-'))/len(v) <  maxSeqGap:
        filteredDict[k] = v
  
  print(f"Filtered out {len(msaDict) - len(filteredDict)} records in msa")
  print(f"Remainging records: {len(filteredDict)} ")
  return(filteredDict) 


# dictionary of species ID and N records assigned to species {speciesID: N sequences}
# sanity check and counts look good
def makeSpeciesDictionary(a3m):
   
   orgIDs = collections.defaultdict()
   for k,v in a3m.items():
      # regex for OX= in the seq header and take taxon ID
      orgID = re.search('OX[=][0-9]+', k)
      if orgID:
        org = re.sub("OX=", "", orgID.group())
        if org not in orgIDs:
          # initialise key and count
          orgIDs[(org)] = 1
        else:
          orgIDs[(org)] += 1
   return(orgIDs)


# deduplicate MSA within species level.. dont want replicate species to artificially inflate..
def deduplicateMSA(msaDict):
   
    cleanMSA = collections.defaultdict(str)
    # track dups within species; so we create a set to capture unique sequences in each
    speciesSet = collections.defaultdict(set)
 
    for k,v in msaDict.items():
      orgID = re.search('OX[=][0-9]+', k)
      if orgID: #get species idx
        org = re.sub("OX=", "", orgID.group())
        
        if v not in speciesSet[org]:
          cleanMSA[k] = v
          speciesSet[org].add(v)

    print(f"{len(cleanMSA)} sequences remaining after deduplication")                  
    return cleanMSA

# get the number of paired records in the paired MSA, to estimate the number of remaining records
# how does AF performing the species pairing? one to many or 1-1? For safety and first test, t
# take as the min size of either set, then use remaing space (from ~16k) to populate the MSA
def getPairedMSADepth(a3m1_path, a3m2_path):
   
   a3m1 = readA3M(a3m1_path)
   a3m2 = readA3M(a3m2_path)
   
   print(f"Removing duplicate MSA records at the species level...")
   clean1 = deduplicateMSA(a3m1)
   clean2 = deduplicateMSA(a3m2)

   sp1Dict = makeSpeciesDictionary(clean1)
   sp2Dict = makeSpeciesDictionary(clean2)

   # converting dicitonary to set keeps names only, and set has intersection function to find overlap
   shared_species = set(sp1Dict).intersection(set(sp2Dict))
   print(f"found {len(shared_species)} shared species...")
   
   #populate subdicts
   print(f"Subsetting MSAs to shared species")
   sp1sub = {}; sp2sub = {}
   for k in shared_species:
      sp1sub[k] = sp1Dict[k]
      sp2sub[k] = sp2Dict[k]

   if sum(sp1sub.values()) < sum(sp2sub.values()):
      print(f'{os.path.basename(a3m1_path)} has the shallowest MSA after deduplication and species pairing:\n{sum(sp1sub.values())}')
      return(sum(sp1sub.values()))
   else:
      print(f'{os.path.basename(a3m2_path)} has the shallowest MSA after deduplication and species pairing:\n{sum(sp2sub.values())}')
      return(sum(sp2sub.values()))


# padding added to A3M MSA rows; I would prefer to add the padding before and after to make the padding clear, but how
# skip records determines how mant MSA records to skip before inserting dummy rows
def padA3MMSARows(a3m_obj, skip_records=1):
   
   outDict = {}
   count = 0 

   # get length of the query seq
   query_key = list(a3m_obj.keys())[0]
   query_len = len(a3m_obj[query_key])

   for k,v in a3m_obj.items():
    outDict[k] = v
    count += 1
    if count >= skip_records:
        outDict[f">dummySeq_{count}"] = "-" * query_len

   print(f"MSA length without padding: {len(a3m_obj.values())}")  
   print(f"MSA length with padding: {len(outDict.values())}")

   outA3M = "\n".join(f"{k}\n{v}" for k,v in outDict.items())
   #result = "\\n".join(f"{k}\\n{v}" for k,v in outDict.items()) # check if the above creates the correct format first..
   return(outA3M)

   # not needed; better to stream
   #with open(out_a3m, mode='w', encoding="utf-8") as outFile:
   #   for k,v in outDict.items():
   #      outFile.write(f'{k}\n{v}\n')

# count the number MSA records (do we need to remove sequences with only padding?)
def getMSADepth(a3mObj):
  a3m_headers = a3mObj.keys()
  return(len(a3m_headers))

def getMSALength(a3mObj):
  a3m_headers = a3mObj.values()
  seq_size = [len(l) for l in a3m_headers]
  print(f"sequence size ranges: {min(seq_size)} - {max(seq_size)} ")


def updateJsonwithPaddedMSA(jsonPath, 
                            outDir,
                            chainIDs):
   
  print(f"Loading {jsonPath}...") 
  with open(jsonPath, mode="r", encoding="utf-8") as json_file:
     json_data = json.load(json_file)

  for idx,id in enumerate(chainIDs):
     
     print(f"Adding padding to unpairedMsa for chain {id} in {jsonPath}...") 
     # might need to simplify this.. change the a3m_obj input to a json file format? Or else different way of running this...
     json_data['sequences'][idx]['protein']['unpairedMsa'] = padA3MMSARows(a3m_obj, skip_records=idx+1)



def updateJsonwithPaddingMSA(jsonPath, 
                             outDir,
                             msa1Path, 
                             msa2Path):

  with open(jsonPath, mode="r", encoding="utf-8") as json_file:
     json_data = json.load(json_file)

  # read the MSAs and assign to the file
  with open(msa1Path, 'r') as a3m_1:
     a3m_1 = a3m_1.read()

  with open(msa2Path, 'r') as a3m_2:
     a3m_2 = a3m_2.read()      

  json_data['sequences'][0]['protein']['unpairedMsa'] = a3m_1
  json_data['sequences'][1]['protein']['unpairedMsa'] = a3m_2
  json_data['sequences'][0]['protein']['pairedMsa'] = ""
  json_data['sequences'][1]['protein']['pairedMsa'] = ""

  print(f"Writing {outDir}/{json_data['name']}_padded.json...")
  # write output 
  with open(outDir+'/'+json_data['name'] +'_padded.json', mode="w", encoding="utf-8") as write_file:
     json.dump(json_data, write_file)
