
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
import Alphafold3_utils as af3
#import pandas as pd

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
  print(f"Filtering out msa records with >= {maxSeqGap} % gapped ('-') characters...")

  filteredDict = collections.defaultdict(str)
  for k,v in msaDict.items():
     if int(v.count('-'))/len(v) <  maxSeqGap:
        filteredDict[k] = v
  
  print(f"Filtered out {len(msaDict) - len(filteredDict)} records in msa")
  print(f"Remainging records: {len(filteredDict)} ")
  return(filteredDict) 


# dictionary of species ID and N records assigned to species {speciesID: N sequences}
# change to  tuple with N sequences and list of sequences?
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
    # track duplicates within species 
    speciesSet = collections.defaultdict(set)
 
    for k,v in msaDict.items():
      orgID = re.search('OX[=][0-9]+', k)
      if orgID: #get species idx
        org = re.sub("OX=", "", orgID.group())
        
        if v not in speciesSet[org]:
          cleanMSA[k] = v
          speciesSet[org].add(v) # append to the species set

    print(f"{len(cleanMSA)} sequences remaining after deduplication")                  
    return cleanMSA

def mapSequence(msa, option='map'):
    '''Convert MSA seq to/from numeric numpy array'''

    assert option %in% ['map', 'backmap']

    backmapDict = { 1:'A', 2:'C', 3:'D', 4:'E', 5:'F',6:'G' ,7:'H',
                    8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P',14:'Q',
                    15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 21:'-'} #Here all unusual AAs and gaps are set to the same char
    
    mapDict = {value: key for key, value in backmap.items()}

    if option  == 'map':
      mapping = mapDict
    else:
      mapping = backmapDict
    # use this to backmap from 
    seq = ''.join([mapping[ch] for ch in msa])
    return(seq)

# another read_a3m file
def read_a3m(infile,max_gap_fraction=0.9):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    parsed = [] #Save extracted msa
    species = [] #  OX id 
    header = [] # seq header

    seqlen = 0
    lc = 0
    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'): #OX=OrganismIdentifier
                if 'OX=' in line:
                    OX=line.split('OX=')[1]
                    if len(OX)>0:
                        try:
                            species.append(int(OX.split(' ')[0]))
                            header.append(line)
                        except:
                            species.append(0)
                            header.append(line)
                    else:
                        species.append(0)
                        header.append(line)
                elif 'TaxID=' in line:
                    OX= line.split('TaxID=')[1]
                    if len(OX)>0:
                        try:
                            species.append(int(OX.split(' ')[0]))
                            header.append(line)
                        except:
                            species.append(0) 
                            header.append(line)                           
                    else:
                        species.append(0)
                        header.append(line)
                else:
                    species.append(0)
                    header.append(line)
                continue
            line = line.rstrip()
            gap_fraction = line.count('-') / float(len(line))
            # if gapped > limit, drop this species
            if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
                parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()]) #convert msa line to numeric
            else:
                if len(species)>1:
                    species = species[:-1] #Remove the previously stored species
                    header = header[:-1] # and header
                    continue
            #Check that the lengths match
            if len(parsed[-1])!=seqlen and lc>=1:
                parsed = parsed[:-1]
                species = species[:-1]
                header = header[:-1]
                continue
            seqlen = len(parsed[-1])
            lc+=1


    return np.array(parsed, dtype=np.int8, order='F'), np.array(species), np.array(header)

# test script from FoldDock as template for what we need
def pair_MSAs(ox1, ox2, msa1, msa2, header1, header2):
    '''
    Following AF3 MSA pairing, concatentate all seq from same species across chains.
    If sp sequences in one chain more abundant than the other, include padding
    TODO try other methods like the Bryant pairing
    '''
    #Don't remove the zeros (no OX), then the query sequences (first line)
    #will be removed
    matching_ox = np.intersect1d(ox1,ox2)
    print(f'Found {len(matching_ox)-1} shared species in msas')

    paired_msa1, paired_msa2 = [], []
    paired_seq1, paired_seq2 = [], []
    headers1, headers2 = [], []

    #Go through all matching and select the first (top) hit
    for ox in matching_ox:
        idx1 = np.argwhere(ox1 == ox).flatten()  # get all occurrences in ox1 in list
        idx2 = np.argwhere(ox2 == ox).flatten()  # get all occurrences in ox2 in list

        num1, num2 = len(idx1), len(idx2) # get number of sequences matching to that species
        max_num = max(num1, num2) #which has the most

        # padds the shortest seq indices with -1 to flag as missing; these are filled with '-'
        idx1 = np.pad(idx1, (0, max_num - num1), constant_values=-1)
        idx2 = np.pad(idx2, (0, max_num - num2), constant_values=-1)

        for i in range(max_num):
            seq1 = msa1[idx1[i]] if idx1[i] != -1 else np.full(msa1.shape[1], 21)
            seq2 = msa2[idx2[i]] if idx2[i] != -1 else np.full(msa2.shape[1], 21)

            paired_msa1.append(seq1)
            paired_msa2.append(seq2)
            
            headers1.append(header1[idx1[i]] if idx1[i] != -1 else f'>{ox}_padding_{i}')
            headers2.append(header2[idx2[i]] if idx2[i] != -1 else f'>{ox}_padding_{i}')

    # convert numeric code back to msa sequence
    for i in range(max_num):
      paired_msa1[i] = backmapSequence(paired_msa1[i])
      paired_msa2[i] = backmapSequence(paired_msa2[i])

    if len(paired_msa1) != len(paired_msa2):
       print('missmatching number of records in msa... aborting')
       exit
    return list(zip(header1, paired_msa1)), list(zip(header2, paired_msa2)), len(paired_msa1)

# dont need to worry about the length, just need to divide the remainder into equally sized blocks depending on chain N
# then find a way to add these to the input json, but dont want to store!
# divide up the remainder into equally sized blocks and concatenate the sequences

# get list of paths to blockchain MSA
def blockChainMSA(unpairedMSAs=[], depth=16384/2):
  assert(type(unpairedMSAs) is list)

  nChains = len(unpairedMSAs)
  chainDepth = depth/nChains
  print(f"returning {chainDepth} records per MSA..")

  chainID = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  msa_list = []
  ox_list = []
  header_list = []
  msa_width = [] # capture N sequecnes in each MSA
  msa_length = []

  for i,msa in enumerate(unpairedMSAs):
    print(f'getting unpaired MSA for {chainID[i]}')
    msa, ox, header = read_a3m(msa)

    msa_length.append(msa.shape[0])
    msa_width.append(msa.shape[1])
    msa_list.append(msa)
    ox_list.append(ox)
    header_list.append(header)

  merge_msa = np.full((sum(msa_length), sum(msa_width)), 21, dtype=int)
  print(merge_msa.shape)
  print(merge_msa)



  return msa_list, ox_list, header_list
  # I think here we use the same method to find the msaFilepath, and read in using the above 


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

   # use intersection set method to find overlapping names
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


# for now adding 1 after other but want a block diagonal matrix; first need to get depth of paired matrix
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
