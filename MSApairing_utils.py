
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
from Bio import SeqIO
#import pandas as pd


# get the HPI mapping 
# basically, using a DB map between the host and the pathogen and will use this for our sequence pairing,
# similar to AF3, take the unpaired uniprot for both species, if there is info in the DB about the species pairing, pair 
# will make a dictionary of the HPI dt key will be the unique identifier (names joined together) and will reorder based on display name

def getMSAsFromJson(outDir, type='unpaired'):

    jsonPaths = [f for f in glob.glob(outDir+'/*_data.json')]
    
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

def read_a3m_msa(path):
  allSeqs = collections.defaultdict(str)
  namesInOrder = []

  with open(path) as fp:
    currentName = ""
    for line in fp.readlines():
      line = line.strip()
      if line.startswith("#"): #skip header if it contains '#'
        continue
      # start processing
      if line.startswith(">"):
        currentName = line
        namesInOrder.append(currentName)
      else:
        assert currentName != "", "Unexpected format, empty sequence name or sequence data before first >"
        allSeqs[currentName] = allSeqs[currentName] + line
  # return as tuple structure (not sure if the MSA has redundancy so safer maybe...)
  return(allSeqs, namesInOrder)


# padding added to A3M MSA rows
# skip records determines how mant MSA records to skip before inserting dummy rows
def padA3MMSARows(a3m_obj, out_a3m, skip_records=1):
   
   outDict = {}
   count = 0 

   # get length of the query
   query_key = list(a3m_obj[0].keys())[0]
   query_len = len(a3m_obj[0][query_key])

   for k,v in a3m_obj[0].items():
    outDict[k] = v
    count += 1
    if count >= skip_records:
        outDict[f">dummySeq_{count}"] = "-" * query_len

   print(f"MSA length without padding: {len(a3m_obj[0].values())}")  
   print(f"MSA length with padding: {len(outDict.values())}")  
   print(f"Writing output to {out_a3m}...")

   with open(out_a3m, mode='w', encoding="utf-8") as outFile:
      for k,v in outDict.items():
         outFile.write(f'{k}\n{v}\n')


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