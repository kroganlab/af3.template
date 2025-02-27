
# modules
import os
import re
import json
import sys
import collections
import glob
import random
import numpy as np
import argparse

def getMSAsFromJson(outDir, type='unpaired'):

    jsonPath = [f for f in glob.glob(outDir+'/*_data.json')]
    
    if len(jsonPath) == 0 or len(jsonPath) > 1:
        raise RuntimeError(f'Expected 1 json file in {output_dir} but found {len(jsonPath)}:\n{jsonPath}')
    
    with open(jsonPath[0], mode="r", encoding="utf-8") as json_file:
        json_data = json.load(json_file)

    # pull out the relevant info for each of the biomoleucles
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


# deduplicate MSA within species level..
#  dont want replicate sequences from same species to artificially inflate signal..
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

    assert option == 'map' or option == 'backmap'

    backmapDict = { 1:'A', 2:'C', 3:'D', 4:'E', 5:'F',6:'G' ,7:'H',
                    8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P',14:'Q',
                    15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 21:'-'} #Here all unusual AAs and gaps are set to the same char
    
    mapDict = {value: key for key, value in backmapDict.items()}

    if option  == 'map':
      mapping = mapDict
    else:
      mapping = backmapDict

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
            if line.startswith('#'): #skip AP introduced lines
              continue
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
            if line: # ignore empty lines
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
      paired_msa1[i] = mapSequence(paired_msa1[i], option='backmap')
      paired_msa2[i] = mapSequence(paired_msa2[i], option='backmap')

    if len(paired_msa1) != len(paired_msa2):
       print('missmatching number of records in msa... aborting')
       exit
    return list(zip(header1, paired_msa1)), list(zip(header2, paired_msa2)), len(paired_msa1)

## TODO function works, but seems it is not the query that is added to the top.. is this unordered?

def blockChainMSA(unpairedMSAs=[], depth=16384/2, saveMSAs=False, outDir=False):
  '''
  Given a list of MSA in A3M format, generate a merged block diagonal matrix
  Chains are concatenated along the n dimension (column) with missing characters populated with -
  Idea here is to avoid merging unrealted sequences and adding noise to the paired representation
  Depth of unpaired MSA is hardcoded in AF3 pipeline
  Note: moving query to top and stays concatenated to avoid failure due to sequence check
  '''
  assert(type(unpairedMSAs) is list)

  nChains = len(unpairedMSAs)
  chainDepth = depth/nChains
  if type(chainDepth) == float:
       chainDepth = round(depth/nChains)
     
  print(f"extracting {chainDepth} records per MSA..")
  
  chainID = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  msa_list = []
  ox_list = []
  header_list = []
  msa_width = []
  msa_length = []
 
  querySeq = []
  #queryLen = 0
  for i,msa in enumerate(unpairedMSAs):

    print(f'getting unpaired MSA for chain {chainID[i]}')
    msa, ox, header = read_a3m(msa)

    # get the query sequence from each 
    querySeq.append(msa[0])
    #queryLen += len(msa[0])

    # drop the query row
    msa = np.delete(msa,0, axis=0); header = np.delete(header,0, axis=0); ox = np.delete(ox,0, axis=0)

    # dim of merge_MSA
    msa_length.append(msa.shape[0])
    msa_width.append(msa.shape[1])

    msa_list.append(msa)
    ox_list.append(ox)
    header_list.append(header)

  # concatenate hte query sequences
  querySeq = np.hstack(querySeq)

  # create null merged msa and populate with empty characters ('-')
  merge_msa = np.full((chainDepth *nChains, sum(msa_width)), 21, dtype=int)
  newNames=np.array([], dtype='<U263')
  currentDepth=0; currentLen=0; 

  for i,msa in enumerate(unpairedMSAs):
     
     # update the merged MSA with sequences 
     merge_msa[currentDepth:chainDepth*(i+1), currentLen:currentLen+msa_width[i]] = msa_list[i][:chainDepth,:msa_width[i]]
     currentDepth=chainDepth*(i+1)
     currentLen=currentLen + msa_width[i]
     # take names and append the merged msa
     newNames = np.append(newNames, np.char.add(header_list[i][:chainDepth], ' mergedMSA'))

  print('Moving query sequences to first record...')
  # cbind the query, then rbind
  merge_msa = np.vstack((querySeq, merge_msa))
  newNames = np.insert(newNames, 0, '>query mergedMSA')

  outMSAs = []
  start = 0

  for width in msa_width:
    subset_msa = merge_msa[:, start:start + width] # Slice columns; keep all rows
    # backmap to chr
    processed_msa = np.array([mapSequence(row, option='backmap') for row in subset_msa])
    outMSAs.append((newNames, processed_msa))
    start += width

  if saveMSAs and outDir:
    print(f'Saving MSA to {outDir}')
    
    for i, (msa_header, msa) in enumerate(outMSAs):  
        file_path = f"{outDir}/{chainID[i]}.blockChain.af3"
        with open(file_path, 'w') as a3m_out:
          # a little mealy, but basically iterate over msa records and write header, followed by concatenated str
          a3m_out.write("\n".join(f"{name}\n{''.join(map(str, row))}" for name, row in zip(msa_header, msa)) + "\n")

  return outMSAs


def unpackMSA(msaObj):
   '''
   Unpack the MSAs produced by blockChainMSA
   '''
   return "\n".join(f"{name}\n{''.join(map(str, row))}" for name, row in zip(msaObj[0], msaObj[1])) + "\n"

def getMSAIndex(unpairedMSA, nrecords):
   'Given n records, extract the indices of records to keep from the MSA'
   'Select based on OX identifer as we want evolutionary diversity in our AF MSA'
   'return list of indices of the extracted records; use this to subset MSAs'
   
   msa,ox,header = read_a3m(unpairedMSA)

   msa_depth, msa_length = msa.shape
   assert nrecords <= msa_depth, f"Number of requested records is greater than MSA depth: {msa_depth}"

   speciesDict = {}
   indices = []

   for org in ox:
      idx = np.argwhere(org == ox).flatten()  # get all occurrences in org in list
      if org not in speciesDict.keys():
        speciesDict[org] = idx
        # randomly assign a record per species
        indices.append(random.choice(idx))
    
   if len(indices) < nrecords:
      print(f"Number of unique organisms (OX) is less than number of requested records:\n{len(indices)} < {nrecords} ")
      print(f"randomly resampling remaining records to meet the requested MSA depth")

   # return a random subset of records from the remainder  
   remainingRecords = list(set(range(1, msa_depth)) - set(indices))

   indices = indices + random.sample(remainingRecords, nrecords-len(indices))
   return(indices)


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
