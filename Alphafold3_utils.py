#!/usr/bin/python3

# Collection of helper functions used to run the AF3 pipeline
# Mostly adapted from https://github.com/kroganlab/bp_utils/blob/master/af.template.dir/AF_saveMSAS.231.py

# to do; fix the dict creation; if duplicates in fa input it concatenates the sequence. want to replace for safety

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
#import pandas as pd


# extract name from header after the ID component |sp|tr|up, else return complete string
def nameFromHeader (header):
  pipeParts = header[1:].split("|") #skips > and splits on |
  if(pipeParts[0] in ("sp", "tr", "up")):
    return(pipeParts[1])
  else:
    return(pipeParts[0])

def af3_alphaFoldRunOutputDirectory(seq_list, outputPath):
  seq_list = [seq.lower() for seq in seq_list] #mimic 'tidy' format of AF3 output folder
  return os.path.join(outputPath, f"{'__'.join(seq_list)}")

# modified to drop duplicate values
# modify to handle dimers?
def read_fasta(path):
  allSeqs = collections.defaultdict(str)
  namesInOrder = []
  with open(path) as fp:
    currentName = ""
    for line in fp.readlines():
      line = line.strip()
      if line[0:1] == ">":
        currentName = nameFromHeader(line)
        if len(namesInOrder) == 0 or namesInOrder[-1] != currentName:
          namesInOrder.append (currentName)
      else:
        assert currentName != "", "Unexpected format, empty sequence name or sequence data before first >"
        allSeqs[currentName] = allSeqs[currentName] + line
  return allSeqs, namesInOrder


# given a sequence get the hex and check alignment repo for matching folder
# if found, check the input sequence matches the query in the a3m file
# if true, return the collapsed a3m, else return 'null' and let AF construct the MSA run
def AF3_getMSADirectoryForSequence(sequence, alignmentRepo, MSAType='unpaired'):
    # creates a hash of the input sequence
    seqHash = hashlib.md5(sequence.upper().encode('utf-8')).hexdigest()
    # subdir take the first two elements
    subDir = seqHash[0:2]
    dir = os.path.join(alignmentRepo, subDir, seqHash)
    #if directory exists, raise an error if seqeunces mismatch
    try:
        msaPaths = []
        msaPaths = [p for p in os.listdir(dir) if p.split(".")[-2] == MSAType]
    except:
        print(f"No MSA found in {dir}..\nReturning null to enforce MSA creation")
        return None
    if len(msaPaths) > 0:
        print(msaPaths)
        for msa in msaPaths:
            print(msa)
            with open(dir+'/'+msa,'r') as msa_file:
              for line in msa_file:
                if line.startswith('>query'):
                    a3m_seq = msa_file.readline().rstrip('\n') # skip to next line and get query seq...assumption is query is exact match
                    if sequence != a3m_seq:
                        print(f"Sequence Mismatch Error in {dir}:")
                        print(f"Mismatch in \n{sequence}\n{a3m_seq}")
                        print(f"Returning null to enforce fresh MSA creation")
                        return None
                    else:
                        print(f"Extracting A3M from {dir}:")
                        with open(dir+'/'+msa, 'r') as a3m:
                          if 'unpaired' in msa:
                            print('Preparing unpaired MSA')
                            next(a3m)  # Skip A3M header
                            a3m = a3m.read()
                            return(a3m)
                          else: 
                            print('Preparing paired MSA')
                            a3m = a3m.read()
                            return(a3m)
                        
# need to check if null enforces template search or whether to provide an empty string
def AF3_getTemplateJSONDirectoryForSequence(sequence, alignmentRepo):
    # creates a hash of the input sequence
    seqHash = hashlib.md5(sequence.upper().encode('utf-8')).hexdigest()
    # subdir take the first two elements
    subDir = seqHash[0:2]
    dir = os.path.join(alignmentRepo, subDir, seqHash)
    #if directory exists, raise an error if seqeunces mismatch
    try:
        templatePaths = []
        templatePaths = [p for p in os.listdir(dir) if p.split(".")[-2] == 'templates']
    except:
        print(f"No template json found in {dir}..\nReturning null to enforce template search")
        return None
    if len(templatePaths) > 0:
        print(templatePaths)
        for file in templatePaths:
            with open(dir+'/'+file,'r') as json_file:
              json_template = json.load(json_file)
              return(json_template)

def AF3_collapseMSAs(a3m_input):
    msa_file =  open(a3m_input, 'r',  encoding="utf-8")
    # read all lines into a list, collapse to \n seperated 
    all_records = msa_file.readlines()
    collapse_records = ''.join(all_records).replace('\n', '\\n')
    return(collapse_records)

# create protein, dna, rna, ligand empty json structures to populate
# both paired and unpaired msa must be both set or unset; here we want to set pairedMSA so AF takes care of the species pairing for multimer
# recommendaiton is to perform this manually, but for now use preset
# https://github.com/google-deepmind/alphafold3/issues/171 here we set pairedMsa so AF performs the merging..
# see https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
# setting paired unpaired MSA but unsetting templates ensures AF also performs the template search

protein_template = {
  "protein": {
    "id": None,
    "sequence": None,
    "modifications": [],
    "unpairedMsa": None,
    "pairedMsa": None,
    "templates": []
  }
}

rna_template = {
  "rna": {
    "id": None,
    "sequence": None,
    "modifications": None,
    "unpairedMsa": None
  }
}
dna_template = {
  "dna": {
    "id": None,
    "sequence": None,
    "modifications": None,
  }
}

ligand_template = {
  "ligand": {
    "id": None,
    "ccdCodes": None
  }
} 
# seems to work for now to construct AF3 job.. going to test with a handful of sequences
def af3_setupJob(job_id, 
                 jobTable, 
                 master_fasta,
                 output_dir, 
                 setup_job,
                 alignmentRepo, 
                 nSeeds, # more models recommended by AP for more accurate inference
                 json_template):
    
    jobID = job_id
    with open(jobTable) as fp:
        for line in fp.readlines():
            jobTableCols = [word.strip() for word in line.strip().split(",")]
            j = jobTableCols[0]
            seq_list = jobTableCols[1:] #take all other cols
            
            if(int(j) == jobID):
                break
    assert int(j) == jobID, f"job_id {jobID} not found in AlphaFoldJobList at {jobTable}"

    outDir = af3_alphaFoldRunOutputDirectory(seq_list, output_dir)
    print(f'Saving AF3 output to {outDir}\n')

    if (setup_job):  
        # prepare the output directory
        try:
            os.makedirs(outDir)
        except FileExistsError:
            pass  
           
        # read in the empty template 
        print('Preparing json input...\n')
        with open(json_template, mode="r", encoding="utf-8") as read_file:
            af3_json_template = json.load(read_file)
        
        af3_json_template['name'] = f"{'__'.join(seq_list)}"   
        
        # add seed numbers to the matrix
        modelSeeds = [n for n in range(int(nSeeds))]
        af3_json_template['modelSeeds'] = modelSeeds
            
         # function returns tuple of dict (record,seq) and list (record )
        seqs,inOrder = read_fasta(master_fasta)
        
        # add chain ID
        chainID = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        
        # create a af3 test input file
        protein_list = []
        
        # loop through the sequences and populate the relevant protein level info,
        # if I can fix this, should be in fit enough state to test with multiple sequences

        for i,seqID in enumerate(seq_list):
            print(seq_list)
            # creating a deep copy of template
            clean_template = copy.deepcopy(protein_template)
            clean_template['protein']['id'] = chainID[i]
            clean_template['protein']['sequence'] = seqs[seqID]

            # add protein template info...
            clean_template['protein']['templates'] = AF3_getTemplateJSONDirectoryForSequence(sequence=seqs[seqID], alignmentRepo=alignmentRepo)

            # pull in MSA if available; otherwise generate. Details scarce on github and misleading; need to supply both MSA to the algorithm..
            # see https://github.com/google-deepmind/alphafold3/issues/171

            # supplying both paired and unpaired seqs as they are different
            clean_template['protein']['unpairedMsa'] = AF3_getMSADirectoryForSequence(sequence=seqs[seqID], alignmentRepo=alignmentRepo, MSAType='unpaired')
            clean_template['protein']['pairedMsa']   = AF3_getMSADirectoryForSequence(sequence=seqs[seqID], alignmentRepo=alignmentRepo, MSAType='paired')

            protein_list.append(clean_template)

        af3_json_template['sequences'] =  protein_list
        
        with open(outDir+'/'+af3_json_template['name']+'.af3_input.json', mode="w", encoding="utf-8") as write_file:
            json.dump(af3_json_template, write_file)
        
        return(outDir+'/'+af3_json_template['name']+'.af3_input.json')
    
# postProcessing work
#need to revisit this and completely modify..
# After pipeline completes, recover the MSA from the AF3 output directory and store in alignmentRepo
# check for chainA and chainB and use that naming in the output
def af3_captureMSAs(output_dir, alignmentRepo):

  jsonPath = [f for f in glob.glob(output_dir+'/*_data.json')]

  if len(jsonPath) == 0 or len(jsonPath) > 1:
    raise RuntimeError(f'Expected 1 json file in {output_dir} but found {len(jsonPath)}:\n{jsonPath}')
  
  seqs = [ seq.upper() for seq in os.path.basename(jsonPath[0]).split('_data.')[0].split('__') ]
  chainID = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  chainIDs = [id for id in chainID[0:len(seqs)]]

  with open(jsonPath[0], mode="r", encoding="utf-8") as json_file:
    json_data = json.load(json_file)

    # pull out the relevant info for each of the biomoleucles
    # focus on MSA and templates as longest to generate..
    for biomolecule in json_data.get("sequences", []):

      protein = biomolecule.get("protein", {})
      dna = biomolecule.get("dna", {})
      rna = biomolecule.get("rna", {})
      ligand = biomolecule.get("ligand", {})

      #for bm in [protein, dna, rna, ligand]:
      for bm in [protein]: 
        if len(bm.keys()) != 0: # skip empty entries
          if "sequence" in bm.keys():
            if "\n" in bm["sequence"]:
              print(f"Removing newline character in {bm['sequence']}")
              bm['sequence'] = bm['sequence'].replace("\n", "")

            seqHash = hashlib.md5(bm['sequence'].upper().encode('utf-8')).hexdigest()
            subDir = seqHash[0:2]
            outDir = os.path.join(alignmentRepo, subDir, seqHash)
            print(f"Checking {alignmentRepo} for {bm['sequence']} MSA...")

            if os.path.isdir(outDir):
              print(f"Directory exists: {outDir}\nNot overwriting")
            else:
              os.makedirs(outDir)
              prot_idx = chainIDs.index(bm['id'])
              prot_name = seqs[prot_idx]
              print(f"Saving MSA and template files for {prot_name} in {outDir}...")

              # with open(outDir+'/'+prot_name+'.paired.a3m', mode="w", encoding="utf-8") as a3m_out:
              #   a3m_out.write(bm['pairedMsa'])
              
              # with open(outDir+'/'+prot_name+'.unpaired.a3m', mode="w", encoding="utf-8") as a3m_out:
              #   a3m_out.write(bm['unpairedMsa'])

              if not glob.glob(outDir+'/'+'*.unpaired.a3m'):
                print('unpaired MSA not found... writing from JSON')
                with open(outDir+'/'+prot_name+'.unpaired.a3m', mode="w", encoding="utf-8") as unpairedMSAout:
                  unpairedMSAout.write(bm['unpairedMsa'])
              else:
                  print('unpaired MSA found. Not overwriting...')

              if not glob.glob(outDir+'/'+'*.paired.a3m'):
                print('paired MSA not found... writing from JSON')
                with open(outDir+'/'+prot_name+'.paired.a3m', mode="w", encoding="utf-8") as pairedMSAout:
                  pairedMSAout.write(bm['pairedMsa'])
              else:
                print('paired MSA found. Not overwriting...')

              if not glob.glob(outDir+'/'+'*.templates.json'):
                print('templates not found... writing from JSON')
                with open(outDir+'/'+prot_name+'.templates.json', mode="w", encoding="utf-8") as templatesout:
                  json.dump(bm['templates'], templatesout)

              


# post run score processing; just capture all ptm, iptm and ranking from each of the individual models and write as csv
def af3_captureSummaryScores(output_dir):

  # capture all json scores
  json_paths  = [f for f in glob.glob(output_dir+'/*/summary_confidences.json')]
  cif_paths  =  [f for f in glob.glob(output_dir+'/*/model.cif')]
  outFile = os.path.basename(output_dir)

  mod = []
  ranking = []
  iptm = []
  ptm = []
  pDockq = [] # this is a summary score for interface; complementary(??) to iptm?

  # for cif set 
  #mod = []
  #plddt = [] # want to capture a per residue score

  # iterate through all the files 
  for json_file in json_paths:
    model_name=outFile + '.' + os.path.dirname(json_file).split('/')[-1]
    mod.append(model_name)

    with open(json_file, mode="r", encoding="utf-8") as jf:
      json_data = json.load(jf)
      iptm.append(json_data['iptm'])
      ptm.append(json_data['ptm'])
      ranking.append(json_data['ranking_score'])

    # now calculate the pDockQ score
    #with open()

  with open(output_dir+'/'+outFile+'_summaryScores.csv', mode="w", encoding="utf-8") as outf:
    writer = csv.writer(outf)
    writer.writerow(['model', 'ranking', 'ptm', 'iptm'])
    writer.writerows(zip(mod, ranking, ptm, iptm))

  # alternative: 
  ## pandas option; dont use as need to load additional package
  # creat dictionary and convert to dt
  #outdict = {'model': mod, 'ranking': ranking, 'iptm': iptm, 'ptm': ptm}

  #df = pd.DataFrame(outdict)  
  # saving the dataframe
  #df.to_csv(output_dir+'/'+'summary.scores.csv')
  
# MSA pairing is an important step for multimer models
# by default, AF3 uses taxonIDs to reorder MSAs rowwise and pair species before concatenating
# we would also like to use this pairing between hosts and pathogens
# read in two MSAs; one for host, one for pathogen
#def af3_hostPathogenMSAPairing(hostMSA, pathogenMSA, mappingDB):
   
def af3_createMasterFastafromMSADirectory(inputFile, msaDirectory, outputFile):
  
  print('Creating master fasta file...')
  # read in protein/gene IDs in input file
  with open(inputFile, mode="r", encoding="utf-8") as f:
    proteinIDs = [id.strip() for id in f.readlines()]
    proteinIDs = set(proteinIDs) # convert to set to remove duplicate IDs if present

  print(f"Found {len(proteinIDs)} protein IDs in {inputFile}")

  # now loop over the MSA directory, pull out the query sequence and write to
  with open(outputFile, mode="w", encoding="utf-8") as outFile:
    counter=0
    for i, proteinID in enumerate(proteinIDs):
      counter=i+1
      msaFile = glob.glob(msaDirectory+'/'+proteinID+'.a3m')

      # check this works correctly
      with open(msaFile[0], mode="r", encoding="utf-8") as a3mFile:
        for line in a3mFile:
          if line.startswith('>query'):
            a3m_seq = a3mFile.readline().rstrip('\n') # skip to next line and get query seq.
            outFile.write(f'>{proteinID}\n{a3m_seq}\n')
            break
      
  # count number of fasta records
  print('Number of sequences in master fasta:')
  nSeqs = subprocess.run("grep '>' "+outputFile+" | wc -l", shell = True, executable="/bin/bash")
  
# aim of this function is to extract features from AF3 monomer JSON output and wrtie to a directory for extraction of features for later runs
# for now pulling everything: unpaired, paired and templates
def af3_extractFeaturesfromAF3JSON(jsonInput, msaDirectory):
  
  if os.path.isfile(jsonInput):
    print(f'Searching {jsonInput} for features to extract...')
  else:
    sys.exit(f'File {jsonInput} not found')

  # read in protein/gene IDs in input file
  with open(jsonInput, mode="r", encoding="utf-8") as f:
    json_data = json.load(f)
    name = json_data.get("name", {})

    # returns val for key (which is list) for sequence
    for feature in json_data.get("sequences", []):
      protein = feature.get("protein", {})

      # extract the features we want to keep for reruns
      sequence=protein.get('sequence', {})
      pairedMSA=protein.get('pairedMsa', {})
      unpairedMSA=protein.get('unpairedMsa', {})
      templates=protein.get('templates', {})

      # want to check if the directory already exists for this sequence
      seqHash = hashlib.md5(sequence.upper().encode('utf-8')).hexdigest()
      subDir = seqHash[0:2]
      outDir = os.path.join(msaDirectory, subDir, seqHash)

      if not os.path.exists(outDir):
       print(f'No directory found for sequence {name}\nCreating {outDir} to store MSA/template features...')
       os.makedirs(outDir)
      else:
       print(f'{name} directory exits:\n{outDir}')
       print(f'Checking directory for paired MSA, unpaired MSA and templates...')
      
      # could we make this safer? enforce that the name of the MSA is the same as the name of the JSON file?
      # probably will create redundant files as user is naming convention
      if not glob.glob(outDir+'/'+'*.unpaired.a3m'):
         print('unpaired MSA not found... writing from JSON')
         with open(outDir+'/'+name+'.unpaired.a3m', mode="w", encoding="utf-8") as unpairedMSAout:
            unpairedMSAout.write(unpairedMSA)
      else:
         print('unpaired MSA found. Not overwriting...')

      if not glob.glob(outDir+'/'+'*.paired.a3m'):
         print('paired MSA not found... writing from JSON')
         with open(outDir+'/'+name+'.paired.a3m', mode="w", encoding="utf-8") as pairedMSAout:
            pairedMSAout.write(pairedMSA)
      else:
         print('paired MSA found. Not overwriting...')

      if not glob.glob(outDir+'/'+'*.templates.json'):
         print('templates not found... writing from JSON')
         with open(outDir+'/'+name+'.templates.json', mode="w", encoding="utf-8") as templatesout:
            json.dump(templates, templatesout)
      else:
         print('templates json found. Not overwriting...')


# update all features in AF directory and write to output
def updateAF3MSADirectory(AFJsonDir, outputMSADir):

  jsonFiles = glob.glob(AFJsonDir+"**/**/*_data.json", recursive=True)
  print(f"Found {len(jsonFiles)} in {outputMSADir}")

  print(f'Writing MSAs to {outputMSADir}')

  for json in jsonFiles:
    print('Command:')
    print(f"af3_extractFeaturesfromAF3JSON(jsonInput={json}, msaDirectory={outputMSADir})")
    af3_extractFeaturesfromAF3JSON(jsonInput=json, msaDirectory=outputMSADir)

  print(f'Done...')


# generate MSA and PAE plots for the top scoring model per PPI
def generate_MSAandPAEplots(outDir):

  msa_file =  glob.glob(outDir+'/'+'*msaOut.txt').pop()
  data_json = glob.glob(outDir+'/'+'*_data.json').pop()

  sample_name = os.path.basename(os.path.dirname(data_json))
  confidence_json = glob.glob(outDir+'/'+sample_name+'_confidences.json').pop()

  print(F'Found files:\n{confidence_json}\n{data_json}\n{msa_file}')

  cmd_args = ['Rscript', 'Visualize_PAEandMSA.R', msa_file, data_json, confidence_json]
  cmd = ' '.join(cmd_args)
  print("\n")
  print(cmd)
  print()

  try:
    subprocess.run(cmd,
        stdout = sys.stdout, 
        stderr = sys.stderr,
        shell = True,  # module command is a csh alias on Wynton
        executable = '/bin/csh',
        check = True)
  except Exception as error:
    print('Error executing Visualize_PAEandMSA.R:\n {}'.format(error))

def get_interchainContactsPAE(outDir):

  cif_file =  glob.glob(outDir+'/'+'*_model.cif').pop()
  sample_name = os.path.basename(os.path.dirname(cif_file))
  confidence_json = glob.glob(outDir+'/'+sample_name+'_confidences.json').pop()

  print(F'Found files:\n{confidence_json}\n{cif_file}')

  cmd_args = ['Rscript', 'getContactsPAE.R', cif_file, confidence_json]
  cmd = ' '.join(cmd_args)
  print("\n")
  print(cmd)
  print()

  try:
    subprocess.run(cmd,
        stdout = sys.stdout, 
        stderr = sys.stderr,
        shell = True,  # module command is a csh alias on Wynton
        executable = '/bin/csh',
        check = True)
  except Exception as error:
    print('Error executing getContactsPAE.R:\n {}'.format(error))

def plot_interChainDistances(outDir):

  cif_file =  glob.glob(outDir+'/'+'*_model.cif').pop()
  sample_name = os.path.basename(os.path.dirname(cif_file))
  confidence_json = glob.glob(outDir+'/'+sample_name+'_confidences.json').pop()

  print(F'Found files:\n{confidence_json}\n{cif_file}')

  cmd_args = ['Rscript', 'plotContactsPLDDT.R', cif_file, confidence_json]
  cmd = ' '.join(cmd_args)
  print("\n")
  print(cmd)
  print()

  try:
    subprocess.run(cmd,
        stdout = sys.stdout, 
        stderr = sys.stderr,
        shell = True,  # module command is a csh alias on Wynton
        executable = '/bin/csh',
        check = True)
  except Exception as error:
    print('Error executing plotContactsPLDDT.R:\n {}'.format(error))