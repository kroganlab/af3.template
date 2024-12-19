
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
        