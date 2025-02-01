

# NB modules
import numpy as np
import glob
from collections import defaultdict

# minor adaptions from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py?ref_type=heads

def parse_atm_record(line):
    '''Get the atm record
    '''
    line = line.split()
    record = defaultdict()
    record['name'] = line[0]
    record['atm_no'] = int(line[1])
    record['atm_name'] = line[2].strip()
    record['atm_alt'] = line[3].strip()
    record['res_name'] = line[5].strip()
    record['chain'] = line[6]
    record['chain_no'] = int(line[7])
    record['resid'] = line[8]
    record['x'] = float(line[10])
    record['y'] = float(line[11])
    record['z'] = float(line[12])
    record['occ'] = float(line[13])
    record['B'] = float(line[14])

    return record

def read_cif(ciffile):
  '''Read a cif file predicted with AF and rewritten to contain all chains; extract plddt and coordinate scores
     Adapted from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py?ref_type=heads to work with cif format
  '''
  chain_coords, chain_plddt = {}, {}
  with open(ciffile, 'r') as file:
      for line in file:
          if not line.startswith('ATOM'):
              continue
          record = parse_atm_record(line)
          #Get CB - CA for GLY
          if record['atm_alt']=='CB' or (record['atm_alt']=='CA' and record['res_name']=='GLY'):
              if record['chain'] in [*chain_coords.keys()]:
                  chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                  chain_plddt[record['chain']].append(record['B'])
              else:
                  chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                  chain_plddt[record['chain']] = [record['B']]

    #Convert to arrays
  for chain in chain_coords:
      chain_coords[chain] = np.array(chain_coords[chain])
      chain_plddt[chain] = np.array(chain_plddt[chain])

  return chain_coords, chain_plddt


def calc_pdockq(chain_coords, chain_plddt,t):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    '''

    #Get coords and plddt per chain
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    #Calc 2-norm
    mat = np.append(coords1, coords2, axis=0)
    a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1,l1:] #upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists<=t)

    if contacts.shape[0]<1:
        pdockq=0
        ppv=0
    else:
        #Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
        #Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt*np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

        #PPV
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,
            0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
            0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,
            0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
            0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
            0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
            0.63555449, 0.55890174])

        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
            0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,
            0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,
            0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,
            0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
            0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
            0.06968505, 0.02946438])
        inds = np.argwhere(pdockq_thresholds>=pdockq)
        if len(inds)>0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]

    return pdockq, ppv

def get_model_pDockq(cifPath, contactThreshold=4):
    '''
    Use above function to capture the pDockq score and append to the output summary csv
    '''
    chain_coords,chain_plddt=read_cif(cifPath)
    pdockq,ppv=calc_pdockq(chain_coords, chain_plddt, contactThreshold)
    return(pdockq, ppv)