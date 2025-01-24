library(data.table)
library(magrittr)
library(rjson)
library(bio3d)

pdbFile = character(0)
confidenceJsonFile = character(0)
contactResOutFile = character(0) # not specified picked up from confidence file

args = commandArgs(trailingOnly=TRUE)

pdbFile = grep("cif$|pdb$", args, value = TRUE)
stopifnot(`Pass only one pdb file` = length(pdbFile) <= 1)

confidenceJsonFile <-  grep("_confidences.json$", args, value = TRUE)
stopifnot(`Pass only one confidences json file` = length(confidenceJsonFile) <= 1)

# name the contactOut file 
contactResOutFile <- sprintf("%s.contacts.csv",
                            gsub("\\_confidences.json$", "", confidenceJsonFile))


message("Loading pdb/cif from ", pdbFile)
message("Loading cofidence.json with pae from ", confidenceJsonFile)
message("Writing contacts to ", contactResOutFile)

## functions 

loadPDBAtoms <- function(path){
  
  if (grep('[.]cif', path)) {
    model = bio3d::read.cif(path) 
  } else {
    model = bio3d::read.pdb(path) 
  } 
  atoms <- setDT(model$atom)
  atoms[, idx := .I]
  return (atoms[])
}

interChainContacts <- function (pdbFile){
  atoms <- loadPDBAtoms(pdbFile)
  
  # all by all atom distance
  message (sprintf("All by all atom distance for %d atoms", nrow(atoms)))
  atomDistance <- as.matrix(dist(atoms[, .(x,y,z)]))
  
  # sweep out row and col radii
  message (sprintf("Done with distance, now calculating contacts"))
  
  # radii copied from Jason Nomburg, line 74 at
  # https://github.com/jnoms/SAT/blob/main/sat/scripts/struc_detect_interaction.py
  vdw.radii <- c(H =  1.2, C =  1.7, N =  1.55, O =  1.52, S =  1.8)
  atomDistance <- sweep (atomDistance, 1, vdw.radii[atoms$elesy])
  atomDistance <- sweep (atomDistance, 2, vdw.radii[atoms$elesy])
  
  # if remaining distance is still less than 0.5, we declare it a contact 
  contactIndeces <- which (atomDistance < 0.5, arr.ind = TRUE) |> as.data.table()
  
  # label with chains from idx in atoms table
  contactIndeces[atoms, chainRow := i.chain , on =  c(row = "idx")]
  contactIndeces[atoms, chainCol := i.chain , on =  c(col = "idx")]
  
  # make crosschain only, and only in one direction:
  contactIndeces <- contactIndeces[chainRow < chainCol]
  
  # label with resno from atoms table
  contactIndeces[atoms, resnoRow := i.resno, on = c(row = "idx")]
  contactIndeces[atoms, resnoCol := i.resno, on = c(col = "idx")]
  
  # collapse from atoms to residues and sort
  contactRes <- setorder(unique(contactIndeces[, .(chainRow, resnoRow, chainCol, resnoCol)]))
  
  # translate per-chain resno to the multimer resno based on chain lengths (max(resno)) for all prior chains
  # assumptions!!!
  cl <- atoms[, .(l = max(resno)), by = chain]
  contactRes[, mmerResnoRow := resnoRow + sum(cl[chain < chainRow, sum(l)]), by = chainRow]
  contactRes[, mmerResnoCol := resnoCol + sum(cl[chain < chainCol, sum(l)]), by = chainCol]
  
  return (contactRes[])  
}

contactRes <- interChainContacts(pdbFile)

confidences.json <- fromJSON(file=confidenceJsonFile)
message('Extracting PAE matrix from confidences.json...')
# rbind pae vec to recover the matrix
pae <- do.call(rbind, confidences.json$pae) %>% 
  as.matrix()

contactRes[, pae := pae[mmerResnoRow, mmerResnoCol], by = .(mmerResnoRow, mmerResnoCol)]

message(paste0('Writing contacts file to: ', contactResOutFile))
fwrite(contactRes, file = contactResOutFile)