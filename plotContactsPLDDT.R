library(data.table)
library(magrittr)
library(ggplot2)
library(jsonlite) # convert json to do
library(stringr)
#library(circlize)
library(bio3d) # open pdb
library(usedist)#package for working with distances

pdbFile = character(0)
confidenceJsonFile = character(0)

args = commandArgs(trailingOnly=TRUE)

pdbFile = grep("cif$|pdb$", args, value = TRUE)
stopifnot(`Pass only one pdb file` = length(pdbFile) <= 1)

confidenceJsonFile <-  grep("_confidences.json$", args, value = TRUE)
stopifnot(`Pass only one confidences json file` = length(confidenceJsonFile) <= 1)

# name the contactOut file 
distancesResOutFile <- sprintf("%s.distances.csv",
                             gsub("\\_confidences.json$", "", confidenceJsonFile))

sName <- gsub('\\_model.cif', '', pdbFile)

message("Samplename: ", sName)
message("Loading pdb/cif from ", pdbFile)
message("Loading cofidence.json with pae from ", confidenceJsonFile)
message("Writing inter-chain distances to ", distancesResOutFile)


### Functions

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


# minB = lower bound of atoms to include in distance, from both ends of comparison.
#         Assume plddt, where higher is better
interChainResDistance <- function(atoms, minBStart = 0, minBStop = 25){
  # all by all atom distance
  atomDistance <- dist(atoms[, .(x,y,z)])
  
  # build pairwise table of distances between chains
  # loop over chains, and compare all atoms to that chain, skipping those within that chain
  .oneChain <- function(curChain){
    curChain.idx <- which (atoms$chain == curChain & atoms$b > minBStop)
    atomdist <- atoms[chain != curChain & b > minBStart, 
                      .(chain, resno, eleno, b, otherChain = curChain, distance = min(usedist::dist_get(atomDistance, idx, curChain.idx))),
                      by = .(idx)] 
    residues <- atomdist[, .(distance = min(distance), bfactor = mean(b)), by = .(chain, resno, otherChain)]
    residues[]
  }
  
  resLong <- rbindlist(lapply(unique(atoms$chain), .oneChain))
  resLong[]
}

###

chainNames <- toupper(unlist(strsplit(sName,'__')))
chainIDmapper <- data.table(chainName=chainNames, chainID=LETTERS[1:length(chainNames)])

message('loading PDB file...')
atoms <- loadPDBAtoms(pdbFile)

message('Finding interchain distances...')
interchain.dt <- interChainResDistance(atoms=atoms)

message('mapping gene names to chainIDs...')
print(chainIDmapper)

interchain.dt[chainIDmapper, gene := i.chainName, on = c(chain = 'chainID')]
interchain.dt[chainIDmapper, otherGene := i.chainName, on = c(otherChain = 'chainID')]

message(paste0('Writing distances file to: ', distancesResOutFile))
fwrite(interchain.dt, file = distancesResOutFile)

message('plotting inter-chain distances')

g <- ggplot(interchain.dt, aes(x = resno, y = distance, color = bfactor)) + 
  geom_line(lwd = 1, alpha = 0.4) +
  geom_point(alpha = 1, stroke = NA) +
  #coord_cartesian(xlim = c(0,400)) +
  scale_y_log10(name = "distance (Angstroms)") +
  #coord_cartesian(ylim = c(.1,20)) +
  #scale_y_continuous( ) +
  geom_hline(yintercept = 4.0) + # considered close contact (< 4 angstroms)
  facet_grid(gene~otherGene, scales = "free", space = "free_x") +
  scale_color_gradientn("plDDT", limits = c(0,100), colors = c(red = "#FE0012",
                                                               orange = "#FFA42B",
                                                               yellow = "#FFFD42","#FFFD42",
                                                               palegreen = "palegreen2",
                                                               blue = "lightblue","lightblue",#"#6097E8",
                                                               darkBlue = "#001DF9"),
                        values = c(0, .5, 0.7,0.75,0.8,0.85, 0.9, 1.0)) +
  
  theme_bw() +
  scale_x_continuous(breaks = seq(0, max(interchain.dt$resno, na.rm = TRUE), by = 200)) +
  ggrepel::geom_text_repel(aes(label = resno))

png(paste0(sName,".plddt.distances.png"), width=10,height=8,units="in",res=1200)
g
dev.off()


