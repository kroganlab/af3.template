library(data.table)
library(magrittr)
library(ggplot2)
library(ComplexHeatmap)
library(rjson)
library(stringr)
library(Biostrings)


msaFile = character(0)
dataJsonFile = character(0)
confidenceJsonFile = character(0)

args = commandArgs(trailingOnly=TRUE)

msaFile <-  grep(".msaOut.txt$", args, value = TRUE)
stopifnot(`Pass only one msa output file` = length(msaFile) <= 1)

dataJsonFile <-  grep("_data.json$", args, value = TRUE)
stopifnot(`Pass only one data json file` = length(dataJsonFile) <= 1)

confidenceJsonFile <-  grep("_confidences.json$", args, value = TRUE)
stopifnot(`Pass only one confidences json file` = length(confidenceJsonFile) <= 1)

outDir <-  dirname(msaFile)

message("Loading msa from ", msaFile)
message("Loading input features from ", dataJsonFile)
message("Loading model scores from", confidenceJsonFile)
message("Writing plots to ", outDir)

#######
# functions necessary to run the Script...
# may need to install these pkgs in my local R env
#########

getQualitativePalette <- function(n) {
  
  col.pal <- c(
    "dodgerblue2", 
    "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black",
    "gold1",
    "skyblue2", 
    "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#121111", # lt orange
    "gray70", 
    "khaki2",
    "maroon", 
    "orchid1", 
    "deeppink1",
    "blue1", 
    "steelblue4",
    "darkturquoise",
    "green1", 
    "yellow4", 
    "yellow3",
    "darkorange4",
    "brown"
  )
  return(col.pal[1:n])
}

#https://www.bioinformatics.nl/~berndb/aacolour.html Cinema
resColorsMulti = c(HKR = "#00FFFF",
                   DE  = "#FF0000",
                   STNQ = "#00FF00",
                   AVLIM = "#BBBBBB",
                   FWY = "#FF00FF",
                   PG = "#996000",
                   C ="#FFFF00",
                   BZX = "grey25",
                   "-" = "grey95")


# python residue feature mapping
# AF3 feature encoding aplahfold/src/data/msa_features.py
# seems to be some redundancy in the code...
#   # AF3 _PROTEIN_TO_ID = {
#     'A': 0,
#     'B': 3,  # Same as D.
#     'C': 4,
#     'D': 3,
#     'E': 6,
#     'F': 13,
#     'G': 7,
#     'H': 8,
#     'I': 9,
#     'J': 20,  # Same as unknown (X).
#     'K': 11,
#     'L': 10,
#     'M': 12,
#     'N': 2,
#     'O': 20,  # Same as unknown (X).
#     'P': 14,
#     'Q': 5,
#     'R': 1,
#     'S': 15,
#     'T': 16,
#     'U': 4,  # Same as C.
#     'V': 19,
#     'W': 17,
#     'X': 20,
#     'Y': 18,
#     'Z': 6,  # Same as E.
#     '-': 21,
# }

readA3Mmatrix <-  function(a3mPath, sep=' '){
  
  message('Reading in a3m matrix from', a3mPath)
  a3m.mat <- fread(a3mPath, sep=sep) %>%  
    as.matrix()
  message('matrix dimensions: ' , paste(dim(a3m.mat), collapse=','))
  return(a3m.mat)
}


convertEncodingToResidueMatrix <-  function(a3mMat){
  
  # AF2 alphafold/common/residue_constants.py follows same encoding as AF3 aplahfold/src/data/msa_features.py
  # seems to be some redundancy in the code...
  restypes <-  c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-')
  
  #add  1 as python indexing
  a3mMat <- a3mMat + 1
  message('Converting numeric encoding to aa.residues...')
  
  .rmNA <- function(x)x[!is.na(x)]
  .singleSeq <- function(intCodes){
    paste0(restypes[.rmNA(intCodes)], collapse = "")
  }
  
  #  get chr str of aa sequence
  resStr <- apply(a3mMat, 1, .singleSeq)
  # convert to chr vector
  resVec <- lapply(resStr, function(x) {unlist(strsplit(x, ''))} )
  resMat <- do.call(rbind, resVec) %>% 
    as.matrix()
  
  return(resMat)
}


removeMSAPadding <- function(dataJson, resMat){
  
  message('Removing padding from MSA...\n')
  
  inputSeq <- lapply(seq_along(1:length(prot.names)), function(x){data.json$sequences[[x]]$protein$sequence}) %>% 
    paste0(collapse='')
  
  message('AF3 Query Sequence:\n', inputSeq)
  
  #get the top sequence from the MSA
  msaSeq <- paste0(resMat[1,],collapse='')
  
  match_idx <- matchPattern(inputSeq, msaSeq, max.mismatch = 0) %>% 
    as.data.table()
  
  # get the final indx, and get the residue
  return(resMat[,1:match_idx$end])
}

message('Reading:\n', confidenceJsonFile,'\n',dataJsonFile)
confidences.json <- fromJSON(file=confidenceJsonFile)
data.json <- fromJSON(file=dataJsonFile)


message('Reading:\n', msaFile)
a3m.mat <- readA3Mmatrix(a3mPath=msaFile, sep=' ')

# residue mat
res.mat <- convertEncodingToResidueMatrix(a3m.mat)

message('getting protein chains')
# assuming first element is A, second is B etc...
prot.names <- strsplit(data.json$name, '__')[[1]]

message('Extracting PAE matrix...')
# rbind pae vec to recover the matrix
pae.mat <- do.call(rbind, confidences.json$pae) %>% 
  as.matrix()

message('Generating color palette for PAE..')
chain.col.dt <- data.table(chainID = LETTERS[seq(1, length(prot.names))],
                           color = getQualitativePalette(n=length(prot.names)) )

chain.col.dt <- merge(chain.col.dt, data.table(chainID=unlist(confidences.json['token_chain_ids'])), by='chainID')

chainColVec <- chain.col.dt$color
names(chainColVec) <-  chain.col.dt$chainID

column_ha = HeatmapAnnotation(chainID=chain.col.dt$chainID, 
                              col=list(chainID=chainColVec))
row_ha = rowAnnotation(chainID = chain.col.dt$chainID, 
                       col=list(chainID=chainColVec))


message('Generating PAE plot..')
hm <- Heatmap(pae.mat, 
              cluster_columns = F,
              cluster_rows = F,
              name = data.json$name, 
              top_annotation = column_ha, 
              right_annotation = row_ha)


png(paste0(outDir,"/", data.json$name,".pae.png"),width=10,height=8,units="in",res=1200)
draw(hm, column_title=data.json$name)
dev.off()

message('Removing MSA padding..')
tidy.mat <- removeMSAPadding(dataJson = data.json, resMat = res.mat)

# color palette for aa residues
cdt <- data.table(resids = names(resColorsMulti), color = resColorsMulti)
cdt <- cdt[, .(resid = unlist(strsplit(resids, ""))), by= color]
residue.col <- cdt$color %>% 
  setNames(cdt$resid)


message('Generating MSA plot..')
hm <- Heatmap(tidy.mat, 
              border=T,
              cluster_rows = FALSE, 
              cluster_columns = FALSE, 
              row_title = sprintf('%s MSA records', nrow(tidy.mat)),
              col = residue.col, 
              top_annotation = column_ha, 
              column_split = chain.col.dt$chainID,
              name = "AAtype",
              column_names_side = "bottom")


png(paste0(outDir,"/", data.json$name,".msa.png"), width=10,height=8,units="in",res=1200)
draw(hm, column_title=data.json$name)
dev.off()


message('Done generating MSA and PAE plots..\nExiting...')