library(readr)
library(entropy)
recSiteSubStr <- function(consSite1, startNuc, stopNuc){
  #subset the given area of the consensus recognition site, used to find the largest string
  #without N's
  return(substr(consSite1, (startNuc + 1), (stopNuc - 1)))
}

maxSiteDist <- function(siteVec, consSite){
  #find the max distance between two N's in the consensus recognition site
  maxDistance <- 0
  recSite <- ""
  for (siteInd in seq_len(length(siteVec)-1)){
    distance <- (siteVec[siteInd + 1] - (siteVec[siteInd] + 1)) 
    if (distance > maxDistance){
      maxDistance <- distance
      recSite <- recSiteSubStr(consSite, siteVec[siteInd], siteVec[siteInd + 1])
    }
  }
  return(c(maxDistance, recSite))
}

##distance formula helper function
distForm <- function(distancesVec){
  distsSq <- distancesVec^2
  return(sqrt(sum(distsSq)))
}

pfm <- read_csv("pfm_fungi.csv", col_names = FALSE)
orf1000 <- read_csv("orf_genomic_1000.csv", col_names = FALSE)
TFs <- list()
posWeightMatList <- list()
#consSeqs <- list()
nucOrder <- c("A", "C", "G", "T")
#consRatio <- .6
#minConsSeqLen <- 6
for (rowIndex in 1:nrow(pfm)){
  if (grepl(">", pfm[rowIndex, 1])){
    TFs[length(TFs) + 1] <- pfm[rowIndex, 1]
    posWeightMat <- data.matrix(pfm[(rowIndex+1):(rowIndex+4),])
    emptyCols <- which(apply(posWeightMat, 2, sum)==0)
    if (length(emptyCols) > 0){
      posWeightMatTrunc <- posWeightMat[, 1:(emptyCols[1]-1)]
    }else{posWeightMatTrunc <- posWeightMat}
    posWeightMatList[[length(posWeightMatList)+1]] <- posWeightMatTrunc
  }
}  


#create genome freq matrices
#possibleTFs <- list()
promoterMatList <- list()
geneIDs <- list()
for (geneORFInd in seq_len(nrow(orf1000))){
  if (charmatch(">", orf1000[geneORFInd,1], nomatch = FALSE)){
    promoter <- paste0(orf1000[(geneORFInd+1):(geneORFInd + 18), 1], collapse = "")
    promoter <- substr(promoter, 1, 1001)
    promoterNucIndices <- as.numeric(sapply(strsplit(promoter, split = ""), match, nucOrder))
    geneIDs[length(geneIDs) + 1] <- orf1000[geneORFInd,1]
    promoterMat <- matrix(data = 0, nrow = 4, ncol = length(promoterNucIndices))
    #vectorize???
    for (nucInd in seq_len(length(promoterNucIndices))){
      promoterMat[promoterNucIndices[nucInd], nucInd] <- 1
    }
    promoterMatList[[length(promoterMatList) + 1]] <- promoterMat
  }
}

start.time <- proc.time()
#multicore
#calc KL divergences using sliding windows
allKLDists <- list()
for (promoterMatInd in seq_len(length(promoterMatList))){
  promoterMat1 <- promoterMatList[[promoterMatInd]]
  promMinKLDistList <- list()
  ##parallelize at tf leveltime
  for (tfpwmInd in seq_len(length(posWeightMatList))){
    tfpwm <- posWeightMatList[[tfpwmInd]]
    recSiteLen <- ncol(tfpwm)
    #use parlapply, don't need 2d apply
    #does not give leniency for pyrimidines or whatever
    klDists <- lapply(seq_len((ncol(promoterMat1)-recSiteLen)), function(nucInd1) 
      distForm(mapply(function(x,y) 
        KL.plugin(promoterMat[,x], tfpwm[,y]), seq(from=nucInd1, to=(nucInd1 + recSiteLen)), seq_len(recSiteLen))))
    promMinKLDistList[[tfpwmInd]] <- min(unlist(klDists, recursive = FALSE))
  }
  allKLDists[[promoterMatInd]] <- promMinKLDistList
}

time.elapsed <- proc.time() - start.time