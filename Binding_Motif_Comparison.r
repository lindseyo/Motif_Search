library(readr)
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


pfm <- read_csv("pfm_fungi.csv", col_names = FALSE)
orf1000 <- read_csv("orf_genomic_1000.csv", col_names = FALSE)
TFs <- list()
consSeqs <- list()
nucOrder <- c("A", "C", "G", "T")
consRatio <- .6
minConsSeqLen <- 6
for (rowIndex in 1:nrow(pfm)){
  if (grepl(">", pfm[rowIndex, 1])){
    TFs[length(TFs) + 1] <- pfm[rowIndex, 1]
    consSeq <- ""
    for (colIndex in 1:ncol(pfm)){
      #lapply?? can this be vectorized rather easily?
      for (nucIndex in 1:4){
        if (sum(as.integer(pfm[(rowIndex+1):(rowIndex+4), colIndex])) == 0){
          break()
        }else if ((as.integer(pfm[(rowIndex + nucIndex), colIndex])/sum(as.integer(pfm[(rowIndex+1):(rowIndex+4), colIndex]))) > consRatio){
          consSeq <- paste0(consSeq, nucOrder[nucIndex])
          break()
        }else if (nucIndex == 4){
          consSeq <- paste0(consSeq, "N")
          }
        }
    }
    consSeqs[length(consSeqs) + 1] <- consSeq
    }
}
maxRecSites <- list()
strongRecSites <- list()
strongTFIDs <- list()

for (recSiteIndex in 1:length(consSeqs)){
  recSiteVec <- unlist(strsplit(consSeqs[[recSiteIndex]][1], ""), recursive = FALSE)
  nonspecSites <- which(recSiteVec == "N")
  if (length(nonspecSites) == 0){
    strongRecSites[length(strongRecSites) + 1] <- consSeqs[[recSiteIndex]][1]
    strongTFIDs[length(strongTFIDs) + 1] <- TFs[[recSiteIndex]][1]
    next()
  }
  maxSiteDistRes <- maxSiteDist(nonspecSites, consSeqs[[recSiteIndex]][1])
  maxSpecificity <- maxSiteDistRes[1]
  maxRecSite <- maxSiteDistRes[2]
  if (maxSpecificity >= minConsSeqLen){
    strongRecSites[length(strongRecSites) + 1] <- consSeqs[[recSiteIndex]][1]
    maxRecSites[length(maxRecSites) + 1] <- maxRecSite
    strongTFIDs[length(strongTFIDs) + 1] <- TFs[[recSiteIndex]][1]
  }
}


#match to genome
#vectorize???
possibleTFs <- list()
geneIDs <- list()
for (geneORFInd in seq_len(nrow(orf1000))){
  if (charmatch(">", orf1000[geneORFInd,1], nomatch = FALSE)){
    promoter <- paste0(orf1000[(geneORFInd+1):(geneORFInd + 14), 1], collapse = "")
    geneIDs[length(geneIDs) + 1] <- orf1000[geneORFInd,1]
    geneTFs <- list()
    for (strongRecSiteInd in seq_len(length(maxRecSites))){
      if (charmatch(maxRecSites[[strongRecSiteInd]][1], promoter, nomatch = FALSE)){
        geneTFs[length(geneTFs) + 1] <- strongTFIDs[strongRecSiteInd]
      }
    }
    if (length(geneTFs) > 0){
      print(geneTFs)
      possibleTFs[length(possibleTFs) + 1] <- geneTFs
    }else{possibleTFs[length(possibleTFs) + 1] <- list(NA)}
    
  }
}

genesWithTFs <- possibleTFs[which(!is.na(possibleTFs))]
genesWithTFsIDs <- geneIDs[which(!is.na(possibleTFs))]
