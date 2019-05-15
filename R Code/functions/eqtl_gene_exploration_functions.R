#=====================================================#
# TITLE: eQTL functions from Pratyay Rudra's code     #
#        and other functions for eQTL analysis        #
# AUTHOR: Gordon Kordas                               #
# DATE: 09/10/2018                                    #
#######################################################

#-----------------------------------------------------#
# FUNCTION: return only genes for specfic miRNA 
#-----------------------------------------------------#

genesMirSpecific <- function(getMultimirObj,mirna){
  geneTable <- addmargins(table(getMultimirObj@data[,c(3,6)]))
  
  miRgenes <- geneTable[mirna,]
  miRgenes <- names(miRgenes[miRgenes != 0])
  
  return(miRgenes)
}



#========================================================#
# FUNCTION: extract miRNA-mRNA pairs with 5 or more 
#           databases used in target prediction 
#========================================================#


predPairs <- function(targets.predicted,numDb,exactDb = NULL){
  library(purrr)
  library(dplyr)
  
  if(is.null(exactDb) == TRUE){
    miRmRNApairs <- unique(targets.predicted@data[,c("mature_mirna_id","target_ensembl")])
  }
  
  # for specific database #
  if(is.null(exactDb) == FALSE){
    targetscan_db <- targets.predicted@data %>% filter(database == exactDb)
    miRmRNApairs <- unique(targetscan_db[,c("mature_mirna_id","target_ensembl")])
  }
  
  # loop over number of unique mirna and ensemble_id pairs #
  keepPairs <- list()
  
  for (ii in 1:nrow(miRmRNApairs)){
      
    # count miRNA and mRNA with predictions from 5 or more databases #
    if (is.null(exactDb) == FALSE){
      keepPairs[[ii]] <- c(miRmRNApairs[ii,1],miRmRNApairs[ii,2])
    }
    
    if(numDb > 1){
      
      numDatabases <- targets.predicted@data %>% 
        dplyr::filter(mature_mirna_id == miRmRNApairs[ii,1] & target_ensembl == miRmRNApairs[ii,2]) %>% 
        dplyr::count(database) %>% 
        nrow()
      
      if (is.null(exactDb) == TRUE & numDatabases >= numDb){
        keepPairs[[ii]] <- c(miRmRNApairs[ii,1],miRmRNApairs[ii,2])
      }
    }
    
  }
  
  return(purrr::compact(keepPairs))
}

#===============================================#
# FUNCTION: determine the mRNA qtl that fall 
#           into the miRNA qtl Credible Interval
#===============================================#

# does the mRNA qtl fall in the miRNA Bayes credible interval? #
ininterval <- function(mRNAtarget,miRNAqtl){
  
  library(dplyr)
  
  # make sure positions and intervals are numeric #
  mRNAtarget$qtlChr <- as.character(mRNAtarget$qtlChr)
  miRNAqtl$qtlChr <- as.character(miRNAqtl$qtlChr)
  
  mRNAtarget$qtlPos <- as.numeric(mRNAtarget$qtlPos)
  miRNAqtl$bayesLow <- as.numeric(miRNAqtl$bayesLow)
  miRNAqtl$bayesHi <- as.numeric(miRNAqtl$bayesHi)
  
  # merge if qtl position of mRNA is within bayes credible interval #
  merged <- mRNAtarget %>% left_join(miRNAqtl, by = "qtlChr") %>% 
    filter(qtlPos.x > bayesLow.y, bayesHi.y < qtlPos.x)
    
  return(merged)
}




#--------------------------------------------------------#
# FUNCTION: extract target information for specific miRNA
#           from valCIOverlap or predCIOverlap tables
#--------------------------------------------------------#

targetextract <- function(overlapTbl,mirna){
  invisible(library(dplyr))
  
  targets <- overlapTbl %>% dplyr::filter(miRNA == mirna) %>% dplyr::select(mRNA)
  names(targets) <- mirna
  return(targets)
}

