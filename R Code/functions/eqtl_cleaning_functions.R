#=====================================================#
# TITLE: eQTL functions from Pratyay Rudra's code     #
#        and other functions for eQTL data cleaning   #
# AUTHOR: Gordon Kordas                               #
# DATE: 09/13/2018                                    #
#######################################################




#===========================================================================#
# Function to organize the genotype data in terms of major or minor alleles #
# and generate sdps                                                         #
#===========================================================================#

geno.alleles=function(genomat){  #takes a genotype matrix with alleles A,T,G,C and 
  #gives the number of minor alleles
  # alleles must be coded as characters
  `%nin%` = Negate(`%in%`)
  
  genomat[genomat %nin% c("A","T","G","C")]=NA
  
  for (i in 1:dim(genomat)[1]){
    tabobj=table(genomat[i,])
    
    minorind=which.min(tabobj)
    minor=grep(names(tabobj)[minorind],genomat[i,])
    major=grep(names(tabobj)[3-minorind],genomat[i,])
    
    genomat[i,minor]=2   #MINOR is 2
    genomat[i,major]=0
  }
  
  mode(genomat)="numeric"
  
  return(genomat)
}


#===================================================#
# return minor allele to match with parent strains
#===================================================#

geno.minor <- function(genomat){ #takes a genotype matrix with alleles A,T,G,C and 
  #gives the number of minor alleles
  # alleles must be coded as characters
  `%nin%` = Negate(`%in%`)
  
  genomat[genomat %nin% c("A","T","G","C")]=NA
  
  minor <- length(1:dim(genomat)[1])
  
  for (i in 1:dim(genomat)[1]){
    tabobj=table(genomat[i,])
    
    minor[i]<- names(which.min(tabobj))
  }
  
  return(minor)
}


#======================================================================#
# FUNCTION to average over miRNA replicates and return avg expression  #
#          per strain                                                  #
#======================================================================#

avgStrain <- function(miRdata,columnNames){
  
  mRNAstring <- colnames(miRdata)
  
  #-----------------------------------#
  #FUNCTION: remove characters before
  #          second _
  #-----------------------------------#
  subMiRNAstr <- function(miRNAstring){
    gsub(".*_","",miRNAstring)
  }
  
  newmiRNAnames <- sapply(mRNAstring,FUN = subMiRNAstr)
  names(newmiRNAnames) <- NULL
  
  # unique names #
  uniqueNames <- unique(newmiRNAnames)
  

  # select non unique columns and average over them #
  
  avg.ls <- list()
  
  for (i in 1:length(uniqueNames)){
    avgCol <- rowMeans(select(as.data.frame.matrix(miRdata),contains(uniqueNames[[i]])))
    avg.ls[[i]] <- avgCol
  }
  
  # covert list to data frame #
  miRNAavg.df <- data.frame(Reduce(cbind,avg.ls))
  names(miRNAavg.df) <- uniqueNames # add colnames 
  
  return(miRNAavg.df)
}


#===================================================================#
#FUNCTION: create clean miRNA names (LXS data specific)
#===================================================================#

lxsCleanMirNames <- function(miRcolnames,miRNAavg.df){
  
  #-----------------------------------#
  #FUNCTION: remove characters before
  #          ':'
  #-----------------------------------#
  remBeforeSemiC <- function(miRcolnames){
    gsub(".*:","",miRcolnames)
  }
  
  newMmuNames <- sapply(rownames(miRNAavg.df),FUN = remBeforeSemiC)
  
  
  
  # remove extraneous info on miRNA.df rownames #
  cleanMiRnames <- sapply(as.character(newMmuNames),function(x){sub(".*mmu", "", x)}) 
  names(cleanMiRnames) <- NULL
  
  # fix double novel miRNA #
  cleanMiRnames[644] <- "chr14_37296"
  
  # add back mmu #
  cleanMiRnames2 <- sapply(cleanMiRnames,function(x){sub("-miR", "mmu-miR", x)})
  names(cleanMiRnames2) <- NULL
  
  cleanMiRnames2 <- sapply(cleanMiRnames2,function(x){sub("-let", "mmu-let", x)}) 
  names(cleanMiRnames2) <- NULL
  
  cleanMiRnames2 <- sapply(cleanMiRnames2,function(x){sub("chr", "novel:chr", x)}) 
  names(cleanMiRnames2) <- NULL
  
  # return clean rownames #
  return(cleanMiRnames2)
  
  
}




#==============================================================================#
# FUNCTION: prep miRNA and mRNA data for R qtl package
#==============================================================================#

qtlPackagePrep <- function(expressionData){
  
  # transpose data #
  datat <- as.data.frame(t(expressionData))
  
  # add id column #
  id <- as.character(seq(1:ncol(expressionData)))
  datat <- cbind(datat,id)
  rownames(datat) <- NULL
  
  return(datat)
}




#==============================================================================#
# FUNCTION: prep genotype data for R qtl package
#==============================================================================#


genotypeqtlPackagePrep <- function(sdpData,minorAllele,majorAllele,sdpLocations){
  
  # convert 0,1 to BB,AA #
  genotypeAB.df <- sdpData
  
  genotypeAB.df[sdpData == minorAllele] <- "BB"
  genotypeAB.df[sdpData == majorAllele] <- "AA"
  
  # combine location and SDP #
  genotypeAB.df <- cbind(sdpLocations,genotypeAB.df)
  
  # transpose miRNA data #
  genotypeABt.df <- as.data.frame(t(genotypeAB.df))
  
  # create 'id' column of strains #
  
  id <- as.character(c(NA,NA,seq(1:ncol(sdpData))))
  genotypeABt.df <- cbind(id,genotypeABt.df)
  
  rownames(genotypeABt.df) <- NULL
  
  return(genotypeABt.df)
  
}




####################################### For miRNA locations ######################################

#=================================================================================#
# FUNCTION: identify most likely location for multi position miRNA expression
#=================================================================================#

localQTLidentify <- function(pvalueMatrix,dupNames,dupPositions,sdpPositions){ 
  #--------------------------------------------------------#
  # dupPositions must include miRNA,chromosome,bp colnames #
  #--------------------------------------------------------#
  
  # combine location info to p value matrix #
  pval.mat <- as.data.frame(cbind(sdpPositions,as.matrix(t(pvalueMatrix))))
  
  #-----------------------------------------------------------#
  # FUNCTION: locate sdps with same chromosome and 
  #           within 5 Mb on either side (Nica et. al.) 
  #-----------------------------------------------------------#
  locateSDPs <- function(pvalmatrix,locOfInt,dupName){
    
    sdpOfinterest <- pvalmatrix[pvalmatrix$chromoSDP == locOfInt[1] & pvalmatrix$medianBpSDP < (as.numeric(locOfInt[2])+5E6)
                                & pvalmatrix$medianBpSDP > (as.numeric(locOfInt[2])-5E6),dupName]
    
    return(sdpOfinterest)   
  }
  #-----------------------------------------------------------#
  
  # iterate over all miRNAs with duplicate locations #
  sdp.interest <- list()
  
  for (ii in 1:length(dupNames)){ #change back to length(dupNames)
    
    # return positions of miRNA of interest #
    locsOfinterest <- dupPositions[dupNames[ii] == dupPositions$miRNA,c(2,3)]
    
    rownames(locsOfinterest) <- paste(unlist(unname(locsOfinterest[1])),unlist(unname(locsOfinterest[2])),sep = "-")
    
    # apply locateSDPs function #
    sdp.interest[[ii]] <- apply(locsOfinterest,1,
                                locateSDPs,pvalmatrix = pval.mat,dupName = as.character(dupNames[ii]))
    
  }
  
  return(sdp.interest)
}
#===============================================================================================#





#=================================================================================#
# FUNCTION: calculate minimum p value while dealing with miRNA in close proximity
#           and miRNA with no local QTL via distance to strongest SDP/SNP
#=================================================================================#

minPvalMir <- function(pvalu,duplicatename,miRpvalueMat,sdpLocs,dupPositions){
  
  # add sdp locations to p value matrix #
  pvalsPlusSdpLoc <- cbind(sdpLocs,t(miRpvalueMat)) # have to transpose to get sdps as rows
  rownames(pvalsPlusSdpLoc) <- rownames(sdpLocs)
  
  #------------------------------------------#
  # if matrix 
  #------------------------------------------#
  if (is.matrix(pvalu) == TRUE){
    
    ### many extremely close locations ### 
    if(length(which(apply(pvalu,2,min) == min(apply(pvalu,2,min)))) != 1){
      
      miRlocsOfinterest <- names(which(apply(pvalu,2,min) == min(apply(pvalu,2,min))))
      
      # which chromosome are we examining? #
      chrOfint <- substr(miRlocsOfinterest[[1]],1,1)
      
      # find location of strongest SDP/SNP for that miRNA #
      sdpLocTemp <- which.min(pvalsPlusSdpLoc[pvalsPlusSdpLoc$chromoSDP %in% chrOfint,duplicatename,
                                              drop = F][[1]])
      
      
      sdpName <- rownames(pvalsPlusSdpLoc[pvalsPlusSdpLoc$chromoSDP %in% chrOfint,duplicatename,
                                          drop = F])
      
      # our strongest associated SDP's location #
      strongSDP <- sdpName[sdpLocTemp]
      
      strongSDPloc <- sdpLocs[strongSDP,]
      
      # calculate miRNA location closest to strong SDP location #
      
      # match chromosome of strongest SDP #
      miRlocsOfinterest <- dupPositions[dupPositions$miRNA == duplicatename & 
                                          dupPositions$chromosome == strongSDPloc$chromoSDP,]
      
      # calculate difference in bp position #
      miRlocsOfinterest$bpDiff <- miRlocsOfinterest$bp - strongSDPloc$medianBpSDP
      
      
      rownames(miRlocsOfinterest) <-paste(unlist(unname(miRlocsOfinterest[2])),
                                          unlist(unname(miRlocsOfinterest[3])),sep = "-") 
      
      # return location of miRNA #
      bestLoc <- rownames(miRlocsOfinterest)[[which.min(abs(miRlocsOfinterest$bpDiff))]]
      
      return(bestLoc)
      
    }
    # if not extremely close locations #
    bestLoc <- names(which.min(apply(pvalu,2,min)))
    
    return(bestLoc)
  }
  
  ### if there is non numeric value ###
  if(length(pvalu) == 0){
    
    # which chromosomes are we examining? #
    nonNumericChr <- dupPositions[dupPositions$miRNA==duplicatename,"chromosome"]
    
    # calculate strongest SDP #  
    sdpLocTemp <- which.min(pvalsPlusSdpLoc[pvalsPlusSdpLoc$chromoSDP %in% nonNumericChr,duplicatename,
                                            drop = F][[1]])
    
    
    sdpName <- rownames(pvalsPlusSdpLoc[pvalsPlusSdpLoc$chromoSDP %in% nonNumericChr,duplicatename,
                                        drop = F])
    
    
    # our strongest associated SDP's location #
    strongSDP <- sdpName[sdpLocTemp]
    
    strongSDPloc <- sdpLocs[strongSDP,]
    
    # calculate miRNA location closest to strong SDP location #
    
    # match chromosome of strongest SDP #
    miRlocsOfinterest <- dupPositions[dupPositions$miRNA == duplicatename & 
                                        dupPositions$chromosome == strongSDPloc$chromoSDP,]
    
    # calculate difference in bp position #
    miRlocsOfinterest$bpDiff <- miRlocsOfinterest$bp - strongSDPloc$medianBpSDP
    
    rownames(miRlocsOfinterest) <-paste(unlist(unname(miRlocsOfinterest[2])),
                                        unlist(unname(miRlocsOfinterest[3])),sep = "-") 
    
    # return location of miRNA #
    bestLoc <- rownames(miRlocsOfinterest)[[which.min(abs(miRlocsOfinterest$bpDiff))]]
    
    return(bestLoc)
    
  }
  
  #-------------------------------------------------------------------#
  
  # calculate minimum p values for all locations #
  minLoca <- suppressWarnings(lapply(pvalu,min))
  
  ### many extremely close locations ### 
  
  if(length(which(minLoca == min(unlist(minLoca)))) != 1){
    
    miRlocsOfinterest <- names(which(minLoca == min(unlist(minLoca))))
    
    # which chromosome are we examining? #
    chrOfint <- substr(miRlocsOfinterest[[1]],1,1)
    
    # find location of strongest SDP/SNP for that miRNA #
    sdpLocTemp <- which.min(pvalsPlusSdpLoc[pvalsPlusSdpLoc$chromoSDP %in% chrOfint,duplicatename,
                                            drop = F][[1]])
    
    sdpName <- rownames(pvalsPlusSdpLoc[pvalsPlusSdpLoc$chromoSDP %in% chrOfint,duplicatename,
                                        drop = F])
    
    # our strongest associated SDP's location #
    strongSDP <- sdpName[sdpLocTemp]
    
    strongSDPloc <- sdpLocs[strongSDP,]
    
    # calculate miRNA location closest to strong SDP location #
    
    # match chromosome of strongest SDP #
    miRlocsOfinterest <- dupPositions[dupPositions$miRNA == duplicatename & 
                                        dupPositions$chromosome == strongSDPloc$chromoSDP,]
    
    # calculate difference in bp position #
    miRlocsOfinterest$bpDiff <- miRlocsOfinterest$bp - strongSDPloc$medianBpSDP
    
    rownames(miRlocsOfinterest) <-paste(unlist(unname(miRlocsOfinterest[2])),
                                        unlist(unname(miRlocsOfinterest[3])),sep = "-") 
    
    # return location of miRNA #
    bestLoc <- rownames(miRlocsOfinterest)[[which.min(abs(miRlocsOfinterest$bpDiff))]]
    
    return(bestLoc)
    
  }
  
  else{
    # return location of miRNA #
    bestLoc <- names(minLoca[which.min(minLoca)])
    return(bestLoc)
  }
  
}

#=================================================================================================#

