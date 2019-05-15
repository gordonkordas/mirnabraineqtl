#=====================================================#
# TITLE: eQTL functions from Pratyay Rudra's code     #
#        and other functions for eQTL analysis        #
# AUTHOR: Gordon Kordas                               #
# DATE: 05/07/2019                                    #
#######################################################


#====================================================================#
# FUNCTION: extract results from R qtl objects
#====================================================================#

extractFromRqtl <- function(scanoneAll, permsAll,alpha,numPheno,qtl.obj,pheno,specific = FALSE,specificMir = NULL){ #pheno = miRNA or mRNA
  invisible(library(tidyr))
  invisible(library(qtl))
  
  # extract results for all phenotypes #
  results <- as.data.frame(summary(scanoneAll,perms = permsAll,pvalues = TRUE,alpha = alpha,format = "allpheno"))
  
  ## extract p values ##
  results.p <- results[,grep("pval",colnames(results))]
  
  # if we want all phenos in qtl object #
  if (specific == FALSE){
    colnames(results.p) <- names(qtl.obj$pheno[1:numPheno])
  }
  else{colnames(results.p) <- specificMir}
  
  ## extract LOD scores ##
  results.lod <- results[,-c(1,2,grep("pval",colnames(results)))]
  
  # add sdp column #
  results.p$sdp <- rownames(results.p)
  results.lod$sdp <- rownames(results.lod)
  
  
  # convert to long data format #
  results.long <- gather(results.p,pheno,pvalue,1:(dim(results.p)[2] - 1),factor_key = TRUE) # generalize
  lod.long <- gather(results.lod,pheno,lod,1:(dim(results.p)[2] - 1),factor_key = TRUE) # generalize
  
  
  # add lod scores to long format data #
  results.long <- cbind(results.long,lod.long$lod)
  
  # only p values below alpha #
  results.sig <- results.long[results.long$pvalue < alpha,]
  colnames(results.sig) <- c("sdps","miRNA","pvalue","lod")
  
  # add location back #
  locs <- results[,c("chr","pos")]
  locs$sdps <- rownames(locs)
  results.final <- merge(results.sig,locs,by = "sdps")
  
  return(results.final)
  
}



#==================================================================#
# FUNCTION: calculate bayes intervals for R qtl objects            
#==================================================================#

bayesRqtls <- function(qtl.miR,numPheno,name = FALSE,phenoName = NULL){
  
  invisible(library(qtl))
  
  # loop over every miRNA and chromosome to obtain bayes intervals #
  bayes.ls <- list()
  bayes.chromo <- list()
  
  if(name == FALSE){
    for (ii in 1:numPheno){
      scan.out <- scanone(qtl.miR,pheno.col = ii,method = "mr")
      
      ### if 'Y' chromosome, else...
      for (jj in c(seq(1:19),'X')){
        bayes.temp <- bayesint(scan.out,chr = jj)
        
        bayes.chromo[[jj]] <- data.frame(chr = c(bayes.temp$chr[2]),start = c(bayes.temp$pos[1]),stop = bayes.temp$pos[3])
      }
      
      # combine all chromosomes #
      bayes.full <- Reduce(rbind,bayes.chromo)
      
      bayes.full[bayes.full$chr == 20,"chr"] <- "X"
      
      bayes.ls[[ii]] <- bayes.full
    }
    
    # name list of miRNA credible intervals #
    names(bayes.ls) <- names(qtl.miR$pheno[1:numPheno])
    
    return(bayes.ls)
  }
  else{
    for (ii in phenoName){
      scan.out <- scanone(qtl.miR,pheno.col = ii,method = "mr")
      
      ### if 'Y' chromosome, else...
      for (jj in c(seq(1:19),'X')){
        bayes.temp <- bayesint(scan.out,chr = jj)
        
        bayes.chromo[[jj]] <- data.frame(chr = c(bayes.temp$chr[2]),start = c(bayes.temp$pos[1]),stop = bayes.temp$pos[3])
      }
      
      # combine all chromosomes #
      bayes.full <- Reduce(rbind,bayes.chromo)
      
      bayes.full[bayes.full$chr == 20,"chr"] <- "X"
      
      bayes.ls[[ii]] <- bayes.full
    }
    
    # name list of miRNA credible intervals #
    names(bayes.ls) <- phenoName
    
    return(bayes.ls)
  }

}



#---------------------------------------------------------------#
# FUNCTION: merge qtls and bayes credible intervals 
#---------------------------------------------------------------#

bayesCImerge <- function(qtldata,bayesCIlist,trait,all = FALSE){ # trait is the column of interest for eQTL
  
  if (all == FALSE){
    sig.ci <- bayesCIlist[unique(qtldata[,trait])]
    allci.ls <- list()
    
    
    for (ii in names(sig.ci)){
      
      trait.df <- qtldata[qtldata[,trait] == ii,]
      trait.ci <- bayesCIlist[[ii]]
      
      # create chromosome factor for merge to be possible #
      trait.ci$chr <- factor(trait.ci$chr,levels = levels(qtldata$chr))
      
      merged <- merge(trait.df,trait.ci,by = "chr")
      
      # iterate over all chromosomes within unique miRNA #
      chromoHighLOD.ls <- list()
      
      for (jj in unique(merged$chr)){
        chromoHighLOD.ls[[jj]] <- merged %>% subset(chr == jj) %>% subset(lod == max(lod))
      }
      
      allci.ls[[ii]] <- Reduce(rbind,chromoHighLOD.ls)
      
    }
  }
  else{
    sig.ci <- bayesCIlist
    allci.ls <- list()
    
    for (ii in names(sig.ci)){
      
      trait.df <- qtldata[qtldata[,trait] == ii,]
      trait.ci <- bayesCIlist[[ii]]
      
      # create chromosome factor for merge to be possible #
      #trait.ci$chr <- factor(trait.ci$chr,levels = levels(qtldata$chr))
      
      merged <- merge(trait.df,trait.ci,by = "chr")
      
      allci.ls[[ii]] <- merged
    }
  }

  qtlAndCi.df <- Reduce(rbind,allci.ls)
  
  return(qtlAndCi.df)
}






#-----------------------------------------------------------------#
# FUNCTION: expand bayes credible interval
#           to width of corresponding snps
#
# input: qtlsubset: qtls with zero bayes interval due to large LOD
#        rsnumberLocs: sdp locations with rsnumbers
#        sdpSnpLocMatch: sdp's corresponding snps and locations
#-----------------------------------------------------------------#

expandToSnps <- function(qtlsubset,rsnumberLocs,sdpSnpLocMatch){

  # calculate min and max snp locations corresponding to sdp #
  sdpTosnp <- function(sdp,snplocs,minim){ # sdp to rsnum
    sub.rs <- snplocs[sdp == snplocs$match.keys,]
    
    # calc lower bound
    if(minim == TRUE){
      minimum <- min(sub.rs$bp.vec)
      return(minimum)
    }
    # calc upper bound
    else{
      maximum <- max(sub.rs$bp.vec)
      return(maximum)
    }
  }
  
  # apply sdpTosnp function to all sdps of interest #
  qtlsubset$cilow <- sapply(X = qtlsubset$sdp, FUN = sdpTosnp, snplocs = sdpSnpLocMatch, minim = TRUE)
  qtlsubset$cihi <- sapply(X = qtlsubset$sdp, FUN = sdpTosnp, snplocs = sdpSnpLocMatch, minim = FALSE)
  return(qtlsubset)
}




#=================================#
# FUNCTION: limit to one eQTL per
# miRNA by maximum LOD score
#=================================#

max_lod <- function(outallmir,name){
  
  outallmir <- as.data.frame(outallmir)
  mir_row <- list()
  
  for (ii in colnames(outallmir)[! colnames(outallmir) %in% c("chr","pos")]){
    
    mir_row[[ii]] <- outallmir[outallmir[ii] == max(outallmir[ii]),c("chr","pos",ii)]
    mir_row[[ii]]$mirna <- colnames(mir_row[[ii]])[3]
    mir_row[[ii]]$sdp <- rownames(mir_row[[ii]])
    colnames(mir_row[[ii]]) <- c("chr","pos","lod",name,"sdp")
  }
  
  max_lod_out <- Reduce(rbind,mir_row)
  return(max_lod_out)
}





#=================================#
# FUNCTION: exact p estimate
#=================================#

exact_p <- function(x,nperm){
  (x + 1)/(nperm + 1)
}




#---------------------------------------------#
# FUNCTION: cis vs. trans qtl
#---------------------------------------------#

cisVsTransQtl <- function(mirnachr,qtlchr,mirnaStartpos,qtlpos,bpwindow = 5E6){
  
  # can pass entire column of positions #
  # returns one if cis, zero if trans #
  cis <- ifelse((mirnachr == qtlchr) & 
                  qtlpos - bpwindow < mirnaStartpos &
                  mirnaStartpos < qtlpos + bpwindow, 1, 0)
  
  return(cis)  # vector
}






#-----------------------------------------------------#
# FUNCTION: Intronic vs intergenic
#-----------------------------------------------------#

intronic <- function(intronRegions,phenoChrs,phenobpLoc,phenoNames,pheno){
  
  if (!require(data.table, quietly=T)) {
    install.packages(data.table)
    library(data.table)
  }
  library(data.table)
  
  # create dataset for function #
  phenoDat <- as.data.frame(cbind(phenoNames,phenoChrs,phenobpLoc),stringsAsFactors = FALSE)
  phenoDat$phenoChrs[phenoDat$phenoChrs == "20"] <- "X"

  # function to apply: is location intronic? #
  inIntron <- function(phenosub,intronSub){
    
    tOrF <- sum(data.table::between(phenosub, intronSub$start,intronSub$end)) != 0
    return(tOrF)
  }
  
  # loop through chromosomes and miRNA locations (in qtl dataframe) #
  introns.ls <- list()
  
  for (jj in unique(intronRegions$chr[intronRegions$chr != "Y"])){ # no Y chromosome
    
    pheno.sub <- as.list(phenoDat[phenoDat$phenoChrs == jj,3])  # V2 = phenoChr
  
    intronSubset <- intronRegions[intronRegions$chr == jj,]
    
    introns.ls[[jj]] <- lapply(pheno.sub, FUN = inIntron,
                               intronSub = intronSubset) %>% unlist() %>% as.vector()
    
    if(!is.null(introns.ls[[jj]])){
      introns.ls[[jj]] <- cbind(phenoDat[phenoDat$phenoChrs == jj,1],introns.ls[[jj]])
    }
    
  }

  
  intron.col <- as.data.frame(Reduce(rbind,introns.ls))

  colnames(intron.col) <- c(pheno,"intron")
  intron.col$intron <- ifelse(intron.col$intron == TRUE,1,0)
  
  intron.col <- unique(intron.col) # should return 149

  return(intron.col) # vector
  
}





#===================================================================#
# FUNCTION: Extract all adjusted p values from R qtl package for use 
#           in the hotspot function
#===================================================================#

allPermPvals <- function(permMat,unadjLods){
  
  # function to calculate permutation p values #
  calcPermPval <- function(miRlod,miRperm){
    p.adj <- mean(as.numeric(miRperm) > as.numeric(miRlod))
    return(p.adj)
  }
  
  
  pAdj.ls <- list()
  
  # loop over all phenotypes (i.e. miRNA)
  for (ii in 1:ncol(permMat)) {  # ncol(permMat)
    
    miR.perm <- as.vector(permMat[,ii])
    miR.lod <- unadjLods[,ii]
    
    pAdj.ls[[ii]] <- sapply(miR.lod,FUN = calcPermPval,miRperm = miR.perm)
    
  }
  
  pAdj.mat <- Reduce(rbind,pAdj.ls)
  colnames(pAdj.mat) <- rownames(unadjLods)
  return(pAdj.mat)
}





#==========================================#
# FUNCTION: Boxplots for genotype data 
#==========================================#

boxplotsForDiffExpress <- function(mirna,sdp,mirnaData,sdpData, title = NULL,xlabel = "Parent Allele",spec_color,size){
  
  # extract miRNA expression across all strains #
  miRexpress <- as.numeric(c(mirnaData[mirna,]))
  names(miRexpress) <- NULL
  
  # extract sdp alleles row #
  sdpCol <- sdpData[sdp,]
  names(sdpCol) <- NULL
  
  boxPlotData <- data.frame(miRexpress,sdpCol)
  
  # run boxplot (base or ggplot will do) #
  ggplot(data = boxPlotData, aes(x = factor(sdpCol), y = miRexpress, fill = factor(sdpCol))) + 
    geom_boxplot() + 
    ggtitle(title) + xlab(xlabel) + ylab(paste(mirna,"Expression",sep = " ")) +
    scale_fill_manual(values=spec_color) +
    ylim(c(3,6)) +
    theme_bw() +
    theme(legend.position = 'none',axis.title=element_text(size=size,face="bold"),plot.title = element_text(hjust = 0.5))
  
}

