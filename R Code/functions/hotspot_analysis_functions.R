#########################################################
# TITLE: Functions used in the mi-eQTL hotspot analysis #
# AUTHOR: Gordon Kordas                                 #
# DATE: 05/07/2019                                      #
#########################################################



#=============================================#
# FUNCTION: Similiar to Brem method but 
#           output p-value to allow FDR
#           correction
#=============================================#

hotspots.brem.pval <- function(window = 3,genolocs,permpvals,mult_thresh = 0.05,
                          outputThresh = FALSE, correction = "BH"){  #need eQTL sig thresh and hotspot thresh
  # define window size,chromosome,genomic locatiion,p values, threshold
  
  sdps <- c(1:nrow(genolocs))
  
  # calculate number of windows necessary #
  locmin <- min(sdps)
  locmax <- max(sdps)
  num.window <- ceiling((locmax-locmin)/window)
  
  #-----------------------------------#
  # iterate over all possible windows #
  #-----------------------------------#
  
  num.sig <- numeric(num.window)
  
  for (i in 1:num.window){
    
    # define current window #
    locmin.curr <- locmin+(i-1)*window
    locmax.curr <- min(locmin+i*window, locmax)
    
    sdps.curr <- sdps[sdps >= locmin.curr & sdps < locmax.curr]
    
    # calculate number significant in each window #
    if (length(sdps.curr)>0){
      curr.win <- c(permpvals[sdps.curr]) # define window, used to be as.matrix(permpvals[,sdps.curr])
      
      # count number of significant eQTLs, one eQTL per miRNA
      num.sig[i] <- sum(curr.win)
    }
    else{
      num.sig[i] <- 0
    }
  }
  
  hotspots <- numeric(num.window) # variable to hold whether a hotspot
  hotspots_pvals <- numeric(num.window)
  hotspots_pvals_fdr <- numeric(num.window)
  tot.num.sig <- sum(num.sig) # for poisson distribution
  lambda <- tot.num.sig/num.window  # mean of the poisson distribution (# hotspots by random chance)
  
  for (i in 1:num.window){
    # Is the number of eQTL in a given window greater than we expect by random chance? #
    hotspots_pvals[i] <- ppois(num.sig[i] - 1,lambda = lambda,lower.tail = FALSE) #-1 since ppois returns P(X > x)
  }
  
  # bonferonni adjustment #
  hotspots_pvals_adj <- p.adjust(hotspots_pvals,method = correction)
  
  for (i in 1:num.window){
    hotspots[i] <- ifelse(hotspots_pvals_adj[i] < mult_thresh,1,0)
  }
  
  if(outputThresh == TRUE){
    return(qpois(((1 - threshhot)^(1/num.window)),lambda = lambda,lower.tail = TRUE))
  }
  else{
    # output hotspots and the window they reside in #
    return(list(num.sig = num.sig, hotspot = hotspots, windowstarts = locmin +(1:num.window-1)*window,
           hotspots_pvals = hotspots_pvals, hotspots_pvals_adj = hotspots_pvals_adj))
  }
}





#==========================================#
# FUNCTION: extract mirna in hotspots 
#==========================================#

mirna_in_hotspots <- function(index,hotmirna_df,sdpLocations){
  
  mirna_out <- list()
  for (ii in 1:nrow(hotmirna_df)){
    
    # get list of sdps from sdplocations for 4 sdp window #
    sdp_start <- hotmirna_df$sdpStart[[ii]]
    
    sdp_set <- c(sdp_start,
    sdpLocations$sdp[which(sdpLocations$sdp == sdp_start) + 1],
    sdpLocations$sdp[which(sdpLocations$sdp == sdp_start) + 2],
    sdpLocations$sdp[which(sdpLocations$sdp == sdp_start) + 3])
    
    mirna_out[[ii]] <- index %>% filter(sdp == sdp_set[1] | sdp == sdp_set[2] | sdp == sdp_set[3] | sdp == sdp_set[4]) %>%
      dplyr::select(mirna) %>% unique() %>% unlist %>% unname
  }
  return(mirna_out)
}
