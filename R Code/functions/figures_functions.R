############################################
# TITLE:  Functions for miRNA paper graphics
#         generation
# AUTHOR: Gordon Kordas
# DATE:   02/20/2019
############################################


#==========================================#
# FUNCTION: scale heatmap to physical 
#           distance
#==========================================#

heatmap_scale <- function(sdplocs,mirnalocs,lodmat,chromosomes,size_mb,mbscale){
  
  #---------------------------------#
  # create new empty, scaled matrix #
  #---------------------------------#
  scale_df <- as.data.frame(cbind(chromosomes,size_mb),stringsAsFactors = F)
  scale_df$size_mb <- as.numeric(size_mb)
  scale_df$size_mb1 <- scale_df$size_mb*mbscale
  
  # convert to "long" format #
  new_chr <- c()
  new_mb1 <- c()
  
  for (ii in 1:nrow(scale_df)){
    chr_rep <- rep(scale_df[ii,"chromosomes"],times = scale_df[ii,"size_mb1"])
    mb_rep <- 1:length(chr_rep)
    new_chr <- c(new_chr,chr_rep)
    new_mb1 <- c(new_mb1, mb_rep)
  }
  axis_scale <- as.data.frame(cbind(new_chr,new_mb1),stringsAsFactors = FALSE)
  axis_scale$new_mb1 <- as.numeric(axis_scale$new_mb1)
  
  # create empty matrix with location columns #
  new_axis_mat <- as.data.frame(matrix(data = 0.00001,ncol = nrow(axis_scale), nrow = nrow(axis_scale)))
  new_axis_mat$chr <- axis_scale$new_chr
  new_axis_mat$dist <- axis_scale$new_mb1
  new_axis_mat$uniq_id <- 1:nrow(new_axis_mat)
  
  #--------------------#
  # prep sdp locations #
  #--------------------#
  sdplocs$mbs <- round(sdplocs$medianBpSDP/1000000)
  
  #----------------------#
  # prep mirna locations #
  #----------------------#
  mirnalocs$mbs <- round(mirnalocs$mir.start/1000000)
  
  #----------------------#
  # iterate over mirna
  #----------------------#
  empty_row <- rep(0.0001,nrow(new_axis_mat))
  
  for (ii in 1:2){ #ncol(lodmat)
    mir_col <- lodmat[,ii, drop = FALSE]
    sdp_names_tmp <- rownames(mir_col)
    locs_tmp <- sdplocs[sdplocs$sdp %in% sdp_names_tmp,c("chromoSDP","mbs")]
    indx <- inner_join(locs_tmp,new_axis_mat,by = c("chromoSDP" = "chr","mbs" = "dist")) %>% dplyr::select("uniq_id")
    indx <- c(indx[,1])

    # fill empty row with lods in correct (scaled) place #
    empty_row[indx] <- lodmat[,ii]
    
    return(empty_row)
  }
  
  #----------------------#
  # extend sdp locations #
  #----------------------#
  
  #---------------------------------------------------#
  # merge extend sdp location rows by mirna locations #
  #---------------------------------------------------#
}







