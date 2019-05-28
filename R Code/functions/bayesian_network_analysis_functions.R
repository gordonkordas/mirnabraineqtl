#=====================================================#
# TITLE: functions for Bayesian Network Analysis of   #
#        miRNA expression                             #
# AUTHOR: Gordon Kordas                               #
# DATE: 05/09/2019                                    #
#######################################################



#=================================================#
# FUNCTION: use boot package and own function to 
#           determine indirect effect and proportion
#           mediated
#=================================================#

mediation_boot <- function(overlap, mrna_expr_data, mirna_expr_data, geno_data,
                                 num_boots = 1000,reverse = FALSE,boottype = "perc",
                           confidence = 0.995){
  invisible(library(boot))
  #====================================#
  # indirect effect function
  #====================================#
  
  indirectsaved <-  function (dataset, random) {
    d = dataset[random, ] ##randomize by row
    apath = lm(m ~ x, data = d)
    bpath = lm(y ~ x + m, data = d)
    indirect <- apath$coefficients[2]*bpath$coefficients[3]
    return(indirect)
  }
  
  #====================================#
  # prop mediated function
  #====================================#
  
  propmedsaved <-  function (dataset, random) {
    d = dataset[random, ] ##randomize by row
    apath = lm(m ~ x, data = d)
    bpath = lm(y ~ x + m, data = d)
    propmed <- abs(apath$coefficients[2]*bpath$coefficients[3]) / (abs(apath$coefficients[2]*bpath$coefficients[3]) + abs(bpath$coefficients[2]))
    return(propmed)
  }
  
  
  indirect_vec <- list()
  propmed_vec <- list()
  
  # iterate over all mediation pairs #
  for (ii in 1:nrow(overlap)){
    
    # extract specificied vectors #
    mrna_vec <- mrna_expr_data[rownames(mrna_expr_data) == overlap[ii,"mrna"],] %>% matrix() %>% unlist
    mirna_vec <- mirna_expr_data[rownames(mirna_expr_data) == overlap[ii,"miRNA"],] %>% matrix() %>% unlist
    X <- geno_data[rownames(geno_data) == overlap[ii,"sdp"],] %>% c() %>% unname
    
    if (reverse == FALSE){
      M <- mirna_vec
      Y <- mrna_vec
    }
    else{
      M <- mrna_vec
      Y <- mirna_vec
    }

    # extract beta coefficients, 95% CI, proportion mediated, and p-value #
    
    med_df <- as.data.frame(cbind(M,X,Y))
    colnames(med_df) <- c("m","x","y")
    

    set.seed(3112019)
    indirectresults <- boot(data = med_df,
                       statistic = indirectsaved,
                       R = num_boots)

    indirect_ci <- boot.ci(indirectresults, 
                        conf = confidence,
                        type = "perc")[[4]]

    indirect_vec[[ii]] <- c(unname(indirectresults[[1]]),indirect_ci[[4]],indirect_ci[[5]])
    
    propmedresults <- boot(data = med_df,
                                  statistic = propmedsaved,
                                  R = num_boots)

    propmed_ci <- boot.ci(propmedresults,
                                 conf = confidence,
                                 type = "perc")[[4]]

    propmed_vec[[ii]] <- c(unname(propmedresults[[1]]), propmed_ci[[4]], propmed_ci[[5]])
    
    print(ii)
  }
  
  # return dataframes for indirect effect and proportion mediated #
  indirect_df <- as.data.frame(Reduce(rbind,indirect_vec))
  colnames(indirect_df) <- c("indirect","ci_low","ci_high")
  indirect_df$sig <- ifelse((indirect_df$ci_low < 0 & indirect_df$ci_high < 0) | (indirect_df$ci_low > 0 & indirect_df$ci_high > 0),
                            1, 0)
  
  propmed_df <- as.data.frame(Reduce(rbind, propmed_vec))
  colnames(propmed_df) <- c("propmed","ci_low","ci_high")
  propmed_df$sig <- ifelse((propmed_df$ci_low != 0),1, 0)
  return(list(indirect_df,propmed_df))
  
}



#==================================#
# FUNCTION: BNA hillclimb function
#==================================#

bna_exhaustive <- function(gene_ls,gene_names,mirna,sdp,mirna_expr_data,
                           geno_data,mrna_expr_data, plot=T,
                           num_ind_tests,sdp_names,blacklist,num_boots){
  
  gene_expr <- list()
  
  for (ii in 1:length(gene_ls)){
    gene_expr[[ii]] <- mrna_expr_data[rownames(mrna_expr_data) == gene_ls[[ii]],] %>% matrix() %>% unlist
  }
  
  gene_mat <- do.call(cbind, gene_expr)
  
  mir <- mirna_expr_data[rownames(mirna_expr_data) == mirna,] %>% matrix() %>% unlist
  snp <- geno_data[rownames(geno_data) == sdp,] %>% c() %>% unname
  
  # set up for independent network #
  snp_mat <- matrix(snp, length(snp), num_ind_tests)
  
  # extract beta coefficients, 95% CI, proportion mediated, and p-value #
  
  smalldata <- as.data.frame(cbind(snp,gene_mat,mir))
  colnames(smalldata) <- c("snp",gene_names,"mir")
  
  # obtain BIC for both directions #
  #best_bna <- hc(smalldata,blacklist = blacklist)
  boot_edges <- boot.strength(smalldata,R = num_boots,algorithm = "hc",algorithm.args = list(blacklist = blacklist))
  avg_boot <- averaged.network(boot_edges,threshold = 0.5)

  return(list(avg_boot,boot_edges))
}

