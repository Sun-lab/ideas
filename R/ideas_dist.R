
# per_cell_adjust = c("NB", "both"); quantile = 0.975; empirical_n = 1024

ideas_dist <-
  function(count_input, meta_cell, meta_ind, var_per_cell, var2test, 
           var2test_type = c("binary", "continuous"), 
           d_metric = c("Was", "JSD"), 
           fit_method = c("nb", "zinb", "kde", "dca_direct", "saver_direct"), 
           per_cell_adjust = c("NB", "both"), quantile = 0.975,
           empirical_n = 1024) {
    
    var2test_type   = var2test_type[1]
    d_metric        = d_metric[1]
    fit_method      = fit_method[1]
    
    
    if(!(is.data.frame(meta_cell))){
      stop("meta_cell should be a data.frame\n")
    }
    
    if(is.data.table(meta_cell)){
      meta_cell = as.data.frame(meta_cell)
    }
    
    if(!(is.data.frame(meta_ind))){
      stop("meta_ind should be a data.frame\n")
    }

    if(is.data.table(meta_ind)){
      meta_ind = as.data.frame(meta_ind)
    }

    if (fit_method == "zinb"){
      per_cell_adjust = per_cell_adjust[1]
    }
    
    if (fit_method == "dca_direct"){
      # -----------------------------------------------------------------
      # check the input data of count_input, when fit_method == dca_direct
      # -----------------------------------------------------------------
      if (length(count_input) != 3){
        stop("when fit_method is dca_direct, the count_input should be a 
             list of 3 matrices")
      }
      # -----------------------------------------------------------------
      # check cell_id order of meta_cell, when fit_method == dca_direct
      # -----------------------------------------------------------------
      if(any(meta_cell$cell_id != colnames(count_input[[1]]))){
        stop("cell_id in meta_cell do not match colnames of mean_norm 
             matrix\n")
      }
      
    }else{
      # -----------------------------------------------------------------
      # check the input data of count_input,when fit_method != dca_direct
      # -----------------------------------------------------------------
      count_matrix = count_input
      
      if(! is.matrix(count_matrix)){
        stop("count_matrix is not a matrix\n")
      }
      
      if( fit_method %in% c("nb", "zinb", "kde")){
        check_count <- function(v){any(v != round(v) | v < 0)}
        not_count = apply(count_matrix, 1, check_count)
        
        if(any(not_count)){
          str1 = "count_matrix should only include non-negative integers"
          str1 = sprintf("%s, violation in row %d\n", str1, which(not_count)[1])
          stop(str1)
        }
      }
      
      n_cell = ncol(count_matrix)
      n_gene = nrow(count_matrix)
      
      gene_ids = rownames(count_matrix)
      cell_ids = colnames(count_matrix)
      
      if(is.null(gene_ids)){
        stop("count_matrix should have row names for gene ids\n")
      }
      
      if(is.null(cell_ids)){
        stop("count_matrix should have col names for cell ids\n")
      }
      
      if(length(unique(gene_ids)) != n_gene){
        stop("row names of count_matrix (gene ids) are not unique\n")
      }
      
      if(length(unique(cell_ids)) != n_cell){
        stop("col names of count_matrix (cell ids) are not unique\n")
      }
      
      message(sprintf("the count_matrix includes %d genes in %d cells\n", 
                      n_gene, n_cell))
      
      # -----------------------------------------------------------------
      # check cell_id order of meta_cell, when fit_method != dca_direct
      # -----------------------------------------------------------------
      
      if(any(meta_cell$cell_id != colnames(count_matrix))){
        stop("cell_id in meta_cell do not match colnames of count_matrix\n")
      }
      
    }

    # -----------------------------------------------------------------
    # check other aspects of meta_cell
    # -----------------------------------------------------------------    
    
    columns.meta.cell = c("cell_id", "individual", var_per_cell)
    
    if(! all(columns.meta.cell %in% names(meta_cell))){
      str1 = paste(columns.meta.cell, collapse=", ")
      stop(sprintf("names of meta_cell should contain %s\n", str1))
    }
    
    if(length(unique(meta_cell$cell_id)) != nrow(meta_cell)){
      stop("the cell_id's in meta_cell are not unique\n")
    }
    
    # -----------------------------------------------------------------
    # check the input data of meta_ind
    # -----------------------------------------------------------------
    
    columns.meta.ind = c("individual", var2test)
    
    if(! all(columns.meta.ind %in% names(meta_ind))){
      str1 = paste(columns.meta.ind, collapse=", ")
      stop(sprintf("names of meta_ind should contain %s\n", str1))
    }
    
    if(! setequal(meta_cell$individual, meta_ind$individual)){
      stop("the individual ids in meta_cell and meta_ind do not match\n")
    }
    
    if(length(unique(meta_ind$individual)) != nrow(meta_ind)){
      stop("the individual ids in meta_ind are not unique\n")
    }
    
    
    #if(any(meta_cell[,..var_per_cell] <= 0.0)){
    if(any(meta_cell[,var_per_cell, drop = FALSE] <= 0.0)){
      str1 = "the variables listed in 'var_per_cell' will be log transformed,"
      stop(paste(str1, "so they must be positive."))
    }
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using zinb
    # -----------------------------------------------------------------
    
    extract_zinb_par <- function(zb_fit, cov_value){
      mean_a = exp(t(zb_fit$logmean) %*% c(1, cov_value))
      disp_a = zb_fit$dispersion
      drop_a = exp(t(zb_fit$logitdropout) %*% c(1, cov_value))
      drop_a = drop_a/(1 + drop_a)
      c(mean_a, disp_a, drop_a)
    }
    
    i_g = 0
    
    if(fit_method == "zinb"){
      message("estimating distribution for each gene and each individaul by zinb\n")
      cov_value = apply(log10(meta_cell[,var_per_cell,drop=FALSE]), 2, median)
      #cov_value = apply(log10(meta_cell[, ..var_per_cell]), 2, median)
      zinb_fit=foreach (i_g = 1:n_gene) %dorng% {
        zinb_fit1 = list()
        # collapose counts per individual and estimate ZINB
        for (j in 1:nrow(meta_ind)) {
          cur_ind = meta_ind$individual[j]
          w2use   = which(meta_cell$individual == cur_ind)
          dat_ind = c(count_matrix[i_g, w2use])
          #z       = log10(meta_cell[w2use, ..var_per_cell])
          z       = log10(meta_cell[w2use, var_per_cell, drop=FALSE])
          zinb_fit1[[j]] = fit_zinb(x=dat_ind, z=z, per_cell_adjust)
        }
        names(zinb_fit1) = as.character(meta_ind$individual)
        zinb_fit1
      }
      
      length(zinb_fit)
      length(zinb_fit[[1]])
      zinb_fit[[1]][[1]]
      
      dist_array_list=foreach (i = 1:n_gene) %dorng% {
        gene_i_fit  = zinb_fit[[i]]
        
        dist_array1 = array(NA, dim=rep(nrow(meta_ind), 2))
        rownames(dist_array1) = meta_ind$individual
        colnames(dist_array1) = meta_ind$individual
        diag(dist_array1) = 0
        
        for (j_a in 1:(nrow(meta_ind)-1)) {
          if(is.na(gene_i_fit[[j_a]]$logmean[1])){next }
          par_a = extract_zinb_par(gene_i_fit[[j_a]], cov_value)
          # occasionally, the mean value fit can be very large
          if(par_a[1] > 1000){next }
          
          for (j_b in (j_a+1):nrow(meta_ind)) {
            if(is.na(gene_i_fit[[j_b]]$logmean[1])){next }
            par_b = extract_zinb_par(gene_i_fit[[j_b]], cov_value)
            # occasionally, the mean value fit can be very large
            if(par_b[1] > 1000){next }
            
            dist_array1[j_a, j_b] = tryCatch(
              divergence(par_a, par_b, d_metric = d_metric, 
                         fit_method = fit_method), 
              error = function(e) { NA }
            )
            
            dist_array1[j_b, j_a] = dist_array1[j_a, j_b]
          }
        }
        
        dist_array1
      }
    }
    
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using nb
    # -----------------------------------------------------------------  
    
    extract_nb_par <- function(nb_fit, cov_value){
      mean_a = exp(t(nb_fit$logmean) %*% c(1, cov_value))
      disp_a = nb_fit$dispersion
      c(mean_a, disp_a)
    }
    
    
    if(fit_method == "nb"){
      message("estimating distribution for each gene and each individual by nb\n")
      cov_value = apply(log10(meta_cell[,var_per_cell,drop=FALSE]), 2, median)
      #cov_value = apply(log10(meta_cell[, ..var_per_cell]), 2, median)
      nb_fit=foreach (i_g = 1:n_gene) %dorng% {
        nb_fit1 = list()
        # collapose counts per individual and estimate NB
        for (j in 1:nrow(meta_ind)) {
          cur_ind = meta_ind$individual[j]
          w2use   = which(meta_cell$individual == cur_ind)
          dat_ind = c(count_matrix[i_g, w2use])
          #z       = log10(meta_cell[w2use, ..var_per_cell])
          z       = log10(meta_cell[w2use, var_per_cell, drop=FALSE])
          nb_fit1[[j]] = fit_nb(x=dat_ind, z=z)
        }
        names(nb_fit1) = as.character(meta_ind$individual)
        nb_fit1
      }
      
      length(nb_fit)
      length(nb_fit[[1]])
      nb_fit[[1]][[1]]
      
      dist_array_list=foreach (i = 1:n_gene) %dorng% {
        gene_i_fit  = nb_fit[[i]]
        
        dist_array1 = array(NA, dim=rep(nrow(meta_ind), 2))
        rownames(dist_array1) = meta_ind$individual
        colnames(dist_array1) = meta_ind$individual
        diag(dist_array1) = 0
        
        for (j_a in 1:(nrow(meta_ind)-1)) {
          if(is.na(gene_i_fit[[j_a]]$logmean[1])){next }
          par_a = extract_nb_par(gene_i_fit[[j_a]], cov_value)
          # occasionally, the mean value fit can be very large
          if(par_a[1] > 1000){next }
          
          for (j_b in (j_a+1):nrow(meta_ind)) {
            if(is.na(gene_i_fit[[j_b]]$logmean[1])){next }
            par_b = extract_nb_par(gene_i_fit[[j_b]], cov_value)
            # occasionally, the mean value fit can be very large
            if(par_b[1] > 1000){next }
            
            dist_array1[j_a, j_b] = tryCatch(
              divergence(par_a, par_b, d_metric = d_metric, 
                         fit_method = fit_method), 
              error = function(e) { NA }
            )
            
            dist_array1[j_b, j_a] = dist_array1[j_a, j_b]
          }
        }
        
        dist_array1
      }
    }
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using kde
    # -----------------------------------------------------------------
    if (fit_method == "kde") {
      message("estimating distribution for each gene and each individual by kde\n")
      cov_value = apply(log10(meta_cell[,var_per_cell,drop=FALSE]), 2, median)
      #cov_value = apply(log10(meta_cell[, ..var_per_cell]), 2, median)
      dat_res=foreach (i_g = 1:n_gene) %dorng% {
        res_ig = list()
        for (j in 1:nrow(meta_ind)) {
          ind_j = meta_ind$individual[j]
          w2use = which(meta_cell$individual == ind_j)
          dat_j = c(count_matrix[i_g, w2use])
          #z     = log10(meta_cell[w2use, ..var_per_cell])
          z     = log10(meta_cell[w2use, var_per_cell, drop = FALSE])
          z     = as.data.frame(z)
          str_z = paste(names(z), collapse= "+")
          fmla  = as.formula(paste("log10(dat_j + 0.5) ~ ", str_z))
          lm_j  = lm(fmla, data = z)
          base_j = c(t(lm_j$coefficients) %*% c(1, cov_value))
          res_ig[[j]] = lm_j$resid + base_j
        }
        names(res_ig) = as.character(meta_ind$individual)
        res_ig
      }
      
      dist_array_list=foreach (i_g = 1:n_gene) %dorng% {
        
        res_ig = dat_res[[i_g]]
        dist_array1 = array(NA, dim=rep(nrow(meta_ind), 2))
        rownames(dist_array1) = meta_ind$individual
        colnames(dist_array1) = meta_ind$individual
        diag(dist_array1) = 0
        
        for (j_a in 1:(nrow(meta_ind)-1)) {
          res_a = res_ig[[j_a]]
          
          for (j_b in (j_a+1):nrow(meta_ind)) {
            res_b = res_ig[[j_b]]
            
            dist_array1[j_a, j_b] = tryCatch(
              divergence(res_a, res_b, d_metric = d_metric,  
                         fit_method = fit_method, empirical_n = empirical_n), 
              error = function(e) { NA }
            )
            
            dist_array1[j_b, j_a] = dist_array1[j_a, j_b]
          }
        }
        dist_array1
      }
    }
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using combined dca outputs
    # -----------------------------------------------------------------
    if(fit_method == "dca_direct"){
      message("estimating distribution for each gene and each individual by dca_direct\n")
      dca_mean = count_input[[1]]
      dca_disp = count_input[[2]]
      dca_pi   = count_input[[3]]
      
      gene_ids = rownames(dca_mean)
      n_gene   = nrow(dca_mean)
      
      # %dopar% to %dorng%, may not be necessary
      dca_mass = foreach (i_g = 1:n_gene) %dorng% { 
        dca_mass1 = list()
        
        for (j in 1:nrow(meta_ind)) {
          cur_ind = meta_ind$individual[j]
          w2use   = which(meta_cell$individual == cur_ind)
          mean_ind = c(dca_mean[i_g, w2use])
          disp_ind = c(dca_disp[i_g, w2use])
          pi_ind = c(dca_pi[i_g, w2use])
          # number of cells for the current individual
          ncell_ind = length(w2use) 
          range_list = rep(0,  ncell_ind)
          # get the max for the range of distributions within the current ind
          for (k in 1: ncell_ind){
            range_list[k] = 
              qnbinom(p = max(0, (quantile - pi_ind[k])/(1 - pi_ind[k])),
                      mu = mean_ind[k], size = disp_ind[k])
          }
          range_max = max(range_list)
          mass_matrix = matrix(ncol = (range_max + 1), nrow = ncell_ind)
          for (k in 1:ncell_ind){
            mass_matrix[k, ] = emdbook::dzinbinom(
              0:range_max,
              mu = mean_ind[k],
              size = disp_ind[k],
              zprob = pi_ind[k],
              log = FALSE
            )
          }
          
          dca_mass1[[j]] = apply(mass_matrix, 2, mean)
        }
        names(dca_mass1) = as.character(meta_ind$individual)
        dca_mass1
      }
      
      length(dca_mass)
      length(dca_mass[[1]])
      dca_mass[[1]][[1]]
      
      dist_array_list=foreach (i = 1:n_gene) %dorng% {
        gene_i_mass  = dca_mass[[i]]
        
        dist_array1 = array(NA, dim=rep(nrow(meta_ind), 2))
        rownames(dist_array1) = meta_ind$individual
        colnames(dist_array1) = meta_ind$individual
        diag(dist_array1) = 0
        
        for (j_a in 1:(nrow(meta_ind)-1)) {
          # unscaled
          mass_a = gene_i_mass[[j_a]]
          for (j_b in (j_a+1):nrow(meta_ind)) {
            mass_b = gene_i_mass[[j_b]]
            
            dist_array1[j_a, j_b] = tryCatch(
              divergence(mass_a, mass_b, d_metric = d_metric, 
                         fit_method = fit_method), 
              error = function(e) { NA }
            )
            
            dist_array1[j_b, j_a] = dist_array1[j_a, j_b]
          }
        }
        
        dist_array1
      }
    }
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using combined saver outputs
    # -----------------------------------------------------------------
    if(fit_method == "saver_direct"){
      message("estimating distribution for each gene and each individual by saver_direct\n")
      saver_mean = count_input

      gene_ids = rownames(saver_mean)
      n_gene   = nrow(saver_mean)
      
      # %dopar% to %dorng%, may not be necessary
      saver_mass = foreach (i_g = 1:n_gene) %dorng% { 
        saver_mass1 = list()
        
        for (j in 1:nrow(meta_ind)) {
          cur_ind = meta_ind$individual[j]
          w2use   = which(meta_cell$individual == cur_ind)
          mean_ind = c(saver_mean[i_g, w2use])
          # number of cells for the current individual
          ncell_ind = length(w2use) 
          range_max = qpois(p = max(0, quantile), lambda = max(mean_ind))
          
          mass_matrix = matrix(ncol = (range_max + 1), nrow = ncell_ind)
          for (k in 1:ncell_ind){
            mass_matrix[k, ] = dpois(0:range_max, lambda=mean_ind[k])
          }
          
          saver_mass1[[j]] = apply(mass_matrix, 2, mean)
        }
        names(saver_mass1) = as.character(meta_ind$individual)
        saver_mass1
      }
      
      length(saver_mass)
      length(saver_mass[[1]])
      saver_mass[[1]][[1]]
      
      dist_array_list=foreach (i = 1:n_gene) %dorng% {
        gene_i_mass  = saver_mass[[i]]
        
        dist_array1 = array(NA, dim=rep(nrow(meta_ind), 2))
        rownames(dist_array1) = meta_ind$individual
        colnames(dist_array1) = meta_ind$individual
        diag(dist_array1) = 0
        
        for (j_a in 1:(nrow(meta_ind)-1)) {
          # unscaled
          mass_a = gene_i_mass[[j_a]]
          for (j_b in (j_a+1):nrow(meta_ind)) {
            mass_b = gene_i_mass[[j_b]]
            
            dist_array1[j_a, j_b] = tryCatch(
              divergence(mass_a, mass_b, d_metric = d_metric, 
                         fit_method = fit_method), 
              error = function(e) { NA }
            )
            
            dist_array1[j_b, j_a] = dist_array1[j_a, j_b]
          }
        }
        
        dist_array1
      }
    }
    
    
    length(dist_array_list)
    dim(dist_array_list[[1]])
    dist_array_list[[1]][1:2,1:2]
    
    nNA = sapply(dist_array_list, function(x){sum(is.na(c(x)))})
    table(nNA)
    
    dist_array = array(
      dim = c(
        n_gene,
        nrow(meta_ind),
        nrow(meta_ind)
      ),
      dimnames = list(gene_ids, meta_ind$individual, meta_ind$individual)
    )
    
    dim(dist_array)
    
    for (i in 1:n_gene){
      dist_array[i,,] = dist_array_list[[i]]
    }
    
    dim(dist_array)
    dist_array[1,1:2,1:2]
    
    dist_array
  }
