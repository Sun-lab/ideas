
permanova <-
  function(dist_array, meta_ind, var2test, var2adjust=NULL, 
           var2test_type = c("binary", "continuous"), 
           n_perm = 999, r.seed=2020, 
           residulize.x = FALSE, delta = 0.5){
    
    # -----------------------------------------------------------------
    # check dist_array
    # -----------------------------------------------------------------
    n_genes = dim(dist_array)[1]
    
    wNA = which(apply(dist_array, 1, anyNA))
    if(length(wNA) > 0){
      message(sprintf("skip %d gene(s) with NA in the dist_array\n", length(wNA)))
      dist_array = dist_array[-wNA,,]
    }
    
    # -----------------------------------------------------------------
    # check the input data of meta_ind
    # -----------------------------------------------------------------
    
    if(!(is.data.frame(meta_ind))){
      stop("meta_ind should be a data.frame\n")
    }
    
    columns.meta.ind = c("individual", var2test, var2adjust)
    
    if(! all(columns.meta.ind %in% names(meta_ind))){
      str1 = paste(columns.meta.ind, collapse=", ")
      stop("names of meta_ind should conttain: ",str1)
    }
    
    if(length(unique(meta_ind$individual)) != nrow(meta_ind)){
      stop("the individual ids in meta_ind are not unique\n")
    }
    
    # -----------------------------------------------------------------
    # check the var2test
    # -----------------------------------------------------------------
    
    x = meta_ind[[var2test]]
    
    if(any(is.na(x))){
      stop("variable to test has NA values.\n")
    }
    
    nUniq_x = length(unique(x))
    
    if(nUniq_x == 2){
      if(var2test_type == "continuous"){
        str1 = "the variable to test has only two unique values, "
        warning(paste0(str1, "will treat it as a binary variable\n"))
        var2test_type = "binary"
      }
      if(min(table(x)) < 5){
        stop("permutation test will be inaccurate with n < 5 in each group\n")
      }
    }else{
      if(var2test_type == "binary"){
        str1 = "the variable to test has only more than two unique values, "
        warning(paste0(str1, "will treat it as a continuous variable\n"))
        var2test_type = "continuous"
      }
    }
    
    message(sprintf("testing for '%s', a %s variable\n", var2test, var2test_type))
    if(var2test_type == "binary"){
      # make sure x is a vector of 0 or 1. 
      if(is.character(x)){ x = as.factor(x) }
      if(is.factor(x)){ x = as.numeric(x) }
      if(is.numeric(x)){ x = as.numeric(x == max(x)) }
    }
    
    # -----------------------------------------------------------------
    # start testing
    # -----------------------------------------------------------------
    
    set.seed(r.seed)
    x_perm = matrix(rep(x, times = n_perm), ncol = n_perm)
    x_perm = apply(x_perm, 2, sample, size = length(x))
    
    # if there is no cavariate, use standard permanova
    if(is.null(var2adjust)){
      if(var2test_type == "binary"){
        F_ob   = calc_F_manova(dist_array, label = x)
        F_perm = apply(x_perm, 2, function(perm1){ 
          calc_F_manova(dist_array, label = perm1)
        })
      }else{
        F_ob    = calc_F_permanovaSZ(dist_array, Rs = x)
        F_perm  = calc_F_permanovaSZ(dist_array, Rs = x_perm)
      }
    }else{
      fm1 = as.formula(paste("~", paste(var2adjust, collapse=" + ")))
      z   = model.matrix(fm1, data = meta_ind)
      
      if(residulize.x){
        resid_perm = matrix(NA, nrow(meta_ind), n_perm)
        
        if(var2test_type == "binary"){
          # step1: fit
          m1 = glm(x ~ -1 + z, family=binomial(link="logit"))  # logistic model
          fitted_x = fitted(m1)
          resid_x  = x - fitted_x
          
          # step2: permutation
          ip = 0 # index for usable permutations
          id = 0 # working index
          sd_e = sd(resid_x) # expected sd
          
          while(ip < n_perm){
            id = id + 1
            perm_x = rbinom(length(fitted_x), 1, prob=fitted_x)
            
            m2 = tryCatch(glm(perm_x ~ -1 + z, family=binomial(link="logit")), 
                          warning = function(w) { NULL }, 
                          error = function(e) { NULL},
                          finally = {})
            
            if(! is.null(m2)){
              resid_p_i  = perm_x - fitted(m2) 
              sd_resid_i = sd(resid_p_i)
              
              if(sd_resid_i > (1-delta)*sd_e & sd_resid_i < (1+delta)*sd_e){
                ip = ip + 1
                resid_perm[,ip] = resid_p_i
              }
            }
          }
          id
          ip
        }else{
          
          #step1: fit
          m1 = lm(x ~ -1 + z)  #linear model
          fitted_x = fitted(m1)
          resid_x  = x - fitted_x
          
          #step2: permutation
          for(ip in 1:n_perm){
            perm_x = sample(resid_x, size=length(resid_x)) + fitted_x 
            m2     = lm(perm_x ~ -1 + z)
            resid_perm[,ip] = perm_x - fitted(m2) 
          }
        }
        
        resid_x = matrix(resid_x, ncol=1)
        F_ob    = calc_F_permanovaSZ(dist_array, Rs = resid_x, z = z)
        F_perm  = calc_F_permanovaSZ(dist_array, Rs = resid_perm, z = z)
        
      }else{
        
        F_ob    = calc_F_permanovaSZ(dist_array, Rs = x, z = z)
        F_perm  = calc_F_permanovaSZ(dist_array, Rs = x_perm, z = z)
      }
    }
    
    F_perm0 = apply(F_perm, 2, function(x){return(x-F_ob)} )
    
    if(length(F_ob)==1){ F_perm0 = t(F_perm0) }
    
    pval = rowMeans(cbind(F_perm0, rep(1,nrow(F_perm0))) >= 0)
    
    if(length(wNA) == 0){
      pval_all = pval
    }else{
      pval_all = rep(NA, n_genes)
      pval_all[-wNA] = pval
    }
    
    return(pval_all)
  }
