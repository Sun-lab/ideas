

cal_G <- 
  function(m){# m is a distance matrix
    n = nrow(m)
    centerM = diag(n) - 1/n
    G  = -0.5 * centerM %*% (m*m) %*% centerM
    eG = eigen(G, symmetric = TRUE)
    G  = eG$vector %*% diag(pmax(0, eG$values)) %*% t(eG$vector)
    G
  }

calc_F_manova <- 
  function(dist_array, label){
    
    uLabel = unique(label)
    a      = length(uLabel)
    N      = length(label)
    
    n.genes = dim(dist_array)[1]
    d2      = dist_array*dist_array
    
    # change the dimension of 3-way array to 2-way array
    dim(d2) = c(n.genes, N^2)
    sst     = rowSums(d2)/N
    
    # this is an N x N matrix with 0 for between group entries and 
    # 1/group_size for within group entries
    divid_vector = matrix(0,N,N)
    
    for(ik in 1:a){
      cur_index = which(label==uLabel[ik])
      divid_vector[cur_index,cur_index] = 1/length(cur_index)
    }
    divid_vector2 = c(divid_vector)
    
    ssw = d2 %*% divid_vector2
    
    Fstat=((sst-ssw)*(N-a))/(ssw*(a-1))
    return(Fstat)
  }


calTrace <- 
  function(G,H){
    sum(diag(H%*%G%*%H))
  }

calc_F_permanovaSZ <- 
  function(dist_array, Rs, z = NULL){
    
    n = dim(dist_array)[2]
    G = apply(dist_array, 1, cal_G)
    dim(G) = c(n,n,dim(dist_array)[1])
    
    if(is.vector(Rs)){
      Rs = matrix(Rs, ncol=1)
    } 
    
    i = 0
    F_stats = foreach (i = 1:ncol(Rs), .combine='cbind') %dorng% {
      
      if(is.null(z)){
        xz = Rs[,i]
      }else{
        xz  = unique(cbind(Rs[,i], z), MARGIN=2)
      }
      
      H   = xz %*% solve(crossprod(xz)) %*% t(xz)
      IH  = diag(nrow(H)) - H
      
      t1 = apply(G, 3, function(g){calTrace(g,H)})
      t2 = apply(G, 3, function(g){calTrace(g,IH)})
      
      F1 = t1/t2
    }
    
    return(F_stats)
  }
