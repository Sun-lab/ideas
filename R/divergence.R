
calDens_nbzinb <- 
  function(x_triple, y_triple, range_max = 0, quantile = 0.975) {
    
    if (range_max == 0) {
      range_max = max(
        qnbinom(p = max(0, (quantile - x_triple[3])/(1 - x_triple[3])),
                mu = x_triple[1], size = x_triple[2]),
        qnbinom(p = max(0, (quantile - y_triple[3])/(1 - y_triple[3])), 
                mu = y_triple[1], size = y_triple[2]),
        na.rm = TRUE
      )
    }
    
    x = emdbook::dzinbinom(
      0:range_max,
      mu = x_triple[1],
      size = x_triple[2],
      zprob = x_triple[3],
      log = FALSE
    )
    y = emdbook::dzinbinom(
      0:range_max,
      mu = y_triple[1],
      size = y_triple[2],
      zprob = y_triple[3],
      log = FALSE
    )
    res = list(dens_x = x, dens_y = y)
    return(res)
    
  }




calDens_nb <- 
function(x_triple, y_triple, range_max = 0, quantile = 0.975) {
  
  if (range_max == 0) {
    range_max = max(
      qnbinom(p = quantile,
              mu = x_triple[1], size = x_triple[2]),
      qnbinom(p = quantile, 
              mu = y_triple[1], size = y_triple[2]),
      na.rm = TRUE
    )
  }
  
  x = dnbinom(
    0:range_max,
    mu = x_triple[1],
    size = x_triple[2],
    log = FALSE
  )
  y = dnbinom(
    0:range_max,
    mu = y_triple[1],
    size = y_triple[2],
    log = FALSE
  )
  res = list(dens_x = x, dens_y = y)
  return(res)
  
}


calDensity <- 
function(x, y, n = 1024) {
  
  cur_range = range(c(x, y), na.rm = TRUE, finite = TRUE)
  x_dens = density(x, from = cur_range[1], to = cur_range[2], 
                   n = n, na.rm = TRUE)
  
  y_dens = density(y, from = cur_range[1], to = cur_range[2], 
                   n = n, na.rm = TRUE)
  
  res = list(dens_x = x_dens$y, dens_y = y_dens$y)
  return(res)
  
}


calc_JSD <- 
function(px, qx, lower.bound=1e-10) {
  px[which(px < lower.bound)] = lower.bound
  qx[which(qx < lower.bound)] = lower.bound
  mx = (px + qx) / 2
  - (sum(px * (log(mx / px)), na.rm=TRUE) +
       sum(qx * (log(mx / qx)), na.rm=TRUE)) / 2
}

divergence <- 
function(x, y, d_metric = c("JSD", "Was"), empirical_n = 1024, 
         zinb.range_max = 0,  zinb.quantile = 0.975, 
         fit_method = c("nb", "dca_direct", "saver_direct", "zinb", "kde")) {
  
  if (fit_method == "zinb") {
    res = calDens_nbzinb(x, y, quantile = zinb.quantile, 
                         range_max = zinb.range_max)
    px = res$dens_x
    qx = res$dens_y
    px = px / sum(px)
    qx = qx / sum(qx)
    
    if(d_metric == "JSD"){
      dist.val = calc_JSD(px, qx)
    }else if(d_metric == "Was"){
      Fx = cumsum(px)
      Fy = cumsum(qx)
      dist.val = sum(abs(Fx - Fy))
    }else{
      stop("unknown distance metric\n")
    }
  }else if(fit_method == "nb"){
    res = calDens_nb(x, y, quantile = zinb.quantile, 
                         range_max = zinb.range_max)
    px = res$dens_x
    qx = res$dens_y
    px = px / sum(px)
    qx = qx / sum(qx)
    
    if(d_metric == "JSD"){
      dist.val = calc_JSD(px, qx)
    }else if(d_metric == "Was"){
      Fx = cumsum(px)
      Fy = cumsum(qx)
      dist.val = sum(abs(Fx - Fy))
    }else{
      stop("unknown distance metric\n")
    }    
  }else if(fit_method == "dca_direct" || fit_method == "saver_direct"){
    range_max = max(length(x), length(y)) - 1
    if (length(x) < (range_max + 1) ){
      x = c(x, rep(0, range_max + 1 - length(x)))
    }
    if (length(y) < (range_max + 1) ){
      y = c(y, rep(0, range_max + 1 - length(y)))
    }    
    px = x / sum(x)
    qx = y / sum(y)
    
    if(d_metric == "JSD"){
      dist.val = calc_JSD(px, qx)
    }else if(d_metric == "Was"){
      Fx = cumsum(px)
      Fy = cumsum(qx)
      dist.val = sum(abs(Fx - Fy))
    }else{
      stop("unknown distance metric\n")
    }    
  }else if (fit_method == "kde" && d_metric == "JSD") {
    res = calDensity(x, y, n = empirical_n)
    px  = res$dens_x + 1 / empirical_n
    qx  = res$dens_y + 1 / empirical_n
    
    px = px / sum(px)
    qx = qx / sum(qx)
    dist.val = calc_JSD(px, qx)
  }else if (fit_method == "kde" && d_metric == "Was"){
    dist.val = wasserstein1d(x, y)
  }else{
    stop("unknown combination of fit_method and metric\n")
  }
  
  return(dist.val)
}
