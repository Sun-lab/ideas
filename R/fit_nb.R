fit_nb <- 
function(x, z, max.theta=1e6) {

  MIN_VAL = 1e-8
  
  if(! is.numeric(x)){
    stop("x should be a numerical variable.")
  }
  
  if(any(x < 0)){
    stop("x cannot be negative.")
  }
  
  if(any(is.na(x))){
    stop("x should not contaion NAs.")
  }
  
  if(is.null(z)){
    nz = 0
  }else if(is.vector(z)){
    nz = 1
  }else if(is.data.frame(z) || is.matrix(z)){
    nz = ncol(z)
  }else{
    stop('z should be data.frame or matrix\n')
  }
  
  z = as.data.frame(z)
  str_z = paste(names(z), collapse= "+")
  fmla.nb   = as.formula(paste("x ~ ", str_z))

  fill4z = rep(0, nz)
  
  if ((max(x) - min(x) == 0)) {
    if(all(x==0)){ # var = mu + mu^2/theta
      fit_total = list(c(log(MIN_VAL),fill4z), max.theta)
    }else{
      fit_total = list(c(log(unique(x)),fill4z), max.theta)
    }
  }else{
      fm_nb = NULL
      if (is.null(z)) {
        fm_nb = tryCatch(
          MASS::glm.nb(x ~  1),
          error = function(e) { NULL }
        )
      }else {
        fm_nb = tryCatch(
          MASS::glm.nb(fmla.nb, data=z),
          error = function(e) { NULL }
        )
      }
      if (! is.null(fm_nb)) {
        fit_logmean      = fm_nb$coefficients
        fit_dispersion   = fm_nb$theta
        fit_total = list(fit_logmean, fit_dispersion)
      }else{
        fit_total = NULL
      }
  
  }
  if(is.null(fit_total)){
    fit_total = list(NA, NA)
  }
  names(fit_total) = c("logmean", "dispersion")
  
  return(fit_total)
}
