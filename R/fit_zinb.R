fit_zinb <- 
function(x, z, per_cell_adjust, max.theta=1e6) {

  MIN_VAL = 1e-8
  
  logit <- function(x){log(x/(1-x)) }
  
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
  fmla.NB   = as.formula(paste("x ~ ", str_z, "| 1"))
  fmla.both = as.formula(paste("x ~ ", str_z))

  fill4z = rep(0, nz)
  
  if ((max(x) - min(x) == 0)) {
    if(all(x==0)){ # var = mu + mu^2/theta
      fit_total = list(c(log(MIN_VAL),fill4z), max.theta, 
                       c(logit(1 - MIN_VAL),fill4z))
    }else{
      fit_total = list(c(log(unique(x)),fill4z), max.theta, 
                       c(logit(MIN_VAL),fill4z))
    }
  }else{
    if (min(x) == 0) {
      fm_zinb = NULL
      if (is.null(z)) {
        fm_zinb = tryCatch(
          pscl::zeroinfl(x ~  1, dist = "negbin"),
          error = function(e) { NULL }
        )
      }else{
        if(per_cell_adjust == "NB"){
          fm_zinb = tryCatch(
            pscl::zeroinfl(fmla.NB, data = z, dist = "negbin"),
            error = function(e) { NULL }
          )
        }else{
          fm_zinb = tryCatch(
            pscl::zeroinfl(fmla.both, data = z, dist = "negbin"),
            error = function(e) { NULL }
          )
        }
      }
      
      if (! is.null(fm_zinb)) {
        fit_logmean      = fm_zinb$coefficients$count
        fit_dispersion   = fm_zinb$theta
        if(per_cell_adjust == "NB"){
          fit_logitdropout = c(fm_zinb$coefficients$zero, fill4z)
        }else{
          fit_logitdropout = fm_zinb$coefficients$zero
        }
        
        if(fit_dispersion > max.theta){ fit_dispersion = max.theta }
        
        fit_total = list(fit_logmean, fit_dispersion, fit_logitdropout)
      }else{
        fit_total = NULL
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
          MASS::glm.nb(fmla.both, data=z),
          error = function(e) { NULL }
        )
      }
      if (! is.null(fm_nb)) {
        fit_logitdropout = c(log(MIN_VAL), fill4z)
        fit_logmean      = fm_nb$coefficients
        fit_dispersion   = fm_nb$theta
        fit_total = list(fit_logmean, fit_dispersion, fit_logitdropout)
      }else{
        fit_total = NULL
      }
      
    }
  }
  if(is.null(fit_total)){
    fit_total = list(NA, NA, NA)
  }
  names(fit_total) = c("logmean", "dispersion", "logitdropout")
  
  return(fit_total)
}
