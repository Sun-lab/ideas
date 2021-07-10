\name{ideas_dist}
\alias{ideas_dist}
\title{
  Calculate distance array for all the genes and all pairs of individuals.
}
\description{
  Calculate distance array for Individual level Differential Expression Analysis for Single cells (IDEAS). 
}
\usage{
ideas_dist(count_input, meta_cell, meta_ind, var_per_cell, 
var2test, var2test_type = c("binary", "continuous"), 
d_metric = c("Was", "JSD"), 
fit_method = c("nb", "zinb", "kde", "dca_direct", "saver_direct"), 
per_cell_adjust = c("NB", "both"), quantile = 0.975, 
empirical_n = 1024)
}
\arguments{
  \item{count_input}{
 The input data to be used to calculate distance arrays. If \code{fit_method} is "nb", "zinb", or "kde", the \code{count_input} should be a matrix of RNAseq counts, with rows for genes and columns for cells. Row names (uniqe gene ids) and column names (unique cell ids) are required. If \code{fit_method} is "dca_direct", the \code{count_input} is a list of 3 matrices for mean, over-dispersion, and zero inflation proportion, respectively. If \code{fit_method} is "saver_direct", the \code{count_input} is a matrix of Poisson mean values after adjusting cell level read-depth. 
}
  \item{meta_cell}{
  A data.frame of meta information of all the cells. The rows of meta_cell should be one to one correpondence to the columns of \code{count_input}. Three columns are required: "cell_id" and "individual" for cell id and individual labels, respectively, and a column for cell level read-depth specified by \code{cell_rd_var}. 
}
  \item{meta_ind}{
  A data.frame of meta information of all the individuals. Each row of \code{meta_ind} corresponds to an individual. At least two columns are required: "individual" for individual labels, and a column specified by \code{var2test} for the variable against which to test DE.  
}
  \item{var_per_cell}{
  The variables used for cell level adjustment. They should be inlcuded in \code{meta_cell}. Cell level read-depth should be included, and other varibles such as percentage of ribosom expression or mitochodria expression could be considered too. NOTE: all the variables included will be log-transformed when used as covariates in zinb regression (when \code{fit_method}="zinb") or linear regression (when \code{fit_method}="kde").
}
  \item{var2test}{
   A string specificying the name of the variable against which to test DE. This variable should be inlcuded in \code{meta_ind}. 
}
  \item{var2test_type}{
  The type of the variable against which to test DE. 
}
  \item{d_metric}{
  The metric to be used to calculate distance between two individuals. JSD: Jensen-Shannon divergence, Was: Wasserstein distance. If \code{fit_method}="zinb", the Wasserstein distance is simply sum(abs(F(x) - F(y))), where F(x) and F(y) are cdf of ZINB at the same grids of 0, 1, 2, .... If \code{fit_method}="kde", we use the observed data (after removing effects due to \code{var_per_cell}, without kernel density estimation) to calcualte the Wasserstein distance, by function wasserstein1d of R package transport. 
}
  \item{fit_method}{
  The method used to estimate the distribution of gene expression across all the cells per individual. "nb": negative binomial distribution; "zinb": zero-inflated negative binomial. If there is no 0 in the count, it will automatically shift to negative binomial distribution. "kde": kernel density estimation. 
}
  \item{per_cell_adjust}{
  If \code{fit_method} is "zinb", \code{per_cell_adjust} for cell level adjustment, whether the varaibles included in \code{var_per_cell} should be used for NB regression only or for both NB and logistic regression. For most UMI-based scRNA-seq data, the zero-inflated component is often very low or not needed, and thus performing the logistic regression for the zero-inflation part can be highly unstable. Thus it is recommended to set \code{per_cell_adjust} = "NB" for UMI-based scRNA-seq data. 
}
  \item{quantile}{
  If \code{fit_method} is "dca_direct" or "saver_direct", \code{quantile} specifies the maximum quantile where the density will be evaluated. 
}
  \item{empirical_n}{
  If \code{fit_method} is "kde", \code{empirical_n} specifies the number of points where the density will be evaluated. 
}
}
\details{
  To calculate the distance array, we do not need to use the data in \code{meta_ind}. It is used to make sure the distance array aligns with the \code{meta_ind} for the next step analysis. 
}
\value{
A three dimensional array with first dimension for the number of genes and the next two dimensions for the the number of individauls. For example, if we study 1000 genes and 20 individuals, it is an arrary of dimension 1000 x 20 x 20. 
}
\references{
Mengqi Zhang*, Si Liu*, Zhen Miao, Fang Han, Raphael Gottardo, Wei Sun, Individual Level Differential Expression Analysis for Single Cell RNA-seq data, bioRxiv 2021.05.10.443350; doi: https://doi.org/10.1101/2021.05.10.443350  (* co-first author)
}
\author{
Wei Sun, Si Liu, Mengqi Zhang
}
\keyword{ methods }
