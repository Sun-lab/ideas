\name{permanova}
\alias{permanova}
\title{
Permutational multivariate analysis of variance (PERMANOVA)
}
\description{
PERMANOVA for differential expression analysis of scRNA-seq data. 
}
\usage{
permanova(dist_array, meta_ind, var2test, var2adjust = NULL, 
var2test_type = c("binary", "continuous"), n_perm = 999, 
r.seed = 2020, residulize.x = FALSE, delta = 0.5)
}
\arguments{
  \item{dist_array}{
  A three dimensional array with first dimension for the number of genes and the next two dimensions for the the number of individuals. For example, if we study 1000 genes and 20 individuals, it is an array of dimension 1000 x 20 x 20. It can be calculated by function \code{ideas_dist}. 
}
  \item{meta_ind}{
  A data.frame of meta information of all the individuals. Each row of \code{meta_ind} corresponds to an individual. At least two columns are required: "individual" for individual labels, and a column specified by \code{var2test} for the variable against which to test DE. The rows of \code{meta_ind} should be one to one correspondence to the rows/columns of each distance matrix dist_array[k,,], where k is the index for the k-th gene. 
}
  \item{var2test}{
   A string specifying the name of the variable against which to test DE. This variable should be included in \code{meta_ind}. 
}
  \item{var2adjust}{
   A character vector specifying the name of the variables included as covariates in DE analysis. They should be included in \code{meta_ind}. 
}
  \item{var2test_type}{
  The type of the variable against which to test DE. 
}
  \item{n_perm}{
  The number of permutations. 
}
  \item{r.seed}{
  Random seed.
}
  \item{residulize.x}{
  To permute \code{var2test}, we can simply permute it or simulate it based on its conditional distribution with respect to \code{var2adjust}. If \code{residulize.x} is TRUE, we simulate it from conditional distribution and permute it otherwise. 
}
  \item{delta}{
  When sample size is small, a variable simulated based on conditional distribution with respect to \code{var2adjust} may have very different variance than the observed one. We will discard the permuted variables if its standard deviation is smaller than (1-\code{delta}) or larger than (1+\code{delta}) times of the stand deviation of the un-permuted variable. Usually this only exclude less than one percent of the permutated variables. 
}
}
\details{
  When there is no covariate and the \code{var2test} is a binary variable, we simply run the standard PERMANOVA by calculating the summation of the between and within group distance. Otherwise we use the more general formula of F = tr(HGH)/tr[(I - H) G (I-H)], where G is the centralized similarity matrix and H is the hat matrix for the space spanned by \code{var2test} and \code{var2adjust}. In other words, we calculated the F-statistics to assess the contribution of both \code{var2test} and \code{var2adjust} on the distance matrix. 
}
\value{
A vector of p-values for each gene.
}
\references{
McArdle, BH and Anderson, MJ (2001). Fitting multivariate models to community data: a comment on distance-based redundancy analysis. Ecology 82(1), 290-297. 

Tang, Z. Z., Chen, G., & Alekseyenko, A. V. (2016). PERMANOVA-S: association test for microbial community composition that accommodates confounders and multiple distances. Bioinformatics, 32(17), 2618-2625.
}
\author{
Wei Sun, Mengqi Zhang
}
\seealso{
\code{\link{ideas_dist}}
}
\keyword{ methods }
