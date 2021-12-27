# ideas
 Individual level Differential Expression Analysis for Single cells

This is the R package for differential expression analysis using single cell RNA-seq data  of  multiple individuals. The input are scRNA-seq data  and cell level and/or individual-level covariates.  


## Installation 
 To install this package in R, use 
 
 ```
    library("devtools");
    install_github("Sun-lab/ideas")
 ```

## Usage

Here is the example code to run IDEAS using simulated data. First load libraries and simulated data. Here we took 100 genes for illustration

```
library(ideas)
library(foreach)
library(doRNG)
RNGkind("L'Ecuyer-CMRG")
library(doParallel)
registerDoParallel(cores=6)

simu_data_rds = "sim_data_ncase_10_nctrl_10_ncell_120_fold_mean_1.2_var_1.5.rds"
sim_data      = readRDS(paste0("data/", simu_data_rds))

count_matrix = sim_data$count_matrix[1:100,]
meta_cell    = sim_data$meta_cell
meta_ind     = sim_data$meta_ind

var2test      = "phenotype"
var2adjust    = "RIN"
var2test_type = "binary"
var_per_cell  = "cell_rd"
```

Next we ran the analysis in two steps. First calculate the distance matrix by function ```ideas_dist```, and then evaluate the p-value using function ```permanova```.
```
dist1 = ideas_dist(count_matrix, meta_cell, meta_ind, 
                   var_per_cell, var2test, var2adjust, 
                   var2test_type, d_metric = "Was", 
                   fit_method = "nb")

pval_ideas = permanova(dist1, meta_ind, var2test, var2adjust, 
  var2test_type, n_perm=999, r.seed=903)
```

## Reference

M Zhang, S Liu et al. (2021) Individual-level Differential Expression Analysis for Single Cell RNA-seq data
