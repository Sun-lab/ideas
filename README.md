# ideas
 Individual level Differential Expression Analysis for Single cells

This is the R package for differential expression analysis using single cell RNA-seq data  of  multiple individuals. The inputs are scRNA-seq data  and cell level and/or individual-level covariates and the outputs are p-values for all genes tested. This project is licensed under the terms of the MIT license.


## Installation 
 To install this package in R, use 
 
 ```
    library("devtools");
    install_github("Sun-lab/ideas")
 ```

## Usage

Here is the example code to run IDEAS using simulated data. First load libraries and simulated data. Here we took 100 genes for illustration. A complete code can be found [here](https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R)

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

From the above usage example, we can see the two functions that user need to use are ```ideas_dist```, which calculate distance across all individuals, and ```permanova```, which calculate the testing p-values given the distance matrix. Here we give a brief description of the input and output of these two functions. 

### ```ideas_dist```

The output of  ```ideas_dist``` is a three dimensional array with first dimension for the number of genes and the next two dimensions for the the number of individuals. For example, if we study 1000 genes and 20 individuals, it is an array of dimension 1000 x 20 x 20. Some parameters of ```ideas_dist``` that often need to be set by the users are listed below.

- ```fit_method```: The method used to estimate the distribution of gene expression across all the cells per individual. If the input is UMI count data, we recommend  "nb" that stands for negative binomial. If the input is denoised scRNA-seq data by DCA or SAVER, ```fit_method``` should be set as "dca_direct" or "saver_direct", respectively. 

- ```count_input```: The input data to be used to calculate distance arrays. If fit_method is "nb", "zinb", or "kde", the count_input should be a matrix of RNAseq counts, with rows for genes and columns for cells. Row names (unique gene ids) and column names (unique cell ids) are required. If fit_method is "dca_direct", the count_input is a list of 3 matrices for mean, over-dispersion, and zero inflation proportion, respectively. If fit_method is "saver_direct", the count_input is a matrix of Poisson mean values after adjusting cell level read-depth.

- ```meta_cell```: A data.frame of meta information of all the cells. The rows of meta_cell should be one to one correspondence to the columns of count_input. Three columns are required: "cell_id" and "individual" for cell id and individual labels, respectively, and a column for cell level read-depth specified by cell_rd_var.

- ```meta_ind```: A data.frame of meta information of all the individuals. Each row of meta_ind corresponds to an individual. At least two columns are required: "individual" for individual labels, and a column specified by var2test for the variable against which to test DE.

- ```var_per_cell```: The variables used for cell level adjustment. They should be included in meta_cell. Cell level read-depth should be included, and other variables such as percentage of ribosome expression or mitochondria expression could be considered too. NOTE: all the variables included will be log-transformed when used as covariates in zinb regression (when fit_method="zinb") or linear regression (when fit_method="kde").

- ```var2test```: A string specifying the name of the variable against which to test DE. This variable should be included in meta_ind.

- ```var2test_type```: "binary" or "continuous".

### ```permanova```

```permanova``` take the distance matrix as input and its output is a vector of p-values for each gene. Most other inputs of ```permanova``` are the same as the inputs for ```ideas_dist```, such as information for cells (```meta_cell```) and individuals (```meta_ind```). 

## Reference

M Zhang, S Liu et al. (2021) IDEAS: Individual Level Differential Expression Analysis for Single Cell RNA-seq data
