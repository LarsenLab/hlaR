# hlaR
## install package
* library(devtools)   
* install_github("LarsenLab/hlaR")
## load package
* library(hlaR) 
## call eplet mis-match function and check the results
* re <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
* re$results_count  
* re$results_detail
