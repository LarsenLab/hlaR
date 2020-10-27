# hlaR
## install package
* library(devtools)   
* install_github("LarsenLab/hlaR")
## load package
* library(hlaR) 
## call eplet mis-match function
* re <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
## check the result
* re$results_count  
* re$results_detail
