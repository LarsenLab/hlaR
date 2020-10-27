# hlaR
## install package
* library(devtools)   
* install_github("LarsenLab/hlaR")
## load package
* library(hlaR) 
## eplet mis-match calculation
* re <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
* or 
* re <- CalEpletMHCII(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
* re$results_count  
* re$results_detail
