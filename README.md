hlaR package

library(devtools). 
<!--- set repo to public for installation  ---> 
install_github("LarsenLab/hlaR")  
library(hlaR). 
<!--- mis-match of MHC class I --->
re <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR")). 
<!--- check result --->. 
re$results_count. 
re$results_detail
