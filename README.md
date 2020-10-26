hlaR package

setup : need to turn the repo to public 
library(devtools) 
install_github("LarsenLab/hlaR") 
library(hlaR) 
# mis-match of MHC class I 
re <- CalEpletMHCI(dat_in = "inst/extdata/MHC_I_test.csv")
