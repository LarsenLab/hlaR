## install and load package
* library(devtools)   
* install_github("LarsenLab/hlaR")
* library(hlaR) 

## eplet
### mis-match, MHC class I
* eplet_mm1 <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
* eplet_mm1$count  
* eplet_mm1$detail
### mis-atc, MHC class II
* eplet_mm2 <- CalEpletMHCII(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
* eplet_mm2$count  
* eplet_mm2$detail

## HLA cleaning
### clean
* clean <- read_csv(system.file("extdata", "HLA_Clean_test.csv", package = "hlaR"))
* clean1 <- CleanHla(clean$RECIPIENT_A1, clean$RECIPIENT_A2, locus = "a")
* clean2 <- CleanHla(clean$DONOR_DRB11, clean$DONOR_DRB11, locus = "drb")
### mis-match
* hla_mm_eval <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
* a <- EvalMism(hla_mm_eval, hla_mm_eval$DONOR_A1, hla_mm_eval$DONOR_A2, hla_mm_eval$RECIPIENT_A1, hla_mm_eval$RECIPIENT_A2, "a")
* hla$mism.a1 <- a$mism_1
* hla$mism.a2 <- a$mism_2
* hla$mism.b1 <- unlist(EvalMism(hla_mm_eval, hla_mm_eval$DONOR_B1, hla_mm_eval$DONOR_B2, hla_mm_eval$RECIPIENT_B1, hla_mm_eval$RECIPIENT_B2, "b")[1])
* hla$mism.b2 <- unlist(EvalMism(hla_mm_eval, hla_mm_eval$DONOR_B1, hla_mm_eval$DONOR_B2, hla_mm_eval$RECIPIENT_B1, hla_mm_eval$RECIPIENT_B2, "b")[2])
### count of mis-match
* hla_mm_cnt <- read_csv(system.file("extdata", "HLA_MisMatch_count_test.csv", package = "hlaR"))
* classI <- CountMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
* classII <- CountMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"  ))
