## Installation
* library(devtools)   
* install_github("LarsenLab/hlaR")
* library(hlaR) 

## Usage example
### eplet mis-match
#### MHC class I
* eplet_mm1 <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))
* eplet_mm1$count  
* eplet_mm1$detail
#### MHC class II
* eplet_mm2 <- CalEpletMHCII(system.file("extdata", "MHC_II_test.csv", package = "hlaR"))
* eplet_mm2$count  
* eplet_mm2$detail

### HLA clean and mis-match
#### clean
* library(readr)
* clean <- read_csv(system.file("extdata", "HLA_Clean_test.csv", package = "hlaR"))
* clean1 <- CleanHla(clean$RECIPIENT_A1, clean$RECIPIENT_A2, locus = "a")
* clean2 <- CleanHla(clean$DONOR_DRB11, clean$DONOR_DRB11, locus = "drb")
#### mis-match
* hla_mm_eval <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
* a <- EvalMism(hla_mm_eval, hla_mm_eval$donor.a1, hla_mm_eval$donor.a2, hla_mm_eval$recipient.a1, hla_mm_eval$recipient.a2, "a")
* hla_mm_eval$mism.a1 <- a$mism_1
* hla_mm_eval$mism.a2 <- a$mism_2
* hla_mm_eval$mism.b1 <- unlist(EvalMism(hla_mm_eval, hla_mm_eval$donor.b1, hla_mm_eval$donor.b2, hla_mm_eval$recipient.b1,hla_mm_eval$recipient.b2, "b")[1])
* hla_mm_eval$mism.b2 <- unlist(EvalMism(hla_mm_eval, hla_mm_eval$donor.b1, hla_mm_eval$donor.b2, hla_mm_eval$recipient.b1,hla_mm_eval$recipient.b2, "b")[2])
#### count of mis-match
* hla_mm_cnt <- read_csv(system.file("extdata", "HLA_MisMatch_count_test.csv", package = "hlaR"))
* classI <- CountMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
* classII <- CountMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"  ))
#### most frequent alleles
* dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
* nms <- c("recipient.a1", "recipient.a2", "donor.a1","donor.a2")
* result <- CountFreq(dat_in = dat, names_in = nms, top_n = 2)
#### most frequent donor allele that is mis-match to recipients
* dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
* don <- c("donor.a1", "donor.a2")
* rcpt <- c("recipient.a1", "recipient.a2")
* result <- CalMismFreq(dat_in = dat, names_don = don, names_rcpt = rcpt) 

