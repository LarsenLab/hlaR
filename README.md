## Installation
library(devtools)<br>
install_github("LarsenLab/hlaR")<br>
library(hlaR)<br> 

## Usage example
### Allele clean and mis-match
#### - clean
library(readr)<br>
clean <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))<br>
clean1 <- CleanAllele(clean$RECIPIENT_A1, clean$RECIPIENT_A2)<br>
clean2 <- CleanAllele(clean$DONOR_DRB11, clean$DONOR_DRB12)<br>

#### - mis-match
dat <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))<br>
a <- EvalAlleleMism(dat$DONOR_A1, dat$DONOR_A2, dat$RECIPIENT_A1, dat$RECIPIENT_A2)<br>
a<br>
b <- EvalAlleleMism(dat$DONOR_B1, dat$DONOR_B2, dat$RECIPIENT_B1, dat$RECIPIENT_B2)<br>
b<br>

### haplotype
dat <- read.csv(system.file("extdata/example", "Haplotype_test.csv", package = "hlaR"))<br>
re <- CompHaploTbl(dat_in = dat[1:5,])<br>

### single molecule level eplet mis-match
#### - MHC class I
dat <- read.csv(system.file("extdata/example", "MHC_I_test.csv", package = "hlaR"), sep = ",", header = TRUE)<br>
eplet_mm1_v2 <- CalEpletMHCI_s(dat, ver = 2)<br>
head(eplet_mm1_v2)<br>
eplet_mm1_v3 <- CalEpletMHCI_s(dat, ver = 3)<br>
(or simply eplet_mm1_v3 <- CalEpletMHCI(dat) )<br>
head(eplet_mm1_v3)
#### - MHC class II
dat <- read.csv(system.file("extdata/example", "MHC_II_test.csv", package = "hlaR"), sep = ",", header = TRUE)<br>
eplet_mm2_v2 <- CalEpletMHCII_s(dat, ver = 2)<br>
head(eplet_mm2_v2)<br>
eplet_mm2_v3 <- CalEpletMHCII_s(dat, ver = 3)<br>
(or simply eplet_mm2_v3 <- CalEpletMHCII(dat) )<br>
head(eplet_mm2_v3)

### other functionalities
#### - count of mis-match
hla_mm_cnt <- read.csv(system.file("extdata/example", "HLA_MisMatch_count_test.csv", package = "hlaR"))<br>
classI <- CountAlleleMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))<br>
classII <- CountAlleleMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"))<br>
#### - topN most frequent recipient/donor alleles 
dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
result <- CalAlleleTopN(dat_in = dat, nms_don = don, nms_rcpt = rcpt, top_n = 2)<br>
result<br>
#### - frequency(freq count > 1) of donor mis-match alleles to recipients
dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
result <- CalAlleleMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)<br> 
result



