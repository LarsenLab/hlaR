## Installation
library(devtools)<br>
install_github("LarsenLab/hlaR")<br>
library(hlaR)<br> 

## Usage example
### eplet mis-match
#### - MHC class I
dat <- read.csv(system.file("extdata", "MHC_I_test.csv", package = "hlaR"), sep = ",", header = TRUE)<br>
eplet_mm1 <- CalEpletMHCI(dat)<br>
eplet_mm1$count<br>
eplet_mm1$detail<br>
#### - MHC class II
dat <- read.csv(system.file("extdata", "MHC_II_test.csv", package = "hlaR"), sep = ",", header = TRUE)<br>
eplet_mm2 <- CalEpletMHCII(dat)<br>
eplet_mm2$count<br> 
eplet_mm2$detail<br>

### Allele clean and mis-match
#### - clean
library(readr)<br>
clean <- read_csv(system.file("extdata", "HLA_Clean_test.csv", package = "hlaR"))<br>
clean1 <- CleanAllele(clean$RECIPIENT_A1, clean$RECIPIENT_A2, locus = "a")<br>
clean2 <- CleanAllele(clean$DONOR_DRB11, clean$DONOR_DRB11, locus = "drb")<br>
#### - mis-match
hla_mm_eval <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
a <- EvalAlleleMism(hla_mm_eval, hla_mm_eval$donor.a1, hla_mm_eval$donor.a2, hla_mm_eval$recipient.a1, hla_mm_eval$recipient.a2, "a")<br>
hla_mm_eval$mism.a1 <- a$mism_1<br>
hla_mm_eval$mism.a2 <- a$mism_2<br>
hla_mm_eval$mism.b1 <- unlist(EvalAlleleMism(hla_mm_eval, hla_mm_eval$donor.b1, hla_mm_eval$donor.b2, hla_mm_eval$recipient.b1,hla_mm_eval$recipient.b2, "b")[1])<br>
hla_mm_eval$mism.b2 <- unlist(EvalAlleleMism(hla_mm_eval, hla_mm_eval$donor.b1, hla_mm_eval$donor.b2, hla_mm_eval$recipient.b1,hla_mm_eval$recipient.b2, "b")[2])<br>
#### - count of mis-match
hla_mm_cnt <- read_csv(system.file("extdata", "HLA_MisMatch_count_test.csv", package = "hlaR"))<br>
classI <- CountAlleleMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))<br>
classII <- CountAlleleMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"))<br>
#### - topN most frequent recipient/donor alleles 
dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
result <- CalAlleleTopN(dat_in = dat, nms_don = don, nms_rcpt = rcpt, top_n = 2)<br>
result<br>
#### - frequency(freq count > 1) of donor mis-match alleles to recipients
dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
result <- CalAlleleMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)<br> 

### - haplotype
dat <- read_csv(system.file("extdata", "Haplotype_test.csv", package = "hlaR"))<br>
re <- CompHaploTbl(dat_in = dat, cut_p = 0.0001, cut_r = 10)<br>
check recipient and donor of subject id 116 <br>
re$rcpt_116<br>
re$don_116<br>

