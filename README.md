## Installation
library(devtools)<br>
install_github("LarsenLab/hlaR")<br>
library(hlaR)<br> 

## Usage example
### eplet mis-match
#### - MHC class I
eplet_mm1 <- CalEpletMHCI(system.file("extdata", "MHC_I_test.csv", package = "hlaR"))<br>
eplet_mm1$count<br>
eplet_mm1$detail<br>
#### - MHC class II
eplet_mm2 <- CalEpletMHCII(system.file("extdata", "MHC_II_test.csv", package = "hlaR"))<br>
eplet_mm2$count<br> 
eplet_mm2$detail<br>

### HLA clean and mis-match
#### - clean
library(readr)<br>
clean <- read_csv(system.file("extdata", "HLA_Clean_test.csv", package = "hlaR"))<br>
clean1 <- CleanHla(clean$RECIPIENT_A1, clean$RECIPIENT_A2, locus = "a")<br>
clean2 <- CleanHla(clean$DONOR_DRB11, clean$DONOR_DRB11, locus = "drb")<br>
#### - mis-match
hla_mm_eval <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
a <- EvalMism(hla_mm_eval, hla_mm_eval$donor.a1, hla_mm_eval$donor.a2, hla_mm_eval$recipient.a1, hla_mm_eval$recipient.a2, "a")<br>
hla_mm_eval$mism.a1 <- a$mism_1<br>
hla_mm_eval$mism.a2 <- a$mism_2<br>
hla_mm_eval$mism.b1 <- unlist(EvalMism(hla_mm_eval, hla_mm_eval$donor.b1, hla_mm_eval$donor.b2, hla_mm_eval$recipient.b1,hla_mm_eval$recipient.b2, "b")[1])<br>
hla_mm_eval$mism.b2 <- unlist(EvalMism(hla_mm_eval, hla_mm_eval$donor.b1, hla_mm_eval$donor.b2, hla_mm_eval$recipient.b1,hla_mm_eval$recipient.b2, "b")[2])<br>
#### - count of mis-match
hla_mm_cnt <- read_csv(system.file("extdata", "HLA_MisMatch_count_test.csv", package = "hlaR"))<br>
classI <- CountMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))<br>
classII <- CountMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"))<br>
#### - topN most frequent recipient/donor alleles 
dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
result <- CalFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt, top_n = 2)<br>
result<br>
#### - frequency(freq count > 1) of donor mis-match alleles to recipients
dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
result <- CalMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)<br> 

