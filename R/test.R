library(devtools)
install_github("LarsenLab/hlaR")

library(hlaR)
library(readr)
library(tidyverse)

#* hla clean *#
dat <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))
clean1 <- CleanAllele(dat$RECIPIENT_A1, dat$RECIPIENT_A2)
clean2 <- CleanAllele(dat$DONOR_DRB11, dat$DONOR_DRB12)
head(clean1)
head(clean2)

rm( dat, clean1, clean2)
#* hla mis-match *#
dat <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))
mism_a <- EvalAlleleMism(dat$DONOR_A1, dat$DONOR_A2, dat$RECIPIENT_A1, dat$RECIPIENT_A2)
head(mism_a)

mism_b <- EvalAlleleMism(dat$DONOR_B1, dat$DONOR_B2, dat$RECIPIENT_B1, dat$RECIPIENT_B2)
head(mism_b)

rm(dat, mism_a, mism_b)

#* haplotype imputation *#
dat <- read.csv(system.file("extdata/example", "Haplotype_test.csv", package = "hlaR"))
re <- CompHaploTbl(dat_in = dat[1:5,])
re

#* eplet MHC-I *#
dat <- read.csv(system.file("extdata/example", "MHC_I_test.csv", package = "hlaR"), sep = ",", header = TRUE)
eplet_mm1_v2 <- CalEpletMHCI(dat, ver = 2)
head(eplet_mm1_v2)
eplet_mm1_v3 <- CalEpletMHCI(dat, ver = 3)
head(eplet_mm1_v3)

#* eplet MHC-II *#
dat <- read.csv(system.file("extdata/example", "MHC_II_test.csv", package = "hlaR"), sep = ",", header = TRUE)
eplet_mm2_v2 <- CalEpletMHCII(dat, ver = 2)
head(eplet_mm2_v2)
eplet_mm2_v3 <- CalEpletMHCII(dat, ver = 3)
head(eplet_mm2_v3)

#* other functions *#
# count of mis-match
hla_mm_cnt <- read.csv(system.file("extdata/example", "HLA_MisMatch_count_test.csv", package = "hlaR"))
classI <- CountAlleleMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
head(classI)
classII <- CountAlleleMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"))
head(classII)

# topN most frequent recipient/donor alleles
dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))
don <- c("donor.a1", "donor.a2")
rcpt <- c("recipient.a1", "recipient.a2")
result <- CalAlleleTopN(dat_in = dat, nms_don = don, nms_rcpt = rcpt, top_n = 2)
result

# frequency(freq count > 1) of donor mis-match alleles to recipients
dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))
don <- c("donor.a1", "donor.a2")
rcpt <- c("recipient.a1", "recipient.a2")
result <- CalAlleleMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)
result

