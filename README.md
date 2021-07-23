## Installation
### CRAN install
install.packages("hlaR")<br>
### github install<br>
library(devtools)<br>
install_github("LarsenLab/hlaR")<br>
library(hlaR)<br> 

## Usage example
### Allele clean and mis-match
#### - clean
library(readr)<br>
clean <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))<br>
clean1 <- CleanAllele(clean$recipient_a1, clean$recipient_a2)<br>
clean2 <- CleanAllele(clean$donor_a1, clean$donor_a2)<br>

#### - mis-match
dat <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))<br>
mm1 <- EvalAlleleMism(dat$donor_a1, dat$donor_a2, dat$recipient_a1, dat$recipient_a2)<br>
mm1 <br>
mm2 <- EvalAlleleMism(dat$donor_b1, dat$donor_b2, dat$recipient_b1, dat$recipient_b2)<br>
mm2<br>

### haplotype
dat <- read.csv(system.file("extdata/example", "Haplotype_test.csv", package = "hlaR"))<br>
re <- ImputeHaplo(dat_in = dat)<br>

### single molecule level eplet mis-match
#### - MHC class I
dat <- read.csv(system.file("extdata/example", "MHC_I_test.csv", package = "hlaR"), sep = ",", header = TRUE)<br>
eplet_mm1_v2 <- CalEpletMHCI(dat, ver = 2)<br>
single_detail <- eplet_mm1_v2$single_detail<br>
overall_count <- eplet_mm1_v2$overall_count<br>

eplet_mm1_v3 <- CalEpletMHCI(dat, ver = 3)<br>
(or simply eplet_mm1_v3 <- CalEpletMHCI(dat) )<br>
single_detail <- eplet_mm1_v3$single_detail<br>
overall_count <- eplet_mm1_v3$overall_count<br>
#### - MHC class II
dat <- read.csv(system.file("extdata/example", "MHC_II_test.csv", package = "hlaR"), sep = ",", header = TRUE)<br>
eplet_mm2_v2 <- CalEpletMHCII(dat, ver = 2)<br>
single_detail <- eplet_mm2_v2$single_detail<br>
overall_count <- eplet_mm2_v2$overall_count<br>
eplet_mm2_v3 <- CalEpletMHCII(dat, ver = 3)<br>
(or simply eplet_mm2_v3 <- CalEpletMHCII_s(dat) )<br>
single_detail <- eplet_mm2_v3$single_detail<br>
overall_count <- eplet_mm2_v3$overall_count<br>

### other functionalities
#### - count of mis-match
hla_mm_cnt <- read.csv(system.file("extdata/example", "HLA_MisMatch_count_test.csv", package = "hlaR"))<br>
classI <- CountAlleleMism(hla_mm_cnt, c("mism_a", "mism_b"))<br>
classII <- CountAlleleMism(hla_mm_cnt, c("mism_drb1", "mism_dqa", "mism_dqb"))<br>
#### - topN most frequent recipient/donor alleles 
dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
re <- CalAlleleTopN(dat_in = dat, nms_don = don, nms_rcpt = rcpt, top_n = 2)<br>
re<br>
#### - frequency(freq count > 1) of donor mis-match alleles to recipients
dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))<br>
don <- c("donor.a1", "donor.a2")<br>
rcpt <- c("recipient.a1", "recipient.a2")<br>
re <- CalAlleleMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)<br> 
re

## ToDo/Discuss CRAN 0.1.1<br>
- ImputeHaplo(): error check on none ":" punctuation (br)correction: only apply to loci columns, skip check on pair_id/subject_type/ethnicity as pair_id may contain some special symbols <br>
~~- CleanAllele(): incorrect logic for v1.2/v2.2 <br>
correction: 1.line#5-57: use letters_only(); adjust statement ifelse(grepl("[^A-Za-z]+$", v1.1), v1.1, "") <br>
code available in e_txki, hlaRLocalFunctions.R, starting from line #9 <br>
2.remove space within string in addionn to around string <br>
- add alleleclean to etxki pipeline for messy or clean data, check hla_mm_May12.Rmd in etx repo <br>
- haplotype final table: add a flag if the max count doesn't reach the number of unique low res antigens <br>
- haplotype final table: if only 1 record imputed, what to do for the other half in the final table? keep this 1 record only, or duplicate it? ( ex: 13982771, 3163) <br>
- haplotype: additional function to deal with if max-count != count-of-low-res, the replace none-low-res with most-commom-low-res
3. fix memory issue given too many NA hlas during imputation calcualtion (code are available in folder "4nextversion")
4. haplotype for NA race
5. discussion: what to do if low-res part of imputed hla?



