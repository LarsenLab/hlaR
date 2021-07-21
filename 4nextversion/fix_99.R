# date: 07/22/2021
# issue:  shiny app "disconnected from the server" on Imputation tab
# cause: - memory issue, too many NA HLAS causes hugh hapotype combination
# solution: filter by topn for hpl_tp_raw(n = top 20% of rank)
# code: function-for_haplotype-99: insert code on line #104, 112, 137, 160, 183, 206
library(hlaR)
library(tidyverse)

# test on pair_id 101, donor

dat_in <- read.csv("~/projects/hlaR/4nextversion/test.csv")

re <- ImputeHaplo99(dat_in)
