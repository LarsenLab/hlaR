#' @keywords internal
# purpose: this script is to construct the MHC-II eplet reference table based on matchmaker v2.0, plus some missing known eplets froom eplet Db (33N, 78Y, 78VS(98ES), 108P, 85V)
# date: 03/10/2020
# author: Larsen Lab

library(tidyverse)
library(here)
library(readxlsb)
library(openxlsx)
library(qdapRegex)
library(janitor)

# notes:
# 1. eplet names:
#    - all converted to upper cases (DQA1*04:02 and DQA1*04:04 were 187a, DQA1*05:01 and DQA1*01:02 were 187A. all of them are 187A in current version)
# 2.

#* function to pull out eplets *#
PullEplet <- function(filename, sheetname){
  tbl_ref <-  read_xlsb(filename,
                        sheet = sheetname, skip = 0, trim_ws = TRUE, col_names = TRUE) %>%
              mutate_all(., .funs = toupper) %>% # convert all eplet names to upper cases
              select(-c(names(.)[str_detect(names(.), "X.")])) %>%
              replace(., is.na(.), "") %>%
              mutate(row_na = rowSums(. == "" )) %>% # delete the row if it have all NAs
              filter(row_na < dim(.)[2] - 1) %>%
              select(-row_na) %>%
              setNames(paste(names(.), .[1,], sep = "_")) %>%
              setNames(rm_between(names(.), ".", "_")) %>%
              setNames(str_remove(names(.), "_")) %>%
              dplyr::rename(allele = "X") %>%
              filter(allele != "") %>%
              mutate(across(where(is.character), str_trim))

  nm_ex <- tbl_ref %>%
            summarise_all(~list(sum(. == "")))
  nm_ex <- which(unlist(nm_ex) == dim(tbl_ref)[1])
  out <- tbl_ref %>%
          select(-all_of(nm_ex)) %>%
          mutate_if(is.character, str_trim) # remove leading and trailing white spaces
  return(out)
}

#* pull out eplet table for each allele *#
q <- PullEplet(filename = here("inst/extdata/matchmaker/5DRDQDPMatchingVs2.xlsb"),
               sheetname = "Q")

p <- PullEplet(filename = here("inst/extdata/matchmaker/5DRDQDPMatchingVs2.xlsb"),
               sheetname = "P")

r <- PullEplet(filename = here("inst/extdata/matchmaker/5DRDQDPMatchingVs2.xlsb"),
               sheetname = "R") %>%
  mutate(count = rowSums(. == "" )) %>%
  filter(count < ncol(.) - 2) %>%
  select(-count)

#* As *#
# loci -> count of non-na -> filter if all eplets are na -> select allele and eplets
dqa <- q %>%
        filter(str_detect(allele, "DQA")) %>%
        select(grep("allele|QA", names(.))) %>%
        mutate(count = rowSums(. == "" )) %>%
        filter(count < ncol(.) - 2) %>%
        select(-count)

dqa <- as_tibble(t(dqa), rownames = "type") %>%
  row_to_names(row_number = 1) %>%
  mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
  dplyr::rename(type = allele) %>%
  group_by(type) %>%
  mutate(index = row_number()) %>%
  select(index, type, everything())

# dpa
dpa <- p %>%
        filter(str_detect(allele, "DPA")) %>%
        select(grep("allele|PA", names(.))) %>%
        mutate(count = rowSums(. == "" )) %>%
        filter(count < ncol(.) - 2) %>%
        select(-count)

dpa <- as_tibble(t(dpa), rownames = "type") %>%
  row_to_names(row_number = 1) %>%
  mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
  dplyr::rename(type = allele) %>%
  group_by(type) %>%
  mutate(index = row_number()) %>%
  select(index, type, everything())

# dqa and dpa
mmv2_a <- dqa %>%
  full_join(., dpa, by = c("index", "type")) %>%
  replace(., is.na(.), "") %>%
  arrange(type, index)

# write.csv(a, "inst/extdata/matchmaker/check/mhcII_v2_A.csv", row.names = FALSE)

#* Bs *#
# dqb
dqb <- q %>%
  filter(str_detect(allele, "DQB")) %>%
  select(grep("allele|QB", names(.))) %>%
  mutate(count = rowSums(. == "" )) %>%
  filter(count < ncol(.) - 2) %>%
  select(-count)

dqb <- as_tibble(t(dqb), rownames = "type") %>%
  row_to_names(row_number = 1) %>%
  mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
  dplyr::rename(type = allele) %>%
  group_by(type) %>%
  mutate(index = row_number()) %>%
  select(index, type, everything())

# dpb
dpb <- p %>%
  filter(str_detect(allele, "DPB")) %>%
  select(grep("allele|PB", names(.))) %>%
  mutate(count = rowSums(. == "" )) %>%
  filter(count < ncol(.) - 2) %>%
  select(-count)

dpb <- as_tibble(t(dpb), rownames = "type") %>%
  row_to_names(row_number = 1) %>%
  mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
  dplyr::rename(type = allele) %>%
  group_by(type) %>%
  mutate(index = row_number()) %>%
  select(index, type, everything())

# r
r <- as_tibble(t(r), rownames = "type") %>%
  row_to_names(row_number = 1) %>%
  mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
  dplyr::rename(type = allele) %>%
  group_by(type) %>%
  mutate(index = row_number()) %>%
  select(index, type, everything())

# r + dqb + dpb
 mmv2_b <- r %>%
      full_join(., dqb, by = c("index", "type")) %>%
      arrange(type, index) %>%
      full_join(., dpb, by = c("index", "type")) %>%
      select(-`DQB1*X`) %>% # no such a allele, remove it
      replace(., is.na(.), "") %>%
      arrange(type, index)

 # write.csv(b, "inst/extdata/matchmaker/check/mhcII_v2_B.csv", row.names = FALSE)
#* mmv2_a and mmv2_b have been validated, code above are stable *#

#* add missing eplets to the table *#



