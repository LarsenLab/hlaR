#' @title pull out haplotypes based on NMDP frequency table ...
#' @param dat_in
#' dataframe of alleles
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
#' result <- ht(dat_in = dat)
#' }
#' @export
ht <- function(dat_in){
  # ori <- read.csv("~/projects/DEV/haplostats_dev/data/csv/A~C~B~DRB1~DQB1.csv") %>%
  #   select(A, C, B, DRB1, DQB1,
  #          AFA_freq, AFA_rank, API_freq, API_rank,CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank ) %>%
  #   filter_at(vars(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank), all_vars(!is.na(.))) %>%
  #   # remove leading loci letter
  #   mutate(A = str_sub(A, start = 3),
  #          B = str_sub(B, start = 3),
  #          C = str_sub(C, start = 3),
  #          DRB1 = str_sub(DRB1, start = 6),
  #          DQB1 = str_sub(DQB1, start = 6)) %>%
  #   # remove trailing g if exists
  #   mutate(A = ifelse(stri_sub(A, -1) == "g", str_sub(A, end = -2), A),
  #          B = ifelse(stri_sub(B, -1) == "g", str_sub(B, end = -2), B),
  #          C = ifelse(stri_sub(C, -1) == "g", str_sub(C, end = -2), C),
  #          DRB1 = ifelse(stri_sub(DRB1, -1) == "g", str_sub(DRB1, end = -2), DRB1),
  #          DQB1 = ifelse(stri_sub(DQB1, -1) == "g", str_sub(DQB1, end = -2), DQB1))
  #
  #
  # dat <- read.csv("~/projects/DEV/haplostats_dev/data/csv/haplotype_test.csv") %>%
  #   # mutate(ck = ifelse(grepl("\\D", sub(".*\\:", "", A)), "", sub(".*\\:", "", A))) %>%
  #   # if last 2 chars are not number, then remove it so only keep first 2 digits of the allele kept
  #   mutate(A = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", A)), "", sub(".*\\:", "", A)) != "", A, str_sub(A, end = -4)),
  #          B = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", B)), "", sub(".*\\:", "", B)) != "", B, str_sub(B, end = -4)),
  #          C = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", C)), "", sub(".*\\:", "", C)) != "", C, str_sub(C, end = -4)),
  #          DRB1 = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", DRB1)), "", sub(".*\\:", "", DRB1)) != "", DRB1, str_sub(DRB1, end = -4)),
  #          DQB1 = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", DQB1)), "", sub(".*\\:", "", DQB1)) != "", DQB1, str_sub(DQB1, end = -4)))
}

