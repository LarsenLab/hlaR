#' @name CountAlleleMism
#' @title Count number of mismatch alleles.
#' @param dat_in
#' A dataframe with allele's mismatch flag.
#' @param names_in
#' A column names to count mismatch.
#' @return
#' dat_in with count of mismatches appended on the last column. Count is NA if input alleles are NA.
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# hla <- read_csv(system.file("extdata/example", "HLA_MisMatch_count_test.csv", package = "hlaR"))
#' classI <- CountAlleleMism(hla_mm_cnt, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
#' classII <- CountAlleleMism(hla_mm_cnt, c("mism.dqa12", "mism.dqb11", "mism.dqb12"))
#' }

CountAlleleMism <- function(dat_in, names_in){
  names <- syms(names_in)
  len <- length(names) + 1
  dat_out <- dat_in %>%
              select(1, !!!names) %>%
              mutate(num_nas = apply(is.na(.), 1, sum)) %>%
              mutate(mism_total = ifelse(num_nas == length(names), NA, rowSums(.[2:len], na.rm = T))) %>%
              select(-num_nas)

  return(dat_out)
}
