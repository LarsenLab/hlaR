#' @name CountAlleleMism
#' @title Count HLA mismatch at the allele level
#' @description Donor and recipient HLA typing data is compared to determine allele level mismatch. The output of EvalAlleleMism is used as input for this function. Allele level mismatch can be calculated for both high and low resolution data. The generated count will return 'NA' if the input alleles are 'NA.'
#' @param dat_in
#' A dataframe with donor and recipient mismatched alleles as output from EvalAlleleMism function.
#' @param names_in
#' A column containing the names of HLA loci for which to count mismatch.
#' @return
#' The input data, 'dat_in' is returned with the count of mismatches appended as the last column.
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
