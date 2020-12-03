#' count number of mis-match alleles
#' @param dat_in
#' dataframe with mis-match flags of allele
#' @param names_in
#' column names of which to count mis-matches
#' @return
#' dat_in with count of mismatches appended on the last column
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# hla <- read_csv("~/projects/hlaR/inst/extdata/HLA_MisMatch_count_test2.csv")
#' classI <- CountAlleleMism(hla, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
#' classII <- CountAlleleMism(hla, c("mism.dqa12", "mism.dqb11", "mism.dqb12"  ))
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
