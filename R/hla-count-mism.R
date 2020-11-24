
#' count number of mis-match HLAs
#' @param dat_in
#' a dataframe with HLAs mis-match flags
#' @param names_in
#' column names of which we want count mis-matches
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
#' classI <- CountMism(hla, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
#' classII <- CountMism(hla, c("mism.dqa12", "mism.dqb11", "mism.dqb12"  ))
#' }

CountMism <- function(dat_in, names_in){
  names <- syms(names_in)
  len <- length(names) + 1
  dat_out <- dat_in %>%
              select(1, !!!names) %>%
              mutate(num_nas = apply(is.na(.), 1, sum)) %>%
              mutate(mism_total = ifelse(num_nas == length(names), NA, rowSums(.[2:len], na.rm = T))) %>%
              select(-num_nas)

  return(dat_out)
}
