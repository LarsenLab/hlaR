
#' Title
#' function to count number of mis-match HLAs
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
# hla <- read_csv("~/projects/hlaR/inst/extdata/HLA_MisMatch_count_test.csv")
#' classI <- CountMism(hla, c("mism.a1", "mism.a2", "mism.b1", "mism.b2"))
#' classII <- CountMism(hla, c("mism.dqa12", "mism.dqb11", "mism.dqb12"  ))
#' }

CountMism <- function(dat_in, names_in){
  names <- syms(names_in)
  dat_out <- dat_in %>% select(!!!names)
  dat_out[dat_out != 1] <- 0
  dat_out <- dat_out %>% mutate(mism_cnt = rowSums(.))
  return(dat_out)
}
