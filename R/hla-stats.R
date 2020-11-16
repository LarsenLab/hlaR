
#' top N frequency of allele(s)
#' @param dat_in
#' dataframe of clean HLA
#' @param names_in
#' allele name(s) of which we want to count frequency
#' @param top_n
#' names of top N most frequent alleles, default is 5
#' @return
#' dataframe of top_n most frequent alleles
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
#  dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR")
#  names <- c("recipient.a1", "recipient.a2", "donor.a1","donor.a2")
#' result <- CountFreq(dat_in = dat), names_in =  names, top_n = 2)
#' }

CountFreq <- function(dat_in, names_in, top_n = 5){
  names <- syms(names_in)
  dat_out <- dat_in %>%
              select(!!!names) %>%
              mutate_all(as.character) %>%
              gather(., name, allele, c(1:length(names)) ) %>%
              group_by(name, allele) %>%
              summarise(freq = n()) %>%
              filter(!is.na(allele)) %>%
              top_n(top_n) %>%
              ungroup() %>%
              arrange(name, allele, freq)

  return(dat_out)
}


