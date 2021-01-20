#' @name CalAlleleTopN
#' @title calculate topN most frequent recipient/donor alleles
#' @param dat_in
#' dataframe of clean HLA
#' @param nms_don
#' donor's allele name(s) to count for frequency
#' @param nms_rcpt
#' recipient's allele name(s) to count frequency count
#' @param top_n
#' names of top N most frequent alleles, default is 5
#' @return
#' a dataframe of top_n most frequent alleles
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
#  dat <- read_csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))
#  don <- c("donor.a1", "donor.a2")
#  rcpt <- c("recipient.a1", "recipient.a2")
#' result <- CalAlleleTopN(dat_in = dat, nms_don = don, nms_rcpt = rcpt, top_n = 2)
#' }
#' @export

CalAlleleTopN <- function(dat_in, nms_don = c(), nms_rcpt = c(), top_n = 5){
  names <- syms(c(nms_don, nms_rcpt))
  dat_out <- dat_in %>%
              select(!!!names) %>%
              mutate_all(as.character) %>%
              gather(., name, allele, c(1:length(names)) ) %>%
              mutate(part_type = ifelse(name %in% nms_don, "don",
                                        ifelse(name %in% nms_rcpt, "rcpt", "none"))) %>%
              select(-name) %>%
              filter(!is.na(allele)) %>%
              group_by(part_type, allele) %>%
              summarise(freq = n(), .groups = 'drop') %>%
              ungroup() %>%
              group_by(part_type) %>%
              top_n(top_n) %>%
              ungroup() %>%
              arrange(part_type, -freq)

  return(dat_out)
}

