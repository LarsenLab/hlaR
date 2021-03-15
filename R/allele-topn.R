#' @name CalAlleleTopN
#' @title topN most frequent recipient/donor alleles. N can be defined by user, default n is 5.
#' @param dat_in
#' An cleaned HLA dataframe.
#' @param nms_don
#' A vector of donor's allele name(s).
#' @param nms_rcpt
#' A vector of recipient's allele name(s).
#' @param top_n
#' Number of top frequent alleles. Default n = 5.
#' @return
#' A dataframe of top_n most frequent alleles.
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
              filter(!is.na(allele) & allele != "") %>%
              group_by(part_type, allele) %>%
              summarise(freq = n(), .groups = 'drop') %>%
              ungroup() %>%
              group_by(part_type) %>%
              top_n(top_n, wt = freq) %>%
              ungroup() %>%
              arrange(part_type, -freq)

  return(dat_out)
}

