#' many-to-many mis-match calculation, find out best match based on mm-count
#' @param dat_in
#' data table with donor/recipient MHC-I info, each row has unique id
#' @param id_don
#' donor ids
#' @param id_rcpt
#' recipient ids
#' @return
#' list of data tables.
#' count: original input data appended with count mis-matched eplet of each pair
#' detail: detailed mis-match eplet of each subject id, plus count and percentage of mis-match across all pairs
#' @export
#'
#' @import
#' devtools
#' tidyverse
#' dplyr
#' schoolmath
#' tibble
#' stringr
#' purrr
#' tidyr
#'
#' @examples
#' \dontrun{
#' dat_in <- read.csv("~/projects/DEV/hlaR_dev/MHC_I_n2n.csv", sep = ",", header = TRUE)
#' re <- CalEpletMHCI_best(dat_in = dat_in,
#' id_don = c(8), id_rcpt = c(6,1))
#' }
#'

CalEpletMHCI_best <- function(dat_in, id_don, id_rcpt) {

  #* pre-step : check input parameters for overlap donor and recipient id, mis-match donor_type info *#
  don <- dat_in %>% filter(id %in% id_don)
  rcpt <- dat_in %>% filter(id %in% id_rcpt)

  if (id_don %in% id_rcpt)
    stop("donor and recipient ids have overlaps, please check.")

  if(unique(don$donor_type) %in% c("recipient", "recip", "rcpt", "r"))
    stop("donor has type 'recipient', please check you donor id")

  if(unique(rcpt$donor_type) %in% c("donor", "don", "dn","d"))
    stop("recipient has type 'donor', please check your recipient id")
  #* end of pre-step *#

  #* step 1: import raw eplet table *#
  raw_eplet <- read.csv(system.file("extdata", "MHC_I_eplet.csv", package = "hlaR"), check.names = FALSE)

  raw_lookup <- as.data.frame(t(raw_eplet)) %>%
    setNames(paste(raw_eplet$type, raw_eplet$index, sep = "_" )) %>%
    rownames_to_column(var = "locus") %>%
    mutate(locus = ifelse(str_detect(locus, "\\*"), sub("\\*.*", "", locus), locus)) %>%
    filter(!locus %in% c("index", "type") ) %>%
    distinct() %>%
    reshape2::melt(id.vars = "locus") %>%
    filter(value != "" ) %>%
    distinct() %>%
    mutate(index = as.numeric(sub(".*\\_", "", variable)),
           type = sub("\\_.*", "", variable)) %>%
    dplyr::rename(eplet = value) %>%
    select(index, type, eplet) %>%
    distinct()
}







