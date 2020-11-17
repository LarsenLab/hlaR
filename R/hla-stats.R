
#' tcalculate op N frequency of allele(s)
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
#' result <- CalFreq(dat_in = dat), names_in =  names, top_n = 2)
#' }

CalFreq <- function(dat_in, names_in, top_n = 5){
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

#' calculate mis-match frequency of donor's alleles to recipient's
#' @param dat_in
#' dataframe of clean HLA
#' @param names_don
#' column names of donor's alleles, must be length of 2
#' @param names_rcpt
#' column names of recipient's alleles, must be length of 2
#' @return
#' dataframe of donor's mis-match alleles with frequent > 1
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
# don <- c("donor.a1", "donor.a2")
# rcpt <- c("recipient.a1", "recipient.a2")
#' result <- CalMismFreq(dat_in = dat, names_don = don, names_rcpt = rcpt)
#' }
#'
#'
#'
CalMismFreq <- function(dat_in, names_don, names_rcpt){
  #* start of data prep *#
  dat_don <- dat_in %>%
    select(all_of(names_don))
  dat_rcpt <- dat_in %>%
    select(all_of(names_rcpt))

  # tmp[,c(1,2)] are donor's alleles
  # tmp[,c(3,4)] are recipient's alleles
  tmp <- as.data.frame(cbind(dat_don, dat_rcpt)) %>%
    mutate_all(as.character)
  #* end of data prep *#

  #* start of calculating mis-match of donor's hla in recipient *#
  len <- dim(dat_in)[1]
  mis_1 <- numeric(length = len)
  mis_2 <- numeric(length = len)

  for (i in 1:len)
  {
    # if both alleles are NA in donor or recipient, then set mis-match value to NA
    if ((is.na(tmp[i,1]) & is.na(tmp[i,2])) | (is.na(tmp[i,3]) & is.na(tmp[i,4]))) {
      mis_1[i] <- NA
      mis_2[i] <- NA
    }

    # if dornor's hla matches with any of recipient's, then mis-match is 0;  otherwise mis-match is 1
    else{
      if (tmp[i,1] %in% c(tmp[i,3],tmp[i,4]))
        mis_1[i] <- 0
      else mis_1[i] <- 1

      if (tmp[i,2] %in% c(tmp[i,3],tmp[i,4]))
        mis_2[i] <- 0
      else mis_2[i] <- 1
    }
  }

  tmp <- tmp %>%
    mutate(mis_1 = mis_1,
           mis_2 = mis_2) %>%
    filter(!is.na(mis_1) & !is.na(mis_2))
  #* end of calculating mis-match of donor's hla in recipient *#

  #* start of calculate mis-match frequency for donor *#
  dat_out <- tmp %>%
    select(names_don) %>%
    mutate_all(as.character) %>%
    gather(., name, allele) %>%
    group_by(name, allele) %>%
    summarise(freq = n()) %>%
    filter(!is.na(allele) & freq > 1) %>%
    ungroup() %>%
    arrange(name, allele, freq)

  return(dat_out)
}
