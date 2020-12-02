#' @title calculate frequency of donor mis-match alleles to recipient's
#' @param dat_in
#' dataframe of clean HLA
#' @param nms_don
#' column names of donor's alleles, must be length of 2
#' @param nms_rcpt
#' column names of recipient's alleles, must be length of 2
#' @return
#' dataframe of donor's mis-match alleles with frequency > 1
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
# don <- c("donor.a1", "donor.a2")
# rcpt <- c("recipient.a1", "recipient.a2")
#' result <- CalAlleleMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)
#' }
#' @export
#'
CalAlleleMismFreq <- function(dat_in, nms_don = c(), nms_rcpt = c()){
  #* start of data prep *#
  dat_don <- dat_in %>%
    select(all_of(nms_don))
  dat_rcpt <- dat_in %>%
    select(all_of(nms_rcpt))

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
    select(all_of(nms_don)) %>%
    mutate_all(as.character) %>%
    gather(., name, allele) %>%
    filter(!is.na(allele)) %>%
    group_by(allele) %>%
    summarise(freq = n(), .groups = 'drop') %>%
    filter(!is.na(allele) & freq > 1) %>%
    ungroup() %>%
    arrange(-freq)

  return(dat_out)
}

