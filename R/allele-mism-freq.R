#' @name CalAlleleMismFreq
#' @title Evaluate the frequency of specific allele mismatches
#' @description This function evaluates allele level mismatch between donor and recipient and then presents the most commonly mismatched alleles. This function is most effectively used to study the most common mismatches within a transplant population.
#' @param dat_in
#' A data frame of clean HLA typing data.
#' @param nms_don
#' A vector of column names of donor's alleles, must be length of 2.
#' @param nms_rcpt
#' A vector of column names of recipient's alleles, must be length of 2.
#' @return
#' A data frame of donor's mismatched alleles with frequency > 1. No mismatch is calculated if input alleles are NA.
#' @import
#' tidyverse
#' @importFrom
#' stats setNames
#' @examples
#' dat <- read.csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))
#' don <- c("donor.a1", "donor.a2")
#' rcpt <- c("recipient.a1", "recipient.a2")
#' re <- CalAlleleMismFreq(dat_in = dat, nms_don = don, nms_rcpt = rcpt)
#' @export

CalAlleleMismFreq <- function(dat_in, nms_don = c(), nms_rcpt = c()){
  #* step 1: data prep *#
  dat_don <- dat_in %>%
              select(all_of(nms_don))

  dat_rcpt <- dat_in %>%
                select(all_of(nms_rcpt))

  tmp <- as.data.frame(cbind(dat_don, dat_rcpt)) %>%
          mutate_all(as.character)
  #* end of step 1*#

  #* step 2: calculate mis-match of donor's hla in recipient *#
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

    # if donor's hla matches with any of recipient's, then mis-match is 0; otherwise mis-match is 1
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
  #* end of step 2 *#

  #* step 3: calculate mis-match frequency for donor *#
  dat_out <- tmp %>%
              select(-all_of(nms_rcpt)) %>%
              setNames(c("d1", "d2", "m1", "m2")) %>%
              mutate(d1 = ifelse(m1 == 1, d1, 0 ),
                     d2 = ifelse(m2 == 1, d2, 0)) %>%
              select(d1, d2) %>%
              gather(., name, allele) %>%
              filter(allele > 0) %>%
              group_by(allele) %>%
              summarise(freq = n(), .groups = 'drop') %>%
              filter(!is.na(allele)) %>%
              ungroup() %>%
              arrange(-freq)
  #* end of step 3 *#

  return(dat_out)
}

