#' @name EvalAlleleMism
#' @title evaluate mis-match alleles
#' @param dat_in
#' dataframe with allele info
#' @param don_1
#' donor's alpha1 domain
#' @param don_2
#' donor's alpha2 or beta1 domain
#' @param recip_1
#' recipient's alpha1 domain
#' @param recip_2
#' recipient's alpha2 or beta1 domain
#' @param locus
#' locus of domains
#' @return
#' a list of mis-match vectors
#' mism_1 is mis-match flag of locus 1
#' mism_2 is mis-match flag of locus 2
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata/example", "HLA_MisMatch_test.csv", package = "hlaR"))
#' a <- EvalAlleleMism(dat, dat$donor.a1, dat$donor.a2, dat$recipient.a1, dat$recipient.a2, "a")
# dat$mism.a1 <- a$mism_1
# dat$mism.a2 <- a$mism_2
#' }

EvalAlleleMism <- function(dat_in, don_1, don_2, recip_1, recip_2, locus)
{
  # start of calling CleanAllele() to clean hla value #
  don_1_c <- unlist(CleanAllele(don_1, don_2, locus)[1])
  don_2_c <- unlist(CleanAllele(don_1, don_2, locus)[2])
  recip_1_c <- unlist(CleanAllele(recip_1, recip_2, locus)[1])
  recip_2_c <- unlist(CleanAllele(recip_1, recip_2, locus)[2])

  tmp <- as.data.frame(cbind(don_1_c, don_2_c, recip_1_c, recip_2_c)) %>%
    mutate_all(as.character)
  # end of calling CleanAllele() #

  # start of calculating mis-match of donor's hla in recipient #
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
  # end of calculating mis-match #

  rm(tmp)

  return(list(mism_1 = mis_1, mism_2 = mis_2))
}
