
#' function to evaluate mis-match HLAs
#' @param dat_in
#' a dataframe with HLAs
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
#' list of mis-match vectors
#' mism_1 is mis-match flag of alpha1
#' mism_2 is mis-match flag of alpha2 or beta
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# hla <- read_csv("~/projects/hlaclean_test/hla_typing_test.csv") %>%
# rename_all(. %>% tolower %>% gsub("[[:blank:]]|[[:punct:]]", ".", .))
#' a <- EvalMism(hla, hla$donor.a1, hla$donor.a2, hla$recipient.a1, hla$recipient.a2, "a")
# hla$mism.a1 <- a$mism_1
# hla$mism.a2 <- a$mism_2
# hla$mism.b1 <- unlist(MismatchHla(hla, hla$donor.b1, hla$donor.b2,
# hla$recipient.b1, hla$recipient.b2, "b")[1])
# hla$mism.b2 <- unlist(MismatchHla(hla, hla$donor.b1, hla$donor.b2,
# hla$recipient.b1, hla$recipient.b2, "b")[2])
#' }

EvalMism <- function(dat_in, don_1, don_2, recip_1, recip_2, locus)
{
  # start of calling CleanHla() to clean hla value #
  don_1_c <- unlist(CleanHla(don_1, don_2, locus)[1])
  don_2_c <- unlist(CleanHla(don_1, don_2, locus)[2])
  recip_1_c <- unlist(CleanHla(recip_1, recip_2, locus)[1])
  recip_2_c <- unlist(CleanHla(recip_1, recip_2, locus)[2])

  tmp <- as.data.frame(cbind(don_1_c, don_2_c, recip_1_c, recip_2_c)) %>%
    mutate_all(as.character)
  # end of calling CleanHla() #

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
