#' @name EvalAlleleMism
#' @title Evaluate mismatched alleles
#' @description Compare donor and recipient HLA(Human Leukocyte Antigen) typing data to determine mismatched alleles. Input data can be high or low resolution, mismatch is evaluated at the allele level.
#' @param don_1
#' Donor's alpha1 domain.
#' @param don_2
#' Donor's alpha2 or beta1 domain.
#' @param recip_1
#' Recipient's alpha1 domain.
#' @param recip_2
#' Recipient's alpha2 or beta1 domain.
#' @param hmz_cnt
#' Use hmz_cnt to determine how mismatch at homozygous alleles should be handled. By default, a mismatch at a homozygous allele is considered a single mismatch. Set hmz_cnt = 2 to count homozygous mismatches as double.
#' @return
#' A data frame of original input columns followed by mism_cnt of each donor/recipient pair.
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' dat <- read.csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))
#' re <- EvalAlleleMism(dat$donor_a1, dat$donor_a2, dat$recipient_a1, dat$recipient_a2, hmz_cnt = 2)

EvalAlleleMism <- function(don_1, don_2, recip_1, recip_2, hmz_cnt = 1)
{
  # step 1: call CleanAllele() to clean hla value #
  don_1_clean <- CleanAllele(don_1, don_2)$locus1_clean
  don_2_clean <- CleanAllele(don_1, don_2)$locus2_clean
  recip_1_clean <- CleanAllele(recip_1, recip_2)$locus1_clean
  recip_2_clean <- CleanAllele(recip_1, recip_2)$locus2_clean

  tmp <- as.data.frame(cbind(don_1_clean, don_2_clean, recip_1_clean, recip_2_clean)) %>%
    mutate_all(as.character)
  # end of step 1 #

  # step 2: calculate mis-match of donor's hla in recipient #
  len <- length(don_1_clean)
  mis_1 <- numeric(length = len)
  mis_2 <- numeric(length = len)

  for (i in 1:len)
  {
    # if both alleles are NA in donor or recipient, then set mis-match value to NA
    # if((is.na(tmp[i,1]) & is.na(tmp[i,2])) | (is.na(tmp[i,3]) & is.na(tmp[i,4]))){
    if((tmp[i,1] == "" & tmp[i,2] == "") | (tmp[i,3] == "" & tmp[i,4] == "")){
      mis_1[i] <- NA
      mis_2[i] <- NA
    }
    # if all of allele are high resolution
    else if(length(tmp[i,1]) == 5 & length(tmp[i,2]) == 5 & length(tmp[i,3]) == 5 & length(tmp[i,4]) == 5){
      if(tmp[i,1] %in% c(tmp[i,3],tmp[i,4])){
        mis_1[i] <- 0}
      else{
        mis_1[i] <- 1}

      if(tmp[i,2] %in% c(tmp[i,3],tmp[i,4])){
        mis_2[i] <- 0}
      else{
        mis_2[i] <- 1}
    }
    else{
      tmp[i,1] <- str_sub(tmp[i,1],1,2)
      tmp[i,2] <- str_sub(tmp[i,2],1,2)
      tmp[i,3] <- str_sub(tmp[i,3],1,2)
      tmp[i,4] <- str_sub(tmp[i,4],1,2)

      if(tmp[i,1] %in% c(tmp[i,3],tmp[i,4])){
        mis_1[i] <- 0}
      else{
        mis_1[i] <- 1}

      if(tmp[i,2] %in% c(tmp[i,3],tmp[i,4])){
        mis_2[i] <- 0}
      else{
        mis_2[i] <- 1}
    }
  }
  rm(tmp)
  # end of step 3 #

  #* step 4: construct final table *#
  # mismatch total count = 1 if mismatch at homozygous alleles
  result <- data.frame(don_1_clean = don_1_clean, don_2_clean = don_2_clean,
                       recip_1_clean = recip_1_clean, recip_2_clean = recip_2_clean,
                       mism_1 = mis_1, mism_2 = mis_2) %>%
            mutate(cnt = mism_1 + mism_2,
                   mism_cnt = ifelse(hmz_cnt == 1 & don_1_clean == don_2_clean & cnt == 2, cnt - 1, cnt)) %>%
            select(-c(mism_1, mism_2, cnt))

  rownames(result) <- seq(1, dim(result)[1])
 #* end of step 4 *#

  return(result)
}
