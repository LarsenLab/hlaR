#' @name EvalAlleleMism
#' @title evaluate mis-match alleles
#' @param don_1
#' donor's alpha1 domain
#' @param don_2
#' donor's alpha2 or beta1 domain
#' @param recip_1
#' recipient's alpha1 domain
#' @param recip_2
#' recipient's alpha2 or beta1 domain
#' @return
#' a data frame with cleaned donor and recipients, and donor to recipient mis-match flags
#' @export
#'
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# input hlas could be either raw or cleaned
# dat <- read_csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR"))
#' a <- EvalAlleleMism99(dat$DONOR_A1, dat$DONOR_A2, dat$RECIPIENT_A1, dat$RECIPIENT_A2)
# a
#' }

EvalAlleleMism <- function(don_1, don_2, recip_1, recip_2)
{
  # start of calling CleanAllele() to clean hla value #
  don_1_clean <- unlist(CleanAllele(don_1, don_2)[1])
  don_2_clean <- unlist(CleanAllele(don_1, don_2)[2])
  recip_1_clean <- unlist(CleanAllele(recip_1, recip_2)[1])
  recip_2_clean <- unlist(CleanAllele(recip_1, recip_2)[2])

  tmp <- as.data.frame(cbind(don_1_clean, don_2_clean, recip_1_clean, recip_2_clean)) %>%
    mutate_all(as.character)
  # end of calling CleanAllele() #

  # start of calculating mis-match of donor's hla in recipient #
  len <- length(don_1_clean)
  mis_1 <- numeric(length = len)
  mis_2 <- numeric(length = len)

  for (i in 1:len)
  {
    # if both alleles are NA in donor or recipient, then set mis-match value to NA
    if((is.na(tmp[i,1]) & is.na(tmp[i,2])) | (is.na(tmp[i,3]) & is.na(tmp[i,4]))){
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
  # end of calculating mis-match #

  rm(tmp)

  result <- data.frame(don_1_clean = don_1_clean, don_2_clean = don_2_clean,
                       recip_1_clean = recip_1_clean, recip_2_clean = recip_2_clean,
                       mis_1 = mis_1, mis_2 = mis_2)
  rownames(result) <- seq(1, dim(result)[1])

  return(result)
}
