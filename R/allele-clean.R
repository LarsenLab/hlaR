#' @name CleanAllele
#' @title clean messy HLA typing for downstream analysis
#' @param var_1
#' hla on locus 1
#' @param var_2
#' hla on locus 2
#' @return
#' a data frame with cleaned hla of each locus
#' @export
#'
#' @import
#' tidyverse
#' utils
#' readr
#'
#' @examples
#' \dontrun{
# dat <-  read_csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR")) %>% rename_all(. %>% tolower %>% gsub("[[:blank:]]|[[:punct:]]", ".", .))
#' clean1 <- CleanAllele(clean$RECIPIENT_A1, clean$RECIPIENT_A2)
#' clean2 <- CleanAllele(clean$DONOR_DRB11, clean$DONOR_DRB11)
#' }
#'

#' @name CleanAllele
#' @title clean messy HLA typing for downstream analysis
#' @param var_1
#' hla on locus 1
#' @param var_2
#' hla on locus 2
#' @return
#' a data frame with cleaned hla of each locus
#' @export
#'
#' @import
#' tidyverse
#' utils
#' readr
#'
#' @examples
#' \dontrun{
# dat <-  read_csv(system.file("extdata/example", "HLA_Clean_test.csv", package = "hlaR")) %>% rename_all(. %>% tolower %>% gsub("[[:blank:]]|[[:punct:]]", ".", .))
#' clean1 <- CleanAllele(clean$RECIPIENT_A1, clean$RECIPIENT_A2)
#' clean2 <- CleanAllele(clean$DONOR_DRB11, clean$DONOR_DRB11)
#' }
#'

CleanAllele <- function(var_1, var_2) {
  # set to NA for some unidentifiable values
  vec_na <- c("0.796527778", "0.175", "1/1/00 11:12")

  # create a temp data frame to hold input antigens
  tmp <- data.frame(cbind(var_1,var_2)) %>%
    mutate(var_1 = gsub(" ", "", str_trim(var_1, side = "both"), fixed = FALSE),
           var_2 = gsub(" ", "", str_trim(var_2, side = "both"), fixed = FALSE))

  # step 1 : remove BW4s
  bw4 <- c("\\(BW4\\)|\\[BW4\\]|\\[BW4\\}|\\{BW4\\)|\\{BW4\\}|\\{BW4\\}|\\(\\{BW4\\}|\\(BW4\\}\\|BW4\\)")

  tmp$var_1 <- str_replace_all(tmp$var_1, bw4, " ")
  tmp$var_2 <- str_replace_all(tmp$var_2, bw4, " ")

  # step 2 : fill out NA antigen with non-na antigen within a locus to assume homozygosity
  # logics : - if there is an NA antigen on either of the copies, replace it with the antigen from the non-na copy
  #          - if both copies are NA, then keep them as NA
  #          - trim white spaces in the string
  #
  # var_1 <- ifelse((is.na(tmp$var_1) | tmp$var_1 %in% vec_na | str_detect(tmp$var_1, '[A-Za-z]')) & (is.na(tmp$var_2) | tmp$var_2 %in% vec_na | str_detect(tmp$var_2, '[A-Za-z]')), "",
  #                 ifelse(is.na(tmp$var_1) | tmp$var_1 %in% vec_na | str_detect(tmp$var_1, '[A-Za-z]') & tmp$var_2 != "", tmp$var_2, tmp$var_1))
  #
  # var_2 <- ifelse((is.na(tmp$var_2) | tmp$var_2 %in% vec_na | str_detect(tmp$var_2, '[A-Za-z]')) & (is.na(tmp$var_1) | tmp$var_1 %in% vec_na | str_detect(tmp$var_1, '[A-Za-z]')), "",
  #                 ifelse(is.na(tmp$var_2) | tmp$var_2 %in% vec_na | str_detect(tmp$var_2, '[A-Za-z]') & tmp$var_1 != "", tmp$var_1, tmp$var_2))

  var_1 <- ifelse((is.na(tmp$var_1) | tmp$var_1 == "" | tmp$var_1 %in% vec_na) & (is.na(tmp$var_2) | tmp$var_2 == "" | tmp$var_2 %in% vec_na ), "",
                  ifelse(is.na(tmp$var_1) | tmp$var_1 == "" | tmp$var_1 %in% vec_na  & tmp$var_2 != "", tmp$var_2, tmp$var_1))

  var_2 <- ifelse((is.na(tmp$var_2) | tmp$var_2 == "" | tmp$var_2 %in% vec_na ) & (is.na(tmp$var_1) | tmp$var_1 == "" | tmp$var_1 %in% vec_na), "",
                  ifelse(is.na(tmp$var_2) | tmp$var_2 == "" | tmp$var_2 %in% vec_na & tmp$var_1 != "", tmp$var_1, tmp$var_2))

  # step 3 : remove everything within a brackets. numbers within a brackets could be either expert or WHO typing, we keep numbers outside of brackets as low resolution typing to feed into haplotype reference table
  # example : keep 15 for 15{72}, keep 40 for 40[60]
  var_1_c1 <- ifelse(str_detect(var_1, "\\("), str_replace(var_1, "\\s*\\([^\\)]+\\)", ""),
                     ifelse(str_detect(var_1, "\\["), str_replace(var_1, "\\s*\\[[^\\]]+\\]", ""),
                            ifelse(str_detect(var_1, "\\{"), str_replace(var_1, "\\s*\\{[^\\}]+\\}", ""), var_1)))
  var_2_c1 <- ifelse(str_detect(var_2, "\\("), str_replace(var_2, "\\s*\\([^\\)]+\\)", ""),
                     ifelse(str_detect(var_2, "\\["), str_replace(var_2, "\\s*\\[[^\\]]+\\]", ""),
                            ifelse(str_detect(var_2, "\\{"), str_replace(var_2, "\\s*\\{[^\\}]+\\}", ""), var_2)))

  # step 4 : remove other special symbols like *, #, :., ), ], }, (., and {. , then drop non-numeric characters
  var_1_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\(.*|\\{.*|@", '', var_1_c1)
  var_2_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\(.*|\\{.*|@", '', var_2_c1)

  # if the string contains more than 1 ":", then remove everything strating from 2nd ":"
  # ex: 39:02:02 will become 39:02
  var_1_c2 <- ifelse(str_count(var_1_c2, ":") <= 1, var_1_c2, sub("^(([^:]*:){1}[^:]*).*", "\\1", var_1_c2))
  var_2_c2 <- ifelse(str_count(var_2_c2, ":") <= 1, var_2_c2, sub("^(([^:]*:){1}[^:]*).*", "\\1", var_2_c2))

  # step 5 : keep only first two characters
  var_1_c3 <- ifelse(str_count(var_1_c2, ":") == 1, var_1_c2, map_chr(var_1_c2, ~substr(.,1,2)))
  var_2_c3 <- ifelse(str_count(var_2_c2, ":") == 1, var_2_c2, map_chr(var_2_c2, ~substr(.,1,2)))

  # step 7 : if it's low resolution and nchar = 1, then add a leading 0
  # rationale : some of antigen are the same but in different notation
  # example : 2 and 02 are the same
  var_1_c4 <- map_chr(var_1_c3, function(x) ifelse(nchar(x) == 1, paste0("0",x), x))
  var_2_c4 <- map_chr(var_2_c3, function(x) ifelse(nchar(x) == 1, paste0("0",x), x))

  # step 8 : - if string contains letters only, then set to NA
  var_1_c5 <- ifelse(!grepl("[^A-Za-z]", var_1_c4), "", var_1_c4)
  var_2_c5 <- ifelse(!grepl("[^A-Za-z]", var_2_c4), "", var_2_c4)

  # step 9 : - if it's high resolution and nchar < 5, then add a leading 0
  #         - if it's hi resolution and nchar > 5 then keep first 5 chars only
  # var_1_out <- ifelse(str_count(var_1_c4, ":") > 0 & nchar(unlist(str_split(var_1_c4,":"))) == 1, paste("0", var_1_c4, sep=""), var_1_c4)
  # var_2_out <- ifelse(str_count(var_2_c4, ":") > 0 & nchar(unlist(str_split(var_2_c4,":"))) == 1, paste("0", var_2_c4, sep=""), var_2_c4)
  tmp1 <- ifelse(str_count(var_1_c5, ":") > 0 & nchar(var_1_c5) < 5, paste("0", var_1_c5, sep=""), var_1_c5)
  tmp2 <- ifelse(str_count(var_2_c5, ":") > 0 & nchar(var_2_c5) < 5, paste("0", var_2_c5, sep=""), var_2_c5)

  tmp1 <- ifelse(nchar(tmp1) > 5, substr(tmp1, 1, 5), tmp1)
  tmp2 <- ifelse(nchar(tmp2) > 5, substr(tmp2, 1, 5), tmp2)

  tmp1 <- ifelse(str_detect(tmp1, "\\d\\d\\:\\w\\w"), str_sub(tmp1, 1 ,2),
                 ifelse(str_detect(tmp1, "\\ww\\ww\\:\\d\\d"), str_sub(tmp1, 4 ,5), tmp1))

  tmp2 <- ifelse(str_detect(tmp2, "\\d\\d\\:\\w\\w"), str_sub(tmp2, 1 ,2),
                 ifelse(str_detect(tmp2, "\\ww\\ww\\:\\d\\d"), str_sub(tmp2, 4 ,5), tmp2))


  result <- data.frame(cbind(tmp1, tmp2)) %>%
    mutate(locus1_clean = ifelse(tmp1 != "", tmp1,
                                 ifelse(tmp2 != "", tmp2, "")),
           locus2_clean = ifelse(tmp2 != "", tmp2,
                                 ifelse(tmp1 != "", tmp1, ""))) %>%
    select(locus1_clean, locus2_clean)

  # step 10: remove temporary variable holders
  rm(var_1_c1, var_1_c2, var_1_c3, var_1_c4, var_1_c5, var_2_c1, var_2_c2, var_2_c3, var_2_c4, var_2_c5, tmp, tmp1, tmp2)

  return(result)
}

