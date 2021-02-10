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
  # vector of NA definition based on the unique XX pattern from all of the locus in the existing data
  # set to NA for some unidentifiable values
  vec_na <- c("0.796527778", "0.175", "1/1/00 11:12")

  # create a temp data frame to hold input antigens
  tmp <- data.frame(cbind(var_1,var_2))

  # step 1 : remove BW4s
  bw4 <- c("\\(BW4\\)|\\[BW4\\]|\\[BW4\\}|\\{BW4\\)|\\{BW4\\}|\\{BW4\\}|\\(\\{BW4\\}|\\(BW4\\}\\|BW4\\)")

  tmp$var_1 <- str_replace_all(tmp$var_1, bw4, " ")
  tmp$var_2 <- str_replace_all(tmp$var_2, bw4, " ")

  # step 2 : fill out NA antigen with non-na antigen within a locus to assume homozygosity
  # logics : - if there is an NA antigen on either of the copies, replace it with the antigen from the non-na copy
  #          - if both copies are NA, then keep them as NA
  #          - trim white spaces in the string
  #          - if string contains letters only, then set to NA
  var_1 <- ifelse((is.na(tmp$var_1) | tmp$var_1 %in% vec_na | str_detect(tmp$var_1, '[A-Za-z]')) & (is.na(tmp$var_2) | tmp$var_2 %in% vec_na | str_detect(tmp$var_2, '[A-Za-z]')), "",
                  ifelse(is.na(tmp$var_1) | tmp$var_1 %in% vec_na | str_detect(tmp$var_1, '[A-Za-z]') & tmp$var_2 != "", tmp$var_2, tmp$var_1)) %>%
                  str_trim(side = "both") %>%
                  gsub(" ", "", ., fixed = FALSE)

  var_2 <- ifelse((is.na(tmp$var_2) | tmp$var_2 %in% vec_na | str_detect(tmp$var_2, '[A-Za-z]')) & (is.na(tmp$var_1) | tmp$var_1 %in% vec_na | str_detect(tmp$var_1, '[A-Za-z]')), "",
                  ifelse(is.na(tmp$var_2) | tmp$var_2 %in% vec_na | str_detect(tmp$var_2, '[A-Za-z]') & tmp$var_1 != "", tmp$var_1, tmp$var_2)) %>%
                  str_trim(side = "both") %>%
                  gsub(" ", "", ., fixed = FALSE)

  var_1 <- ifelse(!grepl("[^A-Za-z]", var_1), "", var_1)
  var_2 <- ifelse(!grepl("[^A-Za-z]", var_2), "", var_2)

  # step 3 : remove everything within a brackets. numbers within a brackets could be either expert or WHO typing, we keep numbers outside of brackets as low resolution typing to feed into haplotype reference table
  # example : keep 15 for 15{72}, keep 40 for 40[60]
  var_1_c1 <- ifelse(str_detect(var_1, "\\("), str_replace(var_1, "\\s*\\([^\\)]+\\)", ""),
                     ifelse(str_detect(var_1, "\\["), str_replace(var_1, "\\s*\\[[^\\]]+\\]", ""),
                            ifelse(str_detect(var_1, "\\{"), str_replace(var_1, "\\s*\\{[^\\}]+\\}", ""), var_1)))
  var_2_c1 <- ifelse(str_detect(var_2, "\\("), str_replace(var_2, "\\s*\\([^\\)]+\\)", ""),
                     ifelse(str_detect(var_2, "\\["), str_replace(var_2, "\\s*\\[[^\\]]+\\]", ""),
                            ifelse(str_detect(var_2, "\\{"), str_replace(var_2, "\\s*\\{[^\\}]+\\}", ""), var_2)))

  # step 4 : remove other special symbols like *, #, :., ), ], }, (., and {. , then drop non-numeric characters
  #var_1_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\:.*|\\(.*|\\{.*", '', var_1_c1)
  #var_2_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\:.*|\\(.*|\\{.*", '', var_2_c1)

  var_1_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\(.*|\\{.*", '', var_1_c1)
  var_2_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\(.*|\\{.*", '', var_2_c1)

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

  # step 8 : - if it's high resolution and nchar < 5, then add a leading 0
  #         - if it's hi resolution and nchar > 5 then keep first 5 chars only
  # var_1_out <- ifelse(str_count(var_1_c4, ":") > 0 & nchar(unlist(str_split(var_1_c4,":"))) == 1, paste("0", var_1_c4, sep=""), var_1_c4)
  # var_2_out <- ifelse(str_count(var_2_c4, ":") > 0 & nchar(unlist(str_split(var_2_c4,":"))) == 1, paste("0", var_2_c4, sep=""), var_2_c4)
  var_1_out <- ifelse(str_count(var_1_c4, ":") > 0 & nchar(var_1_c4) < 5, paste("0", var_1_c4, sep=""), var_1_c4)
  var_2_out <- ifelse(str_count(var_2_c4, ":") > 0 & nchar(var_2_c4) < 5, paste("0", var_2_c4, sep=""), var_2_c4)

  var_1_out <- ifelse(nchar(var_1_out) > 5, substr(var_1_out, 1, 5), var_1_out)
  var_2_out <- ifelse(nchar(var_2_out) > 5, substr(var_2_out, 1, 5), var_2_out)

  # step 9: remove temporary variable holders
  rm(var_1_c1, var_1_c2, var_1_c3, var_1_c4, var_2_c1, var_2_c2, var_2_c3, var_2_c4, tmp)

  result <- data.frame(cbind(var_1_out,var_2_out)) %>%
    setNames(c("locus1_clean", "locus2_clean"))
  return(result)
}
