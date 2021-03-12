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

  #####################********* 03/11/2021 *********########################
  # set to NA for some unidentifiable values
  vec_na <- c("0.796527778", "0.175", "1/1/00 11:12")

  # define all possible combinations of brackets
  brkt <- combn(c("{", "[", "(", "<", "}", "]", ")", ">"), 2, paste)

  # step 1:
  # create a temp data frame to hold input antigens, and
  # - 1. if string in vec_na then set to NA
  # - 2. remove white spaces around the string
  # - 3. if string contains only letters then set to NA
  # - 4. remove everything between brackets including brackets, it covers BW4 condition
  # - 5. remove **XXX pattern start with* end with X -> \\*{0,}+X{1,} = start with * at least 0 times, end of X at least 1 time. eg: XXX, *X, **XX, *XXX
  # - 6. if there are more than 1 colons, remove everything starting from 2nd colon: 131:01:00 -> 131:01  26:08:00 -> 26:08
  # - 7. remove all punctuations except  ":"
  # - 8. remove if letter
  # - 9. add a leading 0 if 1 char 1:02 -> 01:02
  tmp <- data.frame(cbind(var_1, var_2)) %>%
    mutate(v1 = ifelse(var_1 %in% vec_na, "", var_1), # if string in vec_na then set to NA
           v2 = ifelse(var_2 %in% vec_na, "", var_2)) %>%
    mutate(v1.1 = str_trim(v1), # remove white spaces around the string
           v2.1 = str_trim(v2)) %>%
    mutate(v1.2 = ifelse(grepl("[^A-Za-z]+$", v1.1), v1.1, ""), # if string contains ALL letters then set to NA
           v2.2 = ifelse(grepl("[^A-Za-z]+$", v2.1), v2.1, "")) %>%
    mutate(v1.3 = rm_between(v1.2, brkt[1,], brkt[2,]), # remove everything between brackets including brackets
           v2.3 = rm_between(v2.2, brkt[1,], brkt[2,])) %>%
    mutate(v1.4 = ifelse(str_detect(v1.3, "\\*{0,}+X{1,}"), str_replace(v1.3, "\\*{0,}+X{1,}", ""), v1.3), # remove **XXX pattern star
           v2.4 = ifelse(str_detect(v2.3, "\\*{0,}+X{1,}"), str_replace(v2.3, "\\*{0,}+X{1,}", ""), v2.3)) %>%
    mutate(v1.5 = ifelse(str_count(v1.4, ":") > 1, sub("(:[^:]+):.*", "\\1", v1.4), v1.4), # if there are more than 1 ":", remove everything starting from 2nd
           v2.5 = ifelse(str_count(v2.4, ":") > 1, sub("(:[^:]+):.*", "\\1", v2.4), v2.4)) %>%
    mutate(v1.6 =  sub("([:])|[[:punct:]]", "\\1", v1.5), # remove all punctuations except ":"
           v2.6 =  sub("([:])|[[:punct:]]", "\\1", v2.5)) %>%
    mutate(v1.7 =  str_replace(v1.6, "[A-Za-z]", ""), # remove if letter
           v2.7 =  str_replace(v2.6, "[A-Za-z]", "")) %>%
    select(var_1, v1, v1.2, v1.2, v1.3, v1.4, v1.5, v1.6, v1.7,var_2, v2, v2.1, v2.2, v2.3, v2.4, v2.5, v2.6, v2.7)


  #####################******** end of 03/11/2021**********########################

  # create a temp data frame to hold input antigens
  tmp <- data.frame(cbind(var_1,var_2)) %>%
    mutate(var_1 = gsub(" ", "", str_trim(var_1, side = "both"), fixed = FALSE),
           var_2 = gsub(" ", "", str_trim(var_2, side = "both"), fixed = FALSE))

  # step 1: remove BW4s
  bw4 <- c("\\(BW4\\)|\\[BW4\\]|\\[BW4\\}|\\{BW4\\)|\\{BW4\\}|\\{BW4\\}|\\(\\{BW4\\}|\\(BW4\\}\\|BW4\\)")

  tmp$var_1 <- str_replace_all(tmp$var_1, bw4, " ")
  tmp$var_2 <- str_replace_all(tmp$var_2, bw4, " ")
  #* end of step 1 *#

  # step 2: fill out NA antigen with non-na antigen within a locus to assume homozygosity
  # logics: - if there is an NA antigen on either of the copies, replace it with the antigen from the non-na copy
  #        - if both copies are NA, then keep them as NA
  #        - trim white spaces in the string

  var_1 <- ifelse((is.na(tmp$var_1) | tmp$var_1 == "" | tmp$var_1 %in% vec_na) & (is.na(tmp$var_2) | tmp$var_2 == "" | tmp$var_2 %in% vec_na ), "",
                  ifelse(is.na(tmp$var_1) | tmp$var_1 == "" | tmp$var_1 %in% vec_na  & tmp$var_2 != "", tmp$var_2, tmp$var_1))

  var_2 <- ifelse((is.na(tmp$var_2) | tmp$var_2 == "" | tmp$var_2 %in% vec_na ) & (is.na(tmp$var_1) | tmp$var_1 == "" | tmp$var_1 %in% vec_na), "",
                  ifelse(is.na(tmp$var_2) | tmp$var_2 == "" | tmp$var_2 %in% vec_na & tmp$var_1 != "", tmp$var_1, tmp$var_2))
  #* end of step 2 *#

  # step 3: remove everything within a brackets. numbers within a brackets could be either expert or WHO typing, we keep numbers outside of brackets as low resolution typing to feed into haplotype reference table
  # example: keep 15 for 15{72}, keep 40 for 40[60]
  var_1_c1 <- ifelse(str_detect(var_1, "\\("), str_replace(var_1, "\\s*\\([^\\)]+\\)", ""),
                     ifelse(str_detect(var_1, "\\["), str_replace(var_1, "\\s*\\[[^\\]]+\\]", ""),
                            ifelse(str_detect(var_1, "\\{"), str_replace(var_1, "\\s*\\{[^\\}]+\\}", ""), var_1)))
  var_2_c1 <- ifelse(str_detect(var_2, "\\("), str_replace(var_2, "\\s*\\([^\\)]+\\)", ""),
                     ifelse(str_detect(var_2, "\\["), str_replace(var_2, "\\s*\\[[^\\]]+\\]", ""),
                            ifelse(str_detect(var_2, "\\{"), str_replace(var_2, "\\s*\\{[^\\}]+\\}", ""), var_2)))
  #* end of step 3 *#

  # step 4: remove other special symbols like *, #, :., ), ], }, (., and {. , then drop non-numeric characters
  var_1_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\(.*|\\{.*|@", '', var_1_c1)
  var_2_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\(.*|\\{.*|@", '', var_2_c1)

  # if the string contains more than 1 ":", then remove everything strating from 2nd ":"
  # ex: 39:02:02 will become 39:02
  var_1_c2 <- ifelse(str_count(var_1_c2, ":") <= 1, var_1_c2, sub("^(([^:]*:){1}[^:]*).*", "\\1", var_1_c2))
  var_2_c2 <- ifelse(str_count(var_2_c2, ":") <= 1, var_2_c2, sub("^(([^:]*:){1}[^:]*).*", "\\1", var_2_c2))
  #* end of step 4 *#

  # step 5: keep only first two characters
  var_1_c3 <- ifelse(str_count(var_1_c2, ":") == 1, var_1_c2, map_chr(var_1_c2, ~substr(.,1,2)))
  var_2_c3 <- ifelse(str_count(var_2_c2, ":") == 1, var_2_c2, map_chr(var_2_c2, ~substr(.,1,2)))
  #* end of step 5 *#

  # step 6: if it's low resolution and nchar = 1, then add a leading 0
  # rationale : some of antigen are the same but in different notation
  # example : 2 and 02 are the same
  var_1_c4 <- map_chr(var_1_c3, function(x) ifelse(nchar(x) == 1, paste0("0",x), x))
  var_2_c4 <- map_chr(var_2_c3, function(x) ifelse(nchar(x) == 1, paste0("0",x), x))
  #* end of step 6 *#

  # step 7: - if string contains letters only, then set to NA
  var_1_c5 <- ifelse(!grepl("[^A-Za-z]", var_1_c4), "", var_1_c4)
  var_2_c5 <- ifelse(!grepl("[^A-Za-z]", var_2_c4), "", var_2_c4)
  #* end of step 7 *#

  # step 8: - if it's high resolution and nchar < 5, then add a leading 0
  #         - if it's hi resolution and nchar > 5 then keep first 5 chars only
  tmp1 <- ifelse(str_count(var_1_c5, ":") > 0 & nchar(var_1_c5) < 5, paste("0", var_1_c5, sep=""), var_1_c5)
  tmp2 <- ifelse(str_count(var_2_c5, ":") > 0 & nchar(var_2_c5) < 5, paste("0", var_2_c5, sep=""), var_2_c5)

  tmp1 <- str_replace(tmp1, "[A-Z]|[a-z]", "")
  tmp2 <- str_replace(tmp2, "[A-Z]|[a-z]", "")
  #* end of step 8 *#

  #* step 9: final table *#
  result <- data.frame(cbind(tmp1, tmp2)) %>%
    mutate(locus1_clean = ifelse(tmp1 != "", tmp1,
                                 ifelse(tmp2 != "", tmp2, "")),
           locus2_clean = ifelse(tmp2 != "", tmp2,
                                 ifelse(tmp1 != "", tmp1, ""))) %>%
    select(locus1_clean, locus2_clean)

  rm(var_1_c1, var_1_c2, var_1_c3, var_1_c4, var_1_c5, var_2_c1, var_2_c2, var_2_c3, var_2_c4, var_2_c5, tmp, tmp1, tmp2)
  #* end of step 9 *#

  return(result)
}
