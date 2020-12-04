#' clean messy HLA typing for downstream analysis
#' @param var_1
#' alpha1 domain we want to clean
#' @param var_2
#' alpha2 or beta1 domain we want to clean
#' @param locus
#' locus of which we want to work
#' var_1 is hla of alpha1(or alpha for class II) allele
#' var_2 is hla of alpha2(or beta for class II) allele
#' locus is indicator of which locus we are working on
#' @return
#' list of data cleaned hla of the locus
#' returns alpha1 and alpha2 for HLA class I
#' returns alpha and beta for HLA class II
#' @export
#'
#' @import
#' tidyverse
#' utils
#' readr
#'
#' @examples
#' \dontrun{
# dat <- read_csv("~/projects/hlaclean_test/hla_typing_test.csv") %>% rename_all(. %>% tolower %>% gsub("[[:blank:]]|[[:punct:]]", ".", .))
#' re <- CleanAllele(dat$don.a1.ori, dat$don.a2.ori, locus = "a")
#' re <- CleanAllele(dat$don.drb.ori, dat$don.drb.ori, locus = "drb")
#' }
#'

CleanAllele <- function(var_1, var_2, locus) {
  LookupA <- read.table(text = "
      antigen   type
      1:01:01	  1
      1:02	    1
      2:01	    2
      2:01:01	  2
      2:02	    2
      2:02:01	  2
      2:03	    2
      2:03:01	  2
      2:04	    2
      2:05	    2
      3:02	    3
      11:01:01	11
      23:01	    23
      23:01:03	23
      24:02:00	24
      24:07:00	24
      24:07:01	24
      24:10:00	24
      26:08:00	26
      26:12:00	26
      29:02:00	29
      29:02:01	29
      30:02:00	30
      31:01:00	31
      32:01:00	32
      32:04:00	32
      33:01:00	33
      33:01:01	33
      34:01:00	34
      34:01:01	34
      34:02:00	34
      34:02:01	34
      34:05:00	34
      36:01:00	36
      66:01:00	66
      66:02:00	66
      66:03:00	68
      68:01:00	68
      68:02:00	68
      68:02:01	68
      68:03:00	68
      68:03:01	68
      68:28:00	68
      69:01:00	69
      74:03:00	74
      74:10:00	74
      80:01:00	80
      #*32(BW4}  32
      *36:01	  36
      *68:01	  68", header = TRUE)

  LookupB <- read.table(text = "
      antigen   type
      7:02	    7
      7:12	    7
      7:42	    7
      8:01:01	  8
      8:03	    8
      8:20	    8
      14:01	    64
      14:02:01	65
      14:02:02	65
      14:03	    14
      15:02	    75
      15:10	    71
      15:10:01	71
      15:16	    63
      15:17:01	63
      18:03	    18
      27:03:00	27
      27:06:00	27
      27:07:00	35
      35:02:00	35
      35:02:01	35
      35:05:00	35
      35:08:00	35
      35:09:00	35
      35:12:00	35
      35:12:01	35
      35:16:00	35
      35:24:00	35
      37:01:00	37
      37:01:01	37
      38:01:00	38
      38:01:01	38
      39:02:02	39
      39:05:00	39
      39:06:02	39
      39:08:00	39
      39:09:00	39
      39:10:00	39
      39:31:00	39
      40:05:00	50
      40:06:00	61
      41:01:00	41
      41:02:00	41
      41:02:01	41
      41:03:00	41
      42:01:00	42
      42:01:01	42
      42:02:00	42
      44:02:01	44
      44:03:00	44
      44:03:02	44
      44:04:00	44
      44:05:00	44
      44:05:01	44
      44:07:00	44
      44:10:00	44
      47:01:00	47
      47:01:01	47
      47:03:00	47
      49:01:00	49
      49:01:01	49
      50:01:00	50
      50:01:01	50
      51:01:00	51
      51:01:01	51
      51:02:01	51
      51:05:00	51
      51:08:00	51
      51:13:02	51
      51:22:00	51
      52:01:00	52
      52:01:02	52
      53:01:00	53
      53:01:01	53
      55:01:00	55
      57:01:00	57
      57:01:01	57
      57:02:00	57
      57:03:00	57
      57:03:01	57
      57:04:00	57
      58:02:00	58
      73:01:00	73
      78:01:00	78
      78:01:01	78
      82:01:00	82
      82:02:01	82
      *42:01	  42
      *58:02	  58", header = TRUE)

  LookupDQB <- read.table(text = "
      antigen   type
      2:01      2
      2:01:01   2
      2:02      2
      2:03      2
      3:02      8
      3:02:01   8
      3:02:02   8
      3:03:02   9
      3:04      7
      3:05:01   8
      3:09      7
      3:19      7
      4:01:01   4
      4:02      4
      4:02:01   4
      5:01      5
      5:01:01   5
      5:02      5
      5:02:01   5
      5:03:01   5
      5:04      5
      6:01      6
      6:02      6
      6:02:01   6
      6:03      6
      6:03:01   6
      6:04:01   6
      6:08      6
      6:09      6
      6:09:01   6
      6:10      1
      6:11:01   1
      201 2
      202 2
      203 2
      301 7
      302 8
      303 9
      304 7
      316 7
      401 4
      402 4
      501  5
      502 5
      503 5
      50301 5
      504 5
      601 6
      602 6
      603 6
      604 6
      609 6", header = TRUE)

  LookupDRB1 <- read.table(text = "
      antigen   type
      01(103)   103
      1:01      1
      1:01:01   1
      1:02      1
      1:02:01   1
      1:03      1
      3:01      17
      3:01:01   17
      3:02      18
      3:02:01   18
      3:07      3
      4:01      4
      4:01:01   4
      4:02      4
      4:02:01   4
      4:03      4
      4:03:01   4
      4:04      4
      4:04:01   4
      4:05      4
      4:05:01   4
      4:06      4
      4:07      4
      4:07:01   4
      4:08      4
      4:08:01   4
      4:10      4
      4:11      4
      7:01      7
      7:01:01   7
      8:01      8
      8:02      8
      8:02:01   8
      8:03      8
      8:03:02   8
      8:04      8
      8:04:01   8
      8:06      8
      9:01:02   9
      10:01     10
      10:01:01  10
      11:01     11
      11:01:02  11
      11:02     11
      11:02:01  11
      11:03     11
      11:04     11
      11:04:01  11
      11:12     11
      12:02     12
      12:02:01  12
      13:01     13
      13:01:01  13
      13:02     13
      13:02:01  13
      13:03     13
      13:03:01  13
      13:04     13
      13:05:01  13
      14:02     14
      14:04     14
      14:06:01  14
      14:21     14
      14:24     14
      15:01     15
      15:01:01  15
      15:02     15
      15:02:01  15
      15:02:02  15
      15:03     15
      15:03:01  15
      401 4
      402 4
      403 4
      404 4
      405 4
      407 4
      161 4
      408 4
      410 4
      411 4
      417 4
      438 4
      701 7
      801 8
      802 8
      804 8
      806 8
      811 8
      901 9", header = TRUE)

  # vector of NA definition based on the unique XX pattern from all of the locus in the existing data
  # 2 values (DONOR_B1 0.796527778, DONOR_DQB11 0.175) we can't identify, set then to NA
  vec_na <- c("xx", "(*xx)", "**xx", "*xx", "*xx#", "*xxxx", "xx", "xxxx", "xxxxx",
              "(*XX)", "**XX", "*XX", "*XX#", "*XXXX", "XX", "XXXX", "XXXXX", "NOT TESTED", "unk", "", " ",
              "0.796527778", "0.175", "un", "N/", "n/")

  # create a temp data frame to hold input antigens
  tmp <- data.frame(cbind(var_1,var_2))

  # this function is used within CleanHla only to simply the calculation on different lookup table
  TempFunction <- function(){
    tmp2 <- cbind(tmp,tmp2)
    colnames(tmp2) <- c("ori_1","ori_2","new_1","new_2")
    tmp2 <- data.frame(lapply(tmp2, as.character), stringsAsFactors = FALSE)
    tmp2$var_1 <- ifelse(!is.na(tmp2$new_1), tmp2$new_1, tmp2$ori_1)
    tmp2$var_2 <- ifelse(!is.na(tmp2$new_2), tmp2$new_2, tmp2$ori_2)
    tmp <- tmp2[,colnames(tmp2) %in% c("var_1","var_2")]
    rm(tmp2)
    return(tmp)
  }

  # step 1 : use lookup table to convert dd:dd:dd notation to corresponding antigen type
  # rationale : some of antigen type are different if they are different on specific protein (the 4th and 5th digitd), have to use a crosstable to assign correct type
  # e.g. : 40:05:00	typing is 50, 40:06:00	is 61
  if (locus == "a" | locus == "A")
  {
    tmp2 <-  tmp
    tmp2[] <- LookupA$type[match(unlist(tmp), LookupA$antigen)]
    tmp <- TempFunction()

    # unlike other locus, in A, content within braces are BW4 which are not the actually typing we need
    # remove braceBW4Brace pattern at here to avoid confusion in step 3
    # patt_a <- c("(BW4)", "[BW4]", "[BW4}",	"{BW4)", "{BW4}",	"{BW4}",	"({BW4}",		"(BW4}",		"BW4)")
    patt_a <- c("\\(BW4\\)|\\[BW4\\]|\\[BW4\\}|\\{BW4\\)|\\{BW4\\}|\\{BW4\\}|\\(\\{BW4\\}|\\(BW4\\}\\|BW4\\)")
    tmp$var_1 <- str_replace_all(tmp$var_1, patt_a, " ")
    tmp$var_2 <- str_replace_all(tmp$var_2, patt_a, " ")
  }

  if (locus == "b" | locus == "B")
  {
    tmp2 <-  tmp
    tmp2[] <- LookupB$type[match(unlist(tmp), LookupB$antigen)]
    tmp <- TempFunction()
  }

  if (locus == "dqb" | locus == "DQB")
  {
    tmp2 <-  tmp
    tmp2[] <- LookupDQB$type[match(unlist(tmp), LookupDQB$antigen)]
    tmp <- TempFunction()
  }

  if (locus == "drb1" | locus == "DRB1")
  {
    tmp2 <-  tmp
    tmp2[] <- LookupDRB1$type[match(unlist(tmp), LookupDRB1$antigen)]
    tmp <- TempFunction()
  }

  # step 2 : fill out NA antigen with non-na antigen within a locus to assume homozygosity
  # logic : if there is an NA antigen on either of the copies, replace it with the antigen from the non-na copy
  #         if both copies are NA, then keep them as NA
  #var_1 <- ifelse(is.na(tmp$var_1) | tmp$var_1 %in% vec_na & tmp$var_2 != "", tmp$var_2, tmp$var_1)
  #var_2 <- ifelse(is.na(tmp$var_2) | tmp$var_2 %in% vec_na & tmp$var_1 != "", tmp$var_1, tmp$var_2)

  var_1 <- ifelse((is.na(tmp$var_1) | tmp$var_1 %in% vec_na) & (is.na(tmp$var_2) | tmp$var_2 %in% vec_na), "",
                  ifelse(is.na(tmp$var_1) | tmp$var_1 %in% vec_na & tmp$var_2 != "", tmp$var_2, tmp$var_1))
  var_2 <- ifelse((is.na(tmp$var_2) | tmp$var_2 %in% vec_na) & (is.na(tmp$var_1) | tmp$var_1 %in% vec_na), "",
                  ifelse(is.na(tmp$var_2) | tmp$var_2 %in% vec_na & tmp$var_1 != "", tmp$var_1, tmp$var_2))

  # step 3 : pull out antigen if they show up within (), [], or {}
  # rationale: some of antigen are within (), [], or {}
  # example : 15(72), *15[62] , or 0302{18}
  var_1_c1 <- ifelse(grepl("\\[", var_1), sub(".*\\[", "", var_1),
                     ifelse(grepl("\\(", var_1), sub(".*\\(", "", var_1),
                            ifelse(grepl("\\{", var_1), sub(".*\\{", "", var_1),var_1)))

  var_2_c1 <- ifelse(grepl("\\[", var_2), sub(".*\\[", "", var_2),
                     ifelse(grepl("\\(", var_2), sub(".*\\(", "", var_2),
                            ifelse(grepl("\\{", var_1), sub(".*\\{", "", var_2),var_2)))

  # step 4 : remove other special symbols like *, #, :., ), ], }, (., and {. , then drop non-numeric characters
  var_1_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\:.*|\\(.*|\\{.*", '', var_1_c1)
  var_2_c2 <- gsub("\\)|\\]|\\}|\\*|\\#|\\:.*|\\(.*|\\{.*", '', var_2_c1)

  var_1_c2 <- as.character(parse_number(var_1_c2))
  var_2_c2 <- as.character(parse_number(var_2_c2))

  # step 5 : keep only first two characters
  var_1_c3 <- map_chr(var_1_c2, ~substr(.,1,2))
  var_2_c3 <- map_chr(var_2_c2, ~substr(.,1,2))

  # step 6 : add a leading 0 if nchar is 1
  # rationale : some of antigen are the same but in different notation
  # example : 2 and 02 are the same
  var_1_out <- map_chr(var_1_c3, function(x) ifelse(nchar(x) == 1, paste0("0",x), x))
  var_2_out <- map_chr(var_2_c3, function(x) ifelse(nchar(x) == 1, paste0("0",x), x))

  # step 7: remove temporary variable holders
  rm(var_1_c1,var_1_c2,var_1_c3,var_2_c1,var_2_c2,var_2_c3,tmp)

  if (locus %in% c("A", "a", "B", "b", "C", "c")) {
    return(list(locus1 = var_1_out, locus2 = var_2_out))
  } else{
    return(list(locus1 = var_1_out, locus2 = var_2_out))
  }

}
