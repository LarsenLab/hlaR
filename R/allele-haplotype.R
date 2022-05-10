#' @name ImputeHaplo
#' @title Imputation
#' @description Impute low or mixed resolution HLA(Human Leukocyte Antigen) typing to the most likely high resolution equivalent. Imputation is computationally intensive, so large dataset may encounter delays in processing. This function uses data from the NMDP(National Marrow Donor Program), and is currently limited to HLA A, B, C, and DRB loci.
#' @param dat_in
#' A data frame with low resolution HLA data.
#' @return
#' A data frame with high resolution HLA data pulled from the most likely pair of haplotypes matching the input low resolution data.
#' @import
#' tidyverse
#' utils
#'
#' @examples
#' \donttest{
#' dat <- read.csv(system.file("extdata/example", "Haplotype_test.csv", package = "hlaR"))
#' result <- ImputeHaplo(dat_in = dat[c(1:2), ])
#' }
#' @export
#'
ImputeHaplo <- function(dat_in){
  warning("\nPlease use this imputation function with caution; its accuracy is lower than the current publicly available gold standard (HaploStats) and may produce inaccurate results.\nWork is underway on a collaboration to improve the accuracy of this function.\n")

  #* step 1: import and clean raw haplotype frequency table *#
  p1 <- read.csv("https://raw.githubusercontent.com/LarsenLab/public-data/master/A_C_B_DRB345_DRB1_DQB1_part1.csv", check.names = FALSE)
  p2 <- read.csv("https://raw.githubusercontent.com/LarsenLab/public-data/master/A_C_B_DRB345_DRB1_DQB1_part2.csv", check.names = FALSE)
  p3 <- read.csv("https://raw.githubusercontent.com/LarsenLab/public-data/master/A_C_B_DRB345_DRB1_DQB1_part3.csv", check.names = FALSE)

  tbl_ref <- data.frame(rbind(p1, p2, p3)) %>%
    rename_all(. %>% tolower) %>%
    mutate(idx = as.numeric(rownames(.)),
           drb = str_sub(.$drb345, 1, 4)) %>% # use this syntax to suppress "no symbol named ** in the scope" warning
    # remove trailing g
    mutate(a = ifelse(str_detect(a, "g"), str_replace(a, "g", ""), a),
           b = ifelse(str_detect(b, "g"), str_replace(b, "g", ""), b),
           c = ifelse(str_detect(c, "g"), str_replace(c, "g", ""), c),
           drb1 = ifelse(str_detect(drb1, "g"), str_replace(drb1, "g", ""), drb1),
           dqb1 = ifelse(str_detect(dqb1, "g"), str_replace(dqb1, "g", ""), dqb1),
           drb345 = ifelse(str_detect(drb345, "g"), str_replace(drb345, "g", ""), drb345)) %>%
    # remove leading loci letter
    mutate(a = str_sub(a, start = 3),
           b = str_sub(b, start = 3),
           c = str_sub(c, start = 3),
           drb1 = str_sub(drb1, start = 6),
           dqb1 = str_sub(dqb1, start = 6),
           drb345 = str_sub(drb345, start = 6)) %>%
    # split allele by ":" to compare with input data
    mutate(fst_a = sub("\\:.*", "", a),
           lst_a = sub(".*\\:", "", a),
           fst_b = sub("\\:.*", "", b),
           lst_b = sub(".*\\:", "", b),
           fst_c = sub("\\:.*", "", c),
           lst_c = sub(".*\\:", "", c),
           fst_drb1 = sub("\\:.*", "", drb1),
           lst_drb1 = sub(".*\\:", "", drb1),
           fst_dqb1 = sub("\\:.*", "", dqb1),
           lst_dqb1 = sub(".*\\:", "", dqb1),
           fst_drb345 = sub("\\:.*", "", drb345),
           lst_drb345 = sub(".*\\:", "", drb345))
  #* end of step 1 *#

  #* step 2: some definitions *#
  # column names in the return table
  nms_keep <- c("a", "c", "b", "drb345", "drb1", "dqb1",
                "idx", "drb",
                "fst_a", "lst_a", "fst_b", "lst_b", "fst_c", "lst_c",
                "fst_drb1", "lst_drb1", "fst_dqb1", "lst_dqb1", "fst_drb345", "lst_drb345" )

  # list of ethnicity
  nms_race <- names(tbl_ref)[str_detect(names(tbl_ref), "_rank")] %>% str_remove(., "_rank")
  #* end of step 2 *#

  #* step 3: generate search table for known/unknown ethnicity *#
  tbl_all <- data.frame(matrix(ncol = length(nms_keep) + 3, nrow = 0))
  colnames(tbl_all) <- c(nms_keep, "ethnicity", "rank", "freq")

  for(i in 1:length(nms_race)){
    in_rank <- paste(nms_race[i], "rank", sep = "_")
    in_freq <- paste(nms_race[i], "freq", sep = "_")
    tmp <- tbl_ref %>%
      filter(!is.na(!!as.name(in_rank))) %>%
      mutate(ethnicity = nms_race[i]) %>%
      select(c(all_of(nms_keep), ethnicity, all_of(in_rank), all_of(in_freq))) %>%
      setNames(c(nms_keep,"ethnicity", "rank", "freq"))

    tbl_all <- rbind(tbl_all, tmp)
  }

  tbl4eth <- tbl_all

  tbl_uniq <- tbl_all %>%
    group_by(a, c, b, drb345, drb1, dqb1) %>%
    filter(freq == max(freq)) %>% # if duplicate then keep the record with max freq(min rank)
    ungroup() %>%
    data.frame()

  tbl4NAeth <- tbl_uniq
  rm(tbl_all, tbl_uniq)
  #* end of step 3 *#

  #* step 4: format alleles *#
  # throw error if there are punctuation except colon(:) in the string
  non_punc <- sum(data.frame(sapply(dat_in, str_detect, pattern = "(?!\\:)[[:punct:]]")) %>%
                    mutate(across(everything(), as.numeric)), na.rm=TRUE)

  if(non_punc >= 1){
    stop('there are punctuation marks other than ":" in your data, please check!')
  }
  rm(non_punc)

  # 1 -> 01 , 1:03 -> 01:03, 02:03:06 -> 02:03
  simple_clean <- function(in_char){
    out_char <- ifelse(nchar(in_char) == 1, paste0("0", in_char), # paste a leading 0 if it's low resolution and one digit
                       ifelse(str_count(in_char, ":") > 1, sub("(:[^:]+):.*", "\\1", in_char), # remove after 2nd : if number of : > 1
                              ifelse(str_detect(in_char, ":") & nchar(gsub(":.*", "", in_char)) == 1, paste0("0", in_char), in_char))) # paste a leading 0 if high resolution and first part has only 1 digit
    return(out_char)
  }

  dat_in <- dat_in %>%
    replace(., is.na(.), "") %>%
    # if one of the allele is non-NA and the other one is NA, then use non_na value for both loci
    mutate(a1 = ifelse(a1 == "" & a2 != "", a2, a1),
           a2 = ifelse(a1 != "" & a2 == "" , a1, a2),

           b1 = ifelse(b1 == "" & b2 != "", b2, b1),
           b2 = ifelse(b1 != "" & b2 == "", b1, b2),

           c1 = ifelse(c1 == ""  & c2 != "", c2, c1),
           c2 = ifelse(c1 != "" & c2 == "", c1, c2),

           drb1 = ifelse(drb1 == "" & drb2 != "", drb2, drb1),
           drb2 = ifelse(drb1 != "" & drb2 == "", drb1, drb2),

           dqb1 = ifelse(dqb1 == "" & dqb2 != "", dqb2, dqb1),
           dqb2 = ifelse(dqb1 != "" & dqb2 == "", dqb1, dqb2),

           drb31 = ifelse(drb31 == "" & drb32 != "", drb32, drb31),
           drb32 = ifelse(drb31 != "" & drb32 == "", drb31, drb32),

           drb41 = ifelse(drb41 == "" & drb42 != "", drb42, drb41),
           drb42 = ifelse(drb41 != "" & drb42 == "", drb41, drb42),

           drb51 = ifelse(drb51 == "" & drb52 != "", drb52, drb51),
           drb52 = ifelse(drb51 != "" & drb52 == "", drb51, drb52)) %>%
    mutate_all(as.character) %>%
    arrange(pair_id) %>%
    mutate(a1 = simple_clean(a1),
           a2 = simple_clean(a2),
           b1 = simple_clean(b1),
           b2 = simple_clean(b2),
           c1 = simple_clean(c1),
           c2 = simple_clean(c2),
           drb1.1 = simple_clean(drb1), # use num.num naming convention to match names of imputation result in later step
           drb1.2 = simple_clean(drb2),
           dqb1.1 = simple_clean(dqb1),
           dqb1.2 = simple_clean(dqb2),
           drb3.1 = simple_clean(drb31),
           drb3.2 = simple_clean(drb32),
           drb4.1 = simple_clean(drb41),
           drb4.2 = simple_clean(drb42),
           drb5.1 = simple_clean(drb51),
           drb5.2 = simple_clean(drb52)) %>%
    select(-c(drb1, drb2, dqb1, dqb2, drb31, drb32, drb41, drb42, drb51, drb52))
  #* end of step 4 *#

  #* step 5: separate data into impute-ready and append-ready tables *#
  # count number of NA alleles, total 16 alleles
  hla_nm <- names(dat_in)[!(names(dat_in) %in% c("pair_id", "ethnicity", "subject_type"))]
  dat_interm <- dat_in %>% mutate(count = rowSums(.[,hla_nm] == "" ))

  # append-ready if all NA hlas
  dat_4_app <- dat_interm %>%
    filter(count == 16 ) %>%
    mutate(subj = paste(paste(pair_id, subject_type, sep = "_"), ethnicity, sep = "_"),
           dat_type = paste("raw"),
           a = "",
           b = "",
           c = "",
           drb1 = "",
           dqb1 = "",
           drb345 = "",
           freq = "",
           rank = "",
           lores_all_match = "") %>%
    select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, lores_all_match)

  # impute-ready if data has at least one hla value
  dat_4_imp <- dat_interm %>% filter(count != 16 )
  #* end of step 5 *#

  #* step 6: FuncForCompHaplo() for each subjects in dat_4_imp table *#
  num_subj <- dim(dat_4_imp)[1]
  hpl_tp_raw <- vector(mode = "list", length = num_subj)
  hpl_tp_pairs <- vector(mode = "list", length = num_subj)

  pb <- txtProgressBar(min = 0, max = num_subj, initial = 0, style = 3)
  for (i in 1:num_subj){
    setTxtProgressBar(pb, i)
    # print(paste0("working on subject #", i))
    if(dat_4_imp[i, ]$ethnicity == ""){
      tbl_ref <- tbl4NAeth
    } else{
      tbl_ref <- tbl4eth %>% filter(ethnicity == dat_4_imp[i, ]$ethnicity)
    }
    hpl_tp_pairs[[i]] <- FuncForCompHaplo(tbl_raw = tbl_ref, tbl_in = dat_4_imp[i, ])

    Sys.sleep(time = 1)
  }
  close(pb)
  #* end of step 6 *#

  #* step 7: final table - imputed table + dat_4_app *#
  names(hpl_tp_pairs) <- dat_4_imp$pair_id

  hpl_tp_pairs <- as.data.frame(do.call(rbind, hpl_tp_pairs))
  row.names(hpl_tp_pairs) <- seq(1:dim(hpl_tp_pairs)[1])

  hpl_tp_pairs <- rbind(hpl_tp_pairs, dat_4_app) %>%
    arrange(pair_id, desc(subj)) %>%
    mutate(freq = ifelse(as.numeric(freq) <= 0.0001, 0.0001, round(as.numeric(freq), 4)))  %>%
    replace(., is.na(.), "")
  #* end of step 7 *#

  return(hpl_tp_pairs)
}
