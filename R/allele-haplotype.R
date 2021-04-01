#' @name ImputeHaplo
#' @title Impute low resolution hla to high resolution typing.
#' @param dat_in
#' A data frame with recipient/donor alleles info.
#' @return
#' A data frame with the best pairs of haplotype combination.
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata/example", "Haplotype_test.csv", package = "hlaR"))
#' result <- ImputeHaplo(dat_in = dat)
#' }
#' @export

ImputeHaplo <- function(dat_in){
  #* step 1: import and clean raw haplotype frequency table *#
  raw_hap_tbl <- read.csv(system.file("extdata/ref", "A_C_B_DRB345_DRB1_DQB1.csv", package = "hlaR"), check.names = FALSE) %>%
    rename_all(. %>% tolower) %>%
    select(a, c, b, drb1, dqb1, drb345,
           afa_freq, afa_rank, api_freq, api_rank, cau_freq, cau_rank, his_freq, his_rank, nam_freq, nam_rank) %>%
    mutate(idx = as.numeric(rownames(.)),
           drb = str_sub(drb345, 1, 4)) %>%
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

  #* step 2: format alleles *#
  # throw error if there are none-: punctuations in the data
  ck <- sum(data.frame(sapply(dat_in, str_detect, pattern = "(?!\\:)[[:punct:]]")) %>%
          mutate(across(everything(), as.numeric)), na.rm=TRUE)

  if(ck >= 1){
    stop('there are punctuation marks other than ":" in your data, please check!')
   }
  rm(ck)

  # 1 -> 01 , 1:03 -> 01:03, 02:03:06 -> 02:03
  simple_clean <- function(in_char){
    out_char <- ifelse(nchar(in_char) == 1, paste0("0", in_char), # paste a leading 0 if it's low resolution and one digit
                       ifelse(str_count(in_char, ":") > 1, sub("(:[^:]+):.*", "\\1", in_char), # remove after 2nd : if number of : > 1
                              ifelse(str_detect(in_char, ":") & nchar(gsub(":.*", "", in_char)) == 1, paste0("0", in_char), in_char))) # paste a leading 0 if high resolution and first part has only 1 digit
    return(out_char)
  }

  dat_in <- dat_in %>%
              replace(., is.na(.), "") %>% # unify NA to ""
              arrange(pair_id) %>%
              mutate_all(as.character) %>%
              mutate(a1 = simple_clean(a1),
                     a2 = simple_clean(a2),
                     b1 = simple_clean(b1),
                     b2 = simple_clean(b2),
                     c1 = simple_clean(c1),
                     c2 = simple_clean(c2),
                     drb1 = simple_clean(drb1),
                     drb2 = simple_clean(drb2),
                     dqb1 = simple_clean(dqb1),
                     dqb2 = simple_clean(dqb2),
                     drb31 = simple_clean(drb31),
                     drb32 = simple_clean(drb32),
                     drb41 = simple_clean(drb41),
                     drb42 = simple_clean(drb42),
                     drb51 = simple_clean(drb51),
                     drb52 = simple_clean(drb52))
  #* end of step 2 *#

  #* step 3: seprate data into impute-ready and append-ready tables *#
  # length of allele columns
  len <- length(names(dat_in)[!(names(dat_in) %in% c("pair_id", "ethnicity", "subject_type"))])

  # count number of NA alleles
  dat_interm <- dat_in %>% mutate(count = rowSums(. == "" ))

  # append-ready if data without enthinicity or with NA hlas
  dat_4_app <- dat_interm %>%
                filter(count == len | !(ethnicity %in% c("cau", "afa", "his", "nam", "api"))) %>%
                select(-count) %>%
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
                       cnt_pair = "" ) %>%
                select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, cnt_pair)

  # impute-ready if data with enthinicity and has at least one hla value
  dat_4_imp <- dat_interm %>%
                filter(count != len & ethnicity %in% c("cau", "afa", "his", "nam", "api")) %>%
                select(-count)
  #* end of step 3 *#

  #* step 4: FuncForCompHaplo() for each subjects in dat_4_imp table *#
  num_subj <- dim(dat_4_imp)[1]
  hpl_tp_raw <- vector(mode = "list", length = num_subj)
  hpl_tp_pairs <- vector(mode = "list", length = num_subj)

  for (i in 1:num_subj){
    hpl_tp_pairs[[i]] <- FuncForCompHaplo(tbl_raw = raw_hap_tbl, tbl_in = dat_4_imp[i, ])
  }
  #* end of step 4 *#

  #* step 5: final table - imputed table + dat_4_app *#
  names(hpl_tp_pairs) <- dat_4_imp$pair_id

  hpl_tp_pairs <- as.data.frame(do.call(rbind, hpl_tp_pairs))
  row.names(hpl_tp_pairs) <- seq(1:dim(hpl_tp_pairs)[1])

  hpl_tp_pairs <- rbind(hpl_tp_pairs, dat_4_app) %>%
                  arrange(pair_id, desc(subj)) %>%
                  replace(., is.na(.), "")

  #* end of step 4 *#

  return(hpl_tp_pairs)
}
