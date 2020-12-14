#' @name CompHaploTbl
#' @title compare alleles with NMDP frequency table to get top combination of haplotypes based on user cutoffs
#' @param dat_in
#' dataframe with recipient/donor alleles info
#' @param cut_p
#' cutoff frequence value of the population, it passed to FuncForCompHaplo()
#' @param cut_r
#' cutoff rank value of each count group, it passed to FuncForCompHaplo()
#' @return
#' a list of dataframe of most matched allele combination of each subject
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read.csv("~/projects/DEV/haplostats_dev/data/csv/tx_cohort_clean.csv"))
#' result <- CompHaploTbl(dat_in = dat, cut_p= 0.0001, cut_r = 10)
#' }
#' @export

CompHaploTbl <- function(dat_in, cut_p, cut_r){
  #* step 1: import raw haplotype frequenc table and do a brief cleaning *#
    # raw_hap_tbl <- read.csv("~/projects/DEV/haplostats_dev/data/csv/A_C_B_DRB1_DQB1.csv") %>%
    raw_hap_tbl <- read.csv(system.file("extdata", "A_C_B_DRB1_DQB1.csv", package = "hlaR"), check.names = FALSE) %>%
    rename_all(. %>% tolower) %>%
    select(a, c, b, drb1, dqb1,
           afa_freq, afa_rank, api_freq, api_rank, cau_freq, cau_rank, his_freq, his_rank, nam_freq, nam_rank) %>%
    mutate(idx = as.numeric(rownames(.))) %>%
    # remove trailing g
    mutate(a = ifelse(str_detect(a, "g"), str_replace(a, "g", ""), a),
           b = ifelse(str_detect(b, "g"), str_replace(b, "g", ""), b),
           c = ifelse(str_detect(c, "g"), str_replace(c, "g", ""), c),
           drb1 = ifelse(str_detect(drb1, "g"), str_replace(drb1, "g", ""), drb1),
           dqb1 = ifelse(str_detect(dqb1, "g"), str_replace(dqb1, "g", ""), dqb1)) %>%
    # remove leading loci letter
    mutate(a = str_sub(a, start = 3),
           b = str_sub(b, start = 3),
           c = str_sub(c, start = 3),
           drb1 = str_sub(drb1, start = 6),
           dqb1 = str_sub(dqb1, start = 6)) %>%
    # split allele by : for comparison with input data
    mutate(fst_a = sub("\\:.*", "", a),
           lst_a = sub(".*\\:", "", a),
           fst_b = sub("\\:.*", "", b),
           lst_b = sub(".*\\:", "", b),
           fst_c = sub("\\:.*", "", c),
           lst_c = sub(".*\\:", "", c),
           fst_drb1 = sub("\\:.*", "", drb1),
           lst_drb1 = sub(".*\\:", "", drb1),
           fst_dqb1 = sub("\\:.*", "", dqb1),
           lst_dqb1 = sub(".*\\:", "", dqb1))
  #* end of step 1 *#

  #* step 2: rechape input data table by recipient and donor *#
  rcpt <- dat_in %>%
           rename_all(. %>% tolower) %>%
           select(contains(c("id", "ethnicity","race","rcpt"))) %>%
           rename_at(vars(contains('rcpt_')), list(~sub('rcpt_', '', .))) %>%
           mutate(id = paste("rcpt", rowid, sep = "_"),
                  a1 = ifelse(nchar(a1) == 1, paste0("0", a1), a1)) %>%
           select(id, everything())

  don <- dat_in %>%
          rename_all(. %>% tolower) %>%
          select(contains(c("id", "ethnicity","race","don"))) %>%
          rename_at(vars(contains('don_')), list(~sub('don_', '', .))) %>%
          mutate(id = paste("don", rowid, sep = "_")) %>%
          select(id, everything())

  dat2 <- rbind(rcpt, don) %>% arrange(rowid) %>%
    mutate_all(as.character) %>%
    mutate(a1 = ifelse(nchar(a1) == 1, paste0("0", a1), a1),
           a2 = ifelse(nchar(a2) == 1, paste0("0", a2), a2),
           b1 = ifelse(nchar(b1) == 1, paste0("0", b1), b1),
           b2 = ifelse(nchar(b2) == 1, paste0("0", b2), b2),
           c1 = ifelse(nchar(c1) == 1, paste0("0", c1), c1),
           c2 = ifelse(nchar(c2) == 1, paste0("0", c2), c2),
           drb1 = ifelse(nchar(drb1) == 1, paste0("0", drb1), drb1),
           drb2 = ifelse(nchar(drb2) == 1, paste0("0", drb2), drb2),
           dqb1 = ifelse(nchar(dqb1) == 1, paste0("0", dqb1), dqb1),
           dqb2 = ifelse(nchar(dqb2) == 1, paste0("0", dqb2), dqb2))

  rm(rcpt, don)
  #* end of step 2 *#

  #* step 3: call FuncForCompHaplo() for each subjects and get top combination of alleles *#
  num_subj <- dim(dat2)[1]
  result <- vector(mode = "list", length = num_subj)

  for (i in 1:num_subj){
    result[[i]] <- FuncForCompHaplo(tbl_in1 = raw_hap_tbl, tbl_in2 = dat2[i, ], cut_freq = cut_p, cut_rank = cut_r)
  }
  #* end of step 3 *#
  names(result) <- dat2$id
  return(result)
}

