#' @name FuncForCompHaplo
#' @title extracts combinations from halpo frequency table, it's called within CompHaploTbl()
#' @param tbl_raw
#' raw haplotype frequency table downloaded from NMDP website
#' @param tbl_in
#' data frame with alleles info
#' @param cut_freq
#' cutoff frequence value of the population
#' @param cut_num
#' cutoff number of how many top combinations will be kept in the final table
#' @import
#' tidyverse

FuncForCompHaplo <- function(tbl_raw, tbl_in, cut_freq = 0.0001, cut_num = 10) {
  # index holder
  tmp_indx <- c()

  tmp <- tbl_raw %>% filter(fst_a %in% c(tbl_in$a1, tbl_in$a2)) %>% pull(unique(idx))
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_b %in% c(tbl_in$b1, tbl_in$b2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_c %in% c(tbl_in$c1, tbl_in$c2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_drb1 %in% c(tbl_in$drb1, tbl_in$drb2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_drb345 %in% c(tbl_in$drb31, tbl_in$drb32)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_drb345 %in% c(tbl_in$drb41, tbl_in$drb42)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_raw %>% filter(fst_drb345 %in% c(tbl_in$drb51, tbl_in$drb52)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp_indx <- unique(tmp_indx)

  if(length(tmp_indx) > 0) {
    if(tbl_in$ethnicity == "cau"){
      subdat <- tbl_raw %>%
        mutate(id = tbl_in$id,
               cnt_a = ifelse(fst_a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
               cnt_drb1 = ifelse(fst_drb1 %in% c(tbl_in$drb1, tbl_in$drb2), 1, 0),
               cnt_dqb1 = ifelse(fst_dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2), 1, 0),
               cnt_drb3 = ifelse(fst_drb345 %in% c(tbl_in$drb31, tbl_in$drb32), 1, 0),
               cnt_drb4 = ifelse(fst_drb345 %in% c(tbl_in$drb41, tbl_in$drb42), 1, 0),
               cnt_drb5 = ifelse(fst_drb345 %in% c(tbl_in$drb51, tbl_in$drb52), 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb1 + cnt_dqb1 + cnt_drb3 + cnt_drb4 + cnt_drb5) %>%
        filter(cnt > 1 & !is.na(cau_rank) & cau_freq > cut_freq) %>%
        select(id, idx, a, b, c, drb1, dqb1, drb345, cau_freq, cau_rank, cnt) %>%
        arrange(-cnt, cau_rank) %>%
        slice(1:cut_num) %>%
        select(-c(idx))
    } else if(tbl_in$ethnicity == "afa"){
      subdat <- tbl_raw %>%
        mutate(id = tbl_in$id,
               cnt_a = ifelse(fst_a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
               cnt_drb1 = ifelse(fst_drb1 %in% c(tbl_in$drb1, tbl_in$drb2), 1, 0),
               cnt_dqb1 = ifelse(fst_dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2), 1, 0),
               cnt_drb3 = ifelse(fst_drb345 %in% c(tbl_in$drb31, tbl_in$drb32), 1, 0),
               cnt_drb4 = ifelse(fst_drb345 %in% c(tbl_in$drb41, tbl_in$drb42), 1, 0),
               cnt_drb5 = ifelse(fst_drb345 %in% c(tbl_in$drb51, tbl_in$drb52), 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb1 + cnt_dqb1 + cnt_drb3 + cnt_drb4 + cnt_drb5) %>%
        filter(cnt > 1 & !is.na(afa_rank) & afa_freq > cut_freq) %>%
        select(id, idx, a, b, c, drb1, dqb1, drb345, afa_freq, afa_rank, cnt) %>%
        arrange(-cnt, afa_rank) %>%
        slice(1:cut_num) %>%
        select(-c(idx))
    }
  } else{
    subdat <- data.frame(id = tbl_in$id)
  }

  return(subdat)
}
