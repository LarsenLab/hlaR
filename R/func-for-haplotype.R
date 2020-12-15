#' @name FuncForCompHaplo
#' @title extracts combinations from halpo frequency table, it's called within CompHaploTbl()
#' @param tbl_in1
#' raw haplotype frequency table downloaded from NMDP website
#' @param tbl_in2
#' tbl_in2aframe of alleles
#' @param cut_freq
#' cutoff frequence value of the population
#' @param cut_rank
#' cutoff rank value of each count group
#' @import
#' tidyverse

FuncForCompHaplo <- function(tbl_in1, tbl_in2, cut_freq = 0.0001, cut_rank = 10) {
  # index holder
  tmp_indx <- c()

  tmp <- tbl_in1 %>% filter(fst_a %in% c(tbl_in2$a1, tbl_in2$a2)) %>% pull(unique(idx))
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in1 %>% filter(fst_b %in% c(tbl_in2$b1, tbl_in2$b2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in1 %>% filter(fst_c %in% c(tbl_in2$c1, tbl_in2$c2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in1 %>% filter(fst_drb1 %in% c(tbl_in2$drb1, tbl_in2$drb2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in1 %>% filter(fst_dqb1 %in% c(tbl_in2$dqb1, tbl_in2$dqb2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in1 %>% filter(fst_drb345 %in% c(tbl_in2$drb345_1, tbl_in2$drb345_2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp_indx <- unique(tmp_indx)

  if(length(tmp_indx) > 0) {
    if(tbl_in2$ethnicity == "cau"){
      subdat <- tbl_in1 %>%
        mutate(id = tbl_in2$id,
               cnt_a = ifelse(fst_a %in% c(tbl_in2$a1, tbl_in2$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(tbl_in2$b1, tbl_in2$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(tbl_in2$c1, tbl_in2$c2), 1, 0),
               cnt_drb = ifelse(fst_drb1 %in% c(tbl_in2$drb1, tbl_in2$drb2), 1, 0),
               cnt_dqb = ifelse(fst_dqb1 %in% c(tbl_in2$dqb1, tbl_in2$dqb2), 1, 0),
               cnt_drb345 = ifelse(fst_drb345 %in% c(tbl_in2$drb345_1, tbl_in2$drb345_2), 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb + cnt_dqb + cnt_drb345) %>%
        filter(cnt > 1 & !is.na(cau_rank)) %>%
        select(id, idx, a, b, c, drb1, dqb1, drb345, cau_freq, cau_rank, cnt) %>%
        arrange(cnt, cau_rank) %>%
        group_by(cnt) %>%
        mutate(rank = 1:n()) %>% # rank of each count group
        ungroup() %>%
        filter(cau_freq > cut_freq & rank < 10) %>%
        select(-idx) %>%
        arrange(desc(cnt))
    } else if(tbl_in2$ethnicity == "afa"){
      subdat <- tbl_in1 %>%
        mutate(id = tbl_in2$id,
               cnt_a = ifelse(fst_a %in% c(tbl_in2$a1, tbl_in2$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(tbl_in2$b1, tbl_in2$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(tbl_in2$c1, tbl_in2$c2), 1, 0),
               cnt_drb = ifelse(fst_drb1 %in% c(tbl_in2$drb1, tbl_in2$drb2), 1, 0),
               cnt_dqb = ifelse(fst_dqb1 %in% c(tbl_in2$dqb1, tbl_in2$dqb2), 1, 0),
               cnt_drb345 = ifelse(fst_drb345 %in% c(tbl_in2$drb345_1, tbl_in2$drb345_2), 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb + cnt_dqb + cnt_drb345) %>%
        filter(cnt > 1 & !is.na(afa_rank)) %>%
        select(id, idx, a, b, c, drb1, dqb1, drb345, afa_freq, afa_rank, cnt) %>%
        arrange(cnt, afa_rank) %>%
        group_by(cnt) %>%
        mutate(rank = 1:n()) %>%
        ungroup() %>%
        filter(afa_freq > cut_freq & rank < cut_rank) %>%
        select(-idx) %>%
        arrange(desc(cnt))
    }
  } else{
    subdat <- data.frame(id = tbl_in2$id)
  }

  return(subdat)
}
