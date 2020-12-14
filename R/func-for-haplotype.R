#' @name FuncForCompHaplo
#' @title extracts combinations from halpo frequency table, it's called within CompHaploTbl()
#' @param tbl_in
#' raw haplotype frequency table downloaded from NMDP website
#' @param dat
#' dataframe of alleles
#' @param cut_freq
#' cutoff frequence value of the population
#' @param cut_rank
#' cutoff rank value of each count group
#' @import
#' tidyverse

FuncForCompHaplo <- function(tbl_in, dat, cut_freq = 0.0001, cut_rank = 10) {
  # index holder
  tmp_indx <- c()

  tmp <- tbl_in %>% filter(fst_a %in% c(dat$a1, dat$a2)) %>% pull(unique(idx))
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in %>% filter(fst_b %in% c(dat$b1, dat$b2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in %>% filter(fst_c %in% c(dat$c1, dat$c2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in %>% filter(fst_drb1 %in% c(dat$drb1, dat$drb2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp <- tbl_in %>% filter(fst_dqb1 %in% c(dat$dqb1, dat$dqb2)) %>% pull(idx)
  if (length(tmp > 0)){
    tmp_indx <- c(tmp_indx, tmp)
  }

  tmp_indx <- unique(tmp_indx)

  if(length(tmp_indx) > 0) {
    if(dat$ethnicity == "cau"){
      subdat <- tbl_in %>%
        mutate(id = dat$id,
               cnt_a = ifelse(fst_a %in% c(dat$a1, dat$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(dat$b1, dat$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(dat$c1, dat$c2), 1, 0),
               cnt_drb = ifelse(fst_drb1 %in% c(dat$drb1, dat$drb2), 1, 0),
               cnt_dqb = ifelse(fst_dqb1 %in% c(dat$dqb1, dat$dqb2), 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb + cnt_dqb) %>%
        filter(cnt > 1 & !is.na(cau_rank)) %>%
        select(id, idx, a, b, c, drb1, dqb1, cau_freq, cau_rank, cnt) %>%
        arrange(cnt, cau_rank) %>%
        group_by(cnt) %>%
        mutate(rank = 1:n()) %>% # rank of each count group
        ungroup() %>%
        filter(cau_freq > cut_freq & rank < 10) %>%
        select(-idx) %>%
        arrange(desc(cnt))
    } else if(dat$ethnicity == "afa"){
      subdat <- tbl_in %>%
        mutate(id = dat$id,
               cnt_a = ifelse(fst_a %in% c(dat$a1, dat$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(dat$b1, dat$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(dat$c1, dat$c2), 1, 0),
               cnt_drb = ifelse(fst_drb1 %in% c(dat$drb1, dat$drb2), 1, 0),
               cnt_dqb = ifelse(fst_dqb1 %in% c(dat$dqb1, dat$dqb2), 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb + cnt_dqb) %>%
        filter(cnt > 1 & !is.na(afa_rank)) %>%
        select(id, idx, a, b, c, drb1, dqb1, afa_freq, afa_rank, cnt) %>%
        arrange(cnt, afa_rank) %>%
        group_by(cnt) %>%
        mutate(rank = 1:n()) %>%
        ungroup() %>%
        filter(afa_freq > cut_freq & rank < cut_rank) %>%
        select(-idx) %>%
        arrange(desc(cnt))
    }
  } else{
    subdat <- data.frame(id = dat$id)
  }

  return(subdat)
}
