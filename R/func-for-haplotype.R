#' @name FuncForCompHaplo
#' @title extracts combinations from halpo frequency table, it's called within CompHaploTbl()
#' @param tbl_raw
#' raw haplotype frequency table downloaded from NMDP website
#' @param tbl_in
#' data frame with alleles info
#' @import
#' tidyverse

FuncForCompHaplo <- function(tbl_raw, tbl_in) {
  #* step 1: calculate if a subject with all NA values *#
  sum_na <- tbl_in %>%
            rowwise %>%
            summarise(na_per_row = sum(is.na(.)), .groups = 'drop')

  num_tt <- dim(tbl_in)[2]

  num_na <- num_tt - sum_na$na_per_row
  #* end of step 1 *#

  #* step 2: pull out all of indices which contain low res antigen of each locus *#
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
  #* end of step 2 *#

  #* step 3: pull out high res combinations of max goupr count *#
  if(length(tmp_indx) > 0 & num_na > 3) {
      hpl_tp_raw <- tbl_raw %>%
        mutate(id = tbl_in$id,
               cnt_a = ifelse(fst_a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
               cnt_b = ifelse(fst_b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
               cnt_c = ifelse(fst_c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
               cnt_drb1 = ifelse(fst_drb1 %in% c(tbl_in$drb1, tbl_in$drb2), 1, 0),
               cnt_dqb1 = ifelse(fst_dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2), 1, 0),
               cnt_drb3 = ifelse(fst_drb345 %in% c(tbl_in$drb31, tbl_in$drb32), 1, 0),
               cnt_drb4 = ifelse(fst_drb345 %in% c(tbl_in$drb41, tbl_in$drb42), 1, 0),
               cnt_drb5 = ifelse(fst_drb345 %in% c(tbl_in$drb51, tbl_in$drb52), 1, 0)) %>%
        mutate(cnt_drb345 = ifelse(cnt_drb3 + cnt_drb4 + cnt_drb5 > 1, 1, 0)) %>%
        mutate(cnt = cnt_a + cnt_b + cnt_c + cnt_drb1 + cnt_dqb1 + cnt_drb345) %>%
        filter(cnt == max(cnt))

      if(tbl_in$ethnicity == "cau"){
            hpl_tp_raw <- hpl_tp_raw %>%
                      arrange(cau_rank) %>%
                       select(id, a, b, c, drb1, dqb1, drb345, cau_freq, cau_rank, cnt)

            tmp <- hpl_tp_raw %>%
              mutate(indx = row_number())

            hpl_tp_pairs <- data.frame(t(combn(tmp$indx,2))) %>%
              setNames(c("indx1", "indx2")) %>%
              mutate(pair = row_number()) %>%
              pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
              left_join(., tmp, by = c("value" = "indx")) %>%
              select(-index)


      } else if(tbl_in$ethnicity == "afa"){
      hpl_tp_raw <- hpl_tp_raw %>%
                arrange(afa_rank) %>%
                select(id, a, b, c, drb1, dqb1, drb345, afa_freq, afa_rank, cnt)

      tmp <- hpl_tp_raw %>%
        mutate(indx = row_number())

      hpl_tp_pairs <- data.frame(t(combn(tmp$indx,2))) %>%
        setNames(c("indx1", "indx2")) %>%
        mutate(pair = row_number()) %>%
        pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
        left_join(., tmp, by = c("value" = "indx")) %>%
        select(-index)
      }
  } else{
    hpl_tp_raw <- data.frame(id = tbl_in$id)
    hpl_tp_pairs <- data.frame(id = tbl_in$id)
  }
  #* end of step 3 *#
  return(list(hpl_tp_raw = hpl_tp_raw, hpl_tp_pairs = hpl_tp_pairs))
}
