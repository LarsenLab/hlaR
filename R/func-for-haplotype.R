#' @name FuncForCompHaplo
#' @title extracts combinations from halpo frequency table, it's called within CompHaploTbl()
#' @param tbl_raw
#' raw haplotype frequency table downloaded from NMDP website
#' @param tbl_in
#' data frame with alleles info
#' @import
#' tidyverse

FuncForCompHaplo <- function(tbl_raw, tbl_in) {
  #* step 0: raw data *#
  raw <- tbl_in %>%
    mutate(pair = "raw",
           id = ethnicity,
           a = ifelse(!is.na(a1) & !is.na(a2), paste(a1, a2, sep = "/" ),
                      ifelse(!is.na(a1) & is.na(a2), a1,
                             ifelse(is.na(a1) & !is.na(a2), a2, NA))),

           b = ifelse(!is.na(b1) & !is.na(b2), paste(b1, b2, sep = "/" ),
                      ifelse(!is.na(b1) & is.na(b2), b1,
                             ifelse(is.na(b1) & !is.na(b2), b2, NA))),

           c = ifelse(!is.na(c1) & !is.na(c2), paste(c1, c2, sep = "/" ),
                      ifelse(!is.na(c1) & is.na(c2), c1,
                             ifelse(is.na(c1) & !is.na(c2), c2, NA))),

           drb1 = ifelse(!is.na(drb1) & !is.na(drb2), paste(drb1, drb2, sep = "/" ),
                         ifelse(!is.na(drb1) & is.na(drb2), drb1,
                                ifelse(is.na(drb1) & !is.na(drb2), drb2, NA))),

           dqb1 = ifelse(!is.na(dqb1) & !is.na(dqb2), paste(dqb1, dqb2, sep = "/" ),
                         ifelse(!is.na(dqb1) & is.na(dqb2), dqb1,
                                ifelse(is.na(dqb1) & !is.na(dqb2), dqb2, NA))),
           # will modify to bringing drb3/4/5 together - 01/27/2021
           drb345 = ifelse(!is.na(drb41) & !is.na(drb42), paste(drb41, drb42, sep = "/" ),
                           ifelse(!is.na(drb41) & is.na(drb42), drb41,
                                  ifelse(is.na(drb41) & !is.na(drb42), drb42, NA))),
           freq = NA,
           rank = NA,
           cnt_pair = NA) %>%
    select(pair, id, a, b, c, drb1, dqb1, drb345, freq, rank, cnt_pair)
  #* end of step 0 *#

  #* step 1: calculate if a subject with all NA values *#
  sum_na <- tbl_in %>%
    rowwise %>%
    summarise(na_per_row = sum(is.na(.)), .groups = 'drop')

  num_tt <- dim(tbl_in)[2]

  num_na <- num_tt - sum_na$na_per_row
  #* end of step 1 *#

  #* step 2a : pull out all of indices which contain high res antigen of each locus  *#
  PullIndx <- function(gene, nms){
    names <- syms(gene)
    tmp <- tbl_raw %>% select(!!!names, idx) %>%
      filter(.[,1] %in% nms) %>% pull(unique(idx))
    return(tmp)
  }

  a_hi <- PullIndx(gene = "a", nms = c(tbl_in$a1, tbl_in$a2))
  b_hi <- PullIndx(gene = "b", nms = c(tbl_in$b1, tbl_in$b2))
  c_hi <- PullIndx(gene = "c", nms = c(tbl_in$c1, tbl_in$c2))

  drb1_hi <- PullIndx(gene = "drb1", nms = c(tbl_in$drb1, tbl_in$drb2))
  dqb1_hi <- PullIndx(gene = "dqb1", nms = c(tbl_in$dqb1, tbl_in$dqb1))
  drb345_hi <- PullIndx(gene = "drb345", nms = c(tbl_in$drb31, tbl_in$drb32, tbl_in$drb41, tbl_in$drb42, tbl_in$drb51, tbl_in$drb52))

  indx_hi <- unique(c(a_hi, b_hi, c_hi, drb1_hi, dqb1_hi, drb345_hi))
  rm(a_hi, b_hi, c_hi, drb1_hi, dqb1_hi, drb345_hi)
  #* end of 2a *#

  #* step 2b: pull out all of indices which contain low res antigen of each locus *#
  a_lo <- PullIndx(gene = "fst_a", nms = c(tbl_in$a1, tbl_in$a2))
  b_lo <- PullIndx(gene = "fst_b", nms = c(tbl_in$b1, tbl_in$b2))
  c_lo <- PullIndx(gene = "fst_c", nms = c(tbl_in$c1, tbl_in$c2))

  drb1_lo <- PullIndx(gene = "fst_drb1", nms = c(tbl_in$drb1, tbl_in$drb2))
  dqb1_lo <- PullIndx(gene = "fst_dqb1", nms = c(tbl_in$dqb1, tbl_in$dqb1))
  drb345_lo <- PullIndx(gene = "fst_drb345", nms = c(tbl_in$drb31, tbl_in$drb32, tbl_in$drb41, tbl_in$drb42, tbl_in$drb51, tbl_in$drb52))

  indx_lo <- unique(c(a_lo, b_lo, c_lo, drb1_lo, dqb1_lo, drb345_lo))

  rm(a_lo, b_lo, c_lo, drb1_lo, dqb1_lo, drb345_lo)

  #* end of step 2b *#

  #* step 3: pull out high res combinations of max group count *#
  tmp_indx <- unique(c(indx_hi, indx_lo))

  if(length(tmp_indx) > 0 & num_na > 3) {
    hpl_tp_raw <- tbl_raw %>%
      mutate(id = tbl_in$rowid,
             cnt_a = ifelse(fst_a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
             cnt_b = ifelse(fst_b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
             cnt_c = ifelse(fst_c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
             cnt_drb1 = ifelse(fst_drb1 %in% c(tbl_in$drb1, tbl_in$drb2), 1, 0),
             cnt_dqb1 = ifelse(fst_dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2), 1, 0),
             cnt_drb3 = ifelse(fst_drb345 %in% c(tbl_in$drb31, tbl_in$drb32), 1, 0),
             cnt_drb4 = ifelse(fst_drb345 %in% c(tbl_in$drb41, tbl_in$drb42), 1, 0),
             cnt_drb5 = ifelse(fst_drb345 %in% c(tbl_in$drb51, tbl_in$drb52), 1, 0),

             cnt_a_hi = ifelse(a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
             cnt_b_hi = ifelse(b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
             cnt_c_hi = ifelse(c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
             cnt_drb1_hi = ifelse(drb1 %in% c(tbl_in$drb1, tbl_in$drb2), 1, 0),
             cnt_dqb1_hi = ifelse(dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2), 1, 0),
             cnt_drb3_hi = ifelse(drb345 %in% c(tbl_in$drb31, tbl_in$drb32), 1, 0),
             cnt_drb4_hi = ifelse(drb345 %in% c(tbl_in$drb41, tbl_in$drb42), 1, 0),
             cnt_drb5_hi = ifelse(drb345 %in% c(tbl_in$drb51, tbl_in$drb52), 1, 0)) %>%
      mutate(cnt_drb345 = ifelse(cnt_drb3 + cnt_drb4 + cnt_drb5 > 0, 1, 0),
             cnt_drb345_hi = ifelse(cnt_drb3_hi + cnt_drb4_hi + cnt_drb5_hi > 0, 1, 0)) %>%
      mutate(cnt_lo = cnt_a + cnt_b + cnt_c + cnt_drb1 + cnt_dqb1 + cnt_drb345,
             cnt_hi = cnt_a_hi + cnt_b_hi + cnt_c_hi + cnt_drb1_hi + cnt_dqb1_hi + cnt_drb345_hi,
             cnt = cnt_lo + cnt_hi) %>%
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
    hpl_tp_raw <- data.frame(id = tbl_in$rowid)
    hpl_tp_pairs <- data.frame(id = tbl_in$rowid)
  }
  #* end of step 3 *#

  #* step4: generate paired table *#
  if(dim(hpl_tp_pairs)[1] > 1) {
    hpl_tp_pairs_2 <-  hpl_tp_pairs %>%
      mutate(a = sub("\\:.*", "", a),
             b = sub("\\:.*", "", b),
             c = sub("\\:.*", "", c) ,
             drb1 = sub("\\:.*", "", drb1),
             dqb1 = sub("\\:.*", "", dqb1),
             drb345 = sub("\\:.*", "", drb345)) %>%
      mutate(cnt_a = ifelse(a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
             cnt_b = ifelse(b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
             cnt_c = ifelse(c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
             cnt_drb1 = ifelse(drb1 %in% c(tbl_in$drb1, tbl_in$drb2), 1, 0),
             cnt_dqb1 = ifelse(dqb1 %in% c(tbl_in$dqb1, tbl_in$dqb2), 1, 0),
             cnt_drb3 = ifelse(drb345 %in% c(tbl_in$drb31, tbl_in$drb32), 1, 0),
             cnt_drb4 = ifelse(drb345 %in% c(tbl_in$drb41, tbl_in$drb42), 1, 0),
             cnt_drb5 = ifelse(drb345 %in% c(tbl_in$drb51, tbl_in$drb52), 1, 0)) %>%
      mutate(cnt_drb345 = ifelse(cnt_drb3 + cnt_drb4 + cnt_drb5 > 0, 1, 0)) %>%
      mutate(cnt_sum = cnt_a + cnt_b + cnt_c + cnt_drb1 + cnt_dqb1 + cnt_drb345) %>%
      select(-c(cnt_drb3, cnt_drb4, cnt_drb5))

    cnt_by_pair <- data.frame(pair = numeric(),
                              cnt_pair = numeric())

    for (i in 1:max(hpl_tp_pairs_2$pair)){
      tmp <- hpl_tp_pairs_2 %>% filter(pair == i) %>%
        mutate(dup_a = ifelse(a[1] == a[2] & cnt_a == 1, cnt_a, cnt_a + 1), # identical = 0, non-identitcal = 1
               dup_b = ifelse(b[1] == b[2] & cnt_b == 1, cnt_b, cnt_b + 1),
               dup_c = ifelse(c[1] == c[2] & cnt_c == 1, cnt_c, cnt_c + 1),
               dup_drb1 = ifelse(drb1[1] == drb1[2]  & cnt_drb1 == 1, cnt_drb1, cnt_drb1 + 1),
               dup_dqb1 = ifelse(dqb1[1] == dqb1[2] & cnt_dqb1 == 1, cnt_dqb1, cnt_dqb1 + 1),
               dup_drb345 = ifelse(drb345[1] == drb345[2] & cnt_drb345 == 1, cnt_drb345, cnt_drb345 + 1)) %>%
        mutate(cnt_pair = dup_a + dup_b + dup_c + dup_drb1 + dup_dqb1 + dup_drb345) %>%
        select(pair, cnt_pair) %>%
        distinct()

      # paired id
      cnt_by_pair <- rbind(cnt_by_pair, tmp)
    }

    hpl_tp_pairs_3 <- left_join(hpl_tp_pairs, cnt_by_pair,
                                by = "pair") %>%
      filter(cnt_pair == max(cnt_pair)) %>%
      setNames(gsub("afa_|cau_", "", names(.))) %>%
      select(-cnt)

    hpl_tp_pairs_4 <- hpl_tp_pairs_3 %>%
      group_by(pair) %>%
      summarise(avg = mean(rank), .groups = 'drop') %>%
      ungroup() %>%
      left_join(hpl_tp_pairs_3, ., by  = "pair") %>%
      arrange(avg, pair) %>%
      select(-c(value, avg))

    hpl_tp_pairs <- hpl_tp_pairs_4

    num_pairs <- dim(hpl_tp_pairs)[1]/2
    hpl_tp_pairs$pair <- rep(1:num_pairs, each  = 2)

    hpl_tp_pairs <- hpl_tp_pairs %>% filter(pair %in% c(1,2,3))
  }

  #* end of step 4 *#
  result <- rbind(raw, hpl_tp_pairs)

  return(result)
}
