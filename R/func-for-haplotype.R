#' @name FuncForCompHaplo
#' @title Compare Haplotypes
#' @description This function generates the ranked haplotypes for high resolution imputation. It is called by the ImputeHaplo() function.
#' @param tbl_raw
#' A raw haplotype frequency table downloaded from NMDP website.
#' @param tbl_in
#' A dataframe with alleles info.
#' @import
#' tidyverse

FuncForCompHaplo <- function(tbl_raw, tbl_in) {
  #* step 0: reshape raw data so it can be binded with paired imputed haplotype data  *#
  raw <- tbl_in %>%
    mutate(subj = paste(paste(pair_id, subject_type, sep = "_"), ethnicity, sep = "_"),
           dat_type = paste("raw"),
           freq = NA,
           rank = NA,
           cnt_pair = NA) %>%
    unite(a, a1:a2, sep = "/", na.rm = TRUE) %>%
    unite(b, b1:b2, sep = "/", na.rm = TRUE) %>%
    unite(c, c1:c2, sep = "/", na.rm = TRUE) %>%
    unite(drb1, drb1:drb2, sep = "/", na.rm = TRUE) %>%
    unite(dqb1, dqb1:dqb2, sep = "/", na.rm = TRUE) %>%
    unite(drb345, drb31:drb52, sep = "/", na.rm = TRUE) %>%
    select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, cnt_pair)
  #* end of step 0 *#

  #* step 1: index *#
  #* step 1a : pull out all of indexes which contain high res antigen of each locus  *#
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
  #* end of 1a *#

  #* step 1b: pull out all of indices which contain low res antigen of each locus *#
  a_lo <- PullIndx(gene = "fst_a", nms = c(tbl_in$a1, tbl_in$a2))
  b_lo <- PullIndx(gene = "fst_b", nms = c(tbl_in$b1, tbl_in$b2))
  c_lo <- PullIndx(gene = "fst_c", nms = c(tbl_in$c1, tbl_in$c2))

  drb1_lo <- PullIndx(gene = "fst_drb1", nms = c(tbl_in$drb1, tbl_in$drb2))
  dqb1_lo <- PullIndx(gene = "fst_dqb1", nms = c(tbl_in$dqb1, tbl_in$dqb1))
  drb345_lo <- PullIndx(gene = "fst_drb345", nms = c(tbl_in$drb31, tbl_in$drb32, tbl_in$drb41, tbl_in$drb42, tbl_in$drb51, tbl_in$drb52))

  indx_lo <- unique(c(a_lo, b_lo, c_lo, drb1_lo, dqb1_lo, drb345_lo))

  rm(a_lo, b_lo, c_lo, drb1_lo, dqb1_lo, drb345_lo)
  #* end of step 1b *#
  #* end of step 1 *#

  #* step 2: pull out high res combinations of max group count *#
  tmp_indx <- unique(c(indx_hi, indx_lo))

  hpl_tp_raw <- tbl_raw %>%
    mutate(pair_id = tbl_in$pair_id,
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
  #* end of step 2 *#

  #* step 3: subset by ethnicity specific rank and frequency *#
  if(tbl_in$ethnicity == "cau"){
    hpl_tp_raw <- hpl_tp_raw %>%
      arrange(cau_rank) %>%
      select(pair_id, a, b, c, drb1, dqb1, drb345, cau_freq, cau_rank, cnt, drb)

    tmp <- hpl_tp_raw %>%
      mutate(indx = row_number())

    if (dim(tmp)[1] <= 1){
      hpl_tp_pairs <- tmp %>%
        mutate(pair = 1,
               value = 1) %>%
        select(pair, value, everything())

    } else{
      hpl_tp_pairs <- data.frame(t(combn(tmp$indx, 2))) %>%
        setNames(c("indx1", "indx2")) %>%
        mutate(pair = row_number()) %>%
        pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
        left_join(., tmp, by = c("value" = "indx")) %>%
        select(-index)
    }

  }
  if(tbl_in$ethnicity == "afa"){
    hpl_tp_raw <- hpl_tp_raw %>%
      arrange(afa_rank) %>%
      select(pair_id, a, b, c, drb1, dqb1, drb345, afa_freq, afa_rank, cnt, drb)

    tmp <- hpl_tp_raw %>%
      mutate(indx = row_number())

    if (dim(tmp)[1] <= 1){
      hpl_tp_pairs <- tmp %>%
        mutate(pair = 1,
               value = 1) %>%
        select(pair, value, everything())
    } else{
      hpl_tp_pairs <- data.frame(t(combn(tmp$indx, 2))) %>%
        setNames(c("indx1", "indx2")) %>%
        mutate(pair = row_number()) %>%
        pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
        left_join(., tmp, by = c("value" = "indx")) %>%
        select(-index)
    }
  }
  if(tbl_in$ethnicity == "his"){
    hpl_tp_raw <- hpl_tp_raw %>%
      arrange(his_rank) %>%
      select(pair_id, a, b, c, drb1, dqb1, drb345, his_freq, his_rank, cnt, drb)

    tmp <- hpl_tp_raw %>%
      mutate(indx = row_number())

    if (dim(tmp)[1] <= 1){
      hpl_tp_pairs <- tmp %>%
        mutate(pair = 1,
               value = 1) %>%
        select(pair, value, everything())
    } else{
      hpl_tp_pairs <- data.frame(t(combn(tmp$indx, 2))) %>%
        setNames(c("indx1", "indx2")) %>%
        mutate(pair = row_number()) %>%
        pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
        left_join(., tmp, by = c("value" = "indx")) %>%
        select(-index)
    }
  }
  if(tbl_in$ethnicity == "nam"){
    hpl_tp_raw <- hpl_tp_raw %>%
      arrange(nam_rank) %>%
      select(pair_id, a, b, c, drb1, dqb1, drb345, nam_freq, nam_rank, cnt, drb)

    tmp <- hpl_tp_raw %>%
      mutate(indx = row_number())

    if (dim(tmp)[1] <= 1){
      hpl_tp_pairs <- tmp %>%
        mutate(pair = 1,
               value = 1) %>%
        select(pair, value, everything())
    } else{
      hpl_tp_pairs <- data.frame(t(combn(tmp$indx, 2))) %>%
        setNames(c("indx1", "indx2")) %>%
        mutate(pair = row_number()) %>%
        pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
        left_join(., tmp, by = c("value" = "indx")) %>%
        select(-index)
    }
  }
  if(tbl_in$ethnicity == "api"){
    hpl_tp_raw <- hpl_tp_raw %>%
      arrange(api_rank) %>%
      select(pair_id, a, b, c, drb1, dqb1, drb345, api_freq, api_rank, cnt, drb)

    tmp <- hpl_tp_raw %>%
      mutate(indx = row_number())

    if (dim(tmp)[1] <= 1){
      hpl_tp_pairs <- tmp %>%
        mutate(pair = 1,
               value = 1) %>%
        select(pair, value, everything())
    } else{
      hpl_tp_pairs <- data.frame(t(combn(tmp$indx, 2))) %>%
        setNames(c("indx1", "indx2")) %>%
        mutate(pair = row_number()) %>%
        pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
        left_join(., tmp, by = c("value" = "indx")) %>%
        select(-index)
    }
  }
  #* end of step 3 *#

  #* step 4: combination works *#
  #* step 4a:  keep first 500 pairs lowest average rank within each pair if number of combinations > 1000 *#
  if(dim(hpl_tp_pairs)[1] > 1000){
    hpl_tp_pairs <- hpl_tp_pairs %>%
      setNames(gsub("afa_|cau_|his_|api_|nam_", "", names(.))) %>%
      group_by(pair) %>%
      summarise(avg = mean(rank), .groups = 'drop') %>%
      ungroup() %>%
      left_join(hpl_tp_pairs, ., by  = "pair") %>%
      arrange(avg, pair) %>%
      select(-c(value, avg)) %>%
      filter(row_number() <= 1000)

    num_pair_tmp <- dim(hpl_tp_pairs)[1]/2
    hpl_tp_pairs$pair <- rep(1:num_pair_tmp, each  = 2)
  }
  #* end of 4a *#

  #* step 4b: generate paired table *#
  # if only 1 apir, then directly write to result
  # else count: by pair(a pair of combinations), purpose: count UNIQUE number of raw low/high allele occurs in the imputed pair
  if(dim(hpl_tp_pairs)[1] == 1){
    result <- hpl_tp_pairs %>%
      setNames(gsub("afa_|cau_|his_|nam|api_", "", names(.))) %>%
      mutate(subj = paste(paste(tbl_in$pair_id, tbl_in$subject_type, sep = "_"), tbl_in$ethnicity, sep = "_"),
             dat_type = "imputed",
             cnt_pair = cnt,
             drb345 = paste(drb, drb345, sep = "*")) %>%
      select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, cnt_pair)
  } else{ hpl_tp_pairs_2 <-  hpl_tp_pairs %>%
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
  # loop through all the pairs
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

  hpl_tp_pairs_3 <- left_join(hpl_tp_pairs, cnt_by_pair, by = "pair") %>%
    filter(cnt_pair == max(cnt_pair)) %>%
    setNames(gsub("afa_|cau_|his_|nam_|api_", "", names(.))) %>%
    select(-cnt)

  hpl_tp_pairs_4 <- hpl_tp_pairs_3 %>%
    group_by(pair) %>%
    summarise(avg = mean(rank), .groups = 'drop') %>%
    ungroup() %>%
    left_join(hpl_tp_pairs_3, ., by  = "pair") %>%
    arrange(avg, pair) %>%
    select(-c(avg))

  num_pairs <- dim(hpl_tp_pairs_4)[1]/2

  result <- hpl_tp_pairs_4 %>%
    mutate(pair = rep(1:num_pairs, each  = 2)) %>%
    filter(pair %in% c(1))

  #if(dim(result)[1] > 1 ){# if there are more than 1 imputed combinations (at least 1 pair)
  result <- result %>%
    mutate(subj = paste(paste(tbl_in$pair_id, tbl_in$subject_type, sep = "_"), tbl_in$ethnicity, sep = "_"),
           dat_type = "imputed",
           drb345 = paste(drb, drb345, sep = "*")) %>%
    select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, cnt_pair)
  }
  result <-  rbind(raw, result)

  return(result)
}
