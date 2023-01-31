#' Basic functions
#'
#' GenerateLookup() called in CalEpletMHCII()
#' @param lkup_in data table
#' @param locus_in string
#' CalRiskScore() calculate DR DQ risk score, it's called in CalEpletMHCII()
#' @param dat_in dataframe
#' FuncForCompHaplo() called in ImputeHaplo()
#' @param tbl_raw data frame
#' @param tbl_in data frame
#' na_to_empty_string()
#' @param df dataframe
#' @name utils
#' @import
#' janitor
#' tidyverse
#' utils
#' tidyselect
NULL
#> NULL

#' @rdname utils
GenerateLookup <- function(lkup_in, locus_in){
  dat_out <- lkup_in %>%
    filter(locus %in% locus_in) %>%
    mutate(index = as.numeric(sub(".*\\_", "", variable)),
           type = sub("\\_.*", "", variable)) %>%
    select(index, type, value) %>%
    dplyr::rename(eplet = value) %>%
    filter(eplet != "") %>%
    distinct()

  return(dat_out)
}

#' @rdname utils
CalRiskScore <- function(dat_in) {
  #* start of old dq score: dq = sum(max(dqa1, dqb1), max(dqa2, dqb2)) *#
  # pre_ck <- dat_in %>%
  #   mutate(hla = toupper(hla)) %>%
  #   filter(!(mm_eplets %in% c("Not Found"))) %>%
  #   filter(str_detect(hla, "DRB|DQ")) %>%
  #   mutate(hla = gsub("\\*.*", "", hla),
  #          grp = ifelse(str_detect(hla, "DRB"), "DR",
  #                       ifelse(hla %in% c("DQA1", "DQB1") & haplotype_id %in% c("1"), "DQ1",
  #                              ifelse(hla %in% c("DQA1", "DQB1") & haplotype_id %in% c("2"), "DQ2", hla)))) %>%
  #   filter(grp %in% c("DR", "DQ1", "DQ2"))
  #
  # if(dim(pre_ck)[1] == 0 | !all(c("DR", "DQ1", "DQ2") %in% unique(pre_ck$grp))) {
  #   #stop("DR DQ risk score calculation requires eplet mismatch info on both DR DQ alleles. Please check your data.")
  #   warning("DR and DQ risk score calculation require input of both DR and DQ. No score is calculated at this time.")
  #   risk_scr <- data.frame(pair_id = NA, DQ = NA, DR = NA, risk = NA)
  # } else{
  #   risk_scr <- pre_ck %>%
  #     group_by(grp) %>%
  #     summarise(max1 = suppressWarnings(max(mm_cnt, na.rm=TRUE)),
  #               sum1 = sum(mm_cnt, na.rm=TRUE)) %>%
  #     ungroup()  %>%
  #     mutate(max1 = ifelse(grp %in% c("DQ1", "DQ2"), sum1, max1)) %>%
  #     select(-sum1) %>%
  #     mutate(locus = ifelse(grp %in% c("DR"), "DR", "DQ")) %>%
  #     group_by(locus) %>%
  #     summarise(score = suppressWarnings(max(max1, na.rm=TRUE))) %>%
  #     ungroup() %>%
  #     rownames_to_column %>%
  #     gather(var, value, -rowname) %>%
  #     spread(rowname, value) %>%
  #     janitor::row_to_names(1) %>%
  #     mutate(DQ = as.numeric(DQ),
  #            DR = as.numeric(DR),
  #            risk = ifelse(between(DQ, 15, 31), "high",
  #                          ifelse((DR >= 7 & DQ <= 14) | (DR < 7 & between(DQ, 9, 15)), "interm",
  #                                 ifelse(DR < 7 & DQ < 9, "low", "out of bound"))))
  #}
 #* end of old dq score *#

  #* start of new dq score: dq = sum(max(dqa 1, 2), max(dqb 1, 2)) *#
  pre_ck <- dat_in %>%
    mutate(hla = toupper(hla)) %>%
    filter(!(mm_eplets %in% c("Not Found"))) %>%
    filter(str_detect(hla, "DRB|DQ")) %>%
    mutate(hla = gsub("\\*.*", "", hla),
           grp = ifelse(str_detect(hla, "DRB"), "DR",
                        ifelse(hla %in% c("DQA1", "DQA2"), "DQA",
                               ifelse(hla %in% c("DQB1", "DQB2"), "DQB", hla)))) %>%
    filter(grp %in% c("DR", "DQA", "DQB"))

  if(dim(pre_ck)[1] == 0 | !all(c("DR", "DQA", "DQB") %in% unique(pre_ck$grp))) {
    #stop("DR DQ risk score calculation requires eplet mismatch info on both DR DQ alleles. Please check your data.")
    warning("DR and DQ risk score calculation require input of both DR and DQ. No score is calculated at this time.")
    risk_scr <- data.frame(pair_id = NA, DQ = NA, DR = NA, risk = NA)
  } else{
    risk_scr <- pre_ck %>%
      group_by(grp) %>%
      summarise(scr = suppressWarnings(max(mm_cnt, na.rm=TRUE)) ) %>%
      ungroup()  %>%
      mutate(locus = ifelse(grp %in% c("DR"), "DR", "DQ")) %>%
      group_by(locus) %>%
      summarise(score = suppressWarnings(sum(scr, na.rm=TRUE))) %>%
      ungroup() %>%
      as.data.frame() %>%
      rownames_to_column %>%
      gather(var, value, -rowname) %>%
      spread(rowname, value) %>%
      janitor::row_to_names(1) %>%
      mutate(DQ = as.numeric(DQ),
             DR = as.numeric(DR),
             risk = ifelse(between(DQ, 15, 31), "high",
                           ifelse((DR >= 7 & DQ <= 14) | (DR < 7 & between(DQ, 9, 15)), "interm",
                                  ifelse(DR < 7 & DQ < 9, "low", "out of bound"))))
  }
  return(risk_scr)
}

#' @rdname utils
FuncForCompHaplo <- function(tbl_raw, tbl_in) {
  #* step 0: reshape raw data so it can be bind with paired imputed haplotype data  *#
  raw <- tbl_in %>%
    mutate(subj = paste(paste(pair_id, subject_type, sep = "_"), ethnicity, sep = "_"),
           dat_type = paste("raw"),
           freq = "",
           rank = "",
           lores_all_match = "") %>%
    unite(a, a1:a2, sep = "/", na.rm = TRUE) %>%
    unite(b, b1:b2, sep = "/", na.rm = TRUE) %>%
    unite(c, c1:c2, sep = "/", na.rm = TRUE) %>%
    unite(drb1, drb1.1:drb1.2, sep = "/", na.rm = TRUE) %>%
    unite(dqb1, dqb1.1:dqb1.2, sep = "/", na.rm = TRUE) %>%
    unite(drb345, drb3.1:drb5.2, sep = "/", na.rm = TRUE) %>%
    select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, lores_all_match)
  #* end of step 0 *#

  #* step 1: add count of low-res-match and filter by max count of occurrence *#
  tbl_raw <- tbl_raw %>%
    mutate(pair_id = tbl_in$pair_id,
           cnt_a = ifelse(fst_a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
           cnt_b = ifelse(fst_b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
           cnt_c = ifelse(fst_c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
           cnt_drb1 = ifelse(fst_drb1 %in% c(tbl_in$drb1.1, tbl_in$drb1.2), 1, 0),
           cnt_dqb1 = ifelse(fst_dqb1 %in% c(tbl_in$dqb1.1, tbl_in$dqb1.2), 1, 0),
           cnt_drb3 = ifelse(fst_drb345 %in% c(tbl_in$drb3.1, tbl_in$drb3.2), 1, 0),
           cnt_drb4 = ifelse(fst_drb345 %in% c(tbl_in$drb4.1, tbl_in$drb4.2), 1, 0),
           cnt_drb5 = ifelse(fst_drb345 %in% c(tbl_in$drb5.1, tbl_in$drb5.2), 1, 0),

           cnt_a_hi = ifelse(a %in% c(tbl_in$a1, tbl_in$a2), 1, 0),
           cnt_b_hi = ifelse(b %in% c(tbl_in$b1, tbl_in$b2), 1, 0),
           cnt_c_hi = ifelse(c %in% c(tbl_in$c1, tbl_in$c2), 1, 0),
           cnt_drb1_hi = ifelse(drb1 %in% c(tbl_in$drb1.1, tbl_in$drb1.2), 1, 0),
           cnt_dqb1_hi = ifelse(dqb1 %in% c(tbl_in$dqb1.1, tbl_in$dqb1.2), 1, 0),
           cnt_drb3_hi = ifelse(drb345 %in% c(tbl_in$drb3.1, tbl_in$drb3.2), 1, 0),
           cnt_drb4_hi = ifelse(drb345 %in% c(tbl_in$drb4.1, tbl_in$drb4.2), 1, 0),
           cnt_drb5_hi = ifelse(drb345 %in% c(tbl_in$drb5.1, tbl_in$drb5.2), 1, 0)) %>%
    mutate(cnt_drb345 = ifelse(cnt_drb3 + cnt_drb4 + cnt_drb5 > 0, 1, 0),
           cnt_drb345_hi = ifelse(cnt_drb3_hi + cnt_drb4_hi + cnt_drb5_hi > 0, 1, 0)) %>%
    mutate(cnt_lo = cnt_a + cnt_b + cnt_c + cnt_drb1 + cnt_dqb1 + cnt_drb345,
           cnt_hi = cnt_a_hi + cnt_b_hi + cnt_c_hi + cnt_drb1_hi + cnt_dqb1_hi + cnt_drb345_hi,
           cnt = cnt_lo + cnt_hi) %>%
    arrange(-cnt) %>%
    select(names(.)[!str_detect(names(.), "cnt_")]) %>%
    mutate(idx = as.numeric(rownames(.)))

  # filter data by max count of single a/b/c/drb/dqb combination
  # if max count has more than one pair, then filter by max_count, else filter by 2nd_max_count
  cnt_max <- tbl_raw %>% group_by(cnt) %>% tally() %>% filter(cnt == max(cnt))

  if(cnt_max$n > 1){
    tbl_raw <- tbl_raw %>% filter(cnt == cnt_max$cnt)
  } else{
    n = cnt_max$cnt - 1
    tbl_raw <- tbl_raw %>% filter(cnt >= n)
  }

  rm(cnt_max)
  #* end of step 1 *#

  #* step 2: some flags and counts of low-res between raw input data and ref data *#
  uniq_drb345 <- unique(c(tbl_in$drb3.1, tbl_in$drb3.2, tbl_in$drb4.1, tbl_in$drb4.2, tbl_in$drb5.1, tbl_in$drb5.2))
  uniq_drb345 <- sort(uniq_drb345,decreasing = TRUE)

  ori_hla <- tbl_in %>%
    select(-c(pair_id, ethnicity, subject_type, count)) %>%
    mutate(drb345.1 = uniq_drb345[1],
           drb345.2 = uniq_drb345[2]) %>%
    select(-c(drb3.1, drb3.2, drb4.1, drb4.2, drb5.1, drb5.2))

  ori_hla[ori_hla == ""] <- NA

  # count of unique hla across all loci in the input data
  tmp <- ori_hla %>%  mutate_if(.,
                                is.character,
                                str_replace_all,
                                pattern = ":.*",
                                replacement = "")
  cnt_uniq_ori <- 0
  for(i in seq(1, length(tmp), 2)){
    j <- i + 1
    cnt_uniq_ori <- cnt_uniq_ori + n_distinct(as.character(c(tmp[i], tmp[j])), na.rm = TRUE)
  }

  rm(tmp)

  # flag of whether or not input data has identical hlas on each locus across all loci
  # if input data has identical hla of each locus of all loci(flag == "yes"), then duplicate first impute record
  uniq_all_loci <- as.vector(ori_hla %>% as.vector())
  uniq_all_loci <- uniq_all_loci[!is.na(uniq_all_loci)]

  uniq_loci_flag <- ifelse(length(unique(uniq_all_loci))*2 == length(uniq_all_loci),
                           "yes", "no")

  ori_hla2 <- ori_hla %>% select_if(~sum(!is.na(.)) > 0)

  nms <- ori_hla %>% select_if(~sum(!is.na(.)) > 0)
  nms <- as.vector(names(nms))
  nms4match <- case_when(nms %in% c("a1", "a2") ~ "a",
                         nms %in% c("b1", "b2") ~ "b",
                         nms %in% c("c1", "c2") ~ "c",
                         nms %in% c("drb1.1", "drb1.2") ~ "drb1",
                         nms %in% c("dqb1.1", "dqb1.2") ~ "dqb1",
                         nms %in% c("drb345.1", "drb345.2") ~ "drb345",
                         TRUE~ nms) %>%
    unique()
  #* end of step 2 *#

  #* step 3: keep one combination of each homozygous alleles *#
  uniq1 <- tbl_raw %>%
    select(all_of(nms4match)) %>%
    mutate_if(.,
              is.character,
              str_replace_all,
              pattern = ":.*",
              replacement = "") %>%
    distinct() %>%
    unite(flag, all_of(nms4match)) %>%
    mutate(grp = paste0("group", row.names(.)))

  uniq2 <- tbl_raw %>%
    select(all_of(nms4match)) %>%
    mutate_if(.,
              is.character,
              str_replace_all,
              pattern = ":.*",
              replacement = "") %>%
    unite(flag, all_of(nms4match)) %>%
    mutate(rnum = row.names(.)) %>%
    left_join(.,  uniq1, by = "flag")

  if(length(unique(uniq2$grp)) > 1){
    uniq3 <- tbl_raw %>%
      mutate(rnum = row.names(.)) %>%
      left_join(., uniq2, by = "rnum") %>%
      group_by(grp) %>%
      slice(which.max(freq)) %>%
      ungroup()
  } else{
    uniq3 <- tbl_raw %>%
      mutate(rnum = row.names(.)) %>%
      left_join(., uniq2, by = "rnum")
  }
  #* end of step 3 *#

  #* step 4: search for a pair with max number of low-res match, lowest rank of each combination *#

  # pairs0: working table
  pairs0 <- uniq3 %>%
    select(a, b, c, drb1, dqb1, drb345, freq, rank, drb) %>%
    arrange(rank) %>%
    mutate(idx = as.numeric(rownames(.))) #reset index

  # pairs1: transpose single a/b/c/drb/dqb to all possible 2-pairs' combination
  pairs1 <- data.frame(t(combn(pairs0$idx, 2))) %>%
    setNames(c("indx1", "indx2")) %>%
    mutate(pair = row_number()) %>%
    pivot_longer(cols = c("indx1", "indx2"), names_to = "index") %>%
    left_join(., pairs0, by = c("value" = "idx")) %>%
    select(-c(index, value))

  # cnt_match_by_row #
  cnt_match_by_row <- pairs1 %>%
    select(c(pair, all_of(nms4match))) %>%
    mutate_if(.,
              is.character,
              str_replace_all,
              pattern = ":.*",
              replacement = "")

  ori_hla2 <- ori_hla %>% select_if(~sum(!is.na(.)) > 0) %>%
    mutate_if(.,
              is.character,
              str_replace_all,
              pattern = ":.*",
              replacement = "")

  cnt_match_by_row <- cbind(cnt_match_by_row, ori_hla2, row.names = NULL) %>% data.frame() %>% select(sort(names(.))) %>% select(pair, everything())

  cnt_match_by_row$cnt_match_by_row <- "no"

  c_st <- 2
  c_ed <- dim(cnt_match_by_row)[2] - 1

  for(i in 1:dim(cnt_match_by_row)[1]){
    cnt_tmp <- 0
    for(j in seq(c_st, c_ed, 3)){
      plus1 <- j +1
      plus2 <- j + 2
      cnt_tmp <- ifelse(cnt_match_by_row[i,j] %in% cnt_match_by_row[i, plus1:plus2] == TRUE , cnt_tmp + 1, cnt_tmp)
    }
    cnt_match_by_row[i,]$cnt_match_by_row <- cnt_tmp
  }
  cnt_match_by_row <- cnt_match_by_row %>% mutate(rnum = rownames(.)) %>% select(rnum, cnt_match_by_row) # rnum is used to match pair1 record later
  rm(nms, nms4match, ori_hla2, cnt_tmp, c_st, c_ed, plus1, plus2)

  # pairs2: add unique number of low-res, keep those with max count of unique low-res
  pairs2 <- pairs1 %>% select(pair, a, b, c, drb1, dqb1, drb345) %>%
    mutate(a = gsub(":.*", "", a),
           b = gsub(":.*", "", b),
           c = gsub(":.*", "", c),
           drb1 = gsub(":.*", "", drb1),
           dqb1 = gsub(":.*", "", dqb1),
           drb345 = gsub(":.*", "", drb345)) %>%
    group_by(pair) %>%
    mutate(rn = row_number()) %>%
    pivot_wider(names_from = rn,
                values_from = c("a", "b", "c", "drb1", "dqb1", "drb345")) %>%
    data.frame()

  pairs2 <- cbind(pairs2, ori_hla, row.names = NULL) %>% data.frame() %>% select(sort(names(.))) %>% select(pair, everything()) %>%
    mutate_if(.,
              is.character,
              str_replace_all,
              pattern = ":.*",
              replacement = "")

  # starting and ending index of new columns
  pairs3 <- pairs2 %>%
    mutate(cnt_match_by_pair = 0)

  n <- dim(pairs2)[1]
  c_st <- 2
  c_ed <- dim(pairs2)[2]
  c_lst <- c_ed + 1

  for (i in 1:n){ # loop through all of pairs
    tmp <- 0
    for (j in seq(c_st, c_ed, 4)){
      plus1 <- j + 1
      plus2 <- j + 2
      plus3 <- j + 3

      # if(all(all_of(unique(as.character(pairs3[i, plus2:plus3][!is.na(pairs3[i, plus2:plus3])])) %in% pairs3[i, j:plus1]) %in% "TRUE")){
      tmp_vec <- unique(as.character(pairs3[i, plus2:plus3][!is.na(pairs3[i, plus2:plus3])]))
      # if there 2 unique hlas
      if(length(tmp_vec) == 2){
        if(tmp_vec[1] %in% pairs3[i, j:plus1] & tmp_vec[2] %in% pairs3[i, j:plus1]){
          tmp <- tmp + 2
        } else if(tmp_vec[1] %in% pairs3[i, j:plus1] | tmp_vec[2] %in% pairs3[i, j:plus1]){
          tmp <- tmp + 1
        } else{
          tmp <- tmp
        }
      }
      # if only 1 hla
      else{
        if(tmp_vec[1] %in% pairs3[i, j:plus1]){
          tmp <- tmp + 1
        } else{
          tmp <- tmp
        }
      }
    }
    pairs3[i,c_lst] <- tmp
  }

  pairs4 <- pairs3 %>%
    mutate(cnt_uniq_ori = cnt_uniq_ori,
           uniq_loci_flag = uniq_loci_flag) %>%
    select(pair, cnt_uniq_ori, cnt_match_by_pair, uniq_loci_flag)

  pairs5 <- pairs1 %>%
    mutate(rnum = rownames(.)) %>% # rnum is used to match cnt_match_by_row
    left_join(., pairs4, by = "pair") %>%
    filter(cnt_match_by_pair == max(cnt_match_by_pair)) %>%
    left_join(., cnt_match_by_row, by = c("rnum")) %>%
    mutate(cnt_match_by_row = as.numeric(cnt_match_by_row)) %>%
    filter(pair == min(pair)) %>%
    arrange(-cnt_match_by_row)

  # if input hla are identical of each locus, then duplicate the combination
  if(pairs5[2,]$uniq_loci_flag == "yes" & (pairs5[2,]$cnt_match_by_row < pairs5[1,]$cnt_uniq_ori)){
    pairs5[2,] = pairs5[1,]
  }

  result <- pairs5 %>%
    filter(pair == min(pair)) %>%
    mutate(subj = ifelse(tbl_in$ethnicity == "",
                         paste(tbl_in$pair_id, tbl_in$subject_type, sep = "_"),
                         paste(paste(tbl_in$pair_id, tbl_in$subject_type, sep = "_"), tbl_in$ethnicity, sep = "_")),
           dat_type = "imputed",
           drb345 = paste(drb, drb345, sep = "*"),
           pair_id = tbl_in$pair_id,
           lores_all_match = ifelse(cnt_uniq_ori == cnt_match_by_pair, "yes", "no")) %>%
    select(subj, dat_type, pair_id, a, b, c, drb1, dqb1, drb345, freq, rank, lores_all_match)

  result <- rbind(raw, result)

  return(result)
}

#' @rdname utils
na_to_empty_string <- function(df) {
  # thanks to Hadley Wickham
  mutate(df, across(tidyselect::where(is.character), ~ na_if(.x, "")))
}
