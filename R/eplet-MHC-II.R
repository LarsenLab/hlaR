#' Title
#' function to calculate eplet mis-match of HLA-DRDQDP
#' @param dat_in
#' input data set with or without complete allele info
#' it has 27 columns (first 3 columns are record id, recipient id, donor id; the rest of columns are MHC II alleles for recipient and donor)
#' @return
#' list of data tables.
#' count: original input data appended with count mis-matched eplet of each pair
#' detail: detailed mis-match eplet of each subject id, plus count and percentage of mis-match across all pairs
#' @export
#'
#' @import
#' devtools
#' vroom
#' tidyverse
#' dplyr
#' schoolmath
#' tibble
#' stringr
#' purrr
#' tidyr
#'
#' @examples
#' \dontrun{
#' re <- CalEpletMHCII(dat_in = "YourDataFile")
#' re <- CalEpletMHCII(system.file(“extdata”, “MHC_II_test.csv”, package = “hlaR”))
#' }
#'
# notes:
# 1. raw_eplet is pulled from worksheets "A" and "B" from DRDQDP_Antibody_Analysis_3.1.xlsb (password of protected sheets: hla) from http://www.epitopes.net/index.html
# 2. raw_lookup is generated on-the-fly based on raw_eplet table
# 3. result generated from this function may slightly different with result from DRDQDP_Antibody_Analysis_3.1.xlsb due to a few possible formatting issue in the excel
# 4. check result: re$detail or re$count
# below is an example to use globalVariables() to suppress "no visible global variable" note
# utils::globalVariables(c("value", "locus", "index", "type", "mm"))

CalEpletMHCII <- function(dat_in) {
  ###*** step 1: import raw eplet tables ***###
  raw_eplet_A <- vroom(system.file("extdata", "MHC_II_eplet_A.csv", package = "hlaR"))

  raw_eplet_B <- vroom(system.file("extdata", "MHC_II_eplet_B.csv", package = "hlaR"))
  ###*** end of step 1 ***###

  ###*** step 2: generate base lookup tables for MHC II loci As and Bs **###
  #* 2a: master lookup tables for As and Bs *#
  lkup_a <- as.data.frame(t(raw_eplet_A)) %>%
            setNames(paste(raw_eplet_A$type, raw_eplet_A$index, sep = "_" )) %>%
            rownames_to_column(var = "locus") %>%
            mutate(locus = ifelse(str_detect(locus, "\\*"), sub("\\*.*", "", locus), locus)) %>%
            filter(!locus %in% c("index", "type") ) %>%
            distinct() %>%
            reshape2::melt(id.vars = "locus") %>%
            filter(value != "" ) %>%
            distinct()

  lkup_b <- as.data.frame(t(raw_eplet_B)) %>%
            setNames(paste(raw_eplet_B$type, raw_eplet_B$index, sep = "_" )) %>%
            rownames_to_column(var = "locus") %>%
            mutate(locus = ifelse(str_detect(locus, "\\*"), sub("\\*.*", "", locus), locus)) %>%
            filter(!locus %in% c("index", "type") ) %>%
            distinct() %>%
            reshape2::melt(id.vars = "locus") %>%
            filter(value != "" ) %>%
            distinct()
  #* end of 2a *#

  # 2b: sub-lookup tables, these tables are used for mis-match comparison of each locus #
  GenerateLookup <- function(dat_in, locus_in){
    dat_out <- dat_in %>%
                filter(locus %in% locus_in) %>%
                mutate(index = as.numeric(sub(".*\\_", "", variable)),
                       type = sub("\\_.*", "", variable)) %>%
                select(index, type, value) %>%
                dplyr::rename(eplet = value) %>%
                distinct()

    return(dat_out)
  }

  lkup_dpa <- GenerateLookup(lkup_a, "DPA1")
  lkup_dqa <- GenerateLookup(lkup_a, "DQA1")

  lkup_dpb <- GenerateLookup(lkup_b, "DPB1")
  lkup_dqb <- GenerateLookup(lkup_b, "DQB1")
  lkup_drb <- GenerateLookup(lkup_b, c("DRB1", "DRB3", "DRB4", "DRB5"))

  # compare the result with original excel table, need confirm 4 eplets below
  # DQA
  # DQA1 oth_10 161DI - it's 161D in the excel
  # DQA1 oth_5 75I - it's 75 in the excel
  # DQA1 oth_10 161DI - it's 160D in the excel
  # DPA
  # DPA1 other_1 pq34Q - it's 34Qx in the excel

  # some known issues cause incorrect mis-match result in the excel file
  # 1. extra space: "108M "
  # 2. missing letters:
  # 3. mis-spell of eplet name: "rq70RK/R" supposed to be "rq70RR/K"
  # 4. mis-match column number: "185I" col"FU4" supposed to be "col "FT4"

  #* end of 2b *#
  ###*** end of step 2 ***###

  ###*** step 3: import records with recipient/donor MHC II pairs ***###
  #* 3a. import data *#
  # input data format has to be same as test data
  dat <- vroom(here(dat_in)) %>%
          mutate(across(where(is.character), str_trim))

  subj_num <- dim(dat)[1]
  tmp_names <- names(dat[-c(1:3)])

  ep_a <- raw_eplet_A %>% select(index, type)
  ep_b <- raw_eplet_B %>% select(index, type)
  #* end of 3a*#

  #* 3b: pull out eplet of each locus *#
  for (i in 1:subj_num) {
    allele <- unlist(transpose(dat[i,-c(1:3)]), use.names = F)

    for (j in 1:length(allele)) {
      varname <- paste0(tmp_names[j], ".", sep = i)
      # A loci
      if(j %in% c(7, 8, 11, 12, 19, 20, 23, 24)){
        if (!is.na(allele[j])) {
          tmp <- raw_eplet_A %>%
                 select(index, type, allele[j])
          tmp <- tmp %>%
                 setNames(c("index", "type", varname))
          ep_a <- ep_a %>%
                  left_join(., tmp, by = c("index", "type"))
        } else{
          ep_a <- ep_a %>%
                  mutate(!!varname := NA)
        }
      }
      # B loci
      else{
        if (!is.na(allele[j])) {
          tmp <- raw_eplet_B %>%
                 select(index, type, allele[j])
          tmp <- tmp %>%
                 setNames(c("index", "type", varname))
          ep_b <- ep_b %>%
                  left_join(., tmp, by = c("index", "type"))
        } else{
          ep_b <- ep_b %>%
                  mutate(!!varname := NA)
        }
      }
    }
  }
  #* end of 3b *#
  #*** end of step 3 ***###

  ###*** step 4: mis-match eplet calculation ***###
  #* logics
  #* if locus DQB: compare donor's 1 allele to recipient's 2 DDRDRQB alleles (1:2)
  #* if locus DQA: compare donor's 1 allele to recipient's 2 DQA alleles (1:2)
  #* if locus DRB and DRw: compare donor's 1 allele to recipient's DRB+DRw alleles (1:4)
  #* if locus DPB
  #*             if inter-loci eplets (av >= 22 or other <= 8) then compare donor's 1 allele to recipient's ALL B alleles (1:8)
  #*             else compare donor's 1 allele recipient's DPB alleles (1:2)
  #* if locus DPB: compare donor's 1 allele to recipient's ALL B alleles (1:8)
  #* if locus DPA: compare donor's 1 allele to recipient's 2 DPA alleles (1:2)

  #* 4a: A loci *#
  ep_a_mm <- ep_a %>%
             select(index, type)

  # A alleles : DQA, DPA
  a_len <- subj_num * 8
  # exclude index and type, pulling data starting from the 3rd position
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 7
    pos <- c(st:ed)
    tmp <- ep_a %>%
           select(pos)
    subj_indx <- sub(".*\\.", "", names(tmp)[i])

    colnames(tmp) <- c("rec_dqa1", "rec_dqa2", "rec_dpa1", "rec_dpa2",
                       "don_dqa1", "don_dqa2", "don_dpa1", "don_dpa2")

    tmp <- tmp %>%
            mutate(dqa1_mm = ifelse(don_dqa1 %in% c(rec_dqa1, rec_dqa2), NA, don_dqa1),
                   dqa2_mm = ifelse(don_dqa2 %in% c(rec_dqa1, rec_dqa2), NA, don_dqa2),
                   dpa1_mm = ifelse(don_dpa1 %in% c(rec_dpa1, rec_dpa2), NA, don_dpa1),
                   dpa2_mm = ifelse(don_dpa2 %in% c(rec_dpa1, rec_dpa2), NA, don_dpa2)) %>%
            select(dqa1_mm, dqa2_mm, dpa1_mm, dpa2_mm) %>%
            setNames(c(paste0("dqa1_mm_subj", subj_indx),
                       paste0("dqa2_mm_subj", subj_indx),
                       paste0("dpa1_mm_subj", subj_indx),
                       paste0("dpa2_mm_subj", subj_indx)))

    ep_a_mm <- cbind(ep_a_mm, tmp)
    st <- ed + 1
  }
  #* end of 4a *#

  #* 4b. B loci *#
  ep_b_mm <- ep_b %>%
             select(index, type)
  ep_b_rownames <- paste(ep_b_mm$type, ep_b_mm$index, sep = "")

  rownames(ep_b_mm) <- ep_b_rownames
  rownames(ep_b) <- ep_b_rownames
  # if av <= 21 or oth >= 8, DPB, one don : 2 rec comparison
  # if interlocus eplets(names contains "Int"), DPB, one don : ALL rec comparison
  int_loc_eplet<- ep_b_rownames[str_detect(ep_b_rownames,"Int")]

  # B loci
  # B alleles : DRB, DRw, DQB, DPB
  b_len <- subj_num * 16
  # exclude index and type, pulling data starting from the 3rd position
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 15
    pos <- c(st:ed)
    tmp <- ep_b %>%
           select(pos)

    subj_indx <- sub(".*\\.", "", names(tmp)[i])

    colnames(tmp) <- c("rec_drb1", "rec_drb2", "rec_drw1", "rec_drw2", "rec_dqb1", "rec_dqb2", "rec_dpb1","rec_dpb2",
                       "don_drb1", "don_drb2", "don_drw1", "don_drw2", "don_dqb1", "don_dqb2", "don_dpb1","don_dpb2")
    tmp$rname <- ep_b_rownames

    tmp <- tmp %>%
            mutate(drb1_mm = ifelse(don_drb1 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drb1),
                   drb2_mm = ifelse(don_drb2 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drb2),
                   drw1_mm = ifelse(don_drw1 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drw1),
                   drw2_mm = ifelse(don_drw2 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drw2),
                   dqb1_mm = ifelse(don_dqb1 %in% c(rec_dqb1, rec_dqb2), NA, don_dqb1),
                   dqb2_mm = ifelse(don_dqb2 %in% c(rec_dqb1, rec_dqb2), NA, don_dqb2),
                   dpb1_tmp2 = ifelse(don_dpb1 %in% c(rec_dpb1, rec_dpb2), NA, don_dpb1),
                   dpb1_tmp8 = ifelse(don_dpb1 %in% c(rec_dpb1, rec_dpb2, rec_drb1, rec_drb2, rec_drw1, rec_drw2, rec_dqb1, rec_dqb2), NA, don_dpb1),
                   dpb2_tmp2 = ifelse(don_dpb2 %in% c(rec_dpb1, rec_dpb2), NA, don_dpb2),
                   dpb2_tmp8 = ifelse(don_dpb2 %in% c(rec_dpb1, rec_dpb2, rec_drb1, rec_drb2, rec_drw1, rec_drw2, rec_dqb1, rec_dqb2), NA, don_dpb2)) %>%
            mutate(dpb1_mm =  ifelse(rname %in% int_loc_eplet , dpb1_tmp8, dpb1_tmp2),
                   dpb2_mm = ifelse(rname %in% int_loc_eplet, dpb2_tmp8, dpb2_tmp2)) %>%
            select(drb1_mm, drb2_mm, drw1_mm, drw2_mm, dqb1_mm, dqb2_mm, dpb1_mm,   dpb2_mm ) %>%
            setNames(c(paste0("drb1_mm_subj", subj_indx),
                       paste0("drb2_mm_subj", subj_indx),
                       paste0("drw1_mm_subj", subj_indx),
                       paste0("drw2_mm_subj", subj_indx),
                       paste0("dqb1_mm_subj", subj_indx),
                       paste0("dqb2_mm_subj", subj_indx),
                       paste0("dpb1_mm_subj", subj_indx),
                       paste0("dpb2_mm_subj", subj_indx)))

    ep_b_mm <- cbind(ep_b_mm, tmp)
    st <- ed + 1
  }
  #* end of 4b *#
  ###*** end of step 4 ***###

  ###*** step 5: mis-match calculation of A loci ***###
  #* 5a. generate mm tables *#
  subj_names <- unique(sub(".*\\_", "", names(ep_a_mm)[-c(1:2)]))

  re_dqa <- lkup_dqa
  re_dpa <- lkup_dpa
  st <- 3

  for (i in 1:subj_num) {
    ed <- st + 3
    positions <- c(st:ed)
    tmp <- ep_a_mm %>%
            select(c(1, 2, positions)) %>%
            setNames(c("index", "type", "dqa1_mm", "dqa2_mm", "dpa1_mm", "dpa2_mm"))

    # dqa, join lkup_dqa
    lkup <- lkup_dqa
    tmp_dqa <- tmp %>%
                select(c(index, type, dqa1_mm, dqa2_mm)) %>%
                right_join(., lkup, by = c("index", "type")) %>%
                mutate(var = ifelse(eplet %in% c(dqa1_mm, dqa2_mm), eplet, NA)) %>%
                select(index, type, eplet, var) %>%
                setNames(c("index", "type", "eplet", subj_names[i]))

    # dpa, join lkup_dpa
    lkup <- lkup_dpa
    tmp_dpa <- tmp %>%
                select(c(index, type, dpa1_mm, dpa2_mm)) %>%
                right_join(., lkup, by = c("index", "type")) %>%
                mutate(var = ifelse(eplet %in% c(dpa1_mm, dpa1_mm), eplet, NA)) %>%
                select(index, type, eplet, var) %>%
                setNames(c("index", "type", "eplet", subj_names[i]))

    st <- ed + 1
    re_dqa <- left_join(re_dqa, tmp_dqa, by = c("index", "type", "eplet"))
    re_dpa <- left_join(re_dpa, tmp_dpa, by = c("index", "type", "eplet"))
  }
  #* end of 5a *#

  #* 5b. generate result tables *#
  GenerateResult <- function(dat_in){
    dtl <- dat_in %>%
           mutate(mm_cnt = subj_num - rowSums(is.na(.)),
           mm_pect = paste(round((subj_num - rowSums(is.na(.))) / subj_num * 100, 1), "%", sep = ""))

    cnt <- dtl %>%
            select(-c(index, type, eplet, mm_cnt, mm_pect)) %>%
            summarise_all(funs(sum(!is.na(.)))) %>%
            gather() %>%
            setNames(c("subject", "mm_count")) %>%
            select(mm_count)
    return(list(detail = dtl, count = cnt))
  }

  dpa <- GenerateResult(re_dpa)
  dqa <- GenerateResult(re_dqa)

  #* end of 5b *#
  ###*** end of step 5 ***###

  ###* step 6: mis-match calculation of B loci ***###

  #* 6a. generate mm tables *#
  re_dpb <- lkup_dpb
  re_dqb <- lkup_dqb
  re_drb <- lkup_drb

  st <- 3

  for (i in 1:subj_num) {
    ed <- st + 7
    positions <- c(st:ed)
    tmp <- ep_b_mm %>%
            select(c(1, 2, positions)) %>%
            setNames(c("index", "type", "drb1_mm", "drb2_mm", "drw1_mm", "drw2_mm", "dqb1_mm", "dqb2_mm", "dpb1_mm", "dpb2_mm"))

    # dpb
    lkup <- lkup_dpb
    cols <- c("index", "type", "dpb1_mm", "dpb2_mm")
    tmp_dpb <- tmp %>%
                select(cols) %>%
                right_join(., lkup, by = c("index", "type")) %>%
                mutate(var = ifelse(eplet %in% c(dpb1_mm, dpb2_mm), eplet, NA)) %>%
                select(index, type, eplet, var) %>%
                setNames(c("index", "type", "eplet", subj_names[i]))

    # dqb
    lkup <- lkup_dqb
    cols <- c("index", "type", "dqb1_mm", "dqb2_mm")
    tmp_dqb <- tmp %>%
                select(cols) %>%
                right_join(., lkup, by = c("index", "type")) %>%
                mutate(var = ifelse(eplet %in% c(dqb1_mm, dqb2_mm), eplet, NA)) %>%
                select(index, type, eplet, var) %>%
                setNames(c("index", "type", "eplet", subj_names[i]))

    # drb
    lkup <- lkup_drb
    cols <- c("index", "type", "drb1_mm", "drb2_mm", "drw1_mm", "drw2_mm")
    tmp_drb <- tmp %>%
                select(cols) %>%
                right_join(., lkup, by = c("index", "type")) %>%
                mutate(var = ifelse(eplet %in% c(drb1_mm, drb2_mm, drw1_mm, drw2_mm), eplet, NA)) %>%
                select(index, type, eplet, var) %>%
                setNames(c("index", "type", "eplet", subj_names[i]))

    st <- ed + 1

    re_dpb <- left_join(re_dpb, tmp_dpb, by = c("index", "type", "eplet"))
    re_dqb <- left_join(re_dqb, tmp_dqb, by = c("index", "type", "eplet"))
    re_drb <- left_join(re_drb, tmp_drb, by = c("index", "type", "eplet"))
  }
  #* end of 76a *#

  #* 6b. generate result tables *#
  #*
  dpb <- GenerateResult(re_dpb)
  dqb <- GenerateResult(re_dqb)
  drb <- GenerateResult(re_drb)

  #* end of 6b *#
  ###* end of step 6 ***###

  ###*** step 7: final mm table ***#
  dat_result <- dat %>%
    mutate(mm_drb = drb$count,
           mm_dqb = dqb$count,
           mm_dqa = dqa$count,
           mm_dpb = dpb$count,
           mm_dpa = dpa$count)

  dat_out <- list(count = dat_result,
                  detail = list(drb = drb$detail,
                                    dqb = dqb$detail,
                                    dqa = dqa$detail,
                                    dpb = dpb$detail,
                                    dpa = dpa$detail)
  )

  return(dat_out)
  ###*** end of step 7 ***###
}











