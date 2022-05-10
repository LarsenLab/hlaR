#' @name CalEpletMHCI
#' @title Calculate class I HLA eplet mismatch
#' @description Use high resolution HLA(Human Leukocyte Antigen) class I data to calculate class I eplet mismatch for a population of donors and recipients. Mismatch is calculated using logic from 'HLAMatchMaker', developed by Rene Dusquesnoy. Current reference tables supported are 'HLAMatchMaker' v2 and v3.
#' @param dat_in
#' A dataframe of recipient and donor's high resolution MHC I data. Each recipient and donor pair are linked by are the “pair_id” column and differentiated by the “subject_type” column.
#' @param ver
#' Version number of HLAMatchMaker based eplet reference table to use.
#' @return
#' A list of data tables.
#' - `single_detail`: single molecule class I MHC eplet mismatch table, including mismatched eplet names and the count of eplets mismatched at each allele.
#' - `overall_count`: original input data appended with total count of mismatched eplets.
#' @export
#'
#' @import
#' devtools
#' tidyverse
#' dplyr
#' schoolmath
#' tibble
#' stringr
#' purrr
#' tidyr
#' @importFrom
#' stats setNames
#' @examples
#' dat<-read.csv(system.file("extdata/example","MHC_I_test.csv",package="hlaR"),sep=",",header=TRUE)
#' re <- CalEpletMHCI(dat_in = dat, ver = 3)

CalEpletMHCI <- function(dat_in, ver = 2) {
  #* step 0: check if recipient and donor are paired *#
  num_rcpt <- length(dat_in$subject_type[dat_in$subject_type %in% c("recipient", "recip", "rcpt", "r")])
  num_don  <- length(dat_in$subject_type[dat_in$subject_type %in% c("donor", "don", "dn", "d")])
  if(num_rcpt == num_don){
    rm(num_rcpt, num_don)
  } else{
    stop("please check that every pair_id has both recipient and donor data.")
  }
  #* end of step 0 *#

  #* step 1: import eplet reference table *#
  if(ver == 2){
    tbl_raw_eplet <- read.csv(system.file("extdata/ref", "MHC_I_eplet_v2.csv", package = "hlaR"), check.names = FALSE)
  } else{
    tbl_raw_eplet <- read.csv(system.file("extdata/ref", "MHC_I_eplet_v3.csv", package = "hlaR"), check.names = FALSE)
  }

  tbl_ref <- as.data.frame(t(tbl_raw_eplet)) %>%
    setNames(paste(tbl_raw_eplet$type, tbl_raw_eplet$index, sep = "_" )) %>%
    rownames_to_column(var = "locus") %>%
    mutate(locus = ifelse(str_detect(locus, "\\*"), sub("\\*.*", "", locus), locus)) %>%
    filter(!locus %in% c("index", "type") ) %>%
    distinct() %>%
    reshape2::melt(id.vars = "locus") %>%
    filter(value != "" ) %>%
    distinct() %>%
    mutate(index = as.numeric(sub(".*\\_", "", variable)),
           type = sub("\\_.*", "", variable)) %>%
    dplyr::rename(eplet = value) %>%
    select(index, type, eplet) %>%
    distinct()
  #* end of step 1 *#

  #* step 2: subject table *#
  nm_rec <- c("rec_a1", "rec_a2", "rec_b1", "rec_b2", "rec_c1", "rec_c2")
  nm_don <- c("don_a1", "don_a2", "don_b1", "don_b2", "don_c1", "don_c2")

  tmp_rcpt <- dat_in %>%
    filter(subject_type %in% c("recipient", "recip", "rcpt", "r")) %>%
    select(-subject_type) %>%
    setNames(c("pair_id", nm_rec))

  all_nas_rcpt <- tmp_rcpt %>%
    filter(rec_a1=="" &  rec_a2=="" & rec_b1=="" & rec_b2=="" & rec_c1=="" &  rec_c2=="" ) %>% pull(pair_id)

  if(length(all_nas_rcpt) != 0) {
    warning(paste0("no MHC class I alleles detected for recipient(s) ", toString(all_nas_rcpt), " !"))
    cat("\n")
  }

  tmp_don <- dat_in %>%
    filter(subject_type %in% c("donor", "don", "dn", "d")) %>%
    select(-subject_type) %>%
    setNames(c("pair_id", nm_don))

  all_nas_don <- tmp_don %>%
    filter(don_a1=="" &  don_a2=="" & don_b1=="" & don_b2=="" & don_c1=="" &  don_c2=="" ) %>% pull(pair_id)

  if(length(all_nas_don) != 0) {
    warning(paste0("no MHC class I alleles detected for donor(s) ", toString(all_nas_don), " !"))
    cat("\n")
  }

  tbl_ready <- left_join(tmp_rcpt, tmp_don, by = c("pair_id")) %>%
    mutate(pair_id_ori = pair_id) %>%
    arrange(pair_id_ori) %>%
    mutate(pair_id = dense_rank(pair_id_ori))

  rm(tmp_rcpt, tmp_don, all_nas_rcpt)

  # create id_match table in case pari_id is not sequential in the patient table
  id_match <- tbl_ready %>%
    select(pair_id_ori, pair_id) %>%
    mutate(pair_id = as.character(pair_id))

  tbl_ready <- tbl_ready %>%
    select(-pair_id_ori) %>%
    select(pair_id, everything())

  subj_num <- dim(tbl_ready)[1]
  tmp_names <- c(nm_rec, nm_don)
  #* end of step 2 *#

  #* step 3: initialize a dataframe to hold eplets *#
  tbl_ep <- tbl_raw_eplet %>%
    select(index, type)
  #* end of step 3 *#

  #* step 4: pull out eplet of each allele *#
  for (i in 1:subj_num) {
    allele <- toupper(unlist(transpose(tbl_ready[i,-c(1)]), use.names = F))

    # warning if input alleles are not found in Matchmaker
    miss_allele <- allele[!(allele %in% names(tbl_raw_eplet))]
    miss_allele <- miss_allele[miss_allele != ""]

    if(length(miss_allele) >= 1){
      warning(paste0(miss_allele, " is not found in the refernce table. Please check!\n"))
    }

    # if allele is blank, then set it to NA
    allele <- ifelse(allele %in% c(names(tbl_raw_eplet), miss_allele),
                     allele, NA)

    for (j in 1:length(allele)) {
      varname <- paste0(tmp_names[j], ".", sep = i)

      if (!is.na(allele[j]) & !(allele[j] %in% miss_allele)) {
        tmp <- tbl_raw_eplet %>%
          select(index, type, allele[j])
        tmp <- tmp %>%
          setNames(c("index", "type", varname))
        tbl_ep <- tbl_ep %>%
          left_join(., tmp, by = c("index", "type"))
      } else if (!is.na(allele[j]) & allele[j] %in% miss_allele){
        tbl_ep <- tbl_ep %>%
          mutate(!!varname := "Not Found")
      } else{
        tbl_ep <- tbl_ep %>%
          mutate(!!varname := NA)
      }
    }
  }
  #* end of step 4 *#

  #* step 5: mark mis-matches *#
  tbl_ep_mm <- tbl_ep %>%
    select(index, type)

  tbl_ep_mm2 <- tbl_ep_mm

  #  rcpt_list <- c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2)

  # exclude index and type, pulling data for each allele, starting from 3rd position
  st <- 3

  # for each subject, compare eplets of donor's EACH allele for ALL of recipients'
  for (i in 1:subj_num) {
    ed <- st + 11
    pos <- c(st:ed)
    tmp <- tbl_ep %>%
      select(all_of(pos))

    subj_indx <- sub(".*\\.", "", names(tmp)[1])

    colnames(tmp) <- c(nm_rec, nm_don)

    # important : replace recipient's "Not Found" to NA to avoid incorrect mismatch calculation in later step
    # if there are missing alleles on the same locus in both recpt and donor, this way we can make sure the donor's mm_eplet is "Not Found" instead of blank
    nm <- names(tmp)[str_detect(names(tmp), "rec_")]
    tmp <- tmp %>%
      mutate_at(vars(nm),
                list(~ifelse(. == "Not Found", "", .)))

    # comparison
    tmp <- tmp %>%
      mutate(a1_mm = ifelse(don_a1 %in% miss_allele,
                            "Not Found",
                            ifelse(don_a1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_a1)),
             a2_mm = ifelse(don_a2 %in% miss_allele,
                            "Not Found",
                            ifelse(don_a2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_a2)),
             b1_mm = ifelse(don_b1 %in% miss_allele,
                            "Not Found",
                            ifelse(don_b1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_b1)),
             b2_mm = ifelse(don_b2 %in% miss_allele,
                            "Not Found",
                            ifelse(don_b2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_b2)),
             c1_mm = ifelse(don_c1 %in% miss_allele,
                            "Not Found",
                            ifelse(don_c1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_c1)),
             c2_mm = ifelse(don_c2 %in% miss_allele,
                            "Not Found",
                            ifelse(don_c2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_c2))) %>%
      select(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm) %>%
      setNames(c(paste0("a1_mm_subj", subj_indx),
                 paste0("a2_mm_subj", subj_indx),
                 paste0("b1_mm_subj", subj_indx),
                 paste0("b2_mm_subj", subj_indx),
                 paste0("c1_mm_subj", subj_indx),
                 paste0("c2_mm_subj", subj_indx)))

    tbl_ep_mm <- cbind(tbl_ep_mm, tmp)
    st <- ed + 1
  }
  rm(st)

  # if all of recipient alleles are NA, then the result eplets are blank, need change change to NA for mismatch count in later step
  tbl_ep_mm <- tbl_ep_mm %>%
    mutate_all(na_if, "")
  #* end of step 5 *#

  #* step 6: compare mis-match with tbl_ref table *#
  # reset starting position for another round of loop
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 5
    pos <- c(st:ed)
    tmp <- tbl_ep_mm %>%
      select(c(1, 2, all_of(pos)))

    ori_name <- colnames(tmp)[-c(1,2)]
    tmp2 <- left_join(tbl_ref, tmp, by = c("index", "type")) %>%
      setNames(c("index", "type", "eplet", "a1", "a2", "b1", "b2", "c1", "c2")) %>%
      mutate(mm = ifelse(eplet == a1 | eplet == a2 | eplet == b1 | eplet == b2 | eplet == c1 | eplet == c2 , eplet, "")) %>%
      filter(mm != "") %>%
      left_join(tmp, ., by = c("index", "type")) %>%
      select(-c(all_of(ori_name), "eplet", "mm")) %>%
      setNames(c("index", "type", ori_name)) %>%
      distinct()

    st <- ed + 1

    tbl_ep_mm2 <- tbl_ep_mm2 %>%
      left_join(., tmp2, by = c("index", "type"))
  }
  rm(st)
  #* end of step 6 *#

  #* step 7: final result - single molecule report *#
  # subject names
  subj_names <- unique(sub(".*\\_", "", names(tbl_ep_mm2)[-c(1:2)]))

  re_s <-  data.frame(t(tbl_ep_mm2)) %>%
    unite("mm_eplets", names(.), na.rm = TRUE, sep = ",", remove = FALSE) %>%
    mutate(subject = gsub(".*_", "", rownames(.)),
           pair_id = gsub(".*subj", "", subject),
           gene = gsub("_.*", "", rownames(.)),
           mm_cnt = str_count(mm_eplets, ",")) %>%
    select(subject, pair_id, gene, mm_eplets, mm_cnt) %>%
    filter(!subject %in% c("index", "type")) %>%
    mutate(mm_cnt = ifelse((is.na(mm_eplets) | mm_eplets == "NA" | mm_eplets == "") & mm_cnt == 0, 0,
                           ifelse((!is.na(mm_eplets) | mm_eplets != "NA" | mm_eplets != "") & mm_cnt == 0, 1, mm_cnt + 1))) %>%
    # if it's a missing-allele then mm_eplet = "Not Found" and mm_cnt = 999
    mutate(mm_eplets = ifelse(str_detect(mm_eplets, "Not Found"), "Not Found", mm_eplets)) %>%
    mutate(     mm_cnt = ifelse(str_detect(mm_eplets, "Not Found"), 0, mm_cnt)) %>%
    mutate(mm_eplets = gsub("(,)\\1+", "\\1", mm_eplets))

  don_allele <- tbl_ready %>%
    select(pair_id, don_a1, don_a2, don_b1, don_b2, don_c1, don_c2) %>%
    pivot_longer(cols = starts_with("don_"),
                 names_to = "gene",
                 values_to = "hla") %>%
    mutate(gene = str_replace(gene, "don_", ""),
           pair_id = as.character(pair_id))

  # add hla to the result_single table
  re_s <- re_s %>%
    left_join(., don_allele, by =c("pair_id", "gene") ) %>%
    select(pair_id, subject, hla, gene, mm_eplets, mm_cnt) %>%
    left_join(., id_match, by = "pair_id") %>%
    select(-pair_id) %>%
    rename(pair_id = pair_id_ori,
           haplotype_id = gene) %>%
    mutate(haplotype_id = gsub("[a-zA-Z]", "", haplotype_id),
           mm_eplets = ifelse(str_detect(mm_eplets, "Not Found"), "Not Found", mm_eplets)) %>%
    select(pair_id, everything())
  #* end of step 7 *#

  #* step 8: mismatch count total *#
  re_o <- re_s %>%
    group_by(subject) %>%
    mutate(mm_cnt_tt = sum(mm_cnt)) %>%
    ungroup() %>%
    select(pair_id, subject, mm_cnt_tt) %>%
    distinct() %>%
    arrange(pair_id)
  #* end of step 8 *#

  #*step 9: mismatch count unique total*#
  re_o_uni <- tbl_ref

  st <- 3
  for (i in 1:subj_num) {
    ed <- st + 5
    pos <- c(st:ed)
    tmp <- tbl_ep_mm2 %>%
      select(c(1, 2, all_of(pos))) %>%
      setNames(c("index", "type", "a1_mm", "a2_mm", "b1_mm", "b2_mm", "c1_mm", "c2_mm")) %>%
      right_join(., tbl_ref,by = c("index", "type")) %>%
      mutate(var = ifelse(eplet %in% c(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm), eplet, NA)) %>%
      select(index, type, eplet, var) %>%
      setNames(c("index", "type", "eplet", subj_names[i]))

    st <- ed + 1
    re_o_uni <- left_join(re_o_uni, tmp, by = c("index", "type", "eplet"))
  }
  rm(st)

  re_o_uni <- re_o_uni %>%
    mutate(mm_cnt = subj_num - rowSums(is.na(.)),
           mm_pect = paste(round((subj_num - rowSums(is.na(.))) / subj_num * 100, 1), "%", sep = "")) %>%
    map(., ~sum(!is.na(.))) %>%
    data.frame() %>%
    select(-c(index, type, eplet, mm_cnt, mm_pect)) %>%
    gather() %>%
    setNames(c("subject", "mm_cnt_uniq"))

  re_o <- re_o %>% left_join(., re_o_uni, by = "subject")
  #* end of step 9 *#

  #* step 10: set mm_cnt = NA if allele is not found in matchmaker *#
  re_s <- re_s %>%
    rowwise() %>%
    mutate(mm_cnt = ifelse(mm_eplets == "Not Found", NA, mm_cnt)) %>%
    filter(!(pair_id %in% all_nas_don)) # filter out if all NAs on donor
  #* end of step 10 *#

  return(list(single_detail = re_s,
              overall_count = re_o))
}
