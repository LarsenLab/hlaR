#' @name CalEpletMHCI
#' @title Calculate HLA Class I eplet mismatch. Eplet Reference table and mistmatch calculation logic are based on MatchMaker.
#' @param dat_in
#' A dataframe with subject info(first 2 columns) and MHC I allele info.
#' Each unique participant id has 2 rows associated with it, 1 for recipient, 1 for donor.
#' @param ver
#' Version number of MatchMaker based eplet reference to use.
#' @return
#' A list of data tables.
#' single_detail: single moelcule level eplet mismatch table, including mismatch eplet name and count of each allele
#' overall_count: original input data appended with count of unique mis-matched eplet
#' overall_detail: percentage of mismatch across all subject of each eplet
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
#'
#' @examples
#' \dontrun{
# dat <- read.csv(system.file("extdata/example", "MHC_I_test.csv", package = "hlaR"), sep = ",", header = TRUE)
#' re <- CalEpletMHCI(dat_in = dat, ver = 3)
#' }
#'
# below is an example to use globalVariables() to suppress "no visible global variable" note
# utils::globalVariables(c("value", "locus", "index", "type", "mm"))

CalEpletMHCI <- function(dat_in, ver = 3) {

  #* step 1: import raw eplet table *#
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

  #* step 2: patient table *#
  nm_rec <- c("rec_a1", "rec_a2", "rec_b1", "rec_b2", "rec_c1", "rec_c2")
  nm_don <- c("don_a1", "don_a2", "don_b1", "don_b2", "don_c1", "don_c2")

  tmp_rcpt <- dat_in %>%
    filter(donor_type %in% c("recipient", "recip", "tmp_rcpt", "r")) %>%
    select(-donor_type) %>%
    setNames(c("part_id", nm_rec))

  tmp_don <- dat_in %>%
    filter(donor_type %in% c("donor", "don", "dn", "d")) %>%
    select(-donor_type) %>%
    setNames(c("part_id", nm_don))

  tbl_ready <- left_join(tmp_rcpt, tmp_don, by = c("part_id")) %>%
    mutate(part_id_ori = part_id) %>%
    arrange(part_id_ori) %>%
    mutate(part_id = dense_rank(part_id_ori))

  rm(tmp_rcpt, tmp_don)

  # create id_match table in case pari_id is not sequential in the patient table
  id_match <- tbl_ready %>%
    select(part_id_ori, part_id) %>%
    mutate(part_id = as.character(part_id))

  tbl_ready <- tbl_ready %>%
    select(-part_id_ori) %>%
    select(part_id, everything())

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
    allele <- ifelse(allele %in% names(tbl_raw_eplet), allele, NA)

    for (j in 1:length(allele)) {
      varname <- paste0(tmp_names[j], ".", sep = i)
      if (!is.na(allele[j])) {
        tmp <- tbl_raw_eplet %>%
          select(index, type, allele[j])
        tmp <- tmp %>%
          setNames(c("index", "type", varname))
        tbl_ep <- tbl_ep %>%
          left_join(., tmp, by = c("index", "type"))
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

  tbl_ep_mm2 <- tbl_ep %>%
    select(index, type)

  # exclude index and type, pulling data for each allele, starting from 3rd position
  st <- 3

  # for each subject, compare eplets of donor's EACH allele for ALL of recipients'
  for (i in 1:subj_num) {
    ed <- st + 11
    positions <- c(st:ed)
    tmp <- tbl_ep %>%
      select(all_of(positions))
    subj_indx <- sub(".*\\.", "", names(tmp)[1])

    colnames(tmp) <- c(nm_rec, nm_don)

    # comparison
    tmp <- tmp %>%
      mutate(a1_mm = ifelse(don_a1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_a1),
             a2_mm = ifelse(don_a2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_a2),
             b1_mm = ifelse(don_b1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_b1),
             b2_mm = ifelse(don_b2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_b2),
             c1_mm = ifelse(don_c1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_c1),
             c2_mm = ifelse(don_c2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_c2)) %>%
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
  #* end of step 5 *#

  #* step 6: compare mis-match with tbl_ref table *#
  # reset starting position for another round of loop
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 5
    positions <- c(st:ed)
    tmp <- tbl_ep_mm %>%
      select(c(1, 2, all_of(positions)))

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
  #* end of step 6 *#

  #* step 7: final result - single molecule report *#
  # subject names
  subj_names <- unique(sub(".*\\_", "", names(tbl_ep_mm2)[-c(1:2)]))

  result_single <-  data.frame(t(tbl_ep_mm2)) %>%
    unite("mm_eplets", names(.), na.rm = TRUE, sep = ",", remove = FALSE) %>%
    mutate(subject = gsub(".*_", "", rownames(.)),
           part_id = gsub(".*subj", "", subject),
           gene = gsub("_.*", "", rownames(.)),
           mm_cnt = str_count(mm_eplets, ",")) %>%
    select(subject, part_id, gene, mm_eplets, mm_cnt) %>%
    filter(!subject %in% c("index", "type")) %>%
    mutate(mm_cnt = ifelse((is.na(mm_eplets) | mm_eplets == "NA" | mm_eplets == "") & mm_cnt == 0, 0,
                           ifelse((!is.na(mm_eplets) | mm_eplets != "NA" | mm_eplets != "") & mm_cnt == 0, 1, mm_cnt + 1)))

  don_allele <- tbl_ready %>%
    select(part_id, don_a1, don_a2, don_b1, don_b2, don_c1, don_c2) %>%
    pivot_longer(cols = starts_with("don_"),
                 names_to = "gene",
                 values_to = "hla") %>%
    mutate(gene = str_replace(gene, "don_", ""),
           part_id = as.character(part_id))

  # add hla to the result_single table
  result_single <- result_single %>%
    left_join(., don_allele, by =c("part_id", "gene") ) %>%
    select(part_id, subject, hla, mm_eplets, mm_cnt) %>%
    left_join(., id_match, by = "part_id") %>%
    select(-part_id) %>%
    rename(part_id = part_id_ori) %>%
    select(part_id, everything())
  #* end of step 7 *#

  #* step 8: final result - overall report *#
  result_overall <- tbl_ref

  st <- 3
  for (i in 1:subj_num) {
    ed <- st + 5
    positions <- c(st:ed)
    tmp <- tbl_ep_mm2 %>%
      select(c(1, 2, all_of(positions))) %>%
      setNames(c("index", "type", "a1_mm", "a2_mm", "b1_mm", "b2_mm", "c1_mm", "c2_mm")) %>%
      right_join(., tbl_ref,by = c("index", "type")) %>%
      mutate(var = ifelse(eplet %in% c(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm), eplet, NA)) %>%
      select(index, type, eplet, var) %>%
      setNames(c("index", "type", "eplet", subj_names[i]))

    st <- ed + 1
    result_overall <- left_join(result_overall, tmp, by = c("index", "type", "eplet"))
  }

  result_detail <- result_overall %>%
    mutate(mm_cnt = subj_num - rowSums(is.na(.)),
           mm_pect = paste(round((subj_num - rowSums(is.na(.))) / subj_num * 100, 1), "%", sep = ""))

  # count unique mismatch eplets
  result_count <- result_detail %>%
    map(., ~sum(!is.na(.))) %>%
    data.frame() %>%
    select(-c(index, type, eplet, mm_cnt, mm_pect)) %>%
    gather() %>%
    setNames(c("part_id", "mm_count")) %>%
    mutate(part_id = as.numeric(str_replace(part_id, "subj", ""))) %>%
    right_join(., tbl_ready, by = "part_id")
  #* end of step 8 *#

  return(list(single_detail = result_single,
              overall_count = result_count,
              overall_detail = result_detail))
}
