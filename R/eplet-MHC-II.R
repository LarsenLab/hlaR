#' calculate HLA class II eplet mismatch using Matchmaker algorithm
#' @param dat_in
#' dataframe with subject info (first 3 columns) and MHC II allele info
#' @param ver
#' version number of eplet mis-match table from epitopes.net to use
#  note: interlocus info only available in table v3
#' @return
#' data table with detailed single molecule level mis-match eplet info of each subject
#'
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
#' re <- CalEpletMHCII(dat_in = "YourDataFile", ver = 2)
#' re <- CalEpletMHCII(system.file("extdata", "MHC_II_test.csv", package = "hlaR"), ver = 3)
#' }
#'
# below is an example to use globalVariables() to suppress "no visible global variable" note
# utils::globalVariables(c("value", "locus", "index", "type", "mm"))

CalEpletMHCII <- function(dat_in, ver = 3) {
  ###*** step 1: import raw eplet tables ***###
  if(ver == 2){
  raw_eplet_A <- read.csv(system.file("extdata/ref", "MHC_II_eplet_A_v2.csv", package = "hlaR"), check.names = FALSE)
  raw_eplet_B <- read.csv(system.file("extdata/ref", "MHC_II_eplet_B_v2.csv", package = "hlaR"), check.names = FALSE)
  } else{
    raw_eplet_A <- read.csv(system.file("extdata/ref", "MHC_II_eplet_A_v3.csv", package = "hlaR"), check.names = FALSE)
    raw_eplet_B <- read.csv(system.file("extdata/ref", "MHC_II_eplet_B_v3.csv", package = "hlaR"), check.names = FALSE)
  }
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
    mutate(locus = ifelse(str_detect(locus, "\\."), sub("\\..*", "", locus), locus)) %>%
    filter(!locus %in% c("index", "type") ) %>%
    distinct() %>%
    reshape2::melt(id.vars = "locus") %>%
    filter(value != "")%>%
    distinct() %>%
    filter(value != " " )

  #* end of 2a *#

  # 2b: sub-lookup tables, these tables are used for mis-match comparison of each locus #
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

  lkup_dpa <- GenerateLookup(lkup_a, "DPA1")
  lkup_dqa <- GenerateLookup(lkup_a, "DQA1")

  lkup_dpb <- GenerateLookup(lkup_b, "DPB1")
  lkup_dqb <- GenerateLookup(lkup_b, "DQB1")
  lkup_drb <- GenerateLookup(lkup_b, c("DRB1", "DRB3", "DRB4", "DRB5"))

  #* end of 2b *#
  ###*** end of step 2 ***###

  ###*** step 3: import records with recipient/donor MHC II pairs ***###
  #* 3a. import data *#
  dat <- dat_in

  nm_rec <- c("rec_drb1", "rec_drb2", "rec_drw1", "rec_drw2", "rec_dqb1", "rec_dqb2", "rec_dqa1", "rec_dqa2", "rec_dpb1", "rec_dpb2", "rec_dpa1", "rec_dpa2")
  nm_don <- c("don_drb1", "don_drb2", "don_drw1", "don_drw2", "don_dqb1", "don_dqb2", "don_dqa1", "don_dqa2", "don_dpb1", "don_dpb2", "don_dpa1", "don_dpa2")

  rcpt <- dat %>%
    filter(donor_type %in% c("recipient", "recip", "rcpt", "r")) %>%
    select(-donor_type) %>%
    setNames(c("part_id", "part_type", nm_rec))

  don <- dat %>%
    filter(donor_type %in% c("donor", "don", "dn", "d")) %>%
    select(-donor_type) %>%
    setNames(c("part_id", "part_type", nm_don))

  dat <- left_join(rcpt, don, by = c("part_id", "part_type"))

  subj_num <- dim(dat)[1]
  tmp_names <- c(nm_rec, nm_don)

  ep_a <- raw_eplet_A %>% select(index, type)
  ep_b <- raw_eplet_B %>% select(index, type)
  #* end of 3a*#

  #* 3b: pull out eplet of each locus *#
  for (i in 1:subj_num) {
    allele <- toupper(unlist(transpose(dat[i,-c(1:2)]), use.names = F))
    allele <- ifelse(allele %in% c(names(raw_eplet_A), names(raw_eplet_B)), allele, NA)


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
  #* 4a: A loci *#
  results_a <- data.frame(subject = character(),
                          mm_eplets = character(),
                          mm_cnt = integer())

  # A alleles : DQA, DPA
  a_len <- subj_num * 8
  # exclude index and type, pulling data starting from the 3rd position
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 7
    pos <- c(st:ed)
    tmp <- ep_a %>%
      select(all_of(pos))
    subj_indx <- sub(".*\\.", "", names(tmp)[1])

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

    result_a <-  data.frame(t(tmp)) %>%
      tidyr::unite_(., paste(colnames(.), collapse="_"), colnames(.)) %>%
      setNames("mm_eplets") %>%
      mutate(subject = colnames(tmp),
             mm_eplets = gsub(",NA", "", gsub("_", ",", gsub("NA_", "", mm_eplets)))) %>%
      mutate(mm_cnt = str_count(mm_eplets, ",")) %>%
      mutate(mm_cnt = ifelse(mm_eplets == "NA" & mm_cnt == 0, 0,
                             ifelse(mm_eplets != "NA" & mm_cnt == 0, 1, mm_cnt + 1))) %>%
      filter(!subject %in% c("index", "type")) %>%
      select(subject, mm_eplets, mm_cnt)

    results_a <- rbind(results_a, result_a)
    st <- ed + 1
  }
  #* end of 4a *#

  #* 4b. B loci *#
  results_b <- data.frame(subject = character(),
                          mm_eplets = character(),
                          mm_cnt = integer())

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
  st1 <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st1 + 15
    pos <- c(st1:ed)
    tmp <- ep_b %>%
      select(all_of(pos))

    subj_indx <- sub(".*\\.", "", names(tmp)[1])

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

    result_b <-  data.frame(t(tmp)) %>%
      tidyr::unite_(., paste(colnames(.), collapse="_"), colnames(.)) %>%
      setNames("mm_eplets") %>%
      mutate(subject = colnames(tmp),
             mm_eplets = gsub(",NA", "", gsub("_", ",", gsub("NA_", "", mm_eplets)))) %>%
      mutate(mm_cnt = str_count(mm_eplets, ",")) %>%
      mutate(mm_cnt = ifelse(mm_eplets == "NA" & mm_cnt == 0, 0,
                             ifelse(mm_eplets != "NA" & mm_cnt == 0, 1, mm_cnt + 1))) %>%
      filter(!subject %in% c("index", "type")) %>%
      select(subject, mm_eplets, mm_cnt)

    results_b <- rbind(results_b, result_b)

    st1 <- ed + 1
  }
  #* end of 4b *#
  ###*** end of step 4 ***###

  ###*** step 5: generate output table ***###
  results <- rbind(results_a, results_b) %>%
             mutate(name = gsub(".*\\_","", subject),
                    gene = gsub("\\_.*", "", subject)) %>%
             select(name, gene, mm_eplets, mm_cnt) %>%
             arrange(name, gene)

  return(results)
  ###*** end of step 5 ***###
}

