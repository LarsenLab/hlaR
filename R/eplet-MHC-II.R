#' @name CalEpletMHCII
#' @title Calculate class II HLA eplet mismatch.
#' @description Use high resolution HLA(Human Leukocyte Antigen) class II data to calculate class II eplet mismatch for a population of donors and recipients. Mismatch is calculated using logic from 'HLAMatchMaker', developed by Rene Dusquesnoy. Current reference tables supported are 'HLAMatchMaker' v2 and v3. Note: interlocus info only available in v3 reference tables.
#' @param dat_in
#' A dataframe of recipient and donor's high resolution MHC II data. Each recipient and donor pair are linked by are the “pair_id” column and differentiated by the “subject_type” column.
#' @param ver
#' Version number of HLAMatchMaker based eplet reference table to use.
#' @return
#' A list of data tables.
#' - `single_detail`: single molecule class II MHC eplet mismatch table, including mismatched eplet names and the count of eplets mismatched at each allele.
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
#' dat <- read.csv(system.file("extdata/example","MHC_II_test.csv",package="hlaR"),sep=",",header=TRUE)
#' re <- CalEpletMHCII(dat, ver = 2)

CalEpletMHCII <- function(dat_in, ver = 2) {
  #* step 0: check if recipient and donor are paired *#
  num_rcpt <- length(dat_in$subject_type[dat_in$subject_type %in% c("recipient", "recip", "rcpt", "r")])
  num_don  <- length(dat_in$subject_type[dat_in$subject_type %in% c("donor", "don", "dn", "d")])
  if(num_rcpt == num_don){
    rm(num_rcpt, num_don)
  } else {
    stop("please check that every pair_id has both recipient and donor data.")
  }
  #* end of step 0 *#

  #* step 1: import raw eplet tables *#
  if(ver == 2) {
    tbl_raw_eplet_A <- read.csv(system.file("extdata/ref", "MHC_II_eplet_A_v2.csv", package = "hlaR"), check.names = FALSE)
    tbl_raw_eplet_B <- read.csv(system.file("extdata/ref", "MHC_II_eplet_B_v2.csv", package = "hlaR"), check.names = FALSE)
  } else {
    tbl_raw_eplet_A <- read.csv(system.file("extdata/ref", "MHC_II_eplet_A_v3.csv", package = "hlaR"), check.names = FALSE)
    tbl_raw_eplet_B <- read.csv(system.file("extdata/ref", "MHC_II_eplet_B_v3.csv", package = "hlaR"), check.names = FALSE)
  }
  #* end of step 1 *#

  #* step 2: generate base lookup tables for MHC II loci As and Bs *#
  #* 2a: master lookup tables for As and Bs *#
  tbl_ref_A <- as.data.frame(t(tbl_raw_eplet_A)) %>%
    setNames(paste(tbl_raw_eplet_A$type, tbl_raw_eplet_A$index, sep = "_" )) %>%
    rownames_to_column(var = "locus") %>%
    mutate(locus = ifelse(str_detect(locus, "\\*"), sub("\\*.*", "", locus), locus)) %>%
    filter(!locus %in% c("index", "type") ) %>%
    distinct() %>%
    reshape2::melt(id.vars = "locus") %>%
    filter(value != "" ) %>%
    distinct()

  tbl_ref_B <- as.data.frame(t(tbl_raw_eplet_B)) %>%
    setNames(paste(tbl_raw_eplet_B$type, tbl_raw_eplet_B$index, sep = "_" )) %>%
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
  tbl_ref_dpa <- GenerateLookup(tbl_ref_A, "DPA1")
  tbl_ref_dqa <- GenerateLookup(tbl_ref_A, "DQA1")

  tbl_ref_dpb <- GenerateLookup(tbl_ref_B, "DPB1")
  tbl_ref_dqb <- GenerateLookup(tbl_ref_B, "DQB1")
  tbl_ref_drb <- GenerateLookup(tbl_ref_B, c("DRB1", "DRB3", "DRB4", "DRB5"))
  #* end of 2b *#
  #* end of step 2 *#

  #* step 3: import subject data *#
  #* 3a. import data and some definitions *#
  nm_rec <- c("rec_drb1", "rec_drb2", "rec_drw1", "rec_drw2", "rec_dqb1", "rec_dqb2", "rec_dqa1", "rec_dqa2", "rec_dpb1", "rec_dpb2", "rec_dpa1", "rec_dpa2")
  nm_don <- c("don_drb1", "don_drb2", "don_drw1", "don_drw2", "don_dqb1", "don_dqb2", "don_dqa1", "don_dqa2", "don_dpb1", "don_dpb2", "don_dpa1", "don_dpa2")
  tmp_names <- c(nm_rec, nm_don)

  #* start of recipient *#
  tmp_rcpt <- dat_in %>%
    filter(subject_type %in% c("recipient", "recip", "rcpt", "r")) %>%
    select(-subject_type) %>%
    setNames(c("pair_id", nm_rec))

  # warning: if all recipient's alleles are NA
  all_nas <- tmp_rcpt %>%
    filter(rec_drb1 == "" & rec_drb2 == "" &  rec_drw1 == "" &  rec_drw2 == "" &  rec_dqb1 == "" &
             rec_dqb2 == "" &  rec_dqa1 == "" & rec_dqa2 == "" &  rec_dpb1 == "" &  rec_dpb2 == "" &  rec_dpa1 == "" &  rec_dpa2 == "" ) %>% pull(pair_id)

  if(length(all_nas) != 0) {
    warning("No MHC class II alleles detected for recipient(s) ", toString(all_nas), ".")
    cat("\n")
  }

  # warning: if recipient's missing any of dr, dqa, dqb values
  drdq_nas_rcpt <- tmp_rcpt %>%
    filter((rec_drb1 == "" & rec_drb2 == "") | (rec_drw1 == "" & rec_drw2 == "") | (rec_dqb1 == "" & rec_dqb2 == "") | (rec_dqa1 == "" & rec_dqa2 == "" )) %>%
    pull(pair_id)

  if(length(drdq_nas_rcpt) != 0) {
    warning("No DRB or DQA or DQB detected for recipient(s) ", toString(drdq_nas_rcpt), ".")
    cat("\n")
  }
  #* end of recipient *#

  #* start of donor *#
  tmp_don <- dat_in %>%
    filter(subject_type %in% c("donor", "don", "dn", "d")) %>%
    select(-subject_type) %>%
    setNames(c("pair_id", nm_don))

  # warning: if all donot's alleles are NA
  all_nas <- tmp_don %>%
    filter(don_drb1 == "" & don_drb2 == "" &  don_drw1 == "" &  don_drw2 == "" &  don_dqb1 == "" &
             don_dqb2 == "" &  don_dqa1 == "" & don_dqa2 == "" &  don_dpb1 == "" &  don_dpb2 == "" &  don_dpa1 == "" &  don_dpa2 == "" ) %>%
    pull(pair_id)

  if(length(all_nas) != 0) {
    warning("No MHC class II alleles detected for donor(s) ", toString(all_nas), ".")
    cat("\n")
  }

  # warning: if donor's missing any of dr, dqa, dqb
  drdq_nas_don <- tmp_don %>%
    filter((don_drb1 == "" & don_drb2 == "") | (don_dqb1 == "" & don_dqb2 == "") | (don_dqa1 == "" & don_dqa2 == "" )) %>%
    pull(pair_id)

  if(length(drdq_nas_don) != 0) {
    warning("No DRB or DQA or DQB detected for donor(s) ", toString(drdq_nas_don), ".")
    cat("\n")
  }
  #* end of donor *#

  tbl_ready <- merge(tmp_rcpt, tmp_don, by = c("pair_id")) %>%
    arrange(pair_id)

  rm(tmp_rcpt, tmp_don, all_nas)

  # create id_match table to trace non-sequential pair_id
  match_id <- tbl_ready %>%
    select(pair_id) %>%
    mutate(match_id = dense_rank(pair_id))

  subj_num <- dim(tbl_ready)[1]

  tbl_ep_a <- tbl_raw_eplet_A %>% select(index, type)
  tbl_ep_b <- tbl_raw_eplet_B %>% select(index, type)

  # index of names of As and Bs
  indx_as <- grep("dqa|dpa", tmp_names)
  indx_bs <- grep("drb|drw|dqb|dpb", tmp_names)

  # data frame to hold HLA typing information, used to merge to the single molecule level report
  allele_detail <- data.frame(allele = character(),
                              name = character())
  #* end of 3a*#

  #* 3b: pull out eplet of each locus *#
  for (i in 1:subj_num) {
    # convert alleles to upper case
    allele <- toupper(unlist(transpose(tbl_ready[i, tmp_names]), use.names = F))

    # record allele names that are not in the reference table, and generate a warning message
    miss_allele <- allele[!(allele %in% c(names(tbl_raw_eplet_A), names(tbl_raw_eplet_B)))]
    miss_allele <- miss_allele[miss_allele != ""]

    if(length(miss_allele) >= 1) {
      warning(paste0(miss_allele, " is not found in the refernce table. Please check!\n"),
              immediate. = TRUE)
    }

    # if allele is blank, then set it to NA
    allele <- ifelse(allele %in% c(names(tbl_raw_eplet_A), names(tbl_raw_eplet_B), miss_allele),
                     allele, NA)

    # for each allele, find it's eplets
    for (j in 1:length(allele)) {
      varname <- paste0(tmp_names[j], ".", sep = i)

      # keep track on detailed hlas
      allele_detail <-  rbind(allele_detail, c(paste0(tmp_names[j], ".", sep = i),allele[j]))

      # A loci
      if(j %in% indx_as) {
        if (!is.na(allele[j]) & !(allele[j] %in% miss_allele)) {
          tmp <- tbl_raw_eplet_A %>%
            select(index, type, allele[j]) %>%
            setNames(c("index", "type", varname))
          tbl_ep_a <- tbl_ep_a %>%
            left_join(., tmp, by = c("index", "type"))

        } else if (!is.na(allele[j]) & allele[j] %in% miss_allele){
          tbl_ep_a <- tbl_ep_a %>%
            mutate(!!varname := "Not Found")
        } else {
          tbl_ep_a <- tbl_ep_a %>%
            mutate(!!varname := NA)
        }
      }
      # B loci
      else {
        if (!is.na(allele[j]) & !(allele[j] %in% miss_allele)) {
          tmp <- tbl_raw_eplet_B %>%
            select(index, type, allele[j]) %>%
            setNames(c("index", "type", varname))
          tbl_ep_b <- tbl_ep_b %>%
            left_join(., tmp, by = c("index", "type"))
        } else if (!is.na(allele[j]) & allele[j] %in% miss_allele) {
          tbl_ep_b <- tbl_ep_b %>%
            mutate(!!varname := "Not Found")
        } else {
          tbl_ep_b <- tbl_ep_b %>%
            mutate(!!varname := NA)
        }
      }
    }
  }

  # make subject name, gene, and allele clear
  allele_detail <- allele_detail %>%
    setNames(c("subj", "hla")) %>%
    filter(str_detect(subj, "don_")) %>% # keep donor's hla only
    mutate(name = paste0("subj",gsub(".*\\.", "", subj)),
           gene = gsub(".*_", "",gsub("\\..*", "", subj))) %>%
    select(name, gene, hla)

  # set blanks to NA to suit mismatch count calculation logic
  tbl_ep_a <- tbl_ep_a %>% na_if(., "")
  tbl_ep_b <- tbl_ep_b %>% na_if(., "")

  #* end of 3b *#
  #* end of step 3 *#

  #* step 4: mis-match eplet calculation *#
  #* 4a: alpha loci *#
  result_a_single <- data.frame(subject = character(),
                                mm_eplets = character(),
                                mm_cnt = integer())

  # length of alphas : DQA, DPA
  a_len <- subj_num * 8

  # exclude index and type, pulling data starting from the 3rd position
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 7
    pos <- c(st:ed)
    tmp <- tbl_ep_a %>%
      select(all_of(pos))

    subj_indx <- sub(".*\\.", "", names(tmp)[1])

    colnames(tmp) <- c("rec_dqa1", "rec_dqa2", "rec_dpa1", "rec_dpa2",
                       "don_dqa1", "don_dqa2", "don_dpa1", "don_dpa2")

    # important : replace recipient's "Not Found" to NA to avoid incorrect mismatch calculation in later step
    # if there are missing alleles on the same locus in both recpt and donor, this way we can make sure the donor's mm_eplet is "Not Found" instead of blank
    nm <- names(tmp)[str_detect(names(tmp), "rec_")]
    tmp <- tmp %>%
      mutate_at(vars(all_of(nm)),
                list(~ifelse(. == "Not Found", NA, .)))

    # if donor's allele doesn't show in the reference table, then set it to "Not Found",
    # otherwise, check for mismatch
    tmp <- tmp %>%
      mutate(dqa1_mm = ifelse(don_dqa1 %in% miss_allele,
                              "Not Found",
                              ifelse(don_dqa1 %in% c(rec_dqa1, rec_dqa2), NA, don_dqa1)),
             dqa2_mm = ifelse(don_dqa2 %in% miss_allele,
                              "Not Found",
                              ifelse(don_dqa2 %in% c(rec_dqa1, rec_dqa2), NA, don_dqa2)),
             dpa1_mm = ifelse(don_dpa1 %in% miss_allele,
                              "Not Found",
                              ifelse(don_dpa1 %in% c(rec_dpa1, rec_dpa2), NA, don_dpa1)),
             dpa2_mm = ifelse(don_dpa2 %in% miss_allele,
                              "Not Found",
                              ifelse(don_dpa2 %in% c(rec_dpa1, rec_dpa2), NA, don_dpa2))) %>%
      select(dqa1_mm, dqa2_mm, dpa1_mm, dpa2_mm) %>%
      setNames(c(paste0("dqa1_mm_subj", subj_indx),
                 paste0("dqa2_mm_subj", subj_indx),
                 paste0("dpa1_mm_subj", subj_indx),
                 paste0("dpa2_mm_subj", subj_indx)))

    tmp_a <-  data.frame(t(tmp)) %>%
      unite("mm_eplets", names(.), na.rm = TRUE, sep = ",", remove = FALSE)  %>%
      mutate(subject = gsub(".*_", "", rownames(.)),
             pair_id = gsub(".*subj", "", subject),
             gene = gsub("_.*", "", rownames(.)),
             mm_cnt = str_count(mm_eplets, ",")) %>%
      mutate(mm_cnt = ifelse((is.na(mm_eplets) | mm_eplets == "NA" | mm_eplets == "") & mm_cnt == 0, 0,
                             ifelse((!is.na(mm_eplets) | mm_eplets != "NA" | mm_eplets != "") & mm_cnt == 0, 1, mm_cnt + 1))) %>%
      select(subject, pair_id, gene, mm_eplets, mm_cnt)

    result_a_single <- rbind(result_a_single, tmp_a)
    st <- ed + 1
  }

  # if allele is missing from the ref table, then set mm_cnt to 0 ( will change it to NA in the very last step)
  result_a_single <- result_a_single %>%
    mutate(mm_eplets = ifelse(str_detect(mm_eplets, "Not Found"), "Not Found", mm_eplets),
           mm_cnt = ifelse(str_detect(mm_eplets, "Not Found"), 0, mm_cnt))
  #* end of 4a *#

  #* 4b. beta loci *#
  result_b_single <- data.frame(subject = character(),
                                mm_eplets = character(),
                                mm_cnt = integer())

  tbl_ep_b_mm <- tbl_ep_b %>%
    select(index, type)
  ep_b_rownames <- paste(tbl_ep_b_mm$type, tbl_ep_b_mm$index, sep = "")

  rownames(tbl_ep_b_mm) <- ep_b_rownames
  rownames(tbl_ep_b) <- ep_b_rownames

  # if interlocus eplets(names contains "Int"), DPB, then compare one don to ALL of recipients'
  # this condition only applies for matchmaker v3
  int_loc_eplet<- ep_b_rownames[str_detect(ep_b_rownames, regex('int', ignore_case = T))]

  # B loci
  # B alleles : DRB, DRw, DQB, DPB
  b_len <- subj_num * 16

  # exclude index and type, pulling data starting from the 3rd position
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 15
    pos <- c(st:ed)
    tmp <- tbl_ep_b %>%
      select(all_of(pos))

    subj_indx <- sub(".*\\.", "", names(tmp)[1])

    colnames(tmp) <- c("rec_drb1", "rec_drb2", "rec_drw1", "rec_drw2", "rec_dqb1", "rec_dqb2", "rec_dpb1","rec_dpb2",
                       "don_drb1", "don_drb2", "don_drw1", "don_drw2", "don_dqb1", "don_dqb2", "don_dpb1","don_dpb2")

    tmp$rname <- ep_b_rownames

    # important : replace recipient's "Not Found" to NA to avoid incorrect mismatch calculation in later step
    # if there are missing alleles on the same locus in both recpt and donor, this way we can make sure the donor's mm_eplet is "Not Found" instead of blank
    nm <- names(tmp)[str_detect(names(tmp), "rec_")]
    tmp <- tmp %>%
      mutate_at(vars(all_of(nm)),
                list(~ifelse(. == "Not Found", NA, .)))

    # if donor's allele doesn't show in the reference table, then set it to "Not Found",
    # otherwise, check for mismatch
    tmp <- tmp %>%
      mutate(drb1_mm = ifelse(don_drb1 %in% miss_allele,
                              "Not Found",
                              ifelse(don_drb1 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drb1)),
             drb2_mm = ifelse(don_drb2 %in% miss_allele,
                              "Not Found"
                              ,ifelse(don_drb2 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drb2)),
             drw1_mm = ifelse(don_drw1 %in% miss_allele,
                              "Not Found",
                              ifelse(don_drw1 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drw1)),
             drw2_mm = ifelse(don_drw2 %in% miss_allele,
                              "Not Found",
                              ifelse(don_drw2 %in% c(rec_drb1, rec_drb2, rec_drw1, rec_drw2), NA, don_drw2)),
             dqb1_mm = ifelse(don_dqb1 %in% miss_allele,
                              "Not Found",
                              ifelse(don_dqb1 %in% c(rec_dqb1, rec_dqb2), NA, don_dqb1)),
             dqb2_mm = ifelse(don_dqb2 %in% miss_allele,
                              "Not Found",
                              ifelse(don_dqb2 %in% c(rec_dqb1, rec_dqb2), NA, don_dqb2)),
             dpb1_tmp2 = ifelse(don_dpb1 %in% miss_allele,
                                "Not Found",
                                ifelse(don_dpb1 %in% c(rec_dpb1, rec_dpb2), NA, don_dpb1)),
             dpb1_tmp8 = ifelse(don_dpb1 %in% miss_allele,
                                "Not Found",
                                ifelse(don_dpb1 %in% c(rec_dpb1, rec_dpb2, rec_drb1, rec_drb2, rec_drw1, rec_drw2, rec_dqb1, rec_dqb2), NA, don_dpb1)),
             dpb2_tmp2 = ifelse(don_dpb2 %in% miss_allele,
                                "Not Found",
                                ifelse(don_dpb2 %in% c(rec_dpb1, rec_dpb2), NA, don_dpb2)),
             dpb2_tmp8 = ifelse(don_dpb2 %in% miss_allele,
                                "Not Found",
                                ifelse(don_dpb2 %in% c(rec_dpb1, rec_dpb2, rec_drb1, rec_drb2, rec_drw1, rec_drw2, rec_dqb1, rec_dqb2), NA, don_dpb2))) %>%
      mutate(dpb1_mm =  ifelse(rname %in% int_loc_eplet , dpb1_tmp8, dpb1_tmp2),
             dpb2_mm = ifelse(rname %in% int_loc_eplet, dpb2_tmp8, dpb2_tmp2)) %>%
      select(drb1_mm, drb2_mm, drw1_mm, drw2_mm, dqb1_mm, dqb2_mm, dpb1_mm, dpb2_mm ) %>%
      setNames(c(paste0("drb1_mm_subj", subj_indx),
                 paste0("drb2_mm_subj", subj_indx),
                 paste0("drw1_mm_subj", subj_indx),
                 paste0("drw2_mm_subj", subj_indx),
                 paste0("dqb1_mm_subj", subj_indx),
                 paste0("dqb2_mm_subj", subj_indx),
                 paste0("dpb1_mm_subj", subj_indx),
                 paste0("dpb2_mm_subj", subj_indx)))

    tmp_b <-  data.frame(t(tmp)) %>%
      unite("mm_eplets", names(.), na.rm = TRUE, sep = ",", remove = FALSE) %>%
      mutate(subject = gsub(".*_", "", rownames(.)),
             pair_id = gsub(".*subj", "", subject),
             gene = gsub("_.*", "", rownames(.)),
             mm_cnt = str_count(mm_eplets, ",")) %>%
      mutate(mm_cnt = ifelse((is.na(mm_eplets) | mm_eplets == "NA" | mm_eplets == "") & mm_cnt == 0, 0,
                             ifelse((!is.na(mm_eplets) | mm_eplets != "NA" | mm_eplets != "") & mm_cnt == 0, 1, mm_cnt + 1))) %>%
      select(subject, pair_id, gene, mm_eplets, mm_cnt)

    result_b_single <- rbind(result_b_single, tmp_b)

    st <- ed + 1
  }

  # if allele is missing from the ref table, then set mm_cnt to 0 ( will change it to NA in the very last step)
  result_b_single <- result_b_single %>%
    mutate(mm_eplets = ifelse(str_detect(mm_eplets, "Not Found"), "Not Found", mm_eplets),
           mm_cnt = ifelse(str_detect(mm_eplets, "Not Found"), 0, mm_cnt))
  #* end of 4b *#
  #* end of step 4 *#

  #* step 5: final result - single molecule *#
  re_s <- rbind(result_a_single, result_b_single) %>%
    mutate(name = subject,
           match_id = as.numeric(str_replace(name, "subj", ""))) %>%
    select(name, gene, mm_eplets, mm_cnt, match_id) %>%
    arrange(name, gene) %>%
    left_join(., match_id, by = "match_id") %>%
    select(-match_id) %>%
    select(pair_id, everything()) %>%
    left_join(., allele_detail, by = c("name", "gene")) %>% # join raw hla typing back to the output table
    dplyr::rename(subject = name, haplotype_id = gene) %>%
    mutate(haplotype_id = gsub("[a-zA-Z]", "", haplotype_id)) %>%
    select(pair_id, subject, hla, haplotype_id, mm_eplets, mm_cnt) %>%
    filter(!is.na(hla))

  #* step 6: overall count, all and unique *#
  # count unique mismatch eplets
  re_o <- re_s %>%
    group_by(subject) %>%
    mutate(mm_cnt_tt = sum(mm_cnt)) %>%
    ungroup() %>%
    select(pair_id, subject, mm_cnt_tt) %>%
    distinct() %>%
    arrange(pair_id)

  re_o_uni <- rbind(result_a_single, result_b_single) %>%
    # convert "Not Found" and blanks to NA for unite() in later step %>%
    mutate(mm_eplets = ifelse(mm_eplets %in% c("Not Found", ""), NA, mm_eplets),
           pair_id = as.numeric(str_replace(subject, "subj", ""))) %>%
    select(subject, gene, mm_eplets, mm_cnt, pair_id) %>%
    group_by(subject, pair_id) %>%
    select(pair_id, subject, mm_eplets) %>%
    mutate(rn = row_number()) %>%
    pivot_wider(names_from = rn,
                values_from = c("mm_eplets")) %>%
    ungroup() %>%
    select(names(.)[!is.na(names(.))]) %>%
    unite("mm_eplets", `1`:`12`, na.rm = TRUE, sep = ",", remove = FALSE) %>%
    mutate(mm_eplets = sapply(strsplit(mm_eplets, ","),
                              function(x) x = paste(unique(x), collapse = ","))) %>%
    mutate(mm_cn_uniq = str_count(mm_eplets, ",") + 1) %>%
    select(subject, mm_cn_uniq)

  re_o <- re_o %>% left_join(., re_o_uni, by = "subject")
  #* end of step 6 *#

  #* step 7: set mm_cnt to NA if eplet is missing *#
  re_s <- re_s %>%
    rowwise() %>%
    mutate(mm_cnt = ifelse(mm_eplets == "Not Found", NA, mm_cnt))
  #* end of step 7 *#

  #* step 8: calculate risk score based on DR DQ mismatch counts *#
  dqdr_risk <- list()
  subj_id <- unique(re_s$pair_id)

  # for each subject, call CalRiskScr() for risk score
  if (length(subj_id > 0)) {
    for (i in 1:length(subj_id)){
      dqdr_risk[[i]] <- re_s %>%
        filter(pair_id %in% subj_id[i]) %>%
        CalRiskScr()
      dqdr_risk[[i]]$pair_id <- subj_id[i]
    }

    dqdr_risk <- bind_rows(dqdr_risk) %>%
      # mutate(DQ = ifelse(DQ == 0, NA, DQ),
      #        DR = ifelse(DR == 0, NA, DR)) %>%
      select(pair_id, DQ, DR, risk) %>%
      arrange(pair_id)
  } else {
    dqdr_risk = data.frame(pair_id = NA, DQ = NA, DR = NA, risk = NA)
  }

  #* end of step 8 *#

  return(list(single_detail = re_s,
              overall_count = re_o,
              dqdr_risk = dqdr_risk))
}

