#' @name CalEpletMHCI
#' @title calculate HLA Class I eplet mismatch using MatchMaker reference table and algorithm
#' @param dat_in
#' dataframe with subject info(first 3 columns) and MHC I allele info
#' each unique participant id has 2 rows associated with it, 1 for recipient, 1 for donor
#' @param ver
#' version number of eplet mis-match table from epitopes.net to use
#' @return
#' data table with detailed single molecule level mis-match eplet info of each subject
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
#' re2 <- CalEpletMHCI(dat, ver = 2)
#' re3 <- CalEpletMHCI(dat, ver = 3)
#' }
#'
# below is an example to use globalVariables() to suppress "no visible global variable" note
# utils::globalVariables(c("value", "locus", "index", "type", "mm"))

CalEpletMHCI <- function(dat_in, ver = 3) {

  #* step 1: import raw eplet table *#
  if(ver == 2){
    raw_eplet <- read.csv(system.file("extdata/ref", "MHC_I_eplet_v2.csv", package = "hlaR"), check.names = FALSE)
  } else{
    raw_eplet <- read.csv(system.file("extdata/ref", "MHC_I_eplet_v3.csv", package = "hlaR"), check.names = FALSE)
  }

  raw_lookup <- as.data.frame(t(raw_eplet)) %>%
                setNames(paste(raw_eplet$type, raw_eplet$index, sep = "_" )) %>%
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

  #* step 2: import patient table *#
  nm_rec <- c("rec_a1", "rec_a2", "rec_b1", "rec_b2", "rec_c1", "rec_c2")
  nm_don <- c("don_a1", "don_a2", "don_b1", "don_b2", "don_c1", "don_c2")

  rcpt <- dat_in %>%
            filter(donor_type %in% c("recipient", "recip", "rcpt", "r")) %>%
            select(-donor_type) %>%
            setNames(c("part_id", "part_type", nm_rec))

  don <- dat_in %>%
          filter(donor_type %in% c("donor", "don", "dn", "d")) %>%
          select(-donor_type) %>%
          setNames(c("part_id", "part_type", nm_don))

  dat <- left_join(rcpt, don, by = c("part_id")) %>%
          select(-part_type.y) %>%
          dplyr::rename(part_type = part_type.x,
                        part_id_ori = part_id) %>%
          arrange(part_id_ori) %>%
          mutate(part_id = dense_rank(part_id_ori))

  id_match <- dat %>%
              select(part_id_ori, part_id) %>%
              mutate(part_id = as.character(part_id))

  dat <- dat %>%
          select(-part_id_ori)

  subj_num <- dim(dat)[1]
  tmp_names <- c(nm_rec, nm_don)

  #* step 3: initial eplet data frame *#
  dat_ep <- raw_eplet %>%
            select(index, type)

  #* step 4: pull out eplet of each allele *#
  for (i in 1:subj_num) {
    allele <- toupper(unlist(transpose(dat[i,-c(1:2)]), use.names = F))
    allele <- ifelse(allele %in% names(raw_eplet), allele, NA)

    for (j in 1:length(allele)) {
      varname <- paste0(tmp_names[j], ".", sep = i)
      if (!is.na(allele[j])) {
        tmp <- raw_eplet %>%
                select(index, type, allele[j])
        tmp <- tmp %>%
                setNames(c("index", "type", varname))
        dat_ep <- dat_ep %>%
                  left_join(., tmp, by = c("index", "type"))
      } else{
        dat_ep <- dat_ep %>%
                  mutate(!!varname := NA)
      }
    }
  }

  #* step 5: mark mis-matches *#
  dat_ep_mm <- dat_ep %>%
                select(index, type)

  # exclude index and type, pulling data starting from the 3rd position
  st1 <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st1 + 11
    positions <- c(st1:ed)
    tmp <- dat_ep %>%
            select(all_of(positions))
    subj_indx <- sub(".*\\.", "", names(tmp)[1])

    colnames(tmp) <- c(nm_rec, nm_don)

    # compare each donor's allele for ALL of recipient's alleles
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

    dat_ep_mm <- cbind(dat_ep_mm, tmp)
    st1 <- ed + 1
  }

  #* step 6: compare mis-match with raw_look up table *#
  dat_ep_mm2 <- dat_ep_mm %>%
    select(index, type)

  # exclude index and type, pulling data starting from the 3rd position
  st2 <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st2 + 5
    positions <- c(st2:ed)
    tmp <- dat_ep_mm %>%
            select(c(1, 2, all_of(positions)))

    ori_name <- colnames(tmp)[-c(1,2)]
    if(i == 1){
      tmp2 <- left_join(raw_lookup, tmp, by = c("index", "type")) %>%
              setNames(c("index", "type", "eplet", "a1", "a2", "b1", "b2", "c1", "c2")) %>%
              mutate(mm = ifelse(eplet == a1 | eplet == a2 | eplet == b1 | eplet == b2 | eplet == c1 | eplet == c2 , eplet, "")) %>%
              filter(mm != "") %>%
              left_join(tmp, ., by = c("index", "type")) %>%
              select(-c(all_of(ori_name), "eplet", "mm")) %>%
              setNames(c("index", "type", ori_name)) %>%
              distinct()
    }

    if (i > 1 ) {
      tmp2 <- left_join(raw_lookup, tmp, by = c("index", "type")) %>%
              setNames(c("index", "type", "eplet", "a1", "a2", "b1", "b2", "c1", "c2")) %>%
              mutate(mm = ifelse(eplet == a1 | eplet == a2 | eplet == b1 | eplet == b2 | eplet == c1 | eplet == c2 , eplet, "")) %>%
              filter(mm != "") %>%
              left_join(tmp, ., by = c("index", "type")) %>%
              select(-c(ori_name, "eplet", "mm")) %>%
              setNames(c("index", "type", ori_name)) %>%
              distinct() %>%
              select(-c("index", "type"))
    }
    st2 <- ed + 1

    dat_ep_mm2 <- cbind(dat_ep_mm2, tmp2)
  }

  dat_ep_mm2 <- dat_ep_mm2 %>%
                select(-c(1:2))

  #* step 7: final result *#
  subj_names <- unique(sub(".*\\_", "", names(dat_ep_mm2)[-c(1:2)]))

  result <- raw_lookup
  st3 <- 3
  for (i in 1:subj_num) {
    ed <- st3 + 5
    positions <- c(st3:ed)
    tmp <- dat_ep_mm2 %>%
            select(c(1, 2, positions)) %>%
            setNames(c("index", "type", "a1_mm", "a2_mm", "b1_mm", "b2_mm", "c1_mm", "c2_mm")) %>%
            right_join(., raw_lookup,by = c("index", "type")) %>%
            mutate(var = ifelse(eplet %in% c(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm), eplet, NA)) %>%
            select(index, type, eplet, var) %>%
            setNames(c("index", "type", "eplet", subj_names[i]))

    st3 <- ed + 1
    result <- left_join(result, tmp, by = c("index", "type", "eplet"))
  }

  result <-  data.frame(t(dat_ep_mm2)) %>%
              tidyr::unite_(., paste(colnames(.), collapse="_"), colnames(.)) %>%
              setNames("mm_eplets") %>%
              mutate(subject = colnames(dat_ep_mm2),
                     mm_eplets = gsub(",NA", "", gsub("_", ",", gsub("NA_", "", mm_eplets)))) %>%
              mutate(mm_cnt = str_count(mm_eplets, ",")) %>%
              mutate(mm_cnt = ifelse(mm_eplets == "NA" & mm_cnt == 0, 0,
                                     ifelse(mm_eplets != "NA" & mm_cnt == 0, 1, mm_cnt + 1))) %>%
              filter(!subject %in% c("index", "type")) %>%
              select(subject, mm_eplets, mm_cnt) %>%
              mutate(part_id = gsub(".*subj", "", subject),
                     gene = gsub("\\_.*","",subject) )

  don_allele <- dat %>%
                select(part_id, don_a1, don_a2, don_b1, don_b2, don_c1, don_c2) %>%
                pivot_longer(cols = starts_with("don_"),
                             names_to = "gene",
                             values_to = "don_type") %>%
                mutate(gene = str_replace(gene, "don_", ""),
                       part_id = as.character(part_id))

  result <- result %>%
            left_join(., don_allele, by =c("part_id", "gene") ) %>%
            select(part_id, subject, don_type, mm_eplets, mm_cnt) %>%
            left_join(., id_match, by = "part_id") %>%
            select(-part_id) %>%
            rename(part_id = part_id_ori) %>%
            select(part_id, everything())

  return(result)
}






