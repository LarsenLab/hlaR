#' Title
#' function to calculate eplet mis-match of HLA-ABC
#' @param dat_in folder name.
#' input data set with or without complete allele info
#' dat_in has to have 15 columns (first 3 columns are record id, recipient id, donor id; the rest of columns are A/B/C alleles for recipient and donor)
#' @return
#' list of data tables.
#' results_count: original input data appended with count mis-matched eplet of each pair
#' results_detail: detailed mis-match eplet of each subject id, plus count and percentage of mis-match across all pairs
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
#' re <- CalEpletMHCI(dat_in = "inst/extdata/MHC_I_test.csv")
#' }
#'
# notes:
# 1. raw_eplet is pulled from "Ep" worksheet from ABC_Eplet_Matching_3.1.xlsb (password of protected sheets: hla) from http://www.epitopes.net/index.html
# 2. raw_lookup is generated on-the-fly based on raw_eplet table
# 3. result generated from this function may slightly different with result from ABC_Eplet_Matching_3.1.xlsb due to a few possible formatting issue in the excel
# 4. check result: re$results_count or re$results_count
# Version: V1 - 10/20/2020

CalEpletMHCI <- function(dat_in) {

  #* step 1: import reference tables *#
  # dev
  # raw_eplet <- vroom("inst/extdata/MHC_I_eplet.csv")

  # prod
  raw_eplet <- vroom(system.file("extdata", "MHC_I_eplet.csv", package = "hlaR"))
  #raw_eplet <- vroom("inst/extdata/MHC_I_eplet.csv")

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
  dat <- vroom(dat_in)
  subj_num <- dim(dat)[1]
  tmp_names <- names(dat[-c(1:3)])

  #* step 3: initial eplet data frame *#
  dat_ep <- raw_eplet %>%
            select(index, type)

  #* step 4 : pull out eplet of each allele *#
  for (i in 1:subj_num) {
    allele <- unlist(transpose(dat[i,-c(1:3)]), use.names = F)

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
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 11
    positions <- c(st:ed)
    tmp <- dat_ep %>%
      select(positions)
    subj_indx <- sub(".*\\.", "", names(tmp)[i])

    colnames(tmp) <- c("rec_a1", "rec_a2", "rec_b1", "rec_b2", "rec_c1", "rec_c2",
                       "don_a1", "don_a2", "don_b1", "don_b2", "don_c1", "don_c2")

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
    st <- ed + 1
  }

  #* step 6: compare mis-match with raw_look up table *#
  dat_ep_mm2 <- dat_ep_mm %>%
    select(index, type)

  # exclude index and type, pulling data starting from the 3rd position
  st <- 3

  # for each subject
  for (i in 1:subj_num) {
    ed <- st + 5
    positions <- c(st:ed)
    tmp <- dat_ep_mm %>%
      select(c(1,2,positions) )

    ori_name <- colnames(tmp)[-c(1,2)]
    if(i == 1){
      tmp2 <- left_join(raw_lookup, tmp, by = c("index", "type")) %>%
        setNames(c("index", "type", "eplet", "a1", "a2", "b1", "b2", "c1", "c2")) %>%
        mutate(mm = ifelse(eplet == a1 | eplet == a2 | eplet == b1 | eplet == b2 | eplet == c1 | eplet == c2 , eplet, "")) %>%
        filter(mm != "") %>%
        left_join(tmp, ., by = c("index", "type")) %>%
        select(-c(ori_name, "eplet", "mm")) %>%
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
    st <- ed + 1

    dat_ep_mm2 <- cbind(dat_ep_mm2, tmp2)
  }

  dat_ep_mm2 <- dat_ep_mm2 %>% select(-c(1:2))

  #* step 7: final result *#
  # subject names
  subj_names <- unique(sub(".*\\_", "", names(dat_ep_mm2)[-c(1:2)]))

  results_detail <- raw_lookup
  st <- 3
  for (i in 1:subj_num) {
    ed <- st + 5
    positions <- c(st:ed)
    tmp <- dat_ep_mm2 %>%
      select(c(1, 2, positions)) %>%
      setNames(c("index", "type", "a1_mm", "a2_mm", "b1_mm", "b2_mm", "c1_mm", "c2_mm")) %>%
      right_join(., raw_lookup,by = c("index", "type")) %>%
      mutate(var = ifelse(eplet %in% c(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm), eplet, NA)) %>%
      select(index, type, eplet, var) %>%
      setNames(c("index", "type", "eplet", subj_names[i]))

    st <- ed + 1
    results_detail <- left_join(results_detail, tmp, by = c("index", "type", "eplet"))
  }

  results_detail <- results_detail %>%
    mutate(mm_cnt = subj_num - rowSums(is.na(.)),
           mm_pect = paste(round((subj_num - rowSums(is.na(.))) / subj_num * 100, 1), "%", sep = ""))

  # excel bugs
  # type 1 bugs : columns with trailing space for mis-miss-match
  # example: workSheet "Results", row4 colDQ 77D, 90D row4 colBD 90D, colEG 152A, colDS 77S, colEX 184H
  # type 2 bugs: eplet names don't match between "Ep" and "Restuls" tables
  # example: eplet name of Abv21 is "150AHA" or "151AHA" ?
  # Abv21 in "Ep" table is Abv21 is 150AHA, it's 151AHA iin header of "Results"

  results_count <- results_detail %>%
    select(-c(index, type, eplet, mm_cnt, mm_pect)) %>%
    summarise_all(funs(sum(!is.na(.)))) %>%
    gather() %>%
    setNames(c("subject", "mm_count")) %>%
    select(mm_count)

  results_count <- cbind(dat ,results_count)

  return(list(results_count = results_count,
              results_detail = results_detail))
}






