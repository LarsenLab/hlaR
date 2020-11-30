#' many-to-many mis-match calculation, find out best match based on mm-count
#' @param dat_in
#' data table with donor/recipient MHC-I info, each row has unique id
#' @param id_don
#' donor ids
#' @param id_rcpt
#' recipient ids
#' @return
#' list of data tables.
#' count: original input data appended with count mis-matched eplet of each pair
#' detail: detailed mis-match eplet of each subject id, plus count and percentage of mis-match across all pairs
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
#' dat_in <- read.csv("~/projects/DEV/hlaR_dev/MHC_I_n2n.csv", sep = ",", header = TRUE)
#' re <- CalEpletMHCI_best(dat_in = dat_in,
#' id_don = c(8), id_rcpt = c(6,1))
#' }
#'

# library(hlaR)
# library(readr)
# library(tidyverse)
#
# dat_in <- read.csv("~/projects/DEV/hlaR_dev/MHC_I_n2n.csv", sep = ",", header = TRUE)
# id_don <- c(8)
# id_rcpt <- c(6, 1)
#
# dat_in <- dat_in
# id_don <- id_don
# id_rcpt <- id_rcpt
#
# CalEpletMHCI_best(dat_in = dat_in, id_don = id_don, id_rcpt = id_rcpt)

CalEpletMHCI_best <- function(dat_in, id_don, id_rcpt) {

  #* pre-step : check input parameters for overlap donor and recipient id, mis-match donor_type info *#
  don <- dat_in %>% filter(id %in% id_don)
  rcpt <- dat_in %>% filter(id %in% id_rcpt)

  if (id_don %in% id_rcpt)
    stop("donor and recipient ids have overlaps, please check.")

  if(unique(don$donor_type) %in% c("recipient", "recip", "rcpt", "r"))
    stop("donor has type 'recipient', please check you donor id")

  if(unique(rcpt$donor_type) %in% c("donor", "don", "dn","d"))
    stop("recipient has type 'donor', please check your recipient id")
  #* end of pre-step *#

  #* step 1: import raw eplet table *#
  raw_eplet <- read.csv(system.file("extdata", "MHC_I_eplet.csv", package = "hlaR"), check.names = FALSE)

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
}

#   ###########################
#   ###########################
#   ###########################
#
# rcpt <- dat_in %>%
#   filter(id %in% id_rcpt)
#
# don <- dat_in %>%
#   filter(donor_type %in% c("donor", "don", "dn", "d")) %>%
#   select(-donor_type) %>%
#   setNames(c("part_id", "part_type", nm_don))
#
#
# dat <- left_join(rcpt, don, by = c("part_id", "part_type"))
#
#   ###########################
#   ###########################
#   ###########################
#
#   #* step 2: import patient table *#
#   nm_rec <- c("rec_a1", "rec_a2", "rec_b1", "rec_b2", "rec_c1", "rec_c2")
#   nm_don <- c("don_a1", "don_a2", "don_b1", "don_b2", "don_c1", "don_c2")
#
#   rcpt <- dat_in %>%
#     filter(donor_type %in% c("recipient", "recip", "rcpt", "r")) %>%
#     select(-donor_type) %>%
#     setNames(c("part_id", "part_type", nm_rec))
#
#   don <- dat_in %>%
#     filter(donor_type %in% c("donor", "don", "dn", "d")) %>%
#     select(-donor_type) %>%
#     setNames(c("part_id", "part_type", nm_don))
#
#
#   dat <- left_join(rcpt, don, by = c("part_id", "part_type"))
#
#   subj_num <- dim(dat)[1]
#   tmp_names <- c(nm_rec, nm_don)
#
#   #* step 3: initial eplet data frame *#
#   dat_ep <- raw_eplet %>%
#     select(index, type)
#
#   #* step 4: pull out eplet of each allele *#
#   for (i in 1:subj_num) {
#     allele <- toupper(unlist(transpose(dat[i,-c(1:2)]), use.names = F))
#     allele <- ifelse(allele %in% names(raw_eplet), allele, NA)
#
#     for (j in 1:length(allele)) {
#       varname <- paste0(tmp_names[j], ".", sep = i)
#       if (!is.na(allele[j])) {
#         tmp <- raw_eplet %>%
#           select(index, type, allele[j])
#         tmp <- tmp %>%
#           setNames(c("index", "type", varname))
#         dat_ep <- dat_ep %>%
#           left_join(., tmp, by = c("index", "type"))
#       } else{
#         dat_ep <- dat_ep %>%
#           mutate(!!varname := NA)
#       }
#     }
#   }
#
#   #* step 5: mark mis-matches *#
#   dat_ep_mm <- dat_ep %>%
#     select(index, type)
#
#   # exclude index and type, pulling data starting from the 3rd position
#   st1 <- 3
#
#   # for each subject
#   for (i in 1:subj_num) {
#     ed <- st1 + 11
#     positions <- c(st1:ed)
#     tmp <- dat_ep %>%
#       select(all_of(positions))
#     subj_indx <- sub(".*\\.", "", names(tmp)[1])
#
#     colnames(tmp) <- c(nm_rec, nm_don)
#
#     # compare each donor's allele for ALL of recipient's alleles
#     tmp <- tmp %>%
#       mutate(a1_mm = ifelse(don_a1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_a1),
#              a2_mm = ifelse(don_a2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_a2),
#              b1_mm = ifelse(don_b1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_b1),
#              b2_mm = ifelse(don_b2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_b2),
#              c1_mm = ifelse(don_c1 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_c1),
#              c2_mm = ifelse(don_c2 %in% c(rec_a1, rec_a2, rec_b1, rec_b2, rec_c1, rec_c2), NA, don_c2)) %>%
#       select(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm) %>%
#       setNames(c(paste0("a1_mm_subj", subj_indx),
#                  paste0("a2_mm_subj", subj_indx),
#                  paste0("b1_mm_subj", subj_indx),
#                  paste0("b2_mm_subj", subj_indx),
#                  paste0("c1_mm_subj", subj_indx),
#                  paste0("c2_mm_subj", subj_indx)))
#
#     dat_ep_mm <- cbind(dat_ep_mm, tmp)
#     st1 <- ed + 1
#   }
#
#   #* step 6: compare mis-match with raw_look up table *#
#   dat_ep_mm2 <- dat_ep_mm %>%
#     select(index, type)
#
#   # exclude index and type, pulling data starting from the 3rd position
#   st2 <- 3
#
#   # for each subject
#   for (i in 1:subj_num) {
#     ed <- st2 + 5
#     positions <- c(st2:ed)
#     tmp <- dat_ep_mm %>%
#       select(c(1,2,all_of(positions)))
#
#     ori_name <- colnames(tmp)[-c(1,2)]
#     if(i == 1){
#       tmp2 <- left_join(raw_lookup, tmp, by = c("index", "type")) %>%
#         setNames(c("index", "type", "eplet", "a1", "a2", "b1", "b2", "c1", "c2")) %>%
#         mutate(mm = ifelse(eplet == a1 | eplet == a2 | eplet == b1 | eplet == b2 | eplet == c1 | eplet == c2 , eplet, "")) %>%
#         filter(mm != "") %>%
#         left_join(tmp, ., by = c("index", "type")) %>%
#         select(-c(all_of(ori_name), "eplet", "mm")) %>%
#         setNames(c("index", "type", ori_name)) %>%
#         distinct()
#     }
#
#     if (i > 1 ) {
#       tmp2 <- left_join(raw_lookup, tmp, by = c("index", "type")) %>%
#         setNames(c("index", "type", "eplet", "a1", "a2", "b1", "b2", "c1", "c2")) %>%
#         mutate(mm = ifelse(eplet == a1 | eplet == a2 | eplet == b1 | eplet == b2 | eplet == c1 | eplet == c2 , eplet, "")) %>%
#         filter(mm != "") %>%
#         left_join(tmp, ., by = c("index", "type")) %>%
#         select(-c(ori_name, "eplet", "mm")) %>%
#         setNames(c("index", "type", ori_name)) %>%
#         distinct() %>%
#         select(-c("index", "type"))
#     }
#     st2 <- ed + 1
#
#     dat_ep_mm2 <- cbind(dat_ep_mm2, tmp2)
#   }
#
#   dat_ep_mm2 <- dat_ep_mm2 %>% select(-c(1:2))
#
#   #* step 7: final result *#
#   # subject names
#   subj_names <- unique(sub(".*\\_", "", names(dat_ep_mm2)[-c(1:2)]))
#
#   results_detail <- raw_lookup
#   st3 <- 3
#   for (i in 1:subj_num) {
#     ed <- st3 + 5
#     positions <- c(st3:ed)
#     tmp <- dat_ep_mm2 %>%
#       select(c(1, 2, positions)) %>%
#       setNames(c("index", "type", "a1_mm", "a2_mm", "b1_mm", "b2_mm", "c1_mm", "c2_mm")) %>%
#       right_join(., raw_lookup,by = c("index", "type")) %>%
#       mutate(var = ifelse(eplet %in% c(a1_mm, a2_mm, b1_mm, b2_mm, c1_mm, c2_mm), eplet, NA)) %>%
#       select(index, type, eplet, var) %>%
#       setNames(c("index", "type", "eplet", subj_names[i]))
#
#     st3 <- ed + 1
#     results_detail <- left_join(results_detail, tmp, by = c("index", "type", "eplet"))
#   }
#
#   results_detail <- results_detail %>%
#     mutate(mm_cnt = subj_num - rowSums(is.na(.)),
#            mm_pect = paste(round((subj_num - rowSums(is.na(.))) / subj_num * 100, 1), "%", sep = ""))
#
#   results_count <- results_detail %>%
#     select(-c(index, type, eplet, mm_cnt, mm_pect)) %>%
#     summarise_all(list(~sum(!is.na(.)))) %>%
#     gather() %>%
#     setNames(c("subject", "mm_count")) %>%
#     select(mm_count)
#
#   results_count <- cbind(dat ,results_count)
#
#   return(list(count = results_count,
#               detail = results_detail))
# }






