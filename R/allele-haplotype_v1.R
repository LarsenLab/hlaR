#' @title pull out haplotypes based on NMDP frequency table ...
#' @param dat_in
#' dataframe of alleles
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata", "HLA_MisMatch_test.csv", package = "hlaR"))
#' result <- ht(dat_in = dat)
#' }
#' @export
ht <- function(dat_in){
  # discussion
  # 1. which frequency table we are going to use
  # 2. algorithm of filter
  # 3. output? format and topN?
  #   4. genotype typing resolution score?
  #
  #   my steps:
  #   1. data clean(seperate typing by before/after colon for comparison)
  # 2. for each allele,
  # - if it's 4 digits, then pull out row index of rows have the 4 digits code
  #  - if it's 2 digits, then pull out rows index of rows start with the first 2 digits
  # 3. union the indices, subset the raw frequency table
  # 4. keep
  # - all allele combination if non-missing input alleles are subset of the records
  # - allele combination from step3 but filtered by rank of each population
  #
  # ori <- read.csv("~/projects/DEV/haplostats_dev/data/csv/A~C~B~DRB1~DQB1.csv") %>%
  #   select(A, C, B, DRB1, DQB1,
  #          AFA_freq, AFA_rank, API_freq, API_rank,CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank ) %>%
  #   #filter_at(vars(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank), all_vars(!is.na(.))) %>%
  #   mutate(idx = as.numeric(rownames(.))) %>%
  #   # remove trailing g if exists
  #   mutate(A = ifelse(stri_sub(A, -1) == "g", str_sub(A, end = -2), A),
  #          B = ifelse(stri_sub(B, -1) == "g", str_sub(B, end = -2), B),
  #          C = ifelse(stri_sub(C, -1) == "g", str_sub(C, end = -2), C),
  #          DRB1 = ifelse(stri_sub(DRB1, -1) == "g", str_sub(DRB1, end = -2), DRB1),
  #          DQB1 = ifelse(stri_sub(DQB1, -1) == "g", str_sub(DQB1, end = -2), DQB1)) %>%
  #   # remove leading loci letter
  #   mutate(A = str_sub(A, start = 3),
  #          B = str_sub(B, start = 3),
  #          C = str_sub(C, start = 3),
  #          DRB1 = str_sub(DRB1, start = 6),
  #          DQB1 = str_sub(DQB1, start = 6)) %>%
  #   # split allele by :
  #   mutate(fstA = sub("\\:.*", "", A),
  #          lstA = sub(".*\\:", "", A),
  #          fstB = sub("\\:.*", "", B),
  #          lstB = sub(".*\\:", "", B),
  #          fstC = sub("\\:.*", "", C),
  #          lstC = sub(".*\\:", "", C),
  #          fstDRB1 = sub("\\:.*", "", DRB1),
  #          lstDRB1 = sub(".*\\:", "", DRB1),
  #          fstDQB1 = sub("\\:.*", "", DQB1),
  #          lstDQB1 = sub(".*\\:", "", DQB1))
  #
  # dat <- read.csv("~/projects/DEV/haplostats_dev/data/csv/haplotype_test.csv") %>%
  #   # remove leading loci letter
  #   mutate(A = str_sub(A, start = 3),
  #          B = str_sub(B, start = 3),
  #          C = str_sub(C, start = 3),
  #          DRB1 = str_sub(DRB1, start = 6),
  #          DQB1 = str_sub(DQB1, start = 6)) %>%
  #   # mutate(ck = ifelse(grepl("\\D", sub(".*\\:", "", A)), "", sub(".*\\:", "", A))) %>%
  #   # if last 2 chars are not number, then remove it so only keep first 2 digits of the allele kept
  #   mutate(A = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", A)), "", sub(".*\\:", "", A)) != "", A, str_sub(A, end = -3)),
  #          B = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", B)), "", sub(".*\\:", "", B)) != "", B, str_sub(B, end = -3)),
  #          C = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", C)), "", sub(".*\\:", "", C)) != "", C, str_sub(C, end = -3)),
  #          DRB1 = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", DRB1)), "", sub(".*\\:", "", DRB1)) != "", DRB1, str_sub(DRB1, end = -3)),
  #          DQB1 = ifelse(ifelse(grepl("\\D", sub(".*\\:", "", DQB1)), "", sub(".*\\:", "", DQB1)) != "", DQB1, str_sub(DQB1, end = -3))) %>%
  #   mutate(fstA = sub("\\:.*", "", A),
  #          lstA = sub(".*\\:", "", A),
  #          fstB = sub("\\:.*", "", B),
  #          lstB = sub(".*\\:", "", B),
  #          fstC = sub("\\:.*", "", C),
  #          lstC = sub(".*\\:", "", C),
  #          fstDRB1 = sub("\\:.*", "", DRB1),
  #          lstDRB1 = sub(".*\\:", "", DRB1),
  #          fstDQB1 = sub("\\:.*", "", DQB1),
  #          lstDQB1 = sub(".*\\:", "", DQB1)) %>%
  #   mutate_all(., list(~na_if(.,""))) %>%
  #   mutate(flagA = ifelse(!is.na(fstA) & !is.na(lstA), 1, 0),
  #          flagB = ifelse(!is.na(fstB) & !is.na(lstB), 1, 0),
  #          flagC = ifelse(!is.na(fstC) & !is.na(lstC), 1, 0),
  #          flagDRB = ifelse(!is.na(fstDRB1) & !is.na(lstDRB1), 1, 0),
  #          flagDQB = ifelse(!is.na(fstDQB1) & !is.na(fstDQB1), 1, 0)) %>%
  #   select(-c(A, B, C, DRB1, DQB1))
  #
  # # combination of !is.na inputs - 11/09/2020
  # tmpFunc <- function(dat_in) {
  #   tmp_indx <- c()
  #   if (dat_in$flagA == 1){
  #     tmp <- ori %>% filter(fstA == dat_in$fstA & lstA == dat_in$lstA) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstA == dat_in$fstA) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     tmp_indx <- c(tmp_indx, tmp)
  #
  #   if (dat_in$flagB == 1) {
  #     tmp <- ori %>% filter(fstB == dat_in$fstB & lstB == dat_in$lstB) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstB == dat_in$fstB) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     tmp_indx <- c(tmp_indx, tmp)
  #
  #   if (dat_in$flagC == 1) {
  #     tmp <- ori %>% filter(fstC == dat_in$fstC & lstC == dat_in$lstC) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstC == dat_in$fstC) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     tmp_indx <- c(tmp_indx, tmp)
  #
  #   if (dat_in$flagDRB == 1) {
  #     tmp <- ori %>% filter(fstDRB1 == dat_in$fstDRB1 & lstDRB1 == dat_in$lstDRB1) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstDRB1 == dat_in$fstDRB1) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     tmp_indx <- c(tmp_indx, tmp)
  #
  #   if (dat_in$flagDQB == 1) {
  #     tmp <- ori %>% filter(fstDQB1 == dat_in$fstDQB1 & lstDQB1 == dat_in$lstDQB1) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstDQB1 == dat_in$fstDQB1) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     tmp_indx <- c(tmp_indx, tmp)
  #
  #   tmp_indx <- unique(tmp_indx)
  #
  #   subdat <- ori %>%
  #     filter(idx %in% tmp_indx)  %>%
  #     # rank w/o NA
  #     # filter(rank(AFA_rank) <= 20 | rank(API_rank) <= 20 | rank(CAU_rank) <= 20 | rank(HIS_rank) <= 20  | rank(NAM_rank) <= 20) %>%
  #     select(A, C, B, DRB1, DQB1,
  #            fstA, lstA, fstB, lstB, fstC, lstC, fstDRB1, lstDRB1, fstDQB1, lstDQB1,
  #            AFA_freq, AFA_rank, API_freq, API_rank, CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank)
  #
  #   return(subdat)
  # }
  #
  # n <- dim(dat)[1]
  # l <- vector(mode="list", length=n)
  # for (i in 1:n){
  #   l[[i]] <- tmpFunc(dat[i,])
  # }
  #
  # re1 <- l[[1]]
  # re2 <- l[[2]]
  #
  # final1 <- tempfunc(datin1 = dat[1,], datin2 = re1)
  # final2 <- tempfunc(datin1 = dat[2,], datin2 = re2)
  #
  # final <- rbind(final1, final2)
  #
  # tempfunc <- function(datin1, datin2){
  #
  #   # pull out alleles info if column is not empty
  #   colnms <- datin1 %>%
  #     select(id, fstA, lstA, fstB, lstB, fstC, lstC, fstDRB1, lstDRB1, fstDQB1, lstDQB1) %>%
  #     mutate(across(c(fstA, lstA, fstB, lstB, fstC, lstC, fstDRB1, lstDRB1, fstDQB1, lstDQB1), is.na)) %>%  # replace all NA with TRUE and else FALSE
  #     pivot_longer(-id, names_to = "var") %>%  # pivot longer
  #     filter(value != 1) %>%
  #     pull(var)
  #
  #   re_w_alleles <- datin2 %>% left_join(., datin1, by = colnms) %>%
  #     filter(!is.na(id)) %>%
  #     # filter_at(vars(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank), all_vars(!is.na(.))) %>%
  #     filter_at(vars(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank), any_vars(!is.na(.))) %>%
  #     select(c(A, C, B, DRB1, DQB1, AFA_freq, AFA_rank, API_freq, API_rank, CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank))
  #
  #   datin2 <- datin2 %>%
  #     mutate(flag = ifelse(AFA_rank %in% re_w_alleles$AFA_rank[!is.na(re_w_alleles$AFA_rank)] |
  #                            API_rank %in% re_w_alleles$API_rank[!is.na(re_w_alleles$API_rank)] |
  #                            CAU_rank %in% re_w_alleles$CAU_rank[!is.na(re_w_alleles$CAU_rank)] |
  #                            HIS_rank %in% re_w_alleles$HIS_rank[!is.na(re_w_alleles$HIS_rank)] |
  #                            NAM_rank %in% re_w_alleles$NAM_rank[!is.na(re_w_alleles$NAM_rank)], "yes", "no"))
  #
  #   re_full <- datin2 %>% filter(flag == "yes") %>%
  #     select(c(A, C, B, DRB1, DQB1, AFA_freq, AFA_rank, API_freq, API_rank, CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank)) %>%
  #     filter_at(vars(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank), any_vars(!is.na(.))) %>%
  #     gather(., pop, poprank, c(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank)) %>%
  #     filter(!is.na(poprank)) %>%
  #     arrange(poprank) %>%
  #     group_by(pop) %>%
  #     mutate(rank = order(poprank)) %>%
  #     ungroup() %>%
  #     mutate(note = "full")
  #
  #   re_comb <- datin2 %>% filter(flag == "no") %>%
  #     select(c(A, C, B, DRB1, DQB1, AFA_freq, AFA_rank, API_freq, API_rank, CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank)) %>%
  #     filter_at(vars(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank), any_vars(!is.na(.))) %>%
  #     gather(., pop, poprank, c(AFA_rank, API_rank, CAU_rank, HIS_rank, NAM_rank)) %>%
  #     filter(!is.na(poprank)) %>%
  #     arrange(poprank) %>%
  #     group_by(pop) %>%
  #     mutate(rank = order(poprank)) %>%
  #     ungroup() %>%
  #     filter(rank <= 150) %>%
  #     mutate(note = "comb")
  #
  #   re <- rbind(re_full, re_comb)
  #   return(re)
  # }
}

