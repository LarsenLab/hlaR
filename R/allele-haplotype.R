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
  #   select(-c(id, A, B, C, DRB1, DQB1))
  #
  # # dat %>%
  # #   mutate(across(c(Afst, Alst, Bfst, Blst, Cfst, Clst, DRB1fst, DRB1lst, DQB1fst, DQB1lst), is.na)) %>%  # replace all NA with TRUE and else FALSE
  # #   pivot_longer(-id, names_to = "var") %>%  # pivot longer
  # #   filter(value) %>%   # remove the FALSE rows
  # #   group_by(id) %>%    # group by the ID
  # #   summarise(`Missing Variables` = toString(var)) # convert the variable names to a string column
  #
  # tmpFunc <- function(dat_in) {
  #   idx = c()
  #   if (dat_in$flagA == 1){
  #     tmp <- ori %>% filter(fstA == dat_in$fstA & lstA == dat_in$lstA) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstA == dat_in$fstA) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     idx <- c(idx, tmp)
  #
  #   if (dat_in$flagB == 1) {
  #     tmp <- ori %>% filter(fstB == dat_in$fstB & lstB == dat_in$lstB) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstB == dat_in$fstB) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     idx <- c(idx, tmp)
  #
  #   if (dat_in$flagC == 1) {
  #     tmp <- ori %>% filter(fstC == dat_in$fstC & lstC == dat_in$lstC) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstC == dat_in$fstC) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     idx <- c(idx, tmp)
  #
  #   if (dat_in$flagDRB == 1) {
  #     tmp <- ori %>% filter(fstDRB1 == dat_in$fstDRB1 & lstDRB1 == dat_in$lstDRB1) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstDRB1 == dat_in$fstDRB1) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     idx <- c(idx, tmp)
  #
  #   if (dat_in$flagDQB == 1) {
  #     tmp <- ori %>% filter(fstDQB1 == dat_in$fstDQB1 & lstDQB1 == dat_in$lstDQB1) %>% pull(idx)
  #   } else {
  #     tmp <- ori %>% filter(fstDQB1 == dat_in$fstDQB1) %>% pull(idx)
  #   }
  #
  #   if (length(tmp > 0))
  #     idx <- c(idx, tmp)
  #
  #   idx <- unique(idx)
  #
  #   subdat <- ori %>%
  #     filter(idx %in% idx)  %>%
  #     # rank w/o NA
  #     # filter(rank(AFA_rank) <= 20 | rank(API_rank) <= 20 | rank(CAU_rank) <= 20 | rank(HIS_rank) <= 20  | rank(NAM_rank) <= 20) %>%
  #     select(A, C, B, DRB1, DQB1, AFA_freq, AFA_rank, API_freq, API_rank, CAU_freq, CAU_rank, HIS_freq, HIS_rank, NAM_freq, NAM_rank)
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
  # ck1 <- l[[1]]
  # ck2 <- l[[2]]
}

