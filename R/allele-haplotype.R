#' @name CompHaploTbl
#' @title feed low resolution hla typing to NMDP frequency table to get best pairs of haplotype combination based on max count of unique low res in the combination
#' @param dat_in
#' dataframe with recipient/donor alleles info
#' @return
#' dataframe of best pairs of haplotype combination
#' @import
#' tidyverse
#'
#' @examples
#' \dontrun{
# dat <- read_csv(system.file("extdata/example", "Haplotype_test.csv", package = "hlaR"))
#' result <- CompHaploTbl(dat_in = dat)
#' }
#' @export

CompHaploTbl <- function(dat_in){
  #* step 1: import raw haplotype frequency table and do a brief cleaning *#
  raw_hap_tbl <- read.csv(system.file("extdata/ref", "A_C_B_DRB345_DRB1_DQB1.csv", package = "hlaR"), check.names = FALSE) %>%
    rename_all(. %>% tolower) %>%
    select(a, c, b, drb1, dqb1, drb345,
           afa_freq, afa_rank, api_freq, api_rank, cau_freq, cau_rank, his_freq, his_rank, nam_freq, nam_rank) %>%
    mutate(idx = as.numeric(rownames(.))) %>%
    # remove trailing g
    mutate(a = ifelse(str_detect(a, "g"), str_replace(a, "g", ""), a),
           b = ifelse(str_detect(b, "g"), str_replace(b, "g", ""), b),
           c = ifelse(str_detect(c, "g"), str_replace(c, "g", ""), c),
           drb1 = ifelse(str_detect(drb1, "g"), str_replace(drb1, "g", ""), drb1),
           dqb1 = ifelse(str_detect(dqb1, "g"), str_replace(dqb1, "g", ""), dqb1),
           drb345 = ifelse(str_detect(drb345, "g"), str_replace(drb345, "g", ""), drb345)) %>%
    # remove leading loci letter
    mutate(a = str_sub(a, start = 3),
           b = str_sub(b, start = 3),
           c = str_sub(c, start = 3),
           drb1 = str_sub(drb1, start = 6),
           dqb1 = str_sub(dqb1, start = 6),
           drb345 = str_sub(drb345, start = 6)) %>%
    # split allele by ":" for comparison with input data
    mutate(fst_a = sub("\\:.*", "", a),
           lst_a = sub(".*\\:", "", a),
           fst_b = sub("\\:.*", "", b),
           lst_b = sub(".*\\:", "", b),
           fst_c = sub("\\:.*", "", c),
           lst_c = sub(".*\\:", "", c),
           fst_drb1 = sub("\\:.*", "", drb1),
           lst_drb1 = sub(".*\\:", "", drb1),
           fst_dqb1 = sub("\\:.*", "", dqb1),
           lst_dqb1 = sub(".*\\:", "", dqb1),
           fst_drb345 = sub("\\:.*", "", drb345),
           lst_drb345 = sub(".*\\:", "", drb345))
  #* end of step 1 *#

  #* step 2: reshape input data table by recipient and donor *#
  dat_ready <- dat_in %>% arrange(rowid) %>%
                mutate_all(as.character) %>%
                mutate(a1 = ifelse(nchar(a1) == 1, paste0("0", a1), a1),
                       a2 = ifelse(nchar(a2) == 1, paste0("0", a2), a2),
                       b1 = ifelse(nchar(b1) == 1, paste0("0", b1), b1),
                       b2 = ifelse(nchar(b2) == 1, paste0("0", b2), b2),
                       c1 = ifelse(nchar(c1) == 1, paste0("0", c1), c1),
                       c2 = ifelse(nchar(c2) == 1, paste0("0", c2), c2),
                       drb1 = ifelse(nchar(drb1) == 1, paste0("0", drb1), drb1),
                       drb2 = ifelse(nchar(drb2) == 1, paste0("0", drb2), drb2),
                       dqb1 = ifelse(nchar(dqb1) == 1, paste0("0", dqb1), dqb1),
                       dqb2 = ifelse(nchar(dqb2) == 1, paste0("0", dqb2), dqb2),
                       drb31 = ifelse(nchar(drb31) == 1, paste0("0", drb31), drb31),
                       drb32 = ifelse(nchar(drb32) == 1, paste0("0", drb32), drb32),
                       drb41 = ifelse(nchar(drb41) == 1, paste0("0", drb41), drb41),
                       drb42 = ifelse(nchar(drb42) == 1, paste0("0", drb42), drb42),
                       drb51 = ifelse(nchar(drb51) == 1, paste0("0", drb51), drb51),
                       drb52 = ifelse(nchar(drb52) == 1, paste0("0", drb52), drb52))
  #* end of step 2 *#

  #* step 3: call FuncForCompHaplo() for each subjects and get top pairs haplotype combination *#
  num_subj <- dim(dat_ready)[1]
  hpl_tp_raw <- vector(mode = "list", length = num_subj)
  hpl_tp_pairs <- vector(mode = "list", length = num_subj)

  for (i in 1:num_subj){
    print(paste("subject", i))
    hpl_tp_pairs[[i]] <- FuncForCompHaplo_v2(tbl_raw = raw_hap_tbl, tbl_in = dat_ready[i, ])
  }

  #* end of step 3 *#
  names(hpl_tp_pairs) <- dat_ready$rowid

  hpl_tp_pairs <- as.data.frame(do.call(rbind, hpl_tp_pairs))
  row.names(hpl_tp_pairs) <- seq(1:dim(hpl_tp_pairs)[1])

  return(hpl_tp_pairs)
}

