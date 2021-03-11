#' #' @keywords internal
#'
#' # purpose: this script is to construct the MHC-II eplet reference table based on matchmaker v2.0, plus some missing known eplets froom eplet Db (33N, 78Y, 78VS(98ES), 108P, 85V)
#' # uncomment to run
#' # date: Match, 2021
#'
#' library(tidyverse)
#' library(here)
#' library(readxlsb)
#' library(openxlsx)
#' library(qdapRegex)
#' library(janitor)
#'
#' # notes:
#' # 1. eplet names:
#' #    - all converted to upper cases (DQA1*04:02 and DQA1*04:04 were 187a, DQA1*05:01 and DQA1*01:02 were 187A. all of them are 187A in current version)
#' # 2. add a list of updated eplets from epregistry.com
#' #    - antibody reactive: 78VS(==98ES)
#' #    - non-antibody reactive: 33N, 78Y, 108P, 85V
#' # 3. reason of add 78VS/33N/78Y/108P/85V to the reference table: these 5 eplets are not in any of Q/R/P tables in the matchmaker v2.0,
#' #    but they are in mm2.0 Result table, and also can be found from epregistry.com
#' #    adding them to the reference table to make it more complete.
#'
#' #* function to pull out eplets *#
#' PullEplet <- function(filename, sheetname){
#'   tbl_ref <-  read_xlsb(filename,
#'                         sheet = sheetname, skip = 0, trim_ws = TRUE, col_names = TRUE) %>%
#'               mutate_all(., .funs = toupper) %>% # convert all eplet names to upper cases
#'               select(-c(names(.)[str_detect(names(.), "X.")])) %>%
#'               replace(., is.na(.), "") %>%
#'               mutate(row_na = rowSums(. == "" )) %>% # delete the row if it have all NAs
#'               filter(row_na < dim(.)[2] - 1) %>%
#'               select(-row_na) %>%
#'               setNames(paste(names(.), .[1,], sep = "_")) %>%
#'               setNames(rm_between(names(.), ".", "_")) %>%
#'               setNames(str_remove(names(.), "_")) %>%
#'               dplyr::rename(allele = "X") %>%
#'               filter(allele != "") %>%
#'               mutate(across(where(is.character), str_trim))
#'
#'   nm_ex <- tbl_ref %>%
#'             summarise_all(~list(sum(. == "")))
#'   nm_ex <- which(unlist(nm_ex) == dim(tbl_ref)[1])
#'   out <- tbl_ref %>%
#'           select(-all_of(nm_ex)) %>%
#'           mutate_if(is.character, str_trim) # remove leading and trailing white spaces
#'   return(out)
#' }
#'
#' #* pull out eplet table for each allele *#
#' q <- PullEplet(filename = here("inst/extdata/matchmaker/5DRDQDPMatchingVs2.xlsb"),
#'                sheetname = "Q")
#'
#' p <- PullEplet(filename = here("inst/extdata/matchmaker/5DRDQDPMatchingVs2.xlsb"),
#'                sheetname = "P")
#'
#' r <- PullEplet(filename = here("inst/extdata/matchmaker/5DRDQDPMatchingVs2.xlsb"),
#'                sheetname = "R") %>%
#'   mutate(count = rowSums(. == "" )) %>%
#'   filter(count < ncol(.) - 2) %>%
#'   select(-count)
#'
#' #* As *#
#' # loci -> count of non-na -> filter if all eplets are na -> select allele and eplets
#' dqa <- q %>%
#'         filter(str_detect(allele, "DQA")) %>%
#'         select(grep("allele|QA", names(.))) %>%
#'         mutate(count = rowSums(. == "" )) %>%
#'         filter(count < ncol(.) - 2) %>%
#'         select(-count)
#'
#' dqa <- as_tibble(t(dqa), rownames = "type") %>%
#'   row_to_names(row_number = 1) %>%
#'   mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
#'   dplyr::rename(type = allele) %>%
#'   group_by(type) %>%
#'   mutate(index = row_number()) %>%
#'   select(index, type, everything())
#'
#' # dpa
#' dpa <- p %>%
#'         filter(str_detect(allele, "DPA")) %>%
#'         select(grep("allele|PA", names(.))) %>%
#'         mutate(count = rowSums(. == "" )) %>%
#'         filter(count < ncol(.) - 2) %>%
#'         select(-count)
#'
#' dpa <- as_tibble(t(dpa), rownames = "type") %>%
#'   row_to_names(row_number = 1) %>%
#'   mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
#'   dplyr::rename(type = allele) %>%
#'   group_by(type) %>%
#'   mutate(index = row_number()) %>%
#'   select(index, type, everything())
#'
#' # dqa and dpa
#' mmv2_a <- dqa %>%
#'   full_join(., dpa, by = c("index", "type")) %>%
#'   replace(., is.na(.), "") %>%
#'   arrange(type, index)
#'
#' # write.csv(a, "inst/extdata/matchmaker/check/mhcII_v2_A.csv", row.names = FALSE)
#'
#' #* Bs *#
#' # dqb
#' dqb <- q %>%
#'   filter(str_detect(allele, "DQB")) %>%
#'   select(grep("allele|QB", names(.))) %>%
#'   mutate(count = rowSums(. == "" )) %>%
#'   filter(count < ncol(.) - 2) %>%
#'   select(-count)
#'
#' dqb <- as_tibble(t(dqb), rownames = "type") %>%
#'   row_to_names(row_number = 1) %>%
#'   mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
#'   dplyr::rename(type = allele) %>%
#'   group_by(type) %>%
#'   mutate(index = row_number()) %>%
#'   select(index, type, everything())
#'
#' # dpb
#' dpb <- p %>%
#'   filter(str_detect(allele, "DPB")) %>%
#'   select(grep("allele|PB", names(.))) %>%
#'   mutate(count = rowSums(. == "" )) %>%
#'   filter(count < ncol(.) - 2) %>%
#'   select(-count)
#'
#' dpb <- as_tibble(t(dpb), rownames = "type") %>%
#'   row_to_names(row_number = 1) %>%
#'   mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
#'   dplyr::rename(type = allele) %>%
#'   group_by(type) %>%
#'   mutate(index = row_number()) %>%
#'   select(index, type, everything())
#'
#' # r
#' r <- as_tibble(t(r), rownames = "type") %>%
#'   row_to_names(row_number = 1) %>%
#'   mutate(allele = ifelse(str_detect(allele, "av"), "AbV", "oth")) %>%
#'   dplyr::rename(type = allele) %>%
#'   group_by(type) %>%
#'   mutate(index = row_number()) %>%
#'   select(index, type, everything())
#'
#' # r + dqb + dpb
#'  mmv2_b <- r %>%
#'       full_join(., dqb, by = c("index", "type")) %>%
#'       arrange(type, index) %>%
#'       full_join(., dpb, by = c("index", "type")) %>%
#'       select(-`DQB1*X`) %>% # no such a allele, remove it
#'       replace(., is.na(.), "") %>%
#'       arrange(type, index)
#'
#' # write.csv(mmv2_b, "inst/extdata/matchmaker/check/mhcII_v2_B.csv", row.names = FALSE)
#' #* mmv2_a and mmv2_b have been validated, code above are stable *#
#'
#' #* add missing eplets to the table *#
#' #* eplet lists
#' `33N` <- c("DRB1*01:01", "DRB1*01:02", "DRB1*01:03", "DRB1*03:01", "DRB1*03:02", "DRB1*03:03", "DRB1*07:01", "DRB1*08:01",
#'  "DRB1*08:02", "DRB1*09:01", "DRB1*09:02", "DRB1*10:01", "DRB1*11:01", "DRB1*11:03", "DRB1*11:04", "DRB1*12:01",
#'  "DRB1*12:02", "DRB1*13:01", "DRB1*13:02", "DRB1*13:03", "DRB1*13:05", "DRB1*14:01", "DRB1*14:02", "DRB1*14:03",
#'  "DRB1*14:04", "DRB1*14:54", "DRB1*15:01", "DRB1*15:02", "DRB1*15:03", "DRB1*16:01", "DRB1*16:02", "DRB3*01:01",
#'  "DRB3*02:01", "DRB3*02:02", "DRB3*03:01", "DRB4*01:01", "DRB4*01:03", "DRB5*01:01", "DRB5*02:02")
#'
#' `78Y` <- c("DRB1*01:01", "DRB1*01:02", "DRB1*01:03", "DRB1*03:01", "DRB1*03:02", "DRB1*03:03", "DRB1*04:01", "DRB1*04:02",
#'            "DRB1*04:03", "DRB1*04:04", "DRB1*04:05", "DRB1*08:01", "DRB1*08:02", "DRB1*10:01", "DRB1*11:01", "DRB1*11:03",
#'            "DRB1*11:04", "DRB1*12:01", "DRB1*12:02", "DRB1*13:01", "DRB1*13:02", "DRB1*13:03", "DRB1*13:05", "DRB1*14:01",
#'            "DRB1*14:02", "DRB1*14:03", "DRB1*14:04", "DRB1*14:54", "DRB1*15:01", "DRB1*15:02", "DRB1*15:03", "DRB1*16:01",
#'            "DRB1*16:02", "DRB3*01:01", "DRB3*02:01", "DRB3*02:02", "DRB3*03:01", "DRB4*01:01", "DRB4*01:03", "DRB5*01:01",
#'            "DRB5*02:02")
#'
#' `108P` <- c("DRB1*01:01", "DRB1*01:02", "DRB1*01:03", "DRB1*03:01", "DRB1*03:02", "DRB1*03:03", "DRB1*04:01", "DRB1*04:02",
#'             "DRB1*04:03", "DRB1*04:04", "DRB1*04:05", "DRB1*07:01", "DRB1*08:01", "DRB1*08:02", "DRB1*09:01", "DRB1*09:02",
#'             "DRB1*10:01", "DRB1*11:01", "DRB1*11:03", "DRB1*11:04", "DRB1*12:01", "DRB1*12:02", "DRB1*13:01", "DRB1*13:02",
#'             "DRB1*13:03", "DRB1*13:05", "DRB1*14:01", "DRB1*14:02", "DRB1*14:03", "DRB1*14:04", "DRB1*14:54", "DRB1*15:01",
#'             "DRB1*15:02", "DRB1*15:03", "DRB1*16:01", "DRB1*16:02", "DRB3*01:01", "DRB3*02:01", "DRB3*02:02", "DRB3*03:01",
#'             "DRB4*01:01", "DRB4*01:03")
#'
#' `85V` <- c("DRB1*01:01", "DRB1*01:03", "DRB1*03:01", "DRB1*03:02", "DRB1*03:03", "DRB1*04:01", "DRB1*04:02", "DRB1*04:03",
#'            "DRB1*04:04", "DRB1*04:05", "DRB1*07:01", "DRB1*08:01", "DRB1*08:02", "DRB1*09:01", "DRB1*09:02", "DRB1*10:01",
#'            "DRB1*11:01", "DRB1*11:03", "DRB1*11:04", "DRB1*13:01", "DRB1*13:02", "DRB1*13:03", "DRB1*13:05", "DRB1*14:01",
#'            "DRB1*14:02", "DRB1*14:03", "DRB1*14:04", "DRB1*14:54", "DRB1*15:01", "DRB1*15:02", "DRB1*15:03", "DRB1*16:01",
#'            "DRB1*16:02", "DRB3*01:01", "DRB3*02:01", "DRB3*02:02", "DRB3*03:01", "DRB4*01:01", "DRB4*01:03", "DRB5*01:01")
#'
#' `78VS` <- c("DRB1*07:01", "DRB1*09:01")
#'
#' # 33N - other, non-antibody reactive, index max_oth, type oth
#'  max <- mmv2_b %>% group_by(type) %>% mutate(max = max(index)) %>% ungroup() %>% select(type, max) %>% distinct()
#'  max_abv <- max %>% filter(type == "AbV") %>% pull(max)
#'  max_oth <- max %>% filter(type == "oth") %>% pull(max)
#'
#' tmp <- data.frame(names(mmv2_b)) %>%
#'         setNames("name") %>%
#'         mutate(add = ifelse(name %in% c("index"), max_oth,
#'                             ifelse(name %in% c("type"), "oth",
#'                                    ifelse(name %in% `33N`, "33N", "")))) %>%
#'         t() %>%
#'         as.data.frame() %>%
#'         janitor::row_to_names(1) %>%
#'         mutate(index = as.numeric(index) + 1)
#'
#' mmv2_b <- rbind(mmv2_b, tmp)
#'
#' rm(max, max_abv, max_oth, tmp)
#'
#' # 78Y - other, non-antibody reactive, index max_oth, type oth
#' # find max numbers
#' max <- mmv2_b %>% group_by(type) %>% mutate(max = max(index)) %>% ungroup() %>% select(type, max) %>% distinct()
#' max_abv <- max %>% filter(type == "AbV") %>% pull(max)
#' max_oth <- max %>% filter(type == "oth") %>% pull(max)
#'
#' tmp <- data.frame(names(mmv2_b)) %>%
#'   setNames("name") %>%
#'   mutate(add = ifelse(name %in% c("index"), max_oth,
#'                       ifelse(name %in% c("type"), "oth",
#'                              ifelse(name %in% `78Y`, "78Y", "")))) %>%
#'   t() %>%
#'   as.data.frame() %>%
#'   janitor::row_to_names(1) %>%
#'   mutate(index = as.numeric(index) + 1)
#'
#' mmv2_b <- rbind(mmv2_b, tmp)
#' rm(max, max_abv, max_oth, tmp)
#'
#' # 108P - other, non-antibody reactive, index max_oth, type oth
#' max <- mmv2_b %>% group_by(type) %>% mutate(max = max(index)) %>% ungroup() %>% select(type, max) %>% distinct()
#' max_abv <- max %>% filter(type == "AbV") %>% pull(max)
#' max_oth <- max %>% filter(type == "oth") %>% pull(max)
#'
#' tmp <- data.frame(names(mmv2_b)) %>%
#'   setNames("name") %>%
#'   mutate(add = ifelse(name %in% c("index"), max_oth,
#'                       ifelse(name %in% c("type"), "oth",
#'                              ifelse(name %in% `108P`, "108P", "")))) %>%
#'   t() %>%
#'   as.data.frame() %>%
#'   janitor::row_to_names(1) %>%
#'   mutate(index = as.numeric(index) + 1)
#'
#' mmv2_b <- rbind(mmv2_b, tmp)
#' rm(max, max_abv, max_oth, tmp)
#'
#' # 85V - other, non-antibody reactive, index max_oth, type oth
#' max <- mmv2_b %>% group_by(type) %>% mutate(max = max(index)) %>% ungroup() %>% select(type, max) %>% distinct()
#' max_abv <- max %>% filter(type == "AbV") %>% pull(max)
#' max_oth <- max %>% filter(type == "oth") %>% pull(max)
#'
#' tmp <- data.frame(names(mmv2_b)) %>%
#'   setNames("name") %>%
#'   mutate(add = ifelse(name %in% c("index"), max_oth,
#'                       ifelse(name %in% c("type"), "oth",
#'                              ifelse(name %in% `85V`, "85V", "")))) %>%
#'   t() %>%
#'   as.data.frame() %>%
#'   janitor::row_to_names(1) %>%
#'   mutate(index = as.numeric(index) + 1)
#'
#' mmv2_b <- rbind(mmv2_b, tmp)
#' rm(max, max_abv, max_oth, tmp)
#'
#' # 78VS - antibody reactive, index max_abv, type AbV
#' max <- mmv2_b %>% group_by(type) %>% mutate(max = max(index)) %>% ungroup() %>% select(type, max) %>% distinct()
#' max_abv <- max %>% filter(type == "AbV") %>% pull(max)
#' max_oth <- max %>% filter(type == "oth") %>% pull(max)
#'
#' tmp <- data.frame(names(mmv2_b)) %>%
#'   setNames("name") %>%
#'   mutate(add = ifelse(name %in% c("index"), max_abv,
#'                       ifelse(name %in% c("type"), "AbV",
#'                              ifelse(name %in% `78VS`, "78VS", "")))) %>%
#'   t() %>%
#'   as.data.frame() %>%
#'   janitor::row_to_names(1) %>%
#'   mutate(index = as.numeric(index) + 1)
#'
#' mmv2_b <- rbind(mmv2_b, tmp) %>% arrange(type, index)
#' rm(max, max_abv, max_oth, tmp)
#'
#' # write out ref tables
#' # write.csv(mmv2_a, "~/projects/hlaR/inst/extdata/ref/MHC_II_eplet_A_v2.csv", row.names = FALSE)
#' # write.csv(mmv2_b, "~/projects/hlaR/inst/extdata/ref/MHC_II_eplet_B_v2.csv", row.names = FALSE)
#'
#'
#'
