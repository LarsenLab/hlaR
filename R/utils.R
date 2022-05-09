#' Basic functions
#' GenerateLookup() called in CalEpletMHCII()
#' @param lkup_in data table
#' @param locus_in string
#' CalRiskScr() called in CalEpletMHCII()
#' @param dat_scr dataframe
#' @name utils
#' @import
#' janitor
NULL
#> NULL

#' @rdname utils
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

#' @rdname utils
CalRiskScr <- function(dat_scr) {
  risk_scr <- dat_scr %>%
              filter(mm_eplets != "Not Found") %>%
              mutate(hla = gsub("\\*.*", "", hla),
                     grp = ifelse(hla %in% c("DRB1", "DRB3", "DRB4", "DRB5"), "DR",
                                  ifelse(hla %in% c("DQA1", "DQB1") & haplotype_id %in% c("1"), "DQ1",
                                         ifelse(hla %in% c("DQA1", "DQB1") & haplotype_id %in% c("2"), "DQ2", hla)))) %>%
              filter(grp %in% c("DR", "DQ1", "DQ2") & mm_cnt > 0 ) %>%
              group_by(grp) %>%
              summarise(max1 = suppressWarnings(max(mm_cnt, na.rm=TRUE)),
                        sum1 = sum(mm_cnt, na.rm=TRUE)) %>%
              ungroup()  %>%
              mutate(max1 = ifelse(grp %in% c("DQ1", "DQ2"), sum1, max1)) %>%
              select(-sum1) %>%
              mutate(locus = ifelse(grp %in% c("DR"), "DR", "DQ")) %>%
              group_by(locus) %>%
              summarise(score = suppressWarnings(max(max1, na.rm=TRUE))) %>%
              ungroup() %>%
              as.data.frame() %>%
              rownames_to_column %>%
              gather(var, value, -rowname) %>%
              spread(rowname, value) %>%
              janitor::row_to_names(1) %>%
              as.data.frame()

  # set DQ if no DQ
  if(!("DQ" %in% names(risk_scr))) {
    risk_scr$DQ = 0
  }

  # set DR if no DR
  if(!("DR" %in% names(risk_scr))) {
    risk_scr$DR = 0
  }

  risk_scr <- risk_scr %>%
              # change data type to numeric for risk calculation
              mutate(DQ = as.numeric(DQ),
                     DR = as.numeric(DR)) %>%
              mutate(risk = ifelse(between(DQ, 15, 31), "high",
                                   ifelse((DR >= 7 & DQ <= 14) | (DR < 7 & between(DQ, 9, 15)), "interm",
                                          ifelse(DR < 7 & DQ < 9, "low", "unknown"))))

  return(risk_scr)
}
