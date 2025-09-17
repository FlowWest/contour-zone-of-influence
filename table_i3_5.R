# Script to create Table I.3-5
alt_order = c("EXP1", "EXP3", "NAA","ALT1",
              "ALT2d", "ALT2b", "ALT2c", "ALT2a", "ALT3", "ALT4")
inflow_order = c("lolo", "lomed", "lohi", "medlo", "medmed", "medhi", "hilo", "himed", "hihi", "NA")

# Calsim 3 input data-----------------------------
# These data are used for other analyses and were used to determine bins. However, these do not
# reflect actual sample sizes for zoi analysis (similar but slight differences)

# read in data (Calsim inputs from Steve's binning document) -----------------------
bins <- read_excel("data_raw/calsim/Reclamation_2021LTO_SacR_SJR_OMR_Binning_rev01_20230929_result.xlsx", skip = 5)
# rename column names
colnames(bins) <- c("Date", "Flow_EXP1", "OMR_EXP1", "Flow_EXP3", "OMR_EXP3",
                    "Flow_NAA", "OMR_NAA", "Flow_ALT1", "OMR_ALT1",
                    "Flow_ALT2a", "OMR_ALT2a", "Flow_ALT2b", "OMR_ALT2b",
                    "Flow_ALT2c", "OMR_ALT2c", "Flow_ALT2d", "OMR_ALT2d",
                    "Flow_ALT3", "OMR_ALT3", "Flow_ALT4", "OMR_ALT4")
bins2 <- bins[-1,]

omr_bins <- read_excel("data_raw/calsim/Reclamation_2021LTO_OMR_rev01_20231010.xlsx", skip = 11)
bins2_flow <- dplyr::select(bins2, contains("Flow"))
bins_upd <- cbind(omr_bins, bins2_flow)

colnames(bins_upd) <- c("Date", "OMR_EXP1",  "OMR_EXP3","OMR_NAA",
                        "OMR_ALT1","OMR_ALT2a",
                        "OMR_ALT2b","OMR_ALT2c", "OMR_ALT2d","OMR_ALT3","OMR_ALT4",
                        "Flow_EXP1",  "Flow_EXP3", "Flow_NAA", "Flow_ALT1",
                        "Flow_ALT2a", "Flow_ALT2b", "Flow_ALT2c", "Flow_ALT2d",
                        "Flow_ALT3", "Flow_ALT4" )

###  make long ---------

# use bins_upd for looking at frequency in the BA
bins_long <-  bins_upd %>%
  pivot_longer(
    -Date,
    cols_vary = "slowest",
    names_to = c(".value", "Alt"),
    names_sep = "_"
  ) %>%
  mutate(Month = month(Date))

### filter Dec-June --------
bins_months <- bins_long %>% filter(Month %in% c(12, 1, 2, 3, 4, 5, 6))

bins_freq_new <- bins_months |>
  mutate(OMR_group = case_when(OMR >= -2500 & OMR <=-1500 ~ "-2000",
                               OMR >= -4000 & OMR <=-3000 ~ "-3500",
                               OMR >= -5500 & OMR <=-4500 ~ "-5000",
                               OMR <=-5500 ~ "less than -5500",
                               OMR > -1500 ~ "positive",
                               TRUE ~ "Other"))
### calculate sample sizes -----------------
bins_summary_months_new <- bins_freq_new %>%
  group_by(Flow, OMR_group, Alt) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(Alt = factor(Alt, levels = alt_order ),
         Flow = factor(Flow, levels = inflow_order),
         OMR = factor(OMR_group, levels = c("positive", "-2000", "-3500", "-5000", "less than -5500", "Other"))) %>%
  arrange(Flow, OMR, Alt)

# fill in groups that have no samples; replace with zero
bins_summary_complete <- bins_summary_months_new %>%
  complete(Flow, OMR,
           nesting(Alt),
           fill = list(n = 0))

## make tables ------------------------------------

# calculate proportion NA
not_included <- bins_summary_complete %>%
  filter(OMR %in% c("Other", "positive")) %>%
  #filter((Flow == "NA" | OMR == "Other") | (Flow == "NA" & OMR == "Other")) %>%
  group_by(Alt, OMR) %>%
  summarize(num =sum(n),
            prop = round(num/700,2))

not_included_wide <- not_included %>%
  dplyr::select(-num) %>%
  pivot_wider(names_from = "Alt", values_from = "prop")
# write_csv(not_included_wide, "data_export/data_not_included_allalts.csv")
