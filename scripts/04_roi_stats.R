##################################
### Average activation in ROIs ###
##################################

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(tidyr)
library(dplyr)
library(stringr)
library(lubridate)

# roi names and coordinates from Ordaz et al. 2013 paper
roi_names <- read.csv('Data/sarah_rois_20200225.csv', 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

# average activation per roi per session - output from AFNI 3dROIstats
roi_vals <- read.table('Data/jen_roi_stats_20200228.txt', 
                       header = TRUE, 
                       stringsAsFactors = FALSE)

# rename AFNI sub-brick labels to be more descriptive
roi_vals <- roi_vals %>%
  mutate(condition = case_when(Sub.brick == '0[anti_corr]' ~ 'anti_corr_v_baseline',
                               Sub.brick == '1[anti_errc]' ~ 'anti_errcorr_v_baseline',
                               Sub.brick == '2[anti_inco]' ~ 'anti_incorr_v_baseline',
                               Sub.brick == '3[anti_corr]' ~ 'anti_corr_v_errcorr',
                               Sub.brick == '4[anti_corr]' ~ 'anti_corr_v_incorr',
                               Sub.brick == '5[anti_corr]' ~ 'anti_corr_v_vgs')) %>%
  mutate(File2 = str_replace_all(File, '^../subjs/...../', '')) %>%
  mutate(File2 = str_replace_all(File2, '_bucket.nii.gz.*', '')) %>%
  separate(col = File2,
           into = c('lunaid', 'date'),
           sep = '_') %>% 
  select(lunaid, date, condition, everything()) %>%
  select(-c('File', 'Sub.brick')) %>%
  mutate(date = as.Date(date, format = '%Y%m%d'))

# rename columns
colnames(roi_vals) <- c('id', 'date', 'condition',
                        'SEF', 'pre_SMA', 'L_FEF', 'R_FEF',
                        'L_putamen', 'R_putamen', 'L_pPC', 'R_pPC',
                        'L_dlPFC', 'R_dlPFC', 'L_vlPFC', 'R_vlPFC', 
                        'dACC')

# pivot dataframe from long to wide format
roi_vals_wide <- roi_vals %>%
  pivot_wider(id_cols = c('id', 'date'),
              names_from = c('condition'),
              values_from = c('SEF', 'pre_SMA', 'L_FEF', 'R_FEF', 
                              'L_putamen', 'R_putamen', 'L_pPC', 'R_pPC',
                              'L_dlPFC', 'R_dlPFC', 'L_vlPFC', 'R_vlPFC', 
                              'dACC'))

roi_vals_wide %>% write.csv('Data/roi_percent_signal_change_20200228.csv')

