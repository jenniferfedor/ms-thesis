#################################
### Compile data for analysis ###
#################################

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(dplyr)
library(tidyr)
library(lubridate)

# read in scan data (316 deconvolved scans among 50 subjects)
scan_data <- read.csv('Data/roi_percent_signal_change_20200228.csv', 
                      stringsAsFactors = FALSE) %>%
  select(-c('X')) %>%
  rename('lunaid' = 'id') %>%
  mutate(date = as.Date(date, format = '%Y-%m-%d'))

# read in scored anti eye data 
eye_data <- read.csv('Data/whoqol_plus_eye_20200228.csv', 
                     stringsAsFactors = FALSE) %>%
  filter(whoqol_age == 'Scan') %>%
  select(-c('X', 'note', 'age_type')) %>%
  mutate(date = as.Date(date, format = '%Y-%m-%d')) %>%
  mutate(dob = as.Date(dob, format = '%Y-%m-%d'))

# read in maternal education level data
mat_edu <- read.csv('Data/maternal_education_20200219.csv', 
                    stringsAsFactors = FALSE) %>%
  select(-c('X'))

# change reference category for maternal education level
mat_edu <- within(mat_edu, 
                  level_edu_mother_cat <- relevel(as.factor(level_edu_mother_cat), 
                                                  ref = 'Completed high school'))

# merge eye and scan data
data_for_gamms <- scan_data %>%
  left_join(eye_data, by = c('lunaid', 'date')) %>%
  mutate(sex = as.factor(sex)) %>%
  mutate(duplicate_row = duplicated(.)) %>%
  filter(duplicate_row == FALSE) %>%
  select(-c('duplicate_row')) %>%
  mutate(age_c = age - mean(age))

data_for_gamms <- data_for_gamms %>%
  left_join(mat_edu, by = c('lunaid')) %>%
  mutate(lunaid = as.factor(lunaid)) 

data_for_gamms %>% write.csv('Data/final_eye_scan_data_for_analysis_20200228.csv') 