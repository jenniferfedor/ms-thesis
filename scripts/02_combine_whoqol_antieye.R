##################################
### Merge WHO-QOL and eye data ###
##################################

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(dplyr)

# whoqol data
whoqol <- read.csv('Data/whoqol_and_scan_ages_20200228.csv', 
                   header = TRUE, 
                   stringsAsFactors = FALSE)

whoqol <- whoqol %>%
  select(-c('X')) %>%
  mutate(dob = as.Date(as.character(dob), format = '%Y-%m-%d')) %>%
  mutate(vtimestamp = as.Date(as.character(vtimestamp), format = '%Y-%m-%d')) %>%
  rename('date' = 'vtimestamp')

# import scored anti task eye data
eye <- read.csv('Data/anti_perf_lat_20200228.csv', 
                header = TRUE, 
                stringsAsFactors = FALSE)

eye <- eye %>%
  select(-c('X')) %>%
  mutate(date = as.Date(as.character(date), format = '%Y-%m-%d')) %>%
  rename('lunaid' = 'id')

# join whoqol and eye data
whoqol_plus_eye <- whoqol %>%
  left_join(eye, by = c('lunaid', 'date')) %>%
  write.csv('Data/whoqol_plus_eye_20200228.csv')


