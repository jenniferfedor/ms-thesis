###########################################
### Score antisaccade eye-tracking data ###
###########################################

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(dplyr)
library(lubridate)

# autoeyescored eye data
eye_data <- read.table('Data/all_trail_et.tsv', 
                       sep = '\t', 
                       header = TRUE, 
                       stringsAsFactors = FALSE)

eye_data <- eye_data %>%
  mutate(date = as.Date(as.character(date), format = '%Y%m%d')) %>%
  mutate(fstCorrect <- as.numeric(fstCorrect)) %>%
  mutate(ErrCorr <- as.numeric(ErrCorr)) %>%
  mutate(Incorrect = ifelse(Count == 0, 1, 0)) %>%
  mutate(Dropped = ifelse(Count == -1, 1, 0))

## accuracy
# proportion of each trial type in each session (considering all 4 anti runs together)
anti_perf <- eye_data %>%
  filter(AS == 'AS') %>%
  group_by(id, date) %>%
  summarise(anti_perc_correct = mean(fstCorrect),
            anti_perc_error_corrected = mean(ErrCorr),
            anti_perc_incorrect = mean(Incorrect),
            anti_perc_dropped_trials = mean(Dropped)) 

# proportion of each trial type in each session among non-dropped trials
anti_perf_nodrop <- eye_data %>%
  filter(AS == 'AS') %>%
  filter(Dropped == 0) %>%
  group_by(id, date) %>%
  summarise(anti_perc_correct_nodrop = mean(fstCorrect),
            anti_perc_error_corrected_nodrop = mean(ErrCorr)) 

## saccade latency
# average antisaccade latency on correct trials only
anti_lat_corr <- eye_data %>%
  filter(AS == 'AS') %>%
  filter(Count == 1) %>%
  group_by(id, date) %>%
  summarise(anti_avg_lat_correct_trials = mean(lat)) 

# average antisaccade latency on error corrected trials
anti_lat_errcorr <- eye_data %>%
  filter(AS == 'AS') %>%
  filter(Count == 2) %>%
  group_by(id, date) %>%
  summarise(anti_avg_lat_error_corrected_trials = mean(lat)) 


# join accuracy and latency data
anti_perf_lat_20200228 <- anti_perf %>%
  full_join(anti_perf_nodrop, by = c('id', 'date')) %>%
  full_join(anti_lat_corr, by = c('id', 'date')) %>%
  full_join(anti_lat_errcorr, by = c('id', 'date')) %>%
  rename('lunaid' = 'id')

write.csv(anti_perf_lat_20200228, 'Data/anti_perf_lat_20200228.csv')
