####################################################
### Create stim time files for AFNI 3dDeconvolve ###
####################################################

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(readxl)
library(dplyr)
library(tidyr)
library(LNCDR)
library(stringr)

# empty dataframe
antistate_timings <- data.frame(FRAME = integer(),
                      `Start Tm` = double(),
                      `End Tm` = double(),
                      Event = character(),
                      Location = integer(),
                      Run = character(),
                      stringsAsFactors = FALSE)

# read in relevant columns from each sheet of Katerina's task timing Excel sheet and append to the dataframe
for (i in 1:12){
  av <- readxl::read_excel('Data/Anti_Mix_Design_Lists_FINAL_minusfix_fixblockrecord_FINMOD.xls', 
                           sheet = sprintf('List%sAV', i), 
                           range = cell_cols('BB:BF'))
  av$Run <- sprintf('%sAV', i)
  va <- readxl::read_excel('Data/Anti_Mix_Design_Lists_FINAL_minusfix_fixblockrecord_FINMOD.xls', 
                           sheet = sprintf('List%sVA', i), 
                           range = cell_cols('BB:BF'))
  va$Run <- sprintf('%sVA', i) 
  av_va <- rbind(av, va)
  antistate_timings <- rbind(antistate_timings, av_va)
}

write.csv(antistate_timings, 'Data/antistate_timings_from_kat.csv')

# subject-date run orders
sub_runs <- readxl::read_excel('Data/Anti-State&MGSEncode_data.xls', 
                               sheet = 'Anti-State_sorted_by_year')
sub_runs$Scan_Date <- as.Date(sub_runs$Scan_Date)

# cleaning subject, date, run data
sub_runs <- sub_runs %>%
  select(id = Oxford_ID, 
         date = Scan_Date, 
         bircid = BIRC_ID, 
         AS1 = `Anti-State_1`, 
         AS2 = `Anti-State_2`, 
         AS3 = `Anti-State_3`, 
         AS4 = `Anti-State_4`) %>%
  pivot_longer(cols = AS1:AS4,
               names_to = 'run',
               values_to = 'run_type') 

sub_runs$run <- str_replace_all(sub_runs$run, "^AS", "")
sub_runs$run <- as.integer(sub_runs$run)
sub_runs$run_type <- str_replace_all(sub_runs$run_type, "^.\\(", "")
sub_runs$run_type <- str_replace_all(sub_runs$run_type, "^. \\(", "")
sub_runs$run_type <- str_replace_all(sub_runs$run_type, "^..\\(", "")
sub_runs$run_type <- str_replace_all(sub_runs$run_type, "\\)$", "")

# scored eye data (per trial, all subjects)
eye_data <- read.csv('Data/all_trail_et.tsv', 
                     sep = '\t', 
                     stringsAsFactors = FALSE)

eye_data <- eye_data %>%
  mutate(date = as.Date(as.character(date), format = '%Y%m%d')) %>%
  mutate(fstCorrect = as.numeric(fstCorrect)) %>%
  mutate(ErrCorr = as.numeric(ErrCorr)) %>%
  mutate(Incorrect = ifelse(Count == 0, 1, 0)) %>%
  mutate(Dropped = ifelse(Count == -1, 1, 0))

eye_data_sub_runs <- eye_data %>%
  left_join(sub_runs, by = c('id', 'date', 'run'))
write.csv(eye_data_sub_runs, 'Data/antistate_scored_eyd_plus_run_numbers.csv')

eye_data_sub_runs  %>%
  filter(is.na(run_type)) %>%
  select(id, date, run, run_type) %>%
  tally()

# subjects for analysis (have whoqol and 3+ scans)
subs_for_analysis <- read.csv('Data/cog_subs_to_preproc_20200217.csv', 
                              stringsAsFactors = FALSE)

subs_for_analysis <- subs_for_analysis %>%
  mutate(vtimestamp = as.Date(vtimestamp, format = '%m/%d/%Y'))
names(subs_for_analysis) <- c('id', 'date', 'age', 'note')

eye_data_sub_runs_analysis <- subs_for_analysis %>%
  left_join(eye_data_sub_runs, by = c('id', 'date'))

# filter trials that are fixation, merge for saving 1D files later
antistate_timings_excl_fix <- antistate_timings %>%
  filter(!grepl('FIX$', Event)) %>% # including 'FIXCUE'
  group_by(run_type) %>%
  filter(Event != lag(Event) | is.na(lag(Event))) %>%
  mutate(trial = cumsum(grepl('(ANTI|VGS)CUE', Event))) 

write.csv(antistate_timings_excl_fix, 'Data/antistate_timings_from_kat_excl_fix.csv')

merged_for_1d <- eye_data_sub_runs_analysis %>%
  full_join(antistate_timings_excl_fix, by = c('trial', 'run_type'))

write.csv(merged_for_1d, 'Data/cog_antistate_timings_merged_for_1d_20200219.csv')

# create new column in merged file to specify left/right
unique(merged_for_1d$Location)
merged_for_1d <- merged_for_1d %>%
  mutate(LeftRight = case_when(Location < 320 ~ 'Left',
                               Location > 320 ~ 'Right'))

# write 1D files for 3Ddeconvolve
oned_return <- function(f, ...){
  save1D(fname = f, ...)
  return(f)
}

merged_for_1d <- merged_for_1d %>%
  filter(!is.na(id))

for(subid in unique(merged_for_1d$id)){
  temp <- merged_for_1d %>%
    filter(id == subid)
  
  for(subdate in unique(temp$date)) {
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Left') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == 1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_left_correct.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Left') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == 2) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_left_errcorr.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Left') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == 0) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_left_incorr.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Left') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == -1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_left_dropped.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Right') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == 1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_right_correct.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Right') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == 2) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_right_errcorr.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Right') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == 0) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_right_incorr.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Right') %>%
      filter(Event == 'ANTICUE') %>%
      filter(Count == -1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_anti_right_dropped.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Left') %>%
      filter(Event == 'VGSCUE') %>%
      filter(Count != -1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_vgs_left.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Left') %>%
      filter(Event == 'VGSCUE') %>%
      filter(Count == -1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_vgs_left_dropped.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Right') %>%
      filter(Event == 'VGSCUE') %>%
      filter(Count != -1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_vgs_right.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
    temp %>%
      rename(block = run) %>%
      filter(LeftRight == 'Right') %>%
      filter(Event == 'VGSCUE') %>%
      filter(Count == -1) %>%
      filter(date == subdate) %>%
      summarize(d = strftime(first(date), format = '%Y%m%d'),
                f = sprintf('1D_files/%s_%s_vgs_right_dropped.1D', first(id), d),
                r = oned_return(f, ., colname = 'Start Tm', nblocks = 4))
    
  }
}


