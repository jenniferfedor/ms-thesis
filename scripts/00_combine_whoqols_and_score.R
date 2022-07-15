###################################
### Score WHO-QOL questionnaire ###
###################################

## required packages
library(dplyr)
library(lubridate)

## import data collected at behavioral visits
whoqol_beh <- read.csv('CogLong_WHOQOL.csv', 
                       header = TRUE, 
                       stringsAsFactors = FALSE)
whoqol_beh <- whoqol_beh[1:90, 1:29]
names(whoqol_beh)[1] <- 'lunaid'
whoqol_beh$date <- as.POSIXct(mdy(whoqol_beh$date))

## import data collected via R03/online survey
# WHOQOL plus demographics version
whoqol_ro3_a <- read.csv('R03+Research+Follow-up+Battery_November+14,+2019_11.43.csv', 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
whoqol_ro3_a <- whoqol_ro3_a[3:nrow(whoqol_ro3_a),]
whoqol_ro3_a <- filter(whoqol_ro3_a, Status != 'Survey Preview')

id_a <- match('ExternalReference', names(whoqol_ro3_a))
date_a <- match('EndDate', names(whoqol_ro3_a))
start_a <- match('Q60_1', names(whoqol_ro3_a))
stop_a <- match ('Q68', names(whoqol_ro3_a))

whoqol_ro3_a <- whoqol_ro3_a[, c(id_a, date_a, start_a:stop_a)]
names(whoqol_ro3_a) <- names(whoqol_beh)[1:29]
whoqol_ro3_a$date <- as.POSIXct(whoqol_ro3_a$date)

# WHOQOL only version
whoqol_ro3_b <- read.csv('R03+WHOQOL-BREF_November+15,+2019_08.48.csv', 
                         header = TRUE, 
                         stringsAsFactors = FALSE)
whoqol_ro3_b <- whoqol_ro3_b[3:nrow(whoqol_ro3_b),]
whoqol_ro3_b <- filter(whoqol_ro3_b, Status != 'Survey Preview')

id_b <- match('ExternalReference', names(whoqol_ro3_b))
date_b <- match('EndDate', names(whoqol_ro3_b))
start_b <- match('Q3_1', names(whoqol_ro3_b))
stop_b <- match ('Q15', names(whoqol_ro3_b))

whoqol_ro3_b <- whoqol_ro3_b[, c(id_b, date_b, start_b:stop_b)]
names(whoqol_ro3_b) <- names(whoqol_beh)[1:29]
whoqol_ro3_b$date <- as.POSIXct(whoqol_ro3_b$date, format = "%m/%d/%Y %H:%M", tz = "")

# merge two online survey versions
whoqol_ro3 <- rbind(whoqol_ro3_a, whoqol_ro3_b, stringsAsFactors = FALSE)

# recode online survey items from text input to numeric values
whoqol_ro3 <- whoqol_ro3 %>% 
  mutate(X1 = recode(X1, 
                     'Very Poor' = 1, 
                     'Poor' = 2, 
                     'Neither Poor nor Good' = 3, 
                     'Good' = 4, 
                     'Very Good' = 5)) %>%
  mutate(X2 = recode(X2, 
                     'Very Dissatisfied' = 1, 
                     'Dissatisfied' = 2, 
                     'Neither Dissatisfied nor Satisfied' = 3, 
                     'Satisfied' = 4, 
                     'Very Satisfied' = 5)) %>%
  mutate(X3 = recode(X3, 
                     'Not at all' = 5, 
                     'A little' = 4, 
                     'A moderate amount' = 3, 
                     'Very much' = 2, 
                     'An extreme amount' = 1)) %>%
  mutate(X4 = recode(X4, 
                     'Not at all' = 5, 
                     'A little' = 4, 
                     'A moderate amount' = 3, 
                     'Very much' = 2, 
                     'An extreme amount' = 1)) %>%
  mutate(X5 = recode(X5, 
                     'Not at all' = 1, 
                     'A little' = 2, 
                     'A moderate amount' = 3, 
                     'Very much' = 4, 
                     'An extreme amount' = 5)) %>%
  mutate(X6 = recode(X6, 
                     'Not at all' = 1, 
                     'A little' = 2, 
                     'A moderate amount' = 3, 
                     'Very much' = 4, 
                     'An extreme amount' = 5)) %>%
  mutate(X7 = recode(X7, 
                     'Not at all' = 1, 
                     'A little' = 2, 
                     'A moderate amount' = 3, 
                     'Very much' = 4, 
                     'Extremely' = 5)) %>%
  mutate(X8 = recode(X8, 
                     'Not at all' = 1, 
                     'A little' = 2, 
                     'A moderate amount' = 3, 
                     'Very much' = 4, 
                     'Extremely' = 5)) %>%
  mutate(X9 = recode(X9, 
                     'Not at all' = 1, 
                     'A little' = 2, 
                     'A moderate amount' = 3, 
                     'Very much' = 4, 
                     'Extremely' = 5)) %>%
  mutate(X10 = recode(X10, 
                      'Not at all' = 1, 
                      'A little' = 2, 
                      'Moderately' = 3, 
                      'Mostly' = 4, 
                      'Completely' = 5)) %>%
  mutate(X11 = recode(X11, 
                      'Not at all' = 1, 
                      'A little' = 2, 
                      'Moderately' = 3, 
                      'Mostly' = 4, 
                      'Completely' = 5)) %>%
  mutate(X12 = recode(X12, 
                      'Not at all' = 1, 
                      'A little' = 2, 
                      'Moderately' = 3, 
                      'Mostly' = 4, 
                      'Completely' = 5)) %>%
  mutate(X13 = recode(X13, 
                      'Not at all' = 1, 
                      'A little' = 2, 
                      'Moderately' = 3, 
                      'Mostly' = 4, 
                      'Completely' = 5)) %>%
  mutate(X14 = recode(X14, 
                      'Not at all' = 1, 
                      'A little' = 2, 
                      'Moderately' = 3, 
                      'Mostly' = 4, 
                      'Completely' = 5)) %>%
  mutate(X15 = recode(X15, 
                      'Very Poor' = 1, 
                      'Poor' = 2, 
                      'Neither poor nor good' = 3, 
                      'Good' = 4, 
                      'Very good' = 5)) %>%
  mutate(X16 = recode(X16, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X17 = recode(X17, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X18 = recode(X18, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X19 = recode(X19, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X20 = recode(X20, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X21 = recode(X21, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X22 = recode(X22, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X23 = recode(X23, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X24 = recode(X24, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X25 = recode(X25, 
                      'Very Dissatisfied' = 1, 
                      'Dissatisfied' = 2, 
                      'Neither dissatisfied nor satisfied' = 3, 
                      'Satisfied' = 4, 
                      'Very Satisfied' = 5)) %>%
  mutate(X26 = recode(X26, 
                      'Never' = 5, 
                      'Seldom' = 4, 
                      'Quite Often' = 3, 
                      'Very Often' = 2, 
                      'Always' = 1))

## replace missing value with mean of other domain items per scoring guidelines
whoqol_ro3[19, c('X15')] <- mean(unlist(whoqol_ro3[19, c('X3', 'X4', 'X10', 'X16', 'X17', 'X18')]))

## merge online survey data with behavioral visit data
whoqol <- rbind(whoqol_ro3, whoqol_beh, stringsAsFactors = FALSE)

## score whoqol per scoring guidelines
whoqol <- whoqol %>% 
  mutate(D_1_Raw = X3 + X4 + X10 + X15 + X16 + X17 + X18) %>% # physical
  mutate(D_2_Raw = X5 + X6 + X7 + X11 + X19 + X26) %>% # psychological
  mutate(D_3_Raw = X20 + X21 + X22) %>% # social
  mutate(D_4_Raw = X8 + X9 + X12 + X13 + X14 + X23 + X24 + X25) # environment
   
## functions to transform scores to 4-20 scale or 0-100 scale per guidelines
scale_20_to_100 <- function(x) {
  ifelse(x == 4, 0,
   ifelse(x == 5, 6,
    ifelse(x == 6, 13,
     ifelse(x == 7, 19,
      ifelse(x == 8, 25,
       ifelse(x ==9, 31,
        ifelse(x == 10, 38,
         ifelse(x == 11, 44,
          ifelse(x == 12, 50,
           ifelse(x == 13, 56,
            ifelse(x == 14, 63,
             ifelse(x == 15, 69,
              ifelse(x == 16, 75,
               ifelse(x == 17, 81,
                ifelse(x == 18, 88,
                 ifelse(x == 19, 94, 100))))))))))))))))
}

scale_domain_1 <- function(x) {
  ifelse(x == 7, 4, 
   ifelse(x == 8 | x == 9, 5,
    ifelse(x == 10 | x == 11, 6,
     ifelse(x == 12 | x == 13, 7,
      ifelse(x == 14, 8,
       ifelse(x == 15 | x == 16, 9,
        ifelse(x == 17 | x == 18, 10,
         ifelse(x == 19 | x == 20, 11,
          ifelse(x == 21, 12,
           ifelse(x == 22 | x == 23, 13,
            ifelse(x == 24 | x == 25, 14,
             ifelse(x == 26 | x == 27, 15,
              ifelse(x == 28, 16,
               ifelse(x == 29 | x == 30, 17,
                ifelse(x == 31 | x == 32, 18,
                 ifelse(x == 33 | x == 34, 19, 20))))))))))))))))
}

scale_domain_2 <- function(x) {
  ifelse(x == 6, 4, 
   ifelse(x == 7 | x == 8, 5,
    ifelse(x == 9, 6,
     ifelse(x == 10 | x == 11, 7,
      ifelse(x == 12, 8,
       ifelse(x == 13 | x == 14, 9,
        ifelse(x == 15, 10,
         ifelse(x == 16 | x == 17, 11,
          ifelse(x == 18, 12,
           ifelse(x == 19 | x == 20, 13,
            ifelse(x == 21, 14,
             ifelse(x == 22 | x == 23, 15,
              ifelse(x == 24, 16,
               ifelse(x == 25 | x == 26, 17,
                ifelse(x == 27, 18,
                 ifelse(x == 28 | x == 29, 19, 20))))))))))))))))
}

scale_domain_3 <- function(x) {
  ifelse(x == 3, 4, 
   ifelse(x == 4, 5,
    ifelse(x == 5, 7,
     ifelse(x == 6, 8,
      ifelse(x == 7, 9,
       ifelse(x == 8, 11,
        ifelse(x == 9, 12,
         ifelse(x == 10, 13,
          ifelse(x == 11, 15,
           ifelse(x == 12, 16,
            ifelse(x == 13, 17,
             ifelse(x == 14, 19, 20))))))))))))
}

scale_domain_4 <- function(x) {
  ifelse(x == 8, 4, 
   ifelse(x == 9 | x == 10, 5,
    ifelse(x == 11 | x == 12, 6,
     ifelse(x == 13 | x == 14, 7,
      ifelse(x == 15 | x == 16, 8,
       ifelse(x == 17 | x == 18, 9,
        ifelse(x == 19 | x == 20, 10,
         ifelse(x == 21 | x == 22, 11,
          ifelse(x == 23 | x == 24, 12,
           ifelse(x == 25 | x == 26, 13,
            ifelse(x == 27 | x == 28, 14,
             ifelse(x == 29 | x == 30, 15,
              ifelse(x == 31 | x == 32, 16,
               ifelse(x == 33 | x == 34, 17,
                ifelse(x == 35 | x == 36, 18,
                 ifelse(x == 37 | x == 38, 19, 20))))))))))))))))
}

## transform the raw scores to scaled scores
whoqol <- whoqol %>% 
  mutate(D_1_Trans_4_20 = scale_domain_1(D_1_Raw)) %>%
  mutate(D_2_Trans_4_20 = scale_domain_2(D_2_Raw)) %>%
  mutate(D_3_Trans_4_20 = scale_domain_3(D_3_Raw)) %>%
  mutate(D_4_Trans_4_20 = scale_domain_4(D_4_Raw)) %>%
  mutate(D_1_Trans_0_100 = scale_20_to_100(D_1_Trans_4_20)) %>%
  mutate(D_2_Trans_0_100 = scale_20_to_100(D_2_Trans_4_20)) %>%
  mutate(D_3_Trans_0_100 = scale_20_to_100(D_3_Trans_4_20)) %>%
  mutate(D_4_Trans_0_100 = scale_20_to_100(D_4_Trans_4_20))

## create counter variable to keep track of repeated administrations
whoqol <- whoqol %>% 
  group_by(lunaid) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(counter = row_number())

## select scores from subject's most recent administration of the whoqol
whoqol_most_recent <- whoqol %>%
  group_by(lunaid) %>%
  filter(counter == max(counter))

## histograms to check if there is sufficient variability in scores 
domains <- c('Domain 1', 'Domain 2', 'Domain 3', 'Domain 4')
types <- c('Raw', 'Scaled_4_20', 'Scaled_0_100')

whoqol_variability <- data.frame(matrix(NA, nrow = 1, ncol = 6))
whoqol_variability <- whoqol_variability[-1, , drop = FALSE]
names(whoqol_variability) <- c('Measure', 'Mean', 'SD', 'Median', 'Min', 'Max')

columns <- names(whoqol_most_recent)[30:41]

for (i in seq(1, length(columns), by = 1)){
  whoqol_variability[i, 'Measure'] <- columns[i]
  whoqol_variability[i, 'Mean'] <- mean(as.numeric(unlist(whoqol_most_recent[, columns[i]])))
  whoqol_variability[i, 'SD'] <- sd(as.numeric(unlist(whoqol_most_recent[, columns[i]])))
  whoqol_variability[i, 'Median'] <- median(as.numeric(unlist(whoqol_most_recent[, columns[i]])))
  whoqol_variability[i, 'Min'] <- range(as.numeric(unlist(whoqol_most_recent[, columns[i]])))[1]
  whoqol_variability[i, 'Max'] <- range(as.numeric(unlist(whoqol_most_recent[, columns[i]])))[2]
}

par(mfrow = c(3, 4))
for (i in seq(1, length(columns), by = 1)){
  hist(as.numeric(unlist(whoqol_most_recent[, columns[i]])), 
       main = columns[i], 
       xlab = columns[i])
}

# proportion of subjects selecting each option for overall QOL
table(as.factor(whoqol_most_recent$X1))
prop.table(table(as.factor(whoqol_most_recent$X1)))

## number of subjects with i administrations
for (i in seq(1, 5, by = 1)){
  print(nrow(whoqol %>%
  group_by(lunaid) %>%
  filter(max(counter) == i) %>%
  filter(counter == max(counter)))) 
}

## histograms for raw domain scores only
raw_scores <- names(whoqol_most_recent)[30:33]
raw_domain_names <- c('Physical health', 'Psychological health', 'Social relationships', 'Environment')

par(mfrow = c(2, 2))
for (i in seq(1, length(raw_scores), by = 1)){
  hist(as.numeric(unlist(whoqol_most_recent[, raw_scores[i]])), 
       main = raw_domain_names[i], 
       xlab = raw_scores[i])
}

write.csv(whoqol_most_recent, 'whoqol_ages_scores_20200204.csv')
