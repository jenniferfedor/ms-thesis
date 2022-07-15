###########
## GAMMs ##
###########

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(dplyr)
library(ggplot2)
library(gratia)
library(gridExtra)
library(lubridate)
library(mgcv)
library(purrr)
library(stringr)
library(tidyr)
library(xtable)


# function to plot raw data, fitted values and pointwise 95% confidence intervals from a GAMM
plot_gamm <- function(model, y_axis_label = '', sig = TRUE) {
  outcome <- str_split(model$gam$formula, '~')[[2]]
  
  newdat <- data.frame(dat$age_c) # observed centered ages
  names(newdat) <- 'age_c'
  
  pred <- predict.gam(model$gam, 
                      newdata = newdat,
                      type = 'response', 
                      se.fit = TRUE)
  
  plotdat <- data.frame(cbind(dat$age, pred$fit, pred$se))
  plotdat <- plotdat %>% 
    mutate(outcome = 1) 
  names(plotdat) <- c('age', 'fit', 'se', outcome)
  
  if (sig == TRUE){
    plotcolor <- '#476FD1'
  } else {
    plotcolor <- 'grey50'
  }

  gammplot <- dat %>% 
    ggplot(aes(x = age, 
               y = eval(parse(text = outcome)))) + 
    geom_line(alpha = .2, 
              color = plotcolor,
              aes(group = lunaid)) + 
    geom_point(alpha = .2, 
               shape = 16,
               color = plotcolor) +
    geom_ribbon(data = plotdat, 
                alpha = .7, 
                aes(ymin = fit-(1.96*se), ymax = fit+(1.96*se)), 
                show.legend = FALSE, 
                fill = plotcolor) +
    geom_line(data = plotdat, 
              aes(y = fit), 
              show.legend = FALSE, 
              color = plotcolor) +
    labs(x = 'Age (years)', 
         y = y_axis_label) +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = 'black')) +
    theme(axis.title = element_text(size = 11),
          axis.text = element_text(size = 10))
  return(gammplot)
}

# function for post-hoc analysis to identify periods of developmental change
calc_dev_change <- function(model) {
  deriv <- gratia::derivatives(model)
  deriv <- deriv %>%
    mutate(sig = !(0 > lower & 0 < upper)) 
  ages <- list(c(min(deriv$data[deriv$sig == TRUE]) + mean(dat$age), 
                 max(deriv$data[deriv$sig == TRUE]) + mean(dat$age)))
  return(ages[[1]])
}


# read in final data file and format variables
dat <- read.csv('Data/final_eye_scan_data_for_analysis_20200228.csv', 
                stringsAsFactors = FALSE) 
dat <- dat %>%
  select(-c('X')) %>%
  mutate(lunaid = as.factor(lunaid)) %>%
  mutate(date = as.Date(date, format = '%Y-%m-%d')) %>%
  mutate(dob = as.Date(dob, format = '%Y-%m-%d')) %>%
  mutate(sex = as.factor(sex))
dat <- within(dat, level_edu_mother_cat <- relevel(as.factor(level_edu_mother_cat), 
                                                   ref = 'Completed high school'))

### I. Behavioral outcomes 

## 1. Antisaccade proportion correct trials
# fit the model with smooth term for age
model_formula <- as.formula("anti_perc_correct_nodrop ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_anti_perc_correct_nodrop <- gamm(as.formula(model_formula), 
                                      family = 'quasibinomial',
                                      random = list(lunaid=~age_c),
                                      data = dat, 
                                      method = 'REML')
summary(gamm_anti_perc_correct_nodrop$gam)
summary(gamm_anti_perc_correct_nodrop$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_anti_perc_correct_nodrop$gam)

# identify significant periods of developmental change
anti_perc_correct_devchange <- calc_dev_change(gamm_anti_perc_correct_nodrop)
anti_perc_correct_devchange

# plot
anti_perc_correct_plot <- plot_gamm(
  model = gamm_anti_perc_correct_nodrop, 
  y_axis_lab = 'Proportion correct trials',
  sig = TRUE)
anti_perc_correct_plot

# add bar to plot reflecting period of sig dev change
anti_perc_correct_plot <- anti_perc_correct_plot + 
  annotate(geom = "rect", 
           xmin = anti_perc_correct_devchange[1], 
           xmax = anti_perc_correct_devchange[2], 
           ymin = -0.04, 
           ymax = -0.01, 
           fill = '#476FD1')
anti_perc_correct_plot


# fit the model w/ fixed effects of sex 
model_formula <- as.formula("anti_perc_correct_nodrop ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_anti_perc_correct_nodrop_sex <- gamm(as.formula(model_formula), 
                                          family = 'quasibinomial',
                                          random = list(lunaid=~age_c),
                                          data = dat, 
                                          method = 'REML')
summary(gamm_anti_perc_correct_nodrop_sex$gam)
summary(gamm_anti_perc_correct_nodrop_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("anti_perc_correct_nodrop ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_anti_perc_correct_nodrop_medu <- gamm(as.formula(model_formula), 
                                           family = 'quasibinomial',
                                           random = list(lunaid=~age_c),
                                           data = dat, 
                                           method = 'REML')
summary(gamm_anti_perc_correct_nodrop_medu$gam)
summary(gamm_anti_perc_correct_nodrop_medu$lme)


# save random slopes and intercepts from final model
ranef_anti_perc_correct_nodrop <- ranef(gamm_anti_perc_correct_nodrop$lme)$lunaid %>%
  rename(anti_perc_correct_nodrop_rint = `(Intercept)`,
         anti_perc_correct_nodrop_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 2. Antisaccade proportion error-corrected trials
# fit the model with smooth term for age
model_formula <- as.formula("anti_perc_error_corrected_nodrop ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_anti_perc_error_corrected_nodrop <- gamm(as.formula(model_formula), 
                                              family = 'quasibinomial',
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_anti_perc_error_corrected_nodrop$gam)
summary(gamm_anti_perc_error_corrected_nodrop$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_anti_perc_error_corrected_nodrop$gam)

# identify significant periods of developmental change
anti_perc_error_corrected_nodrop_devchange <- calc_dev_change(gamm_anti_perc_error_corrected_nodrop)
anti_perc_error_corrected_nodrop_devchange

# plot
anti_perc_error_corrected_nodrop_plot <- plot_gamm(
  model = gamm_anti_perc_error_corrected_nodrop,
  y_axis_label = 'Proportion error-corrected trials',
  sig = TRUE)
anti_perc_error_corrected_nodrop_plot

anti_perc_error_corrected_nodrop_plot <- anti_perc_error_corrected_nodrop_plot +
  annotate(geom = 'rect', 
           xmin = anti_perc_error_corrected_nodrop_devchange[1], 
           xmax = anti_perc_error_corrected_nodrop_devchange[2], 
           ymin = -0.04, 
           ymax = -0.02, 
           fill = '#476FD1')
anti_perc_error_corrected_nodrop_plot

# fit the model w/ fixed effects of sex 
model_formula <- as.formula("anti_perc_error_corrected_nodrop ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_anti_perc_error_corrected_nodrop_sex <- gamm(as.formula(model_formula), 
                                                  family = 'quasibinomial',
                                                  random = list(lunaid=~age_c),
                                                  data = dat, 
                                                  method = 'REML')
summary(gamm_anti_perc_error_corrected_nodrop_sex$gam)
summary(gamm_anti_perc_error_corrected_nodrop_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("anti_perc_error_corrected_nodrop ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_anti_perc_error_corrected_nodrop_medu <- gamm(as.formula(model_formula), 
                                                   family = 'quasibinomial',
                                                   random = list(lunaid=~age_c),
                                                   data = dat, 
                                                   method = 'REML')
summary(gamm_anti_perc_error_corrected_nodrop_medu$gam)
summary(gamm_anti_perc_error_corrected_nodrop_medu$lme)

# save random slopes and intercepts from final model
ranef_anti_perc_error_corrected_nodrop <- ranef(gamm_anti_perc_error_corrected_nodrop$lme)$lunaid %>%
  rename(anti_perc_error_corrected_nodrop_rint = `(Intercept)`,
         anti_perc_error_corrected_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 3. Average latency on correct trials
# fit the model with smooth term for age
model_formula <- as.formula("anti_avg_lat_correct_trials ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_anti_avg_lat_correct_trials <- gamm(as.formula(model_formula), 
                                         random = list(lunaid=~age_c),
                                         data = dat, 
                                         method = 'REML')
summary(gamm_anti_avg_lat_correct_trials$gam)
summary(gamm_anti_avg_lat_correct_trials$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_anti_avg_lat_correct_trials$gam)

# identify significant periods of developmental change
anti_avg_lat_correct_trials_devchange <- calc_dev_change(gamm_anti_avg_lat_correct_trials)
anti_avg_lat_correct_trials_devchange

# plot
anti_avg_lat_correct_trials_plot <- plot_gamm(
  model = gamm_anti_avg_lat_correct_trials,
  y_axis_label = 'Average latency correct trials',
  sig = TRUE)
anti_avg_lat_correct_trials_plot

anti_avg_lat_correct_trials_plot <- anti_avg_lat_correct_trials_plot +
  annotate(geom = 'rect', 
           xmin = anti_avg_lat_correct_trials_devchange[1], 
           xmax = anti_avg_lat_correct_trials_devchange[2], 
           ymin = 300, 
           ymax = 310, 
           fill = '#476FD1')
anti_avg_lat_correct_trials_plot

# fit the model w/ fixed effects of sex
model_formula <- as.formula("anti_avg_lat_correct_trials ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_anti_avg_lat_correct_trials_sex <- gamm(as.formula(model_formula), 
                                             random = list(lunaid=~age_c),
                                             data = dat, 
                                             method = 'REML')
summary(gamm_anti_avg_lat_correct_trials_sex$gam)
summary(gamm_anti_avg_lat_correct_trials_sex$lme)

# fit the model w/ fixed effects of maternal education level
model_formula <- as.formula("anti_avg_lat_correct_trials ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_anti_avg_lat_correct_trials_medu <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_anti_avg_lat_correct_trials_medu$gam)
summary(gamm_anti_avg_lat_correct_trials_medu$lme)


# save random slopes and intercepts from final model
ranef_anti_avg_lat_correct_trials <- ranef(gamm_anti_avg_lat_correct_trials$lme)$lunaid %>%
  rename(anti_avg_lat_correct_trials_rint = `(Intercept)`,
         anti_avg_lat_correct_trials_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 4. Average latency on error-corrected trials
# fit the model with smooth term for age
model_formula <- as.formula("anti_avg_lat_error_corrected_trials ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_anti_avg_lat_error_corrected_trials <- gamm(model_formula,
                                                 random = list(lunaid=~age_c),
                                                 data = dat,
                                                 method = 'REML')
summary(gamm_anti_avg_lat_error_corrected_trials$gam)
summary(gamm_anti_avg_lat_error_corrected_trials$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_anti_avg_lat_error_corrected_trials$gam)

# plot
anti_avg_lat_error_corrected_trials_plot <- plot_gamm(
  model = gamm_anti_avg_lat_error_corrected_trials,
  y_axis_label = 'Average latency error-corrected trials',
  sig = FALSE)
anti_avg_lat_error_corrected_trials_plot


# fit the model w/ fixed effects of sex 
model_formula <- as.formula("anti_avg_lat_error_corrected_trials ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_anti_avg_lat_error_corrected_trials_sex <- gamm(as.formula(model_formula), 
                                                     random = list(lunaid=~age_c),
                                                     data = dat, 
                                                     method = 'REML')
summary(gamm_anti_avg_lat_error_corrected_trials_sex$gam)
summary(gamm_anti_avg_lat_error_corrected_trials_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("anti_avg_lat_error_corrected_trials ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_anti_avg_lat_error_corrected_trials_medu <- gamm(as.formula(model_formula), 
                                                      random = list(lunaid=~age_c),
                                                      data = dat, 
                                                      method = 'REML')
summary(gamm_anti_avg_lat_error_corrected_trials_medu$gam)
summary(gamm_anti_avg_lat_error_corrected_trials_medu$lme)


# save random slopes and intercepts from final model
ranef_anti_avg_lat_error_corrected_trials <- ranef(gamm_anti_avg_lat_error_corrected_trials$lme)$lunaid %>%
  rename(anti_avg_lat_error_corrected_trials_rint = `(Intercept)`,
         anti_avg_lat_error_corrected_trials_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


### II. Brain function outcomes - antisaccade correct trials vs. baseline 
## 1. SEF
# fit the model with smooth term for age
model_formula <- as.formula("SEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_SEF_anti_corr_v_baseline <- gamm(model_formula,
                                      random = list(lunaid=~age_c),
                                      data = dat,
                                      method = 'REML')
summary(gamm_SEF_anti_corr_v_baseline$gam)
summary(gamm_SEF_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_SEF_anti_corr_v_baseline$gam)

# plot
SEF_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_SEF_anti_corr_v_baseline,
  y_axis_label = 'SEF percent signal change',
  sig = FALSE)
SEF_anti_corr_v_baseline_plot


# fit the model w/ fixed effects of sex 
model_formula <- as.formula("SEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_SEF_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                          random = list(lunaid=~age_c),
                                          data = dat, 
                                          method = 'REML')
summary(gamm_SEF_anti_corr_v_baseline_sex$gam)
summary(gamm_SEF_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("SEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_SEF_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                           random = list(lunaid=~age_c),
                                           data = dat, 
                                           method = 'REML')
summary(gamm_SEF_anti_corr_v_baseline_medu$gam)
summary(gamm_SEF_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts
ranef_SEF_anti_corr_v_baseline <- ranef(gamm_SEF_anti_corr_v_baseline$lme)$lunaid %>%
  rename(SEF_anti_corr_v_baseline_rint = `(Intercept)`,
         SEF_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 2. pre-SMA
# fit the model with smooth term for age
model_formula <- as.formula("pre_SMA_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_pre_SMA_anti_corr_v_baseline <- gamm(model_formula,
                                          random = list(lunaid=~age_c),
                                          data = dat,
                                          method = 'REML')
summary(gamm_pre_SMA_anti_corr_v_baseline$gam)
summary(gamm_pre_SMA_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_pre_SMA_anti_corr_v_baseline$gam)

# plot
pre_SMA_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_pre_SMA_anti_corr_v_baseline,
  y_axis_label = 'Pre-SMA percent signal change',
  sig = FALSE)
pre_SMA_anti_corr_v_baseline_plot


# fit the model w/ fixed effects of sex 
model_formula <- as.formula("pre_SMA_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_pre_SMA_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_pre_SMA_anti_corr_v_baseline_sex$gam)
summary(gamm_pre_SMA_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("pre_SMA_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_pre_SMA_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                               random = list(lunaid=~age_c),
                                               data = dat, 
                                               method = 'REML')
summary(gamm_pre_SMA_anti_corr_v_baseline_medu$gam)
summary(gamm_pre_SMA_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_pre_SMA_anti_corr_v_baseline <- ranef(gamm_pre_SMA_anti_corr_v_baseline$lme)$lunaid %>%
  rename(pre_SMA_anti_corr_v_baseline_rint = `(Intercept)`,
         pre_SMA_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 3. L FEF
# fit the model with smooth term for age
model_formula <- as.formula("L_FEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_L_FEF_anti_corr_v_baseline <- gamm(model_formula,
                                        random = list(lunaid=~age_c),
                                        data = dat,
                                        method = 'REML')
summary(gamm_L_FEF_anti_corr_v_baseline$gam)
summary(gamm_L_FEF_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_L_FEF_anti_corr_v_baseline$gam)

# identify significant periods of developmental change
L_FEF_devchange <- calc_dev_change(gamm_L_FEF_anti_corr_v_baseline)
L_FEF_devchange

# plot
L_FEF_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_L_FEF_anti_corr_v_baseline,
  y_axis_label = 'L FEF percent signal change')
L_FEF_anti_corr_v_baseline_plot

L_FEF_anti_corr_v_baseline_plot <- L_FEF_anti_corr_v_baseline_plot +
  annotate(geom = 'rect', 
           xmin = L_FEF_devchange[1], 
           xmax = L_FEF_devchange[2], 
           ymin =  -5, ymax = -4.5, 
           fill = '#476FD1')
L_FEF_anti_corr_v_baseline_plot


# fit the model w/ fixed effects of sex 
model_formula <- as.formula("L_FEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_L_FEF_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                            random = list(lunaid=~age_c),
                                            data = dat, 
                                            method = 'REML')
summary(gamm_L_FEF_anti_corr_v_baseline_sex$gam)
summary(gamm_L_FEF_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("L_FEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_L_FEF_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                             random = list(lunaid=~age_c),
                                             data = dat, 
                                             method = 'REML')
summary(gamm_L_FEF_anti_corr_v_baseline_medu$gam)
summary(gamm_L_FEF_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_L_FEF_anti_corr_v_baseline <- ranef(gamm_L_FEF_anti_corr_v_baseline$lme)$lunaid %>%
  rename(L_FEF_anti_corr_v_baseline_rint = `(Intercept)`,
         L_FEF_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 4. R FEF
# fit the model with smooth term for age
model_formula <- as.formula("R_FEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_R_FEF_anti_corr_v_baseline <- gamm(model_formula,
                                        random = list(lunaid=~age_c),
                                        data = dat,
                                        method = 'REML',
                                        control = lmeControl(maxIter = 50))
summary(gamm_R_FEF_anti_corr_v_baseline$gam)
summary(gamm_R_FEF_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_R_FEF_anti_corr_v_baseline$gam)

# plot
R_FEF_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_R_FEF_anti_corr_v_baseline,
  y_axis_label = 'R FEF percent signal change',
  sig = FALSE)
R_FEF_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of sex 
model_formula <- as.formula("R_FEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_R_FEF_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                            random = list(lunaid=~age_c),
                                            data = dat, 
                                            method = 'REML',
                                            control = lmeControl(opt = 'optim'))
summary(gamm_R_FEF_anti_corr_v_baseline_sex$gam)
summary(gamm_R_FEF_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level
model_formula <- as.formula("R_FEF_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_R_FEF_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                             random = list(lunaid=~age_c),
                                             data = dat, 
                                             method = 'REML',
                                             control = lmeControl(opt = 'optim'))
summary(gamm_R_FEF_anti_corr_v_baseline_medu$gam)
summary(gamm_R_FEF_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts
ranef_R_FEF_anti_corr_v_baseline <- ranef(gamm_R_FEF_anti_corr_v_baseline$lme)$lunaid %>%
  rename(R_FEF_anti_corr_v_baseline_rint = `(Intercept)`,
         R_FEF_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 5. L putamen
# fit the model with smooth term for age
model_formula <- as.formula("L_putamen_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_L_putamen_anti_corr_v_baseline <- gamm(model_formula,
                                            random = list(lunaid=~age_c),
                                            data = dat,
                                            method = 'REML',
                                            control = lmeControl(maxIter = 50))
summary(gamm_L_putamen_anti_corr_v_baseline$gam)
summary(gamm_L_putamen_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_L_putamen_anti_corr_v_baseline$gam)

# plot
L_putamen_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_L_putamen_anti_corr_v_baseline,
  y_axis_label = 'L putamen percent signal change',
  sig = FALSE)
L_putamen_anti_corr_v_baseline_plot


# fit the model w/ fixed effects of sex 
model_formula <- as.formula("L_putamen_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_L_putamen_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                                random = list(lunaid=~age_c),
                                                data = dat, 
                                                method = 'REML')
summary(gamm_L_putamen_anti_corr_v_baseline_sex$gam)
summary(gamm_L_putamen_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level
model_formula <- as.formula("L_putamen_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_L_putamen_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                                 random = list(lunaid=~age_c),
                                                 data = dat, 
                                                 method = 'REML')
summary(gamm_L_putamen_anti_corr_v_baseline_medu$gam)
summary(gamm_L_putamen_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts
ranef_L_putamen_anti_corr_v_baseline <- ranef(gamm_L_putamen_anti_corr_v_baseline$lme)$lunaid %>%
  rename(L_putamen_anti_corr_v_baseline_rint = `(Intercept)`,
         L_putamen_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 6. R putamen
# fit the model with smooth term for age
model_formula <- as.formula("R_putamen_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_R_putamen_anti_corr_v_baseline <- gamm(model_formula,
                                            random = list(lunaid=~age_c),
                                            data = dat,
                                            method = 'REML',
                                            control = lmeControl(maxIter = 50))
summary(gamm_R_putamen_anti_corr_v_baseline$gam)
summary(gamm_R_putamen_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_R_putamen_anti_corr_v_baseline$gam)

# plot
R_putamen_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_R_putamen_anti_corr_v_baseline,
  y_axis_label = 'R putamen percent signal change',
  sig = FALSE
)
R_putamen_anti_corr_v_baseline_plot


# fit the model w/ fixed effects of sex
model_formula <- as.formula("R_putamen_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_R_putamen_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                                random = list(lunaid=~age_c),
                                                data = dat, 
                                                method = 'REML')
summary(gamm_R_putamen_anti_corr_v_baseline_sex$gam)
summary(gamm_R_putamen_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("R_putamen_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_R_putamen_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                                 random = list(lunaid=~age_c),
                                                 data = dat, 
                                                 method = 'REML')
summary(gamm_R_putamen_anti_corr_v_baseline_medu$gam)
summary(gamm_R_putamen_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_R_putamen_anti_corr_v_baseline <- ranef(gamm_R_putamen_anti_corr_v_baseline$lme)$lunaid %>%
  rename(R_putamen_anti_corr_v_baseline_rint = `(Intercept)`,
         R_putamen_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 7. L pPC
# fit the model with smooth term for age
model_formula <- as.formula("L_pPC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_L_pPC_anti_corr_v_baseline <- gamm(model_formula,
                                        random = list(lunaid=~age_c),
                                        data = dat,
                                        method = 'REML',
                                        control = lmeControl(maxIter = 50))
summary(gamm_L_pPC_anti_corr_v_baseline$gam)
summary(gamm_L_pPC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_L_pPC_anti_corr_v_baseline$gam)

# identify significant periods of developmental change
L_pPC_devchange <- calc_dev_change(gamm_L_pPC_anti_corr_v_baseline)
L_pPC_devchange

# plot
L_pPC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_L_pPC_anti_corr_v_baseline,
  y_axis_label = 'L pPC percent signal change',
  sig = TRUE)
L_pPC_anti_corr_v_baseline_plot

L_pPC_anti_corr_v_baseline_plot <- L_pPC_anti_corr_v_baseline_plot +
  annotate(geom = 'rect', 
           xmin = L_pPC_devchange[1], 
           xmax = L_pPC_devchange[2], 
           ymin = -8, 
           ymax = -7.5, 
           fill = '#476FD1')
L_pPC_anti_corr_v_baseline_plot


# fit the model w/ fixed effects of sex
model_formula <- as.formula("L_pPC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_L_pPC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                            random = list(lunaid=~age_c),
                                            data = dat, 
                                            method = 'REML')
summary(gamm_L_pPC_anti_corr_v_baseline_sex$gam)
summary(gamm_L_pPC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("L_pPC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_L_pPC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                             random = list(lunaid=~age_c),
                                             data = dat, 
                                             method = 'REML')
summary(gamm_L_pPC_anti_corr_v_baseline_medu$gam)
summary(gamm_L_pPC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_L_pPC_anti_corr_v_baseline <- ranef(gamm_L_pPC_anti_corr_v_baseline$lme)$lunaid %>%
  rename(L_pPC_anti_corr_v_baseline_rint = `(Intercept)`,
         L_pPC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 8. R pPC
# fit the model with smooth term for age
model_formula <- as.formula("R_pPC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_R_pPC_anti_corr_v_baseline <- gamm(model_formula,
                                        random = list(lunaid=~age_c),
                                        data = dat,
                                        method = 'REML',
                                        control = lmeControl(maxIter = 50))
summary(gamm_R_pPC_anti_corr_v_baseline$gam)
summary(gamm_R_pPC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_R_pPC_anti_corr_v_baseline$gam)

# identify significant periods of developmental change
R_pPC_devchange <- calc_dev_change(gamm_R_pPC_anti_corr_v_baseline)
R_pPC_devchange

# plot
R_pPC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_R_pPC_anti_corr_v_baseline,
  y_axis_label = 'R pPC percent signal change',
  sig = TRUE)
R_pPC_anti_corr_v_baseline_plot

R_pPC_anti_corr_v_baseline_plot <- R_pPC_anti_corr_v_baseline_plot +
  annotate(geom = 'rect', 
           xmin = R_pPC_devchange[1], 
           xmax = R_pPC_devchange[2], 
           ymin = -8.5, 
           ymax = -8, 
           fill = '#476FD1')
R_pPC_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of sex
model_formula <- as.formula("R_pPC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_R_pPC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                            random = list(lunaid=~age_c),
                                            data = dat, 
                                            method = 'REML')
summary(gamm_R_pPC_anti_corr_v_baseline_sex$gam)
summary(gamm_R_pPC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("R_pPC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_R_pPC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                             random = list(lunaid=~age_c),
                                             data = dat, 
                                             method = 'REML')
summary(gamm_R_pPC_anti_corr_v_baseline_medu$gam)
summary(gamm_R_pPC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts
ranef_R_pPC_anti_corr_v_baseline <- ranef(gamm_R_pPC_anti_corr_v_baseline$lme)$lunaid %>%
  rename(R_pPC_anti_corr_v_baseline_rint = `(Intercept)`,
         R_pPC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 9. L dlPFC
# fit the model with smooth term for age
model_formula <- as.formula("L_dlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_L_dlPFC_anti_corr_v_baseline <- gamm(model_formula,
                                          random = list(lunaid=~age_c),
                                          data = dat,
                                          method = 'REML',
                                          control = lmeControl(maxIter = 50))
summary(gamm_L_dlPFC_anti_corr_v_baseline$gam)
summary(gamm_L_dlPFC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_L_dlPFC_anti_corr_v_baseline$gam)

# plot
L_dlPFC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_L_dlPFC_anti_corr_v_baseline,
  y_axis_label = 'L dlPFC percent signal change',
  sig = FALSE)
L_dlPFC_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of sex
model_formula <- as.formula("L_dlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_L_dlPFC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_L_dlPFC_anti_corr_v_baseline_sex$gam)
summary(gamm_L_dlPFC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("L_dlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_L_dlPFC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                               random = list(lunaid=~age_c),
                                               data = dat, 
                                               method = 'REML')
summary(gamm_L_dlPFC_anti_corr_v_baseline_medu$gam)
summary(gamm_L_dlPFC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts
ranef_L_dlPFC_anti_corr_v_baseline <- ranef(gamm_L_dlPFC_anti_corr_v_baseline$lme)$lunaid %>%
  rename(L_dlPFC_anti_corr_v_baseline_rint = `(Intercept)`,
         L_dlPFC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 10. R dlPFC
# fit the model with smooth term for age
model_formula <- as.formula("R_dlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_R_dlPFC_anti_corr_v_baseline <- gamm(model_formula,
                                          random = list(lunaid=~age_c),
                                          data = dat,
                                          method = 'REML',
                                          control = lmeControl(maxIter = 50))
summary(gamm_R_dlPFC_anti_corr_v_baseline$gam)
summary(gamm_R_dlPFC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_R_dlPFC_anti_corr_v_baseline$gam)

# identify significant periods of developmental change
R_dlPFC_devchange <- calc_dev_change(gamm_R_dlPFC_anti_corr_v_baseline)
R_dlPFC_devchange

# plot
R_dlPFC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_R_dlPFC_anti_corr_v_baseline,
  y_axis_label = 'R dlPFC percent signal change',
  sig = TRUE)
R_dlPFC_anti_corr_v_baseline_plot

R_dlPFC_anti_corr_v_baseline_plot <- R_dlPFC_anti_corr_v_baseline_plot + 
  annotate(geom = 'rect', 
           xmin = R_dlPFC_devchange[1], 
           xmax = R_dlPFC_devchange[2], 
           ymin = -8.5, ymax = -8, 
           fill = '#476FD1')
R_dlPFC_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of sex 
model_formula <- as.formula("R_dlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_R_dlPFC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_R_dlPFC_anti_corr_v_baseline_sex$gam)
summary(gamm_R_dlPFC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("R_dlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_R_dlPFC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                               random = list(lunaid=~age_c),
                                               data = dat, 
                                               method = 'REML')
summary(gamm_R_dlPFC_anti_corr_v_baseline_medu$gam)
summary(gamm_R_dlPFC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_R_dlPFC_anti_corr_v_baseline <- ranef(gamm_R_dlPFC_anti_corr_v_baseline$lme)$lunaid %>%
  rename(R_dlPFC_anti_corr_v_baseline_rint = `(Intercept)`,
         R_dlPFC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 11. L vlPFC
# fit the model with smooth term for age
model_formula <- as.formula("L_vlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_L_vlPFC_anti_corr_v_baseline <- gamm(model_formula,
                                          random = list(lunaid=~age_c),
                                          data = dat,
                                          method = 'REML',
                                          control = lmeControl(maxIter = 100, msMaxIter = 100, niterEM = 50, opt = 'optim'))
summary(gamm_L_vlPFC_anti_corr_v_baseline$gam)
summary(gamm_L_vlPFC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_L_vlPFC_anti_corr_v_baseline$gam)

# plot
L_vlPFC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_L_vlPFC_anti_corr_v_baseline,
  y_axis_label = 'L vlPFC percent signal change',
  sig = FALSE)
L_vlPFC_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of sex 
model_formula <- as.formula("L_vlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_L_vlPFC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML',
                                              control = lmeControl(opt = 'optim'))
summary(gamm_L_vlPFC_anti_corr_v_baseline_sex$gam)
summary(gamm_L_vlPFC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("L_vlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_L_vlPFC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                               random = list(lunaid=~age_c),
                                               data = dat, 
                                               method = 'REML',
                                               control = lmeControl(opt = 'optim'))
summary(gamm_L_vlPFC_anti_corr_v_baseline_medu$gam)
summary(gamm_L_vlPFC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_L_vlPFC_anti_corr_v_baseline <- ranef(gamm_L_vlPFC_anti_corr_v_baseline$lme)$lunaid %>%
  rename(L_vlPFC_anti_corr_v_baseline_rint = `(Intercept)`,
         L_vlPFC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 12. R vlPFC
# fit the model with smooth term for age
model_formula <- as.formula("R_vlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_R_vlPFC_anti_corr_v_baseline <- gamm(model_formula,
                                          random = list(lunaid=~age_c),
                                          data = dat,
                                          method = 'REML',
                                          control = lmeControl(maxIter = 50))
summary(gamm_R_vlPFC_anti_corr_v_baseline$gam)
summary(gamm_R_vlPFC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_R_vlPFC_anti_corr_v_baseline$gam)

# plot
R_vlPFC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_R_vlPFC_anti_corr_v_baseline,
  y_axis_label = 'R vlPFC percent signal change',
  sig = FALSE)
R_vlPFC_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of sex 
model_formula <- as.formula("R_vlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_R_vlPFC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_R_vlPFC_anti_corr_v_baseline_sex$gam)
summary(gamm_R_vlPFC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("R_vlPFC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_R_vlPFC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                               random = list(lunaid=~age_c),
                                               data = dat, 
                                               method = 'REML')
summary(gamm_R_vlPFC_anti_corr_v_baseline_medu$gam)
summary(gamm_R_vlPFC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts from final model
ranef_R_vlPFC_anti_corr_v_baseline <- ranef(gamm_R_vlPFC_anti_corr_v_baseline$lme)$lunaid %>%
  rename(R_vlPFC_anti_corr_v_baseline_rint = `(Intercept)`,
         R_vlPFC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## 13. dACC
# fit the model with smooth term for age
model_formula <- as.formula("dACC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_dACC_anti_corr_v_baseline <- gamm(model_formula,
                                       random = list(lunaid=~age_c),
                                       data = dat,
                                       method = 'REML',
                                       control = lmeControl(maxIter = 50))
summary(gamm_dACC_anti_corr_v_baseline$gam)
summary(gamm_dACC_anti_corr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_dACC_anti_corr_v_baseline$gam)

# plot
dACC_anti_corr_v_baseline_plot <- plot_gamm(
  model = gamm_dACC_anti_corr_v_baseline,
  y_axis_label = 'dACC percent signal change correct trials',
  sig = FALSE)
dACC_anti_corr_v_baseline_plot

# fit the model w/ fixed effects of gender 
model_formula <- as.formula("dACC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_dACC_anti_corr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                           random = list(lunaid=~age_c),
                                           data = dat, 
                                           method = 'REML')
summary(gamm_dACC_anti_corr_v_baseline_sex$gam)
summary(gamm_dACC_anti_corr_v_baseline_sex$lme)

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("dACC_anti_corr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_dACC_anti_corr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                            random = list(lunaid=~age_c),
                                            data = dat, 
                                            method = 'REML')
summary(gamm_dACC_anti_corr_v_baseline_medu$gam)
summary(gamm_dACC_anti_corr_v_baseline_medu$lme)

# save random slopes and intercepts, adjusting for maternal education
ranef_dACC_anti_corr_v_baseline <- ranef(gamm_dACC_anti_corr_v_baseline_medu$lme)$lunaid %>%
  rename(dACC_anti_corr_v_baseline_rint = `(Intercept)`,
         dACC_anti_corr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


### III. Brain function outcomes - antisaccade error-corrected trials vs. baseline
## 1. dACC
# fit the model with smooth term for age
model_formula <- as.formula("dACC_anti_errcorr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp')")
gamm_dACC_anti_errcorr_v_baseline <- gamm(model_formula,
                                          random = list(lunaid=~age_c),
                                          data = dat,
                                          method = 'REML',
                                          control = lmeControl(opt = 'optim'))
summary(gamm_dACC_anti_errcorr_v_baseline$gam)
summary(gamm_dACC_anti_errcorr_v_baseline$lme)

# model diagnostics
par(mfrow = c(2,2))
gam.check(gamm_dACC_anti_errcorr_v_baseline$gam)

# identify significant periods of developmental change
dACC_errcorr_devchange <- calc_dev_change(gamm_dACC_anti_errcorr_v_baseline)
dACC_errcorr_devchange

# plot
dACC_anti_errcorr_v_baseline_plot <- plot_gamm(
  model = gamm_dACC_anti_errcorr_v_baseline,
  y_axis_label = 'dACC percent signal change error-corrected trials',
  sig = TRUE)
dACC_anti_errcorr_v_baseline_plot

# fit the model w/ fixed effects of sex 
model_formula <- as.formula("dACC_anti_errcorr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + sex")
gamm_dACC_anti_errcorr_v_baseline_sex <- gamm(as.formula(model_formula), 
                                              random = list(lunaid=~age_c),
                                              data = dat, 
                                              method = 'REML')
summary(gamm_dACC_anti_errcorr_v_baseline_sex$gam)
summary(gamm_dACC_anti_errcorr_v_baseline_sex$lme)

# identify significant periods of developmental change, adjusting for sex
dACC_errcorr_sex_devchange <- calc_dev_change(gamm_dACC_anti_errcorr_v_baseline_sex)
dACC_errcorr_sex_devchange

dACC_anti_errcorr_v_baseline_plot <- dACC_anti_errcorr_v_baseline_plot +
  annotate(geom = 'rect', 
           xmin = dACC_errcorr_sex_devchange[1],
           xmax = dACC_errcorr_sex_devchange[2], 
           ymin =  -25, 
           ymax = -24, 
           fill = '#476FD1')
dACC_anti_errcorr_v_baseline_plot

# fit the model w/ fixed effects of maternal education level 
model_formula <- as.formula("dACC_anti_errcorr_v_baseline ~ s(age_c, k = 10, fx = FALSE, bs = 'tp') + level_edu_mother_cat")
gamm_dACC_anti_errcorr_v_baseline_medu <- gamm(as.formula(model_formula), 
                                               random = list(lunaid=~age_c),
                                               data = dat, 
                                               method = 'REML')
summary(gamm_dACC_anti_errcorr_v_baseline_medu$gam)
summary(gamm_dACC_anti_errcorr_v_baseline_medu$lme)

# save random slopes and intercepts, adjusting for sex
ranef_dACC_anti_errcorr_v_baseline <- ranef(gamm_dACC_anti_errcorr_v_baseline_sex$lme)$lunaid %>%
  rename(dACC_anti_errcorr_v_baseline_rint = `(Intercept)`,
         dACC_anti_errcorr_v_baseline_rslope = age_c) %>%
  tibble::rownames_to_column(., 'lunaid') %>%
  mutate(lunaid = as.factor(gsub('1/', '', lunaid)))


## test of global effect of maternal education level for each model
mods <- c('gamm_anti_perc_correct_nodrop',
          'gamm_anti_perc_error_corrected_nodrop',
          'gamm_anti_avg_lat_correct_trials',
          "gamm_anti_avg_lat_error_corrected_trials",
          "gamm_SEF_anti_corr_v_baseline",
          "gamm_pre_SMA_anti_corr_v_baseline",
          "gamm_L_FEF_anti_corr_v_baseline",
          "gamm_R_FEF_anti_corr_v_baseline",    
          "gamm_L_putamen_anti_corr_v_baseline",
          "gamm_R_putamen_anti_corr_v_baseline",
          "gamm_L_pPC_anti_corr_v_baseline",
          "gamm_R_pPC_anti_corr_v_baseline",    
          "gamm_L_dlPFC_anti_corr_v_baseline",
          "gamm_R_dlPFC_anti_corr_v_baseline",
          "gamm_L_vlPFC_anti_corr_v_baseline",
          "gamm_R_vlPFC_anti_corr_v_baseline",  
          "gamm_dACC_anti_corr_v_baseline",
          "gamm_dACC_anti_errcorr_v_baseline")

global_test_medu <- matrix(nrow = length(mods), ncol = 4)
for(i in 1:length(mods)){
  mod_medu <- paste0(mods[i], '_medu')
  global_test_medu[i, 1] <- mod_medu
  global_test_medu[i, 2:4] <- anova(get(mod_medu)$gam)$pTerms.table[1:3]    
}
colnames(global_test_medu) <- c('model', 'df', 'F', 'p-value for global medu')
global_test_medu
print(xtable(global_test_medu, type = "latex"), 
      file = "gamm_test_global_sig_medu_terms.tex")

## tests of effect of sex for each model
test_sex <- matrix(nrow = length(mods), ncol = 4)
for (i in 1:length(mods)){
  mod_sex <- paste0(mods[i], '_sex')
  test_sex[i, 1] <- mod_sex
  test_sex[i, 2:4] <- summary(get(mod_sex)$gam)$p.table[2, c(1, 3:4)]
}
colnames(test_sex) <- c('model', 'estimate', 't', 'p-value for sex')
test_sex
print(xtable(test_sex, type = "latex"), 
      file = "gamm_test_sig_sex_terms.tex")


### IV. Compile results and plots
# summary of significance of smooth terms for age for final models
mods <- c('gamm_anti_perc_correct_nodrop',
          'gamm_anti_perc_error_corrected_nodrop',
          'gamm_anti_avg_lat_correct_trials',
          "gamm_anti_avg_lat_error_corrected_trials",
          "gamm_SEF_anti_corr_v_baseline",
          "gamm_pre_SMA_anti_corr_v_baseline",
          "gamm_L_FEF_anti_corr_v_baseline",
          "gamm_R_FEF_anti_corr_v_baseline",    
          "gamm_L_putamen_anti_corr_v_baseline",
          "gamm_R_putamen_anti_corr_v_baseline",
          "gamm_L_pPC_anti_corr_v_baseline",
          "gamm_R_pPC_anti_corr_v_baseline",    
          "gamm_L_dlPFC_anti_corr_v_baseline",
          "gamm_R_dlPFC_anti_corr_v_baseline",
          "gamm_L_vlPFC_anti_corr_v_baseline",
          "gamm_R_vlPFC_anti_corr_v_baseline",  
          "gamm_dACC_anti_corr_v_baseline",
          "gamm_dACC_anti_errcorr_v_baseline_sex")
sumtab <- sapply(mods, function (x) summary(get(x)$gam)$s.table)
rownames(sumtab) <- c('edf', 'Ref.df', 'F', 'p-value')
sumtab <- as.data.frame(t(sumtab)) %>%
  tibble::rownames_to_column(., 'model') %>%
  select(-c('Ref.df'))
sumtab <- cbind(sumtab, p.adjust(sumtab$`p-value`, method = 'fdr'))
names(sumtab)[ncol(sumtab)] <- 'q-value'
sumtab <- sumtab %>% 
  mutate(survive_multcomp = case_when(`q-value` < 0.05 ~ 1, 
                                      `q-value` >= 0.05 ~ 0))
print(xtable(sumtab %>% select(-c(survive_multcomp)), 
             digits = c(0,0,2,2,7,7), 
             type = "latex"), 
      file = "sumtab_age_terms.tex")

# merge random effects and save for further analysis
random_effects_list <- do.call('list', mget(grep('ranef', names(.GlobalEnv), value = TRUE)))
random_effects <- random_effects_list %>%
  purrr::reduce(left_join, by = 'lunaid')
write.csv(random_effects, 'Data/random_effects_from_gamms_20200312.csv')

## plots
# behavioral (accuracy and latency) plot
behav_plot <- grid.arrange(
  anti_perc_correct_plot, 
  anti_avg_lat_correct_trials_plot, 
  anti_perc_error_corrected_nodrop_plot, 
  anti_avg_lat_error_corrected_trials_plot, 
  ncol = 2)
behav_plot
ggsave(filename = 'Plots/behav_plot.png', 
       plot = behav_plot,
       width = 8, height = 7, units = 'in')

# motor response ROIs
motor_resp_rois_plot <- grid.arrange(
  SEF_anti_corr_v_baseline_plot,
  pre_SMA_anti_corr_v_baseline_plot,
  L_FEF_anti_corr_v_baseline_plot,
  R_FEF_anti_corr_v_baseline_plot, 
  L_pPC_anti_corr_v_baseline_plot, 
  R_pPC_anti_corr_v_baseline_plot,
  L_putamen_anti_corr_v_baseline_plot, 
  R_putamen_anti_corr_v_baseline_plot, 
  ncol = 2)
motor_resp_rois_plot
ggsave(filename = 'Plots/motor_rois.png', 
       plot = motor_resp_rois_plot,
       width = 8, height = 10.5, units = 'in')

# executive control rois
exec_rois_plot <- grid.arrange(
  L_vlPFC_anti_corr_v_baseline_plot,            
  R_vlPFC_anti_corr_v_baseline_plot, 
  L_dlPFC_anti_corr_v_baseline_plot,
  R_dlPFC_anti_corr_v_baseline_plot, 
  ncol = 2)
exec_rois_plot
ggsave(filename = 'Plots/exec_rois.png', 
       plot = exec_rois_plot,
       width = 7, height = 6, units = 'in')

# dACC roi
dACC_plot <- grid.arrange(
  dACC_anti_corr_v_baseline_plot,
  dACC_anti_errcorr_v_baseline_plot,  
  ncol = 2)
dACC_plot
ggsave(filename = 'Plots/dACC.png', 
       plot = dACC_plot,
       width = 8, height = 4, units = 'in')
























