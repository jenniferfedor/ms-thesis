######################################
### Bootstrap-enhanced elastic net ###
######################################

setwd('/Users/jenniferfedor/Documents/Biostats MS/Spring 2020/Thesis')

library(caret)
library(dplyr)
library(ensr)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(lubridate)
library(purrr)
library(tidyr)
library(xtable)

# response measure (whoqol score)
whoqol <- read.csv('Data/whoqol_ages_scores_20200204.csv', 
                   stringsAsFactors = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  select(-c('X'))

# predictors (random effects from GAMM models)
ranefs <- read.csv('Data/random_effects_from_gamms_20200312.csv', 
                   stringsAsFactors = FALSE) %>%
  select(-c(X))

# scale predictors and response variable
ranefs <- ranefs %>%
  mutate_at(vars(matches('anti')), scale)
head(ranefs)

# join response and predictors
dat <- whoqol %>%
  right_join(ranefs, by = 'lunaid') %>%
  # z-score WHO-QOL domain scores
  mutate(D_1_Raw_z = scale(D_1_Raw)) %>%
  mutate(D_2_Raw_z = scale(D_2_Raw)) %>%
  mutate(D_3_Raw_z = scale(D_3_Raw)) %>%
  mutate(D_4_Raw_z = scale(D_4_Raw)) %>%
  # for each subject, average their four z-scored domain scores
  rowwise() %>%
  mutate(whoqol_avg_score = mean(c(D_1_Raw_z, D_2_Raw_z, D_3_Raw_z, D_4_Raw_z))) 

whoqol_avg_score_z <- scale(dat$whoqol_avg_score)[1:50]

dat <- cbind(dat, whoqol_avg_score_z) %>%
  select(-(date:counter)) %>%
  select(-(D_1_Raw_z:whoqol_avg_score))

write.csv(dat, 'Data/data_for_elastic_net_20200317.csv')

# function to plot histogram of whoqol scores
plot_whoqol <- function(x_var, x_axis_label, binw, centering) {
  p <- ggplot(dat2, aes_string(x = x_var)) + 
    geom_histogram(color = 'black', 
                   fill = 'white', 
                   binwidth = binw, 
                   center = centering) +
    labs(x = x_axis_label, 
         y = 'Frequency') +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black')) +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 11))
  return(p)
}

# histograms of individual domain scores and composite scores
dom1 <- plot_whoqol(x_var = 'D_1_Raw', 
                    x_axis_label = 'Physical health', 
                    binw = 2, centering = 1)
dom2 <- plot_whoqol(x_var = 'D_2_Raw', 
                    x_axis_label = 'Psychological health', 
                    binw = 2, centering = 0)
dom3 <- plot_whoqol(x_var = 'D_3_Raw', 
                    x_axis_label = 'Social', 
                    binw = 1, centering = 1) + 
  scale_y_continuous(limits = c(0, 13))
dom4 <- plot_whoqol(x_var = 'D_4_Raw', 
                    x_axis_label = 'Environment', 
                    binw = 2, centering = 0) + 
  scale_y_continuous(limits = c(0, 11))
composite <- plot_whoqol(x_var = 'whoqol_avg_score_z', 
                         x_axis_label = 'Standardized composite score', 
                         binw = 0.5, centering = 0.75) + 
  scale_y_continuous(limits = c(0, 13))

# arrange histograms in single plot
layoutmatrix <- rbind(c(1,2,3),
                      c(1,2,4),
                      c(5,6,4),
                      c(5,6,7))
whoqol_hists <- gridExtra::grid.arrange(dom1, dom2, grid::nullGrob(), 
                                        composite,
                                        dom3, dom4, grid::nullGrob(),
                                        layout_matrix = layoutmatrix, 
                                        ncol = 3)
whoqol_hists
ggsave('whoqol_hists.png', whoqol_hists, width = 10, height = 7, units = 'in')


## prep data for elastic net
# split data into response y and predictors X
set.seed(100)
X <- as.matrix(dat %>% select(-c(lunaid, whoqol_avg_score_z)))
y <- as.matrix(dat %>% select(c(whoqol_avg_score_z)))

# predictors for which there were significant age effects (when controlling for FDR) in GAMMs
X_sig <- as.matrix(dat %>% 
                     select(-c(lunaid, whoqol_avg_score_z)) %>%
                     select(c(anti_perc_correct_nodrop_rint,
                              anti_perc_correct_nodrop_rslope,
                              anti_perc_error_corrected_nodrop_rint,
                              anti_perc_error_corrected_rslope,
                              anti_avg_lat_correct_trials_rint,
                              anti_avg_lat_correct_trials_rslope,
                              L_FEF_anti_corr_v_baseline_rint,
                              L_FEF_anti_corr_v_baseline_rslope,
                              L_pPC_anti_corr_v_baseline_rint,
                              L_pPC_anti_corr_v_baseline_rslope,
                              R_pPC_anti_corr_v_baseline_rint,
                              R_pPC_anti_corr_v_baseline_rslope,
                              R_dlPFC_anti_corr_v_baseline_rint,
                              R_dlPFC_anti_corr_v_baseline_rslope,
                              dACC_anti_errcorr_v_baseline_rint,
                              dACC_anti_errcorr_v_baseline_rslope)))

## fit the inital elastic net model
# sequence of alpha parameter values to try 
# does not include 0 or 1 because those values correspond to ridge and lasso
alphas <- seq(from = 0.1, to = 0.9, length.out = 17)

# create k = 10 folds for cross-validation
# explicitly setting folds to allow for reproducibility
set.seed(100)
folds <- createFolds(y = y, k = 10, list = FALSE, returnTrain = FALSE)

# simultaneous tuning of alpha and lambda
fit <- ensr(x = X_sig, 
            y = y, 
            alphas = alphas, 
            foldid = folds,
            standardize = FALSE)

# alpha and lambda values that minimized CV-MSE
summary(fit)[cvm == min(cvm)]
summary(fit)[cvm == min(cvm)]$alpha
summary(fit)[cvm == min(cvm)]$lambda

# fit the model using the alpha and lambda values that minimized CV-MSE
glmnetfit <- glmnet(x = X_sig,
                    y = y,
                    alpha = summary(fit)[cvm == min(cvm)]$alpha,
                    lambda = summary(fit)[cvm == min(cvm)]$lambda,
                    family = 'gaussian',
                    standardize = FALSE)

coef(glmnetfit) # estimated beta coefficients
glmnetfit$dev.ratio # model R^2


# organize coefficients into dataframe
coefs_glmnet_fit <- as.data.frame(as.matrix(coef(glmnetfit)))
colnames(coefs_glmnet_fit) <- 'beta_hat'
coefs_glmnet_fit <- tibble::rownames_to_column(coefs_glmnet_fit, 'predictor')
coefs_glmnet_fit


## use bootstrap-enhanced procedure to:
## (1) derive confidence intervals and 
## (2) calculate variable inclusion probabilities
# number of bootstrap samples and subjects
B <- 5000
subjs <- 1:nrow(data_for_boot) 

coefs_boot <- data.frame(Predictor = predictors, stringsAsFactors = FALSE)

# perform bootstrap procedure
for (i in 1:B){
  # resample the data with replacement
  set.seed(i*10)
  boot_sample <- sample(subjs, size = 50, replace = TRUE)
  y_boot <- y[boot_sample]
  X_sig_boot <- X_sig[boot_sample, 1:ncol(X_sig)]
  
  # fit the elastic net model on the bootstrap sample
  glmnet_boot <- glmnet(x = X_sig_boot,
                        y = y_boot,
                        alpha = summary(fit)[cvm == min(cvm)]$alpha,
                        lambda = summary(fit)[cvm == min(cvm)]$lambda,
                        family = 'gaussian',
                        standardize = FALSE)
  
  # store the estimated coefficients
  coefs <- as.data.frame(as.matrix(coef(glmnet_boot))) %>%
    tibble::rownames_to_column(var = 'Predictor')
  colnames(coefs_nocv)[2] <- paste0('Beta_boot_', i)
  
  coefs_boot <- full_join(coefs_boot, coefs, by = 'Predictor')
}

coefs_boot

# transpose bootstrap coefficient dataframe from wide to long format
coefs_boot_long_colnames <- coefs_boot_nocv[,1]
coefs_boot_long <- as.data.frame(t(coefs_boot[,-1]), stringsAsFactors = FALSE)
colnames(coefs_boot_long) <- coefs_boot_long_colnames
coefs_boot_long

# 95% bootstrapped CIs
coefs_quantiles <- as.data.frame(t(sapply(coefs_boot_long, 
                                          quantile, 
                                          probs = c(0.025, 0.975), 
                                          names = TRUE)))
coefs_quantiles <- tibble::rownames_to_column(coefs_quantiles, 'predictor')
coefs_quantiles

# variable inclusion probabilities
VIP <- as.data.frame(cbind(coefs_boot[,1], 
                           rowSums(coefs_boot[,-1] != 0)), 
                     stringsAsFactors = FALSE)
colnames(VIP) <- c('predictor', 'times_selected')
VIP$times_selected <- as.numeric(VIP$times_selected)
VIP$prop_selected <- VIP$times_selected/B 
VIP

# join the results
coefs_and_CI <- full_join(coefs_glmnet_fit, coefs_quantiles)
coefs_and_CI_and_VIP <- full_join(coefs_and_CI, VIP) %>%
  select(-c('times_selected')) %>%
  rename(VIP = prop_selected)
coefs_and_CI_and_VIP

# save the results
write.csv(coefs_and_CI_and_VIP, 'Data/bootstrap_enhanced_enet_results_20200329.csv')


## create LaTeX table to include in document
# clean up predictor names
coefs_and_CI_and_VIP$predictor <- gsub(pattern = '_', 
                                       replacement = ' ',
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'anti corr v baseline', 
                                       replacement = 'percent signal change correct trials', 
                                        x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'anti errcorr v baseline', 
                                       replacement = 'percent signal change error-corrected trials', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'rint', 
                                       replacement = 'intercept', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'rslope', 
                                       replacement = 'slope', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'anti perc', 
                                       replacement = 'Antisaccade proportion', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'error corrected', 
                                       replacement = 'error-corrected trials', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'nodrop', 
                                       replacement = '', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'anti avg lat', 
                                       replacement = 'Antisaccade average latency', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = '  ', 
                                       replacement = ' ', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'Antisaccade proportion correct', 
                                       replacement = 'Antisaccade proportion correct trials', 
                                       x = coefs_and_CI_and_VIP$predictor)
coefs_and_CI_and_VIP$predictor <- gsub(pattern = 'percent signal change', 
                                       replacement = 'activation', 
                                       x = coefs_and_CI_and_VIP$predictor)

# save latex table code
print(xtable(coefs_and_CI_and_VIP, 
             digits = c(0,0,4,4,4,4), 
             type = 'latex'), 
      file = 'bootstrap_enhanced_enet_20200329.tex')
