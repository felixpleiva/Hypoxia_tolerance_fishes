################################################################################
# Variance decomposition method for brms
# Author: Kerstin Pierick (kerstin.pierick@uni-goettingen.de). Thanks Kerstin!!
# Modifications: Félix P Leiva (felixpleiva@gmail.com)
# Date: 20210824
################################################################################
rm(list=ls()) #clear the work environment
today<-format(Sys.Date(),"%Y%m%d") #setting the date
################################################################################
setwd("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs") #Felix's lab PC
getwd()#to check
################################################################################
#Libraries
library(UpSetR)

# models
best_model_final <- readRDS("best_model_final.rds")

mod_F <- brm(
  pcrit_kpa~temp_test * c_value_log +
            temp_test * mass_log +
            sal_test + res_metab_rate + 
            (1 | gr(animal, cov = B)) + (1 | species),
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

mod_M <- brm(
  pcrit_kpa~1 + 
    (1 | gr(animal, cov = B)) + (1 | species),
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
#  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

mod_P <- brm(
  pcrit_kpa~temp_test * c_value_log +
    temp_test * mass_log +
    sal_test + res_metab_rate + 
    #(1 | gr(animal, cov = C)) +
    (1 | species),
  data = d, 
  family = gaussian(), 
  data2 = list(C = C),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

mod_PM <- brm(
  pcrit_kpa~1 + 
    #(1 | gr(animal, cov = C)) +
    (1 | species),
  data = d, 
  family = gaussian(), 
  data2 = list(C = C),
#  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

mod_PSp <- brm(
  pcrit_kpa~temp_test * c_value_log +
    temp_test * mass_log +
    sal_test + res_metab_rate,
    #(1 | gr(animal, cov = C)) +
#    (1 | species),
  data = d, 
  family = gaussian(), 
  data2 = list(C = C),
#  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

mod_Sp <- brm(
  pcrit_kpa~temp_test * c_value_log +
    temp_test * mass_log +
    sal_test + res_metab_rate + 
    (1 | gr(animal, cov = B)),
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

mod_SpM <- brm(
  pcrit_kpa~1 + 
    (1 | gr(animal, cov = B)),
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
 # prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

#-------------------------------------------------------------------------------
#variance partitioning
##-------------------------------------------------------------------------------
# get observed values
obs <- d$pcrit

# Step 1: marginal ##########

best_model_final <- readRDS("best_model_final.rds")
best_model_final<-mod_P
best_model_final<-mod_M
best_model_final<-mod_F
best_model_final<-mod_Sp
best_model_final<-mod_PM
best_model_final<-mod_PSp

# marginal predictions
pred_mar <- posterior_linpred(best_model_final, re_formula = NA) 
pred_mar
# variance of marginal predictions
var_pred_mar <- pred_mar %>% apply(1, var)

# variance of residuals of marginal predictions
var_res_mar <- pred_mar %>% 
  sweep(2, obs, "-") %>%
  apply(1, var)

# Marginal R2
r2_mar <- var_pred_mar / (var_pred_mar + var_res_mar)
mar <- mean(r2_mar) #0.2865047

# Step 2: one RE ##########

# phylogeny
pred_phy <- posterior_linpred(best_model_final, re_formula = ~ (1|animal)) 
var_pred_phy <- pred_phy %>% apply(1, var)
var_res_phy <- pred_phy %>% 
  sweep(2, obs, "-") %>%
  apply(1, var)
r2_phy <- var_pred_phy / (var_pred_phy + var_res_phy)
phy <- mean(r2_phy)#0.5937011

# species
pred_spe <- posterior_linpred(best_model_final, re_formula = ~ (1|species)) 
var_pred_spe <- pred_spe %>% apply(1, var)
var_res_spe <- pred_spe %>% 
  sweep(2, obs, "-") %>%
  apply(1, var)
r2_spe <- var_pred_spe / (var_pred_spe + var_res_spe)
spe <- mean(r2_spe)#0.5347925

# Step 3: two REs ##########

# phylogeny + species
pred_phy_spe <- posterior_linpred(best_model_final, re_formula = ~ (1|animal) + 
                                    (1|species)) 
var_pred_phy_spe <- pred_phy_spe %>% apply(1, var)
var_res_phy_spe <- pred_phy_spe %>% 
  sweep(2, obs, "-") %>%
  apply(1, var)
r2_phy_spe <- var_pred_phy_spe / (var_pred_phy_spe + var_res_phy_spe)
phy_spe <- mean(r2_phy_spe)#0.80699

# Step 4: conditional ##########

pred_con <- posterior_linpred(best_model_final) 
var_pred_con <- pred_con %>% apply(1, var)
var_res_con <- pred_con %>% 
  sweep(2, obs, "-") %>%
  apply(1, var)
r2_con <- var_pred_con / (var_pred_con + var_res_con)
con <- mean(r2_con)#0.80699

# calculate proportions
vd <- tibble(
  ord = "phy_spe",
  mar = mean(r2_mar), 
  phy = mean(r2_phy - r2_mar), 
  spe = mean(r2_phy_spe - r2_phy), 
  res = mean(1 - r2_con))

vd <- vd %>%
  add_row(
    ord = "spe_phy",
    mar = mean(r2_mar),
    spe = mean(r2_spe - r2_mar),
    phy = mean(r2_phy_spe - r2_spe),
    res = mean(1 - r2_con)
  ) %>%
  
  mutate(sum = mar + phy + spe + res) %>%
  summarise_if(is.numeric, mean)

vd
#    mar   phy   spe   res     sum
#    0.287 0.290 0.231 0.193     1


#Marginal
1-0.195-(1-0.276)

#Phylogeny
1-0.195-(1-0.195)

#Species
1-0.195-(1-0.197)

#Marginal & Phylogeny
(0.277-0.195)-(1-0.195-(1-0.276))

#Marginal & Species
(0.280-0.195)-(1-0.195-(1-0.276))-(1-0.195-(1-0.197))

#Species & Phylogeny
(0.491-0.195)-(1-0.195-(1-0.197))

#species & Phylogeny & Marginal
(1-0.195)-(1-0.195-(1-0.276))-(1-0.195-(1-0.195))-(1-0.195-(1-0.197))-((0.277-0.195)-(1-0.195-(1-0.276)))-((0.280-0.195)-(1-0.195-(1-0.276))-(1-0.195-(1-0.197)))-(0.491-0.195)-(1-0.195-(1-0.197))

# Dataset
input <- c(
  Marginal = (1-0.195-(1-0.276)),
  Phylogeny = (1-0.195-(1-0.195)),
  Species = (1-0.195-(1-0.197)),
  "Marginal&Phylogeny" = ((0.277-0.195)-(1-0.195-(1-0.276))),
  "Marginal&Species" = ((0.280-0.195)-(1-0.195-(1-0.276))-(1-0.195-(1-0.197))),
  "Phylogeny&Species" = ((0.491-0.195)-(1-0.195-(1-0.197))),
  "Marginal&Phylogeny&Species" = ((1-0.195)-(1-0.195-(1-0.276))-(1-0.195-(1-0.195))-(1-0.195-(1-0.197))-((0.277-0.195)-(1-0.195-(1-0.276)))-((0.280-0.195)-(1-0.195-(1-0.276))-(1-0.195-(1-0.197)))-(0.491-0.195)-(1-0.195-(1-0.197)))
)


# Plot
upset(fromExpression(input*100), 
      nintersects = 6, 
      intersections = list(list("Marginal","Phylogeny","Species"),list("Marginal","Phylogeny"),list("Marginal","Species"),list("Phylogeny","Species"),list("Species"),list("Phylogeny"),list("Marginal")),
      nsets = 3, 
      order.by = "freq", 
      mainbar.y.label = "Explained variation (%)",
      decreasing = T, 
      mb.ratio = c(0.7, 0.3),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      cutoff = 0,
      sets.x.label = "Cum. expl. variation",
      line.size = 1)

