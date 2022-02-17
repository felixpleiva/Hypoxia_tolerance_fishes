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
library(brms)
library(ape)
library(dplyr)
library(phytools)
library(purrr)
################################################################################
# General STAN specifications
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(6955)

#load data
data <- read.csv("Pcrit data for 195 fish species ( 20210823 ).csv")
head(data)
names(data)

# Models are saved so that they don't need to be rerun for every session.
best_model_final <- readRDS("best_model_final.rds")

# load phylogenetic tree
tree <- read.tree("Phylogenetic tree for 195 species of fish with Pcrit data.tre")
tree$tip.label#195 species
plotTree(tree,fsize=0.5,ftype="i",type="fan")

# The tree object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

A <- ape::vcv.phylo(tree)
A
B <- ape::vcv.phylo(compute.brlen(tree, power=0.4))
#-------------------------------------------------------------------------------
# apply complete.clases (exclude NAs) for the columns of interest
d<-data[complete.cases(data[,c("Estimated_C.Value",
                               "RESPIROMETRY_TYPE",
                               "logpercentcalc",
                               "Temp",
                               "Pcrit..kPa.",
                               "Salt",
                               "animal",
                               "Taccl",
                               "residMR")]),]


# explore whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

setdiff(d$animal, as.character(tree$tip.label))#listed in our database but not in the tree
setdiff(as.character(tree$tip.label),d$animal)# listed in the tree but not in our database

# There are 36 species not included in the data frame. Lets now drop
# these species from the phylogenetic tree.
tree<-drop.tip(tree, setdiff(tree$tip.label, d$animal))

#check again
setdiff(d$animal, as.character(tree$tip.label))#listed in our database but not in the tree
setdiff(as.character(tree$tip.label),d$animal)# listed in the tree but not in our database
#perfectirijillo

#count number of species to be included in the models
length(unique(d$animal))

# create new columns
d$species<-d$animal
d$c_value_log<-log10(d$Estimated_C.Value)
d$rel_mass_log_sq<-(d$logpercentcalc)^2

#use simplest names
names(d)
d<-rename(d,resp_type=RESPIROMETRY_TYPE)
d<-rename(d,rel_mass_log=logpercentcalc)
d<-rename(d,temp_test=Temp)
d<-rename(d,pcrit_kpa=Pcrit..kPa.)
d<-rename(d,sal_test=Salt)
d<-rename(d,temp_acclim=Taccl)
d<-rename(d,res_metab_rate=residMR)
d<-rename(d,mass_log=logMass)
d$mass_log_sq<-(d$mass_log)^2
names(d)
#-------------------------------------------------------------------------------
#variance partitioning
##-------------------------------------------------------------------------------
# get observed values
obs <- d$pcrit

# Step 1: marginal ##########

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
#-------------------------------------------------------------------------------
# saving session information with all packages versions for reproducibility purposes
sink("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs/variance_partitioning_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################