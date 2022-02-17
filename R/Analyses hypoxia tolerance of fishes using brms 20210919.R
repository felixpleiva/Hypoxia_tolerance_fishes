################################################################################
# Script to fit phylogenetic multilevel model using brms
# Author: Félix P Leiva (felixpleiva@gmail.com)
# Date: 20210824
# Modifications: Wilco CEP Verberk (YYYYMMDD); Félix P Leiva (20210920 )
# ------------------------------------------------------------------------------
# Cite as:

#Verberk WCEP, Leiva FP, Sandkler JF, van de Pol I, Urbina M, Wilson R,
#McKenzie DJ (2021). Draft version of paper data and code of manuscript: Hypoxia
#tolerance in fish varies with body size and genome size in a
#temperature-dependent manner (). Zenodo ().
# ------------------------------------------------------------------------------
rm(list=ls()) #clear the work environment
today<-format(Sys.Date(),"%Y%m%d") #setting the date
# ------------------------------------------------------------------------------
setwd("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs") #Felix's lab PC
#setwd("C:/Users/wcepv/surfdrive/papers in progress/pcrit fish")
getwd()#to check
# ------------------------------------------------------------------------------
#Libraries
library(brms)
library(ape)
library(dplyr)
library(phytools)
library(tidybayes)
library(bayestestR)
library(ggtree)
library(treeio)
library(ggplot2)
# ------------------------------------------------------------------------------
set.seed(6955)# we need this to replicate the results

# General STAN specifications
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Models are saved so that they don't need to be rerun for every session.
best_model_final <- readRDS("best_model_final.rds")

#load data
data <- read.csv("Pcrit data for 195 fish species ( 20210823 ).csv")
head(data)
names(data)

# load phylogenetic tree
tree <- read.tree("Phylogenetic tree for 195 species of fish with Pcrit data.tre")
tree$tip.label
plotTree(tree,fsize=0.5,ftype="i",type="fan")

# The tree object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

A <- ape::vcv.phylo(tree)
A
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
#Define priors
priors1=c(
  prior(normal(0,10), "b"),
  prior(normal(0,50), "Intercept"),
  prior(student_t(3,0,20), "sd"),
  prior(student_t(3,0,20), "sigma")
)
#-------------------------------------------------------------------------------
# fit phylogenetic multilevel models with multiple observations per species
# following vignette by Paul Bürkner
# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html

#Prior predictive checks
intra1_prior <- brm(
  pcrit_kpa~temp_test * mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate + 
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

intra1_prior
#This warning can safely be ignored and will be removed in the next release.
# Warning message:
# In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
# '-E' not found.
saveRDS(intra1_prior, file = "intra1_prior.rds")

#Build model
intra1 <- brm(
  pcrit_kpa~temp_test * mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate + 
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra1 <- add_criterion(intra1, "loo", reloo = TRUE) 
intra1$criteria$loo
saveRDS(intra1, file = "intra1.rds")

#Posterior predictive checks                                    
pp_check(intra1, type="scatter_avg", nsamples=100)                      
pp_check(intra1, type = "loo_pit_overlay",nsamples = NULL)                                     
pp_check(intra1, type = "dens_overlay",nsamples = 99)          
pp_check(intra1, type = "hist", nsamples = 11, binwidth = 10)  
#-------------------------------------------------------------------------------
#Prior predictive checks
intra2_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

intra2_prior
saveRDS(intra2_prior, file = "intra2_prior.rds")

#Built model
intra2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra2 <- add_criterion(intra2, "loo", reloo = TRUE) 
intra2$criteria$loo
saveRDS(intra2, file = "intra2.rds")  

#Posterior predictive checks                                    
pp_check(intra2, type="scatter_avg", nsamples=100)                      
pp_check(intra2, type = "loo_pit_overlay",nsamples = NULL)                                     
pp_check(intra2, type = "dens_overlay",nsamples = 99)          
pp_check(intra2, type = "hist", nsamples = 11, binwidth = 10) 
#-------------------------------------------------------------------------------
#Prior predictive checks
intra3_prior <- brm(
  pcrit_kpa ~ temp_test * rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate + 
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

intra3_prior
saveRDS(intra3_prior, file = "intra3_prior.rds")

#Built model
intra3 <- brm(
  pcrit_kpa~temp_test * rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate + 
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra3 <- add_criterion(intra3, "loo", reloo = TRUE)
intra3$criteria$loo
saveRDS(intra3, file = "intra3.rds") 

#Posterior predictive checks                                    
pp_check(intra3, type="scatter_avg", nsamples=100)        
pp_check(intra3, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra3, type = "dens_overlay",nsamples = 99)
pp_check(intra3, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
# Comparing the weights of the models using LOO
model_weights<-model_weights(intra1, intra2, intra3, weights = "loo") %>% 
  round(digits = 2)

model_weights
# intra1 intra2 intra3 
# 0      1      0 
#-------------------------------------------------------------------------------
#adding a 2ndtemp*mass interaction 

#Prior predictive checks
intra1.1_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

intra1.1_prior
saveRDS(intra1.1_prior, file = "intra1.1_prior.rds")

#Built the model
intra1.1 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra1.1 <- add_criterion(intra1.1, "loo", reloo = TRUE)
intra1.1$criteria$loo
saveRDS(intra1.1, file = "intra1.1.rds")

#Posterior predictive checks                                    
pp_check(intra1.1, type="scatter_avg", nsamples=100)        
pp_check(intra1.1, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra1.1, type = "dens_overlay",nsamples = 99)
pp_check(intra1.1, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
#Prior predictive checks
intra3.2_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * rel_mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

intra3.2_prior
saveRDS(intra3.2_prior, file = "intra3.2_prior.rds")

#Built model
intra3.2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * rel_mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra3.2 <- add_criterion(intra3.2, "loo", reloo = TRUE)
intra3.2$criteria$loo
saveRDS(intra3.2, file = "intra3.2.rds")

#Posterior predictive checks                                    
pp_check(intra3.2, type="scatter_avg", nsamples=100)        
pp_check(intra3.2, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra3.2, type = "dens_overlay",nsamples = 99)
pp_check(intra3.2, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
# Comparing the weights of the models using LOO
model_weights<-model_weights(intra2, intra1.1,intra3.2, weights = "loo") %>% 
  round(digits = 2)

model_weights#best model is"intra1.1"
#-------------------------------------------------------------------------------
#Prior predictive checks
intra1.1.1_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

saveRDS(intra1.1.1_prior, file = "intra1.1.1_prior.rds")

#Built model (#adding the last mass metric)
intra1.1.1 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra1.1.1 <- add_criterion(intra1.1.1, "loo", reloo = TRUE)
intra1.1.1$criteria$loo
saveRDS(intra1.1.1, file = "intra1.1.1.rds")

#Posterior predictive checks                                    
pp_check(intra1.1.1, type="scatter_avg", nsamples=100)        
pp_check(intra1.1.1, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra1.1.1, type = "dens_overlay",nsamples = 99)
pp_check(intra1.1.1, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
#Prior predictive checks
intra1.1.2_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    temp_test * rel_mass_log_sq +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

saveRDS(intra1.1.2_prior, file = "intra1.1.2_prior.rds")

#Built model
intra1.1.2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    temp_test * rel_mass_log_sq +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra1.1.2 <- add_criterion(intra1.1.2, "loo", reloo = TRUE)
intra1.1.2$criteria$loo
saveRDS(intra1.1.2, file = "intra1.1.2.rds")

#Posterior predictive checks                                    
pp_check(intra1.1.2, type="scatter_avg", nsamples=100)        
pp_check(intra1.1.2, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra1.1.2, type = "dens_overlay",nsamples = 99)
pp_check(intra1.1.2, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
# Comparing the weights of the models using LOO
model_weights<-model_weights(intra1.1,intra1.1.1,intra1.1.2, weights = "loo") %>% 
  round(digits = 2)
model_weights#best model is"intra1.1.1"

#-------------------------------------------------------------------------------
#Prior predictive checks
intra1.1.1.b_prior <- brm(
  pcrit_kpa~Lat_abs * c_value_log  +
    Lat_abs * mass_log  +
    rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

saveRDS(intra1.1.1.b_prior, file = "intra1.1.1.b_prior.rds")

#Built model (model testing the effect of latitude)
intra1.1.1.b <- brm(
  pcrit_kpa~Lat_abs * c_value_log  +
    Lat_abs * mass_log  +
    rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra1.1.1.b <- add_criterion(intra1.1.1.b, "loo", reloo = TRUE)
intra1.1.1.b$criteria$loo
saveRDS(intra1.1.1.b, file = "intra1.1.1.b.rds")

#Posterior predictive checks                                    
pp_check(intra1.1.1.b, type="scatter_avg", nsamples=100)        
pp_check(intra1.1.1.b, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra1.1.1.b, type = "dens_overlay",nsamples = 99)
pp_check(intra1.1.1.b, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
# Comparing the weights of the models using LOO
model_weights<-model_weights(intra1.1.1,intra1.1.1.b,weights = "loo") %>% 
  round(digits = 2)
model_weights
#-------------------------------------------------------------------------------
#check p values to remove the non significant
p_map(intra1.1.1)
#-------------------------------------------------------------------------------
#Prior predictive checks
intra1.1.1.f_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

saveRDS(intra1.1.1.f_prior, file = "intra1.1.1.f_prior.rds")

# Built model
intra1.1.1.f <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

intra1.1.1.f <- add_criterion(intra1.1.1.f, "loo", reloo = TRUE)
intra1.1.1.f$criteria$loo
saveRDS(intra1.1.1.f, file = "intra1.1.1.f.rds")

#Posterior predictive checks                                    
pp_check(intra1.1.1.f, type="scatter_avg", nsamples=100)        
pp_check(intra1.1.1.f, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(intra1.1.1.f, type = "dens_overlay",nsamples = 99)
pp_check(intra1.1.1.f, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
#INTRASPECIFIC VS INTERSPECIFIC COMPARISONS
#Prior predictive checks
inter1.1.1.f_prior <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = A)), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = "only", 
  cores = 2,
  control = list(adapt_delta = 0.999),
  save_pars = save_pars(all = TRUE))

saveRDS(inter1.1.1.f_prior, file = "inter1.1.1.f_prior.rds")

#Built model
inter1.1.1.f <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = A)), 
    data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

inter1.1.1.f <- add_criterion(inter1.1.1.f, "loo", reloo = TRUE)
inter1.1.1.f$criteria$loo
saveRDS(inter1.1.1.f, file = "inter1.1.1.f.rds")

#Posterior predictive checks                                    
pp_check(inter1.1.1.f, type="scatter_avg", nsamples=100)        
pp_check(inter1.1.1.f, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(inter1.1.1.f, type = "dens_overlay",nsamples = 99)
pp_check(inter1.1.1.f, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
model_weights<-model_weights(intra1.1.1.f,inter1.1.1.f, weights = "loo") %>% 
  round(digits = 2)

model_weights
###############################################################################
## calculating the best value of rho
w<-1:10
plot(NA, xlab = "rho transformation value",ylab = 'Model support relative to best model', 
     xlim = c(0.1, 0.9), ylim = c(0, 1),xpd = NA,las=1)
for (i in 1:10){
p<-0.2-(i*0.002)
B <- ape::vcv.phylo(compute.brlen(tree, power=p))
}
B <- ape::vcv.phylo(compute.brlen(tree, power=0.4))

best_model <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = B)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

best_model <- add_criterion(best_model, "loo", reloo = TRUE)
best_model$criteria$loo
saveRDS(best_model, file = "best_model.rds")

#Posterior predictive checks                                    
pp_check(best_model, type="scatter_avg", ndraws = 100)        
pp_check(best_model, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(best_model, type = "dens_overlay",nsamples = 99)
pp_check(best_model, type = "hist", nsamples = 11, binwidth = 10)
#-------------------------------------------------------------------------------
model_weights<-model_weights(intra1.1.1.f,best_model, weights = "loo") %>% 
  round(digits = 2)
model_weights
#-------------------------------------------------------------------------------
#running time is about 20 minutes
best_model_final <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = A)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 3, cores = 3, 
  iter = 1.5e4, warmup = 1.5e4 / 2,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  save_pars = save_pars(all = TRUE))

#running time is about 4 hours!!!!!

#best_model_final <- add_criterion(best_model_final, "loo", reloo = TRUE)
best_model_final$criteria$loo
saveRDS(best_model_final, file = "best_model_final.rds")

#Posterior predictive checks                                    
pp_check(best_model_final, type="scatter_avg", ndraws = 100)        
pp_check(best_model_final, type = "loo_pit_overlay",nsamples = NULL)                              
pp_check(best_model_final, type = "dens_overlay",nsamples = 99)
pp_check(best_model_final, type = "hist", nsamples = 11, binwidth = 10)

summary(best_model_final,prob = 0.95) %>% print(digits = 2)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: pcrit_kpa ~ temp_test * c_value_log + temp_test * mass_log + sal_test + res_metab_rate + (1 | gr(animal, cov = A)) + (1 | species) 
# Data: d (Number of observations: 590) 
# Draws: 3 chains, each with iter = 15000; warmup = 7500; thin = 1;
# total post-warmup draws = 22500
# 
# Group-Level Effects: 
#   ~animal (Number of levels: 165) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.29      0.54     1.38     3.48 1.00     3969     7622
# 
# ~species (Number of levels: 165) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.50      0.16     1.20     1.81 1.00     6445    11664
# 
# Population-Level Effects: 
#                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                 6.06      1.19     3.71     8.39 1.00     7982    12000
# temp_test                -0.03      0.03    -0.08     0.02 1.00    10255    14668
# c_value_log              -8.18      1.85   -11.82    -4.59 1.00     9035    13626
# mass_log                 -1.37      0.32    -1.99    -0.74 1.00     9889    12713
# sal_test                  0.08      0.01     0.06     0.09 1.00    13967    16166
# res_metab_rate            2.27      0.48     1.34     3.21 1.00    16512    16939
# temp_test:c_value_log     0.46      0.07     0.32     0.61 1.00     9970    14191
# temp_test:mass_log        0.05      0.01     0.03     0.08 1.00    10152    12600
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     1.64      0.06     1.54     1.75 1.00    17155    17049
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


bayes_R2(best_model_final) * 100
#      Estimate Est.Error     Q2.5  Q97.5
# R2   80.699 0.9047613 78.77105 82.313

p_map(best_model_final)
# MAP-based p-value 
# 
# Parameter             | p (MAP)
# -------------------------------
# (Intercept)           |  < .001
# temp_test             |  0.525 
# c_value_log           |  < .001
# mass_log              |  < .001
# sal_test              |  < .001
# res_metab_rate        |  < .001
# temp_test:c_value_log |  < .001
# temp_test:mass_log    |  < .001
#-------------------------------------------------------------------------------

#PHYLOGENETIC SIGNAL
hyp <- "sd_animal__Intercept^2 / (sd_animal__Intercept^2 +sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(best_model_final, hyp, class = NULL))
# Hypothesis Tests for class :
#                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
# 1 (sd_animal__Inter... = 0      0.5      0.12     0.26     0.73          0
#    Post.Prob Star
#    1         0    *
#      ---
#      'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#    '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#    for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#    Posterior probabilities of point hypotheses assume equal prior probabilities.
#-------------------------------------------------------------------------------
# extract fixed effect estimates of the model
fixedEstimates                       <-  brms::fixef(best_model_final, estimate = 'mean')
intercept                            <-  fixedEstimates['Intercept', 'Estimate']
c_value_log_Slope                    <-  fixedEstimates['c_value_log', 'Estimate']
temp_test_Slope                      <-  fixedEstimates['temp_test', 'Estimate']
mass_log_Slope                       <-  fixedEstimates['mass_log', 'Estimate']
sal_test_Slope                       <-  fixedEstimates['sal_test', 'Estimate']
res_metab_rate_Slope                 <-  fixedEstimates['res_metab_rate', 'Estimate']
temp_test_and_c_value_log_Slope      <-  fixedEstimates['temp_test:c_value_log', 'Estimate']
temp_test_and_mass_log_Slope         <-  fixedEstimates['temp_test:mass_log', 'Estimate']

# extract random effect estimates of the model
animal_Random  <-  brms::ranef(best_model_final, estimate = 'mean')$animal[, , 'Intercept'][as.character(d$animal), 'Estimate']
species_Random     <-  brms::ranef(best_model_final, estimate = 'mean')$species[, , 'Intercept'][d$species, 'Estimate']
#-------------------------------------------------------------------------------
# FIGURE 1

summary(lm(animal_Random~species_Random))# random effects are correlated
treetiptrait<-unique(data.frame(label=names(animal_Random+species_Random),
                                traitvalue=(animal_Random+species_Random)+mean(intercept)))
min(treetiptrait$traitvalue)

traitdata<-treetiptrait$traitvalue
names(traitdata)<-treetiptrait$label
fit <- phytools::fastAnc(compute.brlen(tree, power=0.4),traitdata,vars=TRUE,CI=TRUE)


td <- data.frame(node = nodeid(tree, names(traitdata)),trait = traitdata)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
dd <- rbind(td, nd)
dd$node <- as.numeric(dd$node)
fish.tree <- full_join(compute.brlen(tree, power=0.4), dd, by = 'node')

p =
  ggtree(fish.tree, layout='circular', ladderize = FALSE, size=1.5) + 
  geom_tree(aes(color=trait), continuous=T, size=1) + 
  scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
  geom_tiplab(aes(color=trait), hjust = -.1) + 
  xlim(0, 1.2) + 
  theme(legend.position = c(0,.2), plot.margin = margin(2,2,2,2,"cm"))

pdf("Figure 1 Phylogenetic tree.pdf",width = 12,height = 12,useDingbats = FALSE)
print(p)
dev.off()

png("Figure 1 Phylogenetic tree.png",width = 12,height = 12,units = "in",res = 600)
print(p)
dev.off()
#-------------------------------------------------------------------------------
# Calculate predicted response (by Wilco)

d$pcrit_predicted   <-  intercept+(temp_test_Slope * (d$temp_test)) +
  (mass_log_Slope * (d$mass_log)) + 
  (c_value_log_Slope * (d$c_value_log)) + 
  (res_metab_rate_Slope * (d$res_metab_rate)) +
  (sal_test_Slope * mean(d$sal_test)) +
  (temp_test_and_mass_log_Slope * (d$temp_test)* (d$mass_log)) +
  (temp_test_and_c_value_log_Slope * (d$temp_test) * (d$c_value_log)) + 
  species_Random + 
  animal_Random
plot(pcrit_kpa~pcrit_predicted,data=d) # good correlation
plot(I(pcrit_kpa-pcrit_predicted)~pcrit_predicted,data=d) # plotting residuals vs predicted values
abline(a=0,b=1)

min(d$pcrit_kpa-d$pcrit_predicted+intercept)
max(d$pcrit_kpa-d$pcrit_predicted+intercept)
d[which(d$pcrit_kpa-d$pcrit_predicted+intercept<0),]
#-------------------------------------------------------------------------------
# FIGURE 2
{
#pdf("Figure 2 Critical oxygen tension as a function of body mass and genome size.pdf",width = 10,height = 10,useDingbats = FALSE)
#png("Figure 2 Critical oxygen tension as a function of body mass and genome size.png",width = 10,height = 10,units = "in",res = 600)
par(mfrow = c(2, 2),tcl=-0.4, family="Times",omi=c(0,0,0,0), cex = 1, cex.lab = 1.2)
  
#plot interactive effect of temperature and body size
par(mai = c(0.82, 0.82, 0.4, 0.2))
plot(NA, xlab=expression(paste(log[10],"[body mass (g)]")),ylab = 'Critical oxygen tension (kPa)', 
     xlim = c(-4, 4), ylim = c(0, 15),xpd = NA,las=1)

d$Tempcol1<-gsub(" ","",paste("#0000FF",round(70-d$temp_test*2.1,digits=0)))
d$Tempcol1[which(d$temp_test>28)]<-"#0000FF05"
gsub(" ","",paste("#FF0000",round(12+d$temp_test[55]*1.6,digits=0)))
PlotTemps<-c(15,24,28)

meanEffects_TxM  <-  (sal_test_Slope * mean(d$sal_test)) + 
  (c_value_log_Slope * mean(d$c_value_log)) + 
  (res_metab_rate_Slope * mean(d$res_metab_rate)) + 
  (temp_test_and_c_value_log_Slope * mean(d$temp_test) * mean(d$c_value_log))

points((I(pcrit_kpa-pcrit_predicted)+intercept+
          mass_log*mass_log_Slope+
          temp_test*temp_test_Slope+
          temp_test*mass_log*temp_test_and_mass_log_Slope+
          meanEffects_TxM)~mass_log,
       data=d,cex=1.4,pch=16,col=d$Tempcol1) # plotting residuals vs predicted values

iters    <-  brms::posterior_samples(best_model_final, pars = 'b_')
allInts  <-  iters[, 'b_Intercept']
allSlp1  <-  iters[, 'b_mass_log']
allSlp2  <-  iters[, 'b_temp_test']
allSlp3  <-  iters[, 'b_temp_test:mass_log']

xs  <-  seq(-4, 4, length.out = 30)
lw  <-  up  <- med<-  numeric(length = length(xs))

for (i in seq_along(xs)) {
  vals  <-  allInts + allSlp1 * xs[i]+ 
    allSlp2 * PlotTemps[1]+ 
    allSlp3 * xs[i] * PlotTemps[1] + meanEffects_TxM
  lw[i] <-  quantile(vals, probs = 0.025)
  med[i] <-  quantile(vals, probs = 0.5)
  up[i] <-  quantile(vals, probs = 0.975)
}

lines(xs, med, lwd = 2.5, lty = 1, col = c("Blue","Orange","Red")[1])
lines(xs, lw, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[1])
lines(xs, up, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[1])
text(0, 15, "model prediction at 15°C", col="Blue",cex=1.1)
text(-4, 17, '(a)', adj = c(0, 0.5), xpd = NA, font = 3)

par(mai = c(0.82, 0.82, 0.4, 0.2))
plot(NA, xlab=expression(paste(log[10],"[body mass (g)]")),ylab = 'Critical oxygen tension (kPa)', 
     xlim = c(-4, 4), ylim = c(0,15),xpd = NA,las=1)

d$Tempcol2<-gsub(" ","",paste("#FF0000",round(-14+d$temp_test*2.2,digits=0)))
d$Tempcol2[which(d$temp_test<15)]<-"#FF000005"

points((I(pcrit_kpa-pcrit_predicted)+intercept+
          mass_log*mass_log_Slope+
          temp_test*temp_test_Slope+
          temp_test*mass_log*temp_test_and_mass_log_Slope+
          meanEffects_TxM)~mass_log,
       data=d,cex=1.4,pch=16,col=d$Tempcol2) # plotting residuals vs predicted values

PlotTemps<-c(15,24,28)
iters    <-  brms::posterior_samples(best_model_final, pars = 'b_')
allInts  <-  iters[, 'b_Intercept']
allSlp1  <-  iters[, 'b_mass_log']
allSlp2  <-  iters[, 'b_temp_test']
allSlp3  <-  iters[, 'b_temp_test:mass_log']

xs  <-  seq(-4, 4, length.out = 30)
lw  <-  up  <- med<-  numeric(length = length(xs))

for (i in seq_along(xs)) {
  vals  <-  allInts + allSlp1 * xs[i]+ 
    allSlp2 * PlotTemps[3]+ 
    allSlp3 * xs[i] * PlotTemps[3] + meanEffects_TxM
  lw[i] <-  quantile(vals, probs = 0.025)
  med[i] <-  quantile(vals, probs = 0.5)
  up[i] <-  quantile(vals, probs = 0.975)
}

lines(xs, med, lwd = 2.5, lty = 1, col = c("Blue","Orange","Red")[3])
lines(xs, lw, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[3])
lines(xs, up, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[3])
text(0, 15, "model prediction at 28°C", col="Red",cex=1.1)
text(-4, 17, '(b)', adj = c(0, 0.5), xpd = NA, font = 3)

# plot interactive effect of temperature and genome size
par(mai = c(0.82, 0.82, 0.6, 0.2))
plot(NA, xlab=expression(paste(log[10],"[genome size (pg)]")),ylab = 'Critical oxygen tension (kPa)', 
     xlim = c(-0.5, 0.8), ylim = c(0, 15),xpd = NA,las=1)

meanEffects_TxGS  <-  (sal_test_Slope * mean(d$sal_test)) + 
  (res_metab_rate_Slope * mean(d$res_metab_rate)) + 
  (mass_log_Slope * mean(d$mass_log)) + 
  (temp_test_and_mass_log_Slope * mean(d$temp_test) * mean(d$mass_log))

points((I(pcrit_kpa-pcrit_predicted)+intercept+
          c_value_log*c_value_log_Slope+
          temp_test*temp_test_Slope+
          temp_test*c_value_log*temp_test_and_c_value_log_Slope+
          meanEffects_TxGS)~c_value_log,
       data=d,cex=1.4,pch=16,col=d$Tempcol1) # plotting residuals vs predicted values

iters    <-  brms::posterior_samples(best_model_final, pars = 'b_')
allInts  <-  iters[, 'b_Intercept']
allSlp1  <-  iters[, 'b_c_value_log']
allSlp2  <-  iters[, 'b_temp_test']
allSlp3  <-  iters[, 'b_temp_test:c_value_log']

xs  <-  seq(-0.5, 0.8, length.out = 30)
lw  <-  up  <- med<-  numeric(length = length(xs))

for (i in seq_along(xs)) {
  vals  <-  allInts + allSlp1 * xs[i]+ 
    allSlp2 * PlotTemps[1]+ 
    allSlp3 * xs[i] * PlotTemps[1] + meanEffects_TxGS
  lw[i] <-  quantile(vals, probs = 0.025)
  med[i] <-  quantile(vals, probs = 0.5)
  up[i] <-  quantile(vals, probs = 0.975)
}

lines(xs, med, lwd = 2.5, lty = 1, col = c("Blue","Orange","Red")[1])
lines(xs, lw, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[1])
lines(xs, up, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[1])
text(0.15, 15, "model prediction at 15°C", col="Blue",cex=1.1)
text(-0.5, 17, '(c)', adj = c(0, 0.5), xpd = NA, font = 3)

par(mai = c(0.82, 0.82, 0.6, 0.2))
plot(NA, xlab=expression(paste(log[10],"[genome size (pg)]")),ylab = 'Critical oxygen tension (kPa)', 
     xlim = c(-0.5, 0.8), ylim = c(0, 15),xpd = NA,las=1)

points((I(pcrit_kpa-pcrit_predicted)+intercept+
          c_value_log*c_value_log_Slope+
          temp_test*temp_test_Slope+
          temp_test*c_value_log*temp_test_and_c_value_log_Slope+
          meanEffects_TxGS)~c_value_log,
       data=d,cex=1.4,pch=16,col=d$Tempcol2) # plotting residuals vs predicted values

iters    <-  brms::posterior_samples(best_model_final, pars = 'b_')
allInts  <-  iters[, 'b_Intercept']
allSlp1  <-  iters[, 'b_c_value_log']
allSlp2  <-  iters[, 'b_temp_test']
allSlp3  <-  iters[, 'b_temp_test:c_value_log']

xs  <-  seq(-0.5, 0.8, length.out = 30)
lw  <-  up  <- med<-  numeric(length = length(xs))

for (i in seq_along(xs)) {
  vals  <-  allInts + allSlp1 * xs[i]+ 
    allSlp2 * PlotTemps[3]+ 
    allSlp3 * xs[i] * PlotTemps[3] + meanEffects_TxGS
  lw[i] <-  quantile(vals, probs = 0.025)
  med[i] <-  quantile(vals, probs = 0.5)
  up[i] <-  quantile(vals, probs = 0.975)
}

lines(xs, med, lwd = 2.5, lty = 1, col = c("Blue","Orange","Red")[3])
lines(xs, lw, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[3])
lines(xs, up, lwd = 2.5, lty = 2, col = c("Blue","Orange","Red")[3])
text(0.15, 15, "model prediction at 28°C", col="Red",cex=1.1)
text(-0.5, 17, '(d)', adj = c(0, 0.5), xpd = NA, font = 3)
#dev.off()
}
#-------------------------------------------------------------------------------
# FIGURE 3
{
#pdf("Figure 3 Critical oxygen tension as a function of salinity and residual metabolic rate.pdf",width = 10,height = 5,useDingbats = FALSE)
#png("Figure 3 Critical oxygen tension as a function of salinity and residual metabolic rate.png",width = 10,height = 5,units = "in",res = 600)
par(mfrow = c(1, 2),tcl=-0.4, family="Times",omi=c(0,0,0,0), cex = 1, cex.lab = 1.2)

par(mai = c(0.82, 0.82, 0.4, 0.2))
plot(NA, xlab = "Test salinity (PSU)",ylab = 'Critical oxygen tension (kPa)', 
     xlim = c(0, 100), ylim = c(0, 20),xpd = NA,las=1)

meanEffects_Sal  <-  (res_metab_rate_Slope * mean(d$res_metab_rate))+ 
  (c_value_log_Slope * mean(d$c_value_log)) + 
  (temp_test_and_c_value_log_Slope * mean(d$temp_test) * mean(d$c_value_log))+
  (mass_log_Slope * mean(d$mass_log)) + 
  (temp_test_and_mass_log_Slope * mean(d$temp_test) * mean(d$mass_log))

# plotting residuals vs predicted values
points((I(pcrit_kpa-pcrit_predicted)+intercept+sal_test*sal_test_Slope)+meanEffects_Sal~sal_test,data=d,cex=1.4,pch=16,col="#4DAC2670") 

iters    <-  brms::posterior_samples(best_model_final, pars = 'b_')
allInts  <-  iters[, 'b_Intercept']
allSlps  <-  iters[, 'b_sal_test']

xs  <-  seq(min(d$sal_test), max(d$sal_test), length.out = 30)
lw  <-  up  <- med<-  numeric(length = length(xs))
for (i in seq_along(xs)) {
  vals  <-  allInts + allSlps * xs[i] + meanEffects_Sal
  lw[i] <-  quantile(vals, probs = 0.025)
  med[i] <-  quantile(vals, probs = 0.5)
  up[i] <-  quantile(vals, probs = 0.975)
}

lines(xs, med, lwd = 2.5, lty = 1, col = '#4DAC26')
lines(xs, lw, lwd = 2.5, lty = 2, col = '#4DAC26')
lines(xs, up, lwd = 2.5, lty = 2, col = '#4DAC26')
text(65, 1.5, "590 observations", adj = c(0, 0.5),cex=0.9)
text(65, 0, paste0(length(unique(d$species)), ' fish species'), adj = c(0, 0.5),cex=0.9)
text(0, 22, '(a)', adj = c(0, 0.5), xpd = NA, font = 3)

## plot residual metabolic rate vs pcrit using partial residuals by Wilco ##
par(mai = c(0.8, 0.8, 0.4, 0.2))
plot(NA, xlab = "Residuals of metabolic rate",ylab = 'Critical oxygen tension (kPa)', 
     xlim = c(-1.2, 1.2), ylim = c(0, 20),xpd = NA,las=1)

iters    <-  brms::posterior_samples(best_model_final, pars = 'b_')
allInts  <-  iters[, 'b_Intercept']
allSlps  <-  iters[, 'b_res_metab_rate']

meanEffects_resMR  <-  (sal_test_Slope * mean(d$sal_test)) + 
  (c_value_log_Slope * mean(d$c_value_log)) + 
  (temp_test_and_c_value_log_Slope * mean(d$temp_test) * mean(d$c_value_log))+
  (mass_log_Slope * mean(d$mass_log)) + 
  (temp_test_and_mass_log_Slope * mean(d$temp_test) * mean(d$mass_log))

# plotting residuals vs predicted values
points((I(pcrit_kpa-pcrit_predicted)+intercept+res_metab_rate*res_metab_rate_Slope)+meanEffects_resMR~res_metab_rate,data=d,cex=1.4,pch=16,col="#4DAC2670") 

xs  <-  seq(min(d$res_metab_rate), max(d$res_metab_rate), length.out = 30)
lw  <-  up  <- med<-  numeric(length = length(xs))
for (i in seq_along(xs)) {
  vals  <-  allInts + allSlps * xs[i] + meanEffects_resMR
  lw[i] <-  quantile(vals, probs = 0.025)
  med[i] <-  quantile(vals, probs = 0.5)
  up[i] <-  quantile(vals, probs = 0.975)
}

lines(xs, med, lwd = 2.5, lty = 1, col = '#4DAC26')
lines(xs, lw, lwd = 2.5, lty = 2, col = '#4DAC26')
lines(xs, up, lwd = 2.5, lty = 2, col = '#4DAC26')
text(-1.2, 22, '(b)', adj = c(0, 0.5), xpd = NA, font = 3)
#dev.off()
}
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs/phylogenetic_hierarchical_models_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################