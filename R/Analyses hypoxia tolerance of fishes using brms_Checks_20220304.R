################################################################################
# Script to fit phylogenetic multilevel model using brms
# Author: Félix P Leiva (felixpleiva@gmail.com)
# Date: 20220304
# Modifications: Wilco CEP Verberk (20220304); Félix P Leiva (20210920 )
# ------------------------------------------------------------------------------
# Cite as:

#Verberk WCEP, Sandkler JF, van de Pol I, Urbina M, Wilson R,
#McKenzie DJ, Leiva FP (2022). Data and code of manuscript: Hypoxia
#tolerance in fish varies with body size and cell size in a
#temperature-dependent manner (). Zenodo ().
# ------------------------------------------------------------------------------
rm(list=ls()) #clear the work environment
today<-format(Sys.Date(),"%Y%m%d") #setting the date
# ------------------------------------------------------------------------------
setwd("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs") #Felix's lab PC
#setwd("C:/Users/wcepv/surfdrive/papers in progress/pcrit fish") #Wilco's laptop
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
library(corrplot)
# ------------------------------------------------------------------------------
set.seed(6955)# we need this to replicate the results

# General STAN specifications
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load data
data <- read.csv("Pcrit data for 195 fish species ( 2022_02_14).csv")
#write.csv(data,"Pcrit data for 195 fish species ( 2022_02_14).csv")

# load phylogenetic tree
tree <- read.tree("Phylogenetic tree for 195 species of fish with Pcrit data.tre")
plotTree(tree,fsize=0.5,ftype="i",type="fan") # branch lengths are Grafen
#plotTree(compute.brlen(tree, method = "Grafen"),fsize=0.5,ftype="i",type="fan")

# The tree object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

#A <- ape::vcv.phylo(tree)
B <- ape::vcv.phylo(compute.brlen(tree, power=0.4)) # preliminary analyses showed that this tree transformation resulted in better model performance.

#-------------------------------------------------------------------------------

# apply complete.clases (exclude NAs) for the columns of interest
d<-data[complete.cases(data[,c("Estimated_C.Value",
                               "Temp",
                               "logMass",
                               #"delta_bm",
                               #"logpercentcalc",
                               "Pcrit..kPa.",
                               "Salt",
                               "animal",
                               #"Taccl",
                               "residMR")]),]

#write.csv(d,"Pcrit data (selection) for 171 fish species ( 20211127 ).csv")

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

length(unique(d$animal))#count number of species to be included in the models

# create new columns
d$species<-d$animal

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
#d$rel_mass_log_sq<-(d$logpercentcalc)^2
#d$mass_log_sq<-(d$mass_log)^2
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
# Comparing the estimate for salinity effect
m_salt_excl <- brm(
  pcrit_kpa~ c_value_log * temp_test  +
    mass_log * temp_test  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = B)) + 
    (1|species), 
  data = d[which(d$Salt_type_data=="quantatative"),], 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 3, cores = 2, 
  iter = 15000, warmup = 7500,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  save_pars = save_pars(all = TRUE))
saveRDS(m_salt_excl,"m_salt_excl_15_02_22.rds")
best_model_final <- readRDS("best_model_final.rds")

summary(best_model_final)
summary(m_salt_excl)

fixedEstimatesBM                       <-  brms::fixef(best_model_final, estimate = 'mean')
fixedEstimatesBM['sal_test',]#0.075454757 (range 0.056212709 0.094904226)
fixedEstimatesxS                       <-  brms::fixef(m_salt_excl, estimate = 'mean')
fixedEstimatesxS['sal_test',]#0.078265708 (range 0.059004873 0.097599532)

# comparing the explanatory power of genome size vs Maximum mass
e<-d[complete.cases(d[,c("MaxMass_log")]),]

m_MaxMass <- brm(
  pcrit_kpa~ MaxMass_log * temp_test  +
    mass_log * temp_test  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = B)) + 
    (1|species), 
  data = e, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 3, cores = 2, 
  iter = 15000, warmup = 7500,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  save_pars = save_pars(all = TRUE))
m_MaxMass <- readRDS("m_MaxMass_15_02_22.rds")
saveRDS(m_MaxMass,"m_MaxMass_15_02_22.rds")

m_best_model <- brm(
  pcrit_kpa~ c_value_log * temp_test  +
    mass_log * temp_test  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = B)) + 
    (1|species), 
  data = e, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 3, cores = 2, 
  iter = 15000, warmup = 7500,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  save_pars = save_pars(all = TRUE))
saveRDS(m_best_model,"m_best_model_15_02_22.rds")

summary(m_MaxMass)
summary(m_best_model)

bayes_R2(m_MaxMass)
bayes_R2(m_best_model)

p_map(m_MaxMass)
p_map(m_best_model)

loo(m_MaxMass)
loo(m_best_model)
model_weights<-model_weights(m_MaxMass, m_best_model, weights = "loo") %>% 
  round(digits = 4)

# Analysing the relationship between genome size and body size within a phylogenetic context
e<-d[complete.cases(d[,c("MaxMass_log")]),]

gs_mass <- brm(
  mass_log~ c_value_log  +
    (1|gr(animal, cov = B)),
    #(1|species), 
  data = e, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 3, cores = 2, 
  iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  save_pars = save_pars(all = TRUE))

gs_mass <- readRDS("gs_mass_23_02_22.rds")
saveRDS(gs_mass,"gs_mass_23_02_22.rds")
saveRDS(gs_mass,"gs_mass_phylo_24_02_22.rds")
gs_mass <- readRDS("gs_mass_phylo_24_02_22.rds")

gs_mass_max <- brm(
  MaxMass_log~c_value_log +
    (1|gr(animal, cov = B)), 
    #(1|species), 
  data = e,
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))

gs_mass_max <- readRDS("gs_mass_max_23_02_22.rds")
saveRDS(gs_mass_max,"gs_mass_23_02_22.rds")

summary(gs_mass)
summary(gs_mass_max)

p_map(gs_mass)
p_map(gs_mass_max)

par (mfrow=c(1,2))
plot(logMass~c_value_log,data=d,pch=19,col="#00000040",
     xlab="Genome size (log transformed C-value)",
     ylab="log10 individual mass"
     )
plot(MaxMass_log~c_value_log,data=d,pch=19,col="#00000040",
     xlab="Genome size (log transformed C-value)",
     ylab="log10 species maximum mass"
)

cor.test(d$logMass,d$c_value_log,method="spearman")
cor.test(d$c_value_log,d$MaxMass_log,method="spearman")

# Checking the accuracy of predicting the mass using the length mass conversion relation from another closely related species.
BMest<-read.csv("% body mass.csv")
names(BMest)[10]<-"SpeciesA"
names(BMest)[17]<-"SpeciesB"
names(BMest)[7]<-"maxMA"
names(BMest)[20]<-"maxMB"
BMest<-BMest[!is.na(BMest$maxMB),]
plot(log10(maxMA)~log10(maxMB),data=BMest,xlab="Estimated maximum mass Species A",ylab="Estimated maximum mass Species B",pch=19,col="#00000070",cex.lab=1.5,cex=1.5,cex.axis=1.5)
abline(a=0,b=1,lty=2,lwd=3)
m_est<-lm(log10(maxMA)~log10(maxMB)-1,data=BMest)
summary(m_est)     
text(1,4.5,"y=1.02 x; R2 = 0.987")
#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs/phylogenetic_hierarchical_models_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################