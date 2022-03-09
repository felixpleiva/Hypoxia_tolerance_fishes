################################################################################
# Script to fit phylogenetic multilevel model using brms
# Author: Félix P Leiva (felixpleiva@gmail.com)
# Date: 20220304
# Modifications: Félix P Leiva (20210920): Wilco CEP Verberk (20220304); 
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
#save the data
#write.csv(data,"Pcrit data for 195 fish species ( 2022_02_14).csv")

#MR data residuals
data$logMR<-log10(data$MR.µM.h)
data$logMass<-log10(data$Mean.Mass..g.)
res_MR<-lm(logMR~logMass+I(1/(Temp+273))+Salt,data=data)
summary(res_MR)
data$residMR<-0
data$residMR[-c(which(is.na(data$logMass)),which(is.na(data$logMR)),which(is.na(data$Salt)))]<-res_MR$res
hist(data$residMR)
visreg(res_MR,"logMass",by="Temp",overlay=TRUE)
visreg(res_MR,"logMass",by="Salt",overlay=TRUE)

#Percentage body mass
data$Calculatedpercent<-100*data$Mean.Mass..g./data$Calculated.mass..g.
data$Calculatedpercent[which(data$Calculatedpercent>100)]<-100
data$logpercentcalc<-log10(data$Calculatedpercent)

# load phylogenetic tree
tree <- read.tree("Phylogenetic tree for 195 species of fish with Pcrit data.tre")
plotTree(tree,fsize=0.5,ftype="i",type="fan") # branch lengths are Grafen
#plotTree(compute.brlen(tree, method = "Grafen"),fsize=0.5,ftype="i",type="fan")

# The tree object contains information on the relationship between species. 
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).

#A <- ape::vcv.phylo(tree)
B <- ape::vcv.phylo(compute.brlen(tree, power=0.4)) # preliminary analyses showed that this tree transformation resulted in better model performance.

#-------------------------------------------------------------------------------
# correlation plots
data$c_value_log<-log10(data$Estimated_C.Value)
data$MaxMass_log<-log10(data$Calculated.mass..g.)

cor.matrix <- cor(data[,c('Lat_abs',
                          'Temp',
                          'Acclimation.Temperature..OC.',
                          'Taccl',
                          'Salt',
                          'c_value_log',
                          'MaxMass_log',
                          'logpercentcalc',
                          'logMass',
                          'residMR')],
                  use = "pairwise.complete.obs", method = "spearman") # rank correlations
par(mfrow=c(1,1))
rownames(cor.matrix)<-colnames(cor.matrix)<-c("Latitude","Test temperature","Acclimation temperature","AccT - TestT","Salinity","Genome size","Maximum mass","Percentage body mass","Individual mass","Residual Metabolic rate")
corrplot::corrplot(cor.matrix, method = 'number', type = 'lower', number.cex = 0.8)

# apply complete.clases (exclude NAs) for the columns of interest
d<-data[complete.cases(data[,c("Estimated_C.Value",
                               #"Temp",
                               "logMass",
                               #"delta_bm",
                               #"logpercentcalc",
                               "Pcrit..kPa.",
                               "Salt",
                               "animal",
                               #"Taccl",
                               "residMR")]),]

nrow(data[complete.cases(data[,"Estimated_C.Value"]),])   #668 records
nrow(data[complete.cases(data[,"Temp"]),])                #702 records
nrow(data[complete.cases(data[,"logMass"]),])             #685 records
nrow(data[complete.cases(data[,"logpercentcalc"]),])      #679 records
nrow(data[complete.cases(data[,"Pcrit..kPa."]),])         #667 records
nrow(data[complete.cases(data[,"Salt"]),])                #673 records (of which 157 are quantitative)
nrow(data[complete.cases(data[,"Taccl"]),])               #696 records
nrow(data[complete.cases(data[,"residMR"]),])             #574 records (other records are set to 0)
nrow(data[complete.cases(data[,"Acclimation.Temperature..OC."]),]) #696 records
nrow(data[complete.cases(data[,"RESPIROMETRY_TYPE"]),])   #702 records (154 records classified as unkown)

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
d$rel_mass_log_sq<-(d$logpercentcalc)^2
d<-rename(d,rel_mass_log=logpercentcalc)
d<-rename(d,temp_test=Temp)
d<-rename(d,pcrit_kpa=Pcrit..kPa.)
d<-rename(d,sal_test=Salt)
d<-rename(d,temp_acclim=Taccl)
d<-rename(d,res_metab_rate=residMR)
d<-rename(d,mass_log=logMass)
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

#models with random effect per species and phylogeny

A1 <- brm(
  pcrit_kpa~temp_test * mass_log +
    resp_type +
    sal_test +
    temp_acclim +
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
saveRDS(A1,"A1_15_02_22.rds")
A1<-readRDS("A1_15_02_22.rds")

A2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    resp_type+
    sal_test +
    temp_acclim +
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
saveRDS(A2,"A2_15_02_22.rds")
A2<-readRDS("A2_15_02_22.rds")

A3 <- brm(
  pcrit_kpa~temp_test * rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
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
saveRDS(A3,"A3_15_02_22.rds")
A3<-readRDS("A3_15_02_22.rds")

#adding a 2ndtemp*mass interaction 
A2.1 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
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
saveRDS(A2.1,"A2_1_15_02_22.rds")
A2.1<-readRDS("A2_1_15_02_22.rds")

A2.2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * rel_mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
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
saveRDS(A2.2,"A2_2_15_02_22.rds")
A2.2<-readRDS("A2_2_15_02_22.rds")

#adding the last mass metric 
A2.1.1 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
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
saveRDS(A2.1.1,"A2_1_1_15_02_22.rds")
A2.1.1<-readRDS("A2_1_1_15_02_22.rds")

A2.1.2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    resp_type +
    sal_test +
    temp_acclim +
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
saveRDS(A2.1.2,"A2_1_2_15_02_22.rds")
A2.1.2<-readRDS("A2_1_2_15_02_22.rds")

A2.1.3 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    rel_mass_log +
    rel_mass_log_sq +
    resp_type+
    sal_test +
    temp_acclim +
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
saveRDS(A2.1.3,"A2_1_3_23_02_22.rds")
A2.1.3<-readRDS("A2_1_3_23_02_22.rds")

#check p values to revome the non significant
p_map(A2.1)

A4 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
  # resp_type+
    sal_test +
  # temp_acclim +
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
saveRDS(A4,"A4_15_02_22.rds")
A4<-readRDS("A4_15_02_22.rds")

msA1<-loo(A1)$estimates
msA2<-loo(A2)$estimates
msA3<-loo(A3)$estimates
msA2.1<-loo(A2.1)$estimates
msA2.2<-loo(A2.2)$estimates
msA2.1.1<-loo(A2.1.1)$estimates
msA2.1.2<-loo(A2.1.2)$estimates
msA2.1.3<-loo(A2.1.3)$estimates

modellist<-data.frame(model=NA,looic=NA,p_loo=NA,elpd_loo=NA,elpd_loo_SE=NA)
#modellist<-read.csv("modellist.csv")
#modellist<-modellist[,-1]
modellist[1,]<-c("A1",msA1[3,1],msA1[2,1],msA1[1,1],msA1[1,2])
modellist[2,]<-c("A2",msA2[3,1],msA2[2,1],msA2[1,1],msA2[1,2])
modellist[3,]<-c("A3",msA3[3,1],msA3[2,1],msA3[1,1],msA3[1,2])
modellist[4,]<-c("A2.1",msA2.1[3,1],msA2.1[2,1],msA2.1[1,1],msA2.1[1,2])
modellist[5,]<-c("A2.2",msA2.2[3,1],msA2.2[2,1],msA2.2[1,1],msA2.2[1,2])
modellist[6,]<-c("A2.1.1",msA2.1.1[3,1],msA2.1.1[2,1],msA2.1.1[1,1],msA2.1.1[1,2])
modellist[7,]<-c("A2.1.2",msA2.1.2[3,1],msA2.1.2[2,1],msA2.1.2[1,1],msA2.1.2[1,2])
modellist[8,]<-c("A2.1.3",msA2.1.3[3,1],msA2.1.3[2,1],msA2.1.3[1,1],msA2.1.3[1,2])
modellist$elpd_loo_diff<-as.numeric(modellist$elpd_loo)-as.numeric(min(modellist$elpd_loo))
write.csv(modellist,"modellist.csv")

model_weights(A1,A2,A3,A2.1,A2.2,A2.1.1,A2.1.2,A2.1.3, weights = "loo") %>% 
  round(digits = 2)

#models without random effect per species, only phylogeny

E1 <- brm(
  pcrit_kpa~temp_test * mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate + 
    (1|gr(animal, cov = B)), 
    #(1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E1,"A1_18_02_22.rds")
E1<-readRDS("A1_18_02_22.rds")

E2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)),
   # (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E2,"E2_18_02_22.rds")
E2<-readRDS("E2_18_02_22.rds")

E3 <- brm(
  pcrit_kpa~temp_test * rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate + 
    (1|gr(animal, cov = B)), 
  #  (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E3,"E3_18_02_22.rds")
E3<-readRDS("E3_18_02_22.rds")

#adding a 2ndtemp*mass interaction 
E2.1 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)),
    #(1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E2.1,"E2_1_18_02_22.rds")
E2.1<-readRDS("E2_1_18_02_22.rds")

E2.2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * rel_mass_log  +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)),
   # (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E2.2,"E2_2_18_02_22.rds")
E2.2<-readRDS("E2_2_18_02_22.rds")

#adding the last mass metric 
E2.1.1 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    rel_mass_log +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)),
  #  (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E2.1.1,"E2_1_1_18_02_22.rds")
E2.1.1<-readRDS("E2_1_1_18_02_22.rds")

E2.1.2 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    resp_type +
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)), 
 #   (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E2.1.2,"E2_1_2_18_02_22.rds")
E2.1.2<-readRDS("E2_1_2_18_02_22.rds")

E2.1.3 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    rel_mass_log +
    rel_mass_log_sq +
    resp_type+
    sal_test +
    temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)),
 #   (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E2.1.3,"E2_1_3_23_02_22.rds")
E2.1.3<-readRDS("E2_1_3_23_02_22.rds")

#check p values to remove the non significant
p_map(E2.1)

E4 <- brm(
  pcrit_kpa~temp_test * c_value_log  +
    temp_test * mass_log  +
    # resp_type+
    sal_test +
    # temp_acclim +
    res_metab_rate +  
    (1|gr(animal, cov = B)),
#    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 2, cores = 2, 
  iter = 4000, warmup = 1000,
  save_pars = save_pars(all = TRUE))
saveRDS(E4,"E4_18_02_22.rds")
E4<-readRDS("E4_18_02_22.rds")

msE1<-loo(E1)$estimates
msE2<-loo(E2)$estimates
msE3<-loo(E3)$estimates
msE2.1<-loo(E2.1)$estimates
msE2.2<-loo(E2.2)$estimates
msE2.1.1<-loo(E2.1.1)$estimates
msE2.1.2<-loo(E2.1.2)$estimates
msE2.1.3<-loo(E2.1.3)$estimates

modellist<-data.frame(model=NA,looic=NA,p_loo=NA,elpd_loo=NA,elpd_loo_SE=NA)
#modellist<-read.csv("modellist.csv")
#modellist<-modellist[,-1]
modellist[1,]<-c("E1",msE1[3,1],msE1[2,1],msE1[1,1],msE1[1,2])
modellist[2,]<-c("E2",msE2[3,1],msE2[2,1],msE2[1,1],msE2[1,2])
modellist[3,]<-c("E3",msE3[3,1],msE3[2,1],msE3[1,1],msE3[1,2])
modellist[4,]<-c("E2.1",msE2.1[3,1],msE2.1[2,1],msE2.1[1,1],msE2.1[1,2])
modellist[5,]<-c("E2.2",msE2.2[3,1],msE2.2[2,1],msE2.2[1,1],msE2.2[1,2])
modellist[6,]<-c("E2.1.1",msE2.1.1[3,1],msE2.1.1[2,1],msE2.1.1[1,1],msE2.1.1[1,2])
modellist[7,]<-c("E2.1.2",msE2.1.2[3,1],msE2.1.2[2,1],msE2.1.2[1,1],msE2.1.2[1,2])
modellist[8,]<-c("E2.1.3",msE2.1.3[3,1],msE2.1.3[2,1],msE2.1.3[1,1],msE2.1.3[1,2])
modellist$elpd_loo_diff<-as.numeric(modellist$elpd_loo)-as.numeric(min(modellist$elpd_loo))
write.csv(modellist,"modellist.csv")

model_weights(E1,E2,E3,E2.1,E2.2,E2.1.1,E2.1.2,E2.1.3, weights = "loo")



msA4<-loo(A4)$estimates
msE4<-loo(E4)$estimates

modellist<-data.frame(model=NA,looic=NA,p_loo=NA,elpd_loo=NA,elpd_loo_SE=NA)
modellist[1,]<-c("A4",msA4[3,1],msA4[2,1],msA4[1,1],msA4[1,2])
modellist[2,]<-c("E4",msE4[3,1],msE4[2,1],msE4[1,1],msE4[1,2])
modellist$elpd_loo_diff<-as.numeric(modellist$elpd_loo)-as.numeric(min(modellist$elpd_loo))
model_weights(A4,E4, weights = "loo")

write.csv(modellist,"modellist.csv")

best_model_final <- brm(
  pcrit_kpa~ c_value_log * temp_test  +
    mass_log * temp_test  +
    sal_test +
    res_metab_rate +  
    (1|gr(animal, cov = B)) + 
    (1|species), 
  data = d, 
  family = gaussian(), 
  data2 = list(B = B),
  prior = priors1,
  sample_prior = TRUE, 
  chains = 3, cores = 2, 
  iter = 15000, warmup = 7500,
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

plot(conditional_effects(best_model_final), points = TRUE) 
plot(best_model_final, N = 2, ask = FALSE)

summary(best_model_final,prob = 0.95) %>% print(digits = 2)
# Family: gaussian 
# Links: mu = identity; sigma = identity 
# Formula: pcrit_kpa ~ c_value_log * temp_test + mass_log * temp_test + sal_test + res_metab_rate + (1 | gr(animal, cov = B)) + (1 | species) 
# Data: d (Number of observations: 600) 
# Draws: 3 chains, each with iter = 15000; warmup = 7500; thin = 1;
# total post-warmup draws = 22500
# 
# Group-Level Effects: 
#   ~animal (Number of levels: 171) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     2.22      0.45     1.40     3.11 1.00     1471     3121
# 
# ~species (Number of levels: 171) 
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     1.02      0.34     0.18     1.57 1.00     1157     1396
# 
# Population-Level Effects: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                 5.56      1.02     3.56     7.59 1.00     9599    13167
# c_value_log              -7.47      1.85   -11.11    -3.87 1.00    11618    14179
# temp_test                -0.03      0.03    -0.08     0.02 1.00    11894    15327
# mass_log                 -1.37      0.32    -1.99    -0.74 1.00    11312    14283
# sal_test                  0.08      0.01     0.06     0.09 1.00    10384    16807
# res_metab_rate            2.20      0.47     1.27     3.14 1.00    17375    17728
# c_value_log:temp_test     0.44      0.07     0.30     0.58 1.00    10840    15631
# temp_test:mass_log        0.06      0.01     0.03     0.08 1.00    11437    14360
# 
# Family Specific Parameters: 
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     1.64      0.06     1.53     1.75 1.00    17594    15942
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

bayes_R2(best_model_final) * 100
#    Estimate Est.Error     Q2.5    Q97.5
# R2 80.51155 0.9112743 78.56735 82.15607

p_map(best_model_final)
# MAP-based p-value 
# 
# Parameter             | p (MAP)
# -------------------------------
# (Intercept)           |  < .001
# c_value_log           |  < .001
# temp_test             |  0.536 
# mass_log              |  < .001
# sal_test              |  < .001
# res_metab_rate        |  < .001
# temp_test:c_value_log |  < .001
# temp_test:mass_log    |  < .001
#-------------------------------------------------------------------------------

# Models are saved so that they don't need to be rerun for every session.
best_model_final <- readRDS("best_model_final.rds")

#PHYLOGENETIC SIGNAL
hyp <- "sd_animal__Intercept^2 / (sd_animal__Intercept^2 +sd_species__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(best_model_final, hyp, class = NULL))
plot(hyp)
# Hypothesis Tests for class :
#                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio 
# 1 (sd_animal__Inter... = 0     0.56      0.13      0.3     0.78          0  
#   Post.Prob  Star
#          0      1    *
#      ---
#      'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#    '*': For one-sided hypotheses, the posterior probability exceeds 95%;
#    for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#    Posterior probabilities of point hypotheses assume equal prior probabilities.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs/phylogenetic_hierarchical_models_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################