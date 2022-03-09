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
                               "Temp",
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

# Models are saved so that they don't need to be rerun for every session.
best_model_final <- readRDS("best_model_final.rds")
# extract fixed effect estimates of the model
fixedEstimates                       <-  brms::fixef(best_model_final, estimate = 'mean')
intercept                            <-  fixedEstimates['Intercept', 'Estimate']
c_value_log_Slope                    <-  fixedEstimates['c_value_log', 'Estimate']
temp_test_Slope                      <-  fixedEstimates['temp_test', 'Estimate']
mass_log_Slope                       <-  fixedEstimates['mass_log', 'Estimate']
sal_test_Slope                       <-  fixedEstimates['sal_test', 'Estimate']
res_metab_rate_Slope                 <-  fixedEstimates['res_metab_rate', 'Estimate']
temp_test_and_c_value_log_Slope      <-  fixedEstimates['temp_test:c_value_log', 'Estimate']
#temp_test_and_c_value_log_Slope      <-  fixedEstimates['c_value_log:temp_test', 'Estimate']
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
#  scale_color_gradientn(colours=inlmisc::GetColors(length(seq(0,20,5)),scheme='inferno')) +
  geom_tiplab(aes(colour=trait), hjust = -.1) + 
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
#allSlp3  <-  iters[, 'b_c_value_log:temp_test']
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
#allSlp3  <-  iters[, 'b_c_value_log:temp_test']
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