################################################################################
# Script to create a phylogenetic tree from Open Tree of Life
# Original script: Alfredo Sanchez-Tojar; https://github.com/ASanchez-Tojar (Gracias Alfredo!!)
# Modifications: Félix P Leiva; Date: 20210617
################################################################################
rm(list=ls()) #clear your work environment
today<-format(Sys.Date(),"%Y%m%d") #setting the date
################################################################################
#From Windows
setwd("C:/Users/Invunche/Dropbox/collaboration/Jeroen/data") #Felix's personal laptop
getwd()#to check
################################################################################
#Libraries
library(ape) 
library(phytools)
library(geiger)
library(caper)
library(stringr)
library(fishtree)
library(rotl)
################################################################################
#load data
fish<-read.csv("data internship Jeroen.csv",sep = ";",header = TRUE)
names(fish)
head(fish)

#create a new variable (Species2) with the binomial scientific name (e.g., Danio rerio)
fish$Species2<-as.character(paste(fish$Genus,fish$Species,sep = " "))

#Number of species
length(unique(fish$Species2))#217 species

# Removes white spaces. Be care with the "double spaces" in excel!!!!
fish$Species2<-str_squish(fish$Species2)

#Number of species again
length(unique(fish$Species2))#214 species

#correct some names
fish$Species2[fish$Species2 == "Ctenopharyngodon idellus"] <- "Ctenopharyngodon idella"
fish$Species2[fish$Species2 == "Gobiodon erythrophyllus"] <- "Gobiodon erythrospilus"
fish$Species2[fish$Species2 == "Lumpeninae lampretaeformis"] <- "Lumpenus lampretaeformis"
fish$Species2[fish$Species2 == "Fundulus heteroclitus macrolepidotus"] <- "Fundulus heteroclitus heteroclitus"
fish$Species2[fish$Species2 == "Pseudocrenilabrus multicolor victoriae"] <- "Pseudocrenilabrus multicolor"
fish$Species2[fish$Species2 == "Cyprinus Cyprinus carpio var. Jian"] <- "Cyprinus carpio carpio"
fish$Species2[fish$Species2 == "Carassius auratus"] <- "Carassius auratus auratus"
fish$Species2[fish$Species2 == "Carassius Carassius auratus"] <- "Carassius auratus auratus"
fish$Species2[fish$Species2 == "Carassius Carassius carassius"] <- "Carassius carassius"
fish$Species2[fish$Species2 == "Archamia fucata"] <- "Taeniamia fucata"
fish$Species2[fish$Species2 == "Gadus ogac"] <- "Gadus macrocephalus"
fish$Species2[fish$Species2 == "Apogon leptacanthus"] <- "Zoramia leptacantha"
fish$Species2[fish$Species2 == "Fundulus heteroclitus"] <- "Fundulus heteroclitus heteroclitus"
fish$Species2[fish$Species2 == "Apogon doederleini"] <- "Ostorhinchus doederleini"
fish$Species2[fish$Species2 == "Helcogramma medium"] <- "Bellapiscis medius"
fish$Species2[fish$Species2 == "Aristichthys nobilis"] <- "Hypophthalmichthys nobilis"
fish$Species2[fish$Species2 == "Apogon fragilis"] <- "Zoramia fragilis"
fish$Species2[fish$Species2 == "Cyprinus carpio"] <- "Cyprinus carpio carpio"
fish$Species2[fish$Species2 == "Chaetopsylla globiceps"] <- "Clinocottus globiceps" #wrong label in the original paper (Mandic et al 2009)

#Number of species again
length(unique(fish$Species2))#200 species

#exclude some entries
fish<-subset(fish,fish$Species2 != "Hypseleotris Hypseleotris sp.")# genus-level
fish<-subset(fish,fish$Species2 != "Rostaraja eglanteria")# no phylo information
fish<-subset(fish,fish$Species2!="Lepomis hybrid between bluegill & pumpkinseed")# exclude because is an hybrid

################################################################################
#Number of species again
length(unique(fish$Species2))#197 species

# generating list of species
species <- sort(unique(fish$Species2))

# obtaining data frame listing the Open Tree identifiers potentially matching our
# list of species.
taxa <- tnrs_match_names(names = species)

#exclude some entries because "ott3635470" which are the species Chrysiptera
#flavipinnis and Acanthoclinus fuscus has the flag "incertae_sedis_inherited",
#i.e., this taxon is in our taxonomy but does not appear in the latest synthetic
#tree.

fish<-subset(fish,fish$Species2 != "Chrysiptera flavipinnis")# incertae_sedis_inherited
fish<-subset(fish,fish$Species2 != "Acanthoclinus fuscus")# incertae_sedis_inherited

#Number of species again
length(unique(fish$Species2))#195 species

# generating list of species again
species <- sort(unique(fish$Species2))

# obtaining data frame listing the Open Tree identifiers potentially matching our
# list of species.
taxa <- tnrs_match_names(names = species)

# according to the `approximate_match`column, there might be 
# one typos in the species list (well actually 1 with the final subset)
nrow(taxa[taxa$approximate_match==TRUE,])
# 0K
taxa[taxa$approximate_match==TRUE,]
# 0k
 
# In case of typos,follow the example bellow to fix it
#fish$Species2[fish$Species2 == "Cyprinus carpis"] <- "Cyprinus carpio"


# rerun
taxa.c <- tnrs_match_names(names = species)

# exploring which species return more than one match
taxa.c[taxa.c$number_matches != 1,]

#Reasons to make sure we retrieve the correct data.
ott_id_tocheck <- taxa.c[taxa.c$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa.c, ott_id = ott_id_tocheck[i]))
}

# There are some species with synonyms
# fixing those synonyms in the list for the phylo

species[species=="carassius auratus auratus"] <- "Carassius auratus auratus"
species[species=="carassius auratus grandoculis"] <- "Carassius auratus grandoculis"
species[species=="cyprinus carpio carpio"] <- "Cyprinus carpio carpio"
species[species=="fundulus heteroclitus heteroclitus"] <- "Fundulus heteroclitus heteroclitus"
species[species=="gobiodon erythrospilus"] <- "Gobiodon erythrospilus"
species[species=="pagrus auratus"] <- "Pagrus auratus"
species[species=="petrocephalus degeni"] <- "Petrocephalus degeni"
species[species=="trematomus centronotus"] <- "Trematomus pennellii"
##############################################################
# Retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
fish.tree <- tol_induced_subtree(ott_ids = taxa.c[["ott_id"]], label_format = "name")

# Notice that the species names shown in the tree are not exactly 
# the same as the species names that we had in our list. This is 
# because those names had synonyms in the tree of life database, 
# and we are using those names for the plot.
##############################################################
# Dealing with polytomies
##############################################################

# we can check for the existence of polytomies by running the 
# following code. If polytomies exist, the output will be 
# `FALSE`, and vice versa.

is.binary.tree(fish.tree) # there are some polytomies

# to take care of these polytomies, we are going to use a 
# randomization approach
set.seed(111)
tree_random <- multi2di(fish.tree,random=TRUE)
is.binary.tree(tree_random)

# exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))

setdiff(species, as.character(tree_random$tip.label)) #listed in our database but not in the tree
setdiff(as.character(tree_random$tip.label),species) # listed in the tree but not in our database

# they are the same species, the "problem" is that synonyms 
# have been used in the tree (see also `taxa.c`). We are going
# to leave all the names as in Open Tree of Life as it seems
# to be the most updated nomenclature

# Lets fix some evident wrong labels in the tree
tree_random$tip.label
tree_random$tip.label[c(32,181)]<-c("Oncorhynchus mykiss","Gadus morhua")
##############################################################
# Computing branch lengths
##############################################################
# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we need to make sure that tree labels and database
# use the same nomenclature

setdiff(fish$Species2, as.character(tree_random$tip.label))
setdiff(as.character(tree_random$tip.label),fish$Species2)

tree_random.fixed <- tree_random
tree_random.fixed$tip.label <- gsub("Pristiapogon exostigma","Apogon exostigma", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Ostorhinchus compressus","Apogon compressus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Mylopharyngodon piceus","Ctenopharyngodon piceus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Ostorhinchus cyanosoma","Apogon cyanosoma", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Asterropteryx semipunctata","Asterropteryx semipunctatus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Koumansetta rainfordi","Amblygobius rainfordi", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Haplochromis velifer","Astatotilapia velifer", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Cichlasoma urophthalmum","Mayaheros urophthalmus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Scolopsis bilineata","Scolopsis bilineatus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Paragobiodon xanthosoma","Paragobiodon xanthosomus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("mrcaott147479ott199966","Cyprinodon variegatus", tree_random.fixed$tip.label)

# check again
setdiff(fish$Species2, as.character(tree_random.fixed$tip.label))
setdiff(as.character(tree_random.fixed$tip.label),fish$Species2)
# perfectirijillo!

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(tree_random.fixed) # if not run line below
tree_random.fixed<-compute.brlen(tree_random.fixed,method = "Grafen")
is.ultrametric(tree_random.fixed)

#Set new directory to export outputs
setwd("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs/")

# Phylogenetic matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)

# finally, save matrix for future analyses
save(phylo_cor, file = "phylo_cor.Rdata")

#Plot the tree and export as figure
#pdf(paste("Phylogenetic tree for 194 species of fish with Pcrit data (",today,").pdf"),width=10,height=10)
jpeg(paste("Phylogenetic tree for 195 species of fish with Pcrit data (",today,").jpeg"),width=10,height=10,units = "in", res = 300)
plotTree(tree_random.fixed,fsize=0.5,ftype="i",type="fan")
dev.off()

# and as a .tre for phylogenetic corrections
write.tree(tree_random.fixed,"Phylogenetic tree for 195 species of fish with Pcrit data.tre")

#Create a binomial scientific name separated by (_)
fish$species <- as.factor(gsub(" ","_", fish$Species2))

# exporting dataset for analyses
write.csv(fish,paste("Pcrit data for 195 fish species (",today,").csv"),row.names=FALSE)

# saving session information with all packages versions for reproducibility purposes
sink("C:/Users/Invunche/Dropbox/collaboration/Jeroen/Outputs/phylogenetic_tree_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################