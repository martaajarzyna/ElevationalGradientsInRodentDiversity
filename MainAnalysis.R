###Rodent diversity patterns along elevational gradients  - Empirical calculations
##Kohli, Miyajima and Jarzyna, GEB
##last modified by BAK: 25 June 2021


require(data.table)
require(gawdis)
require(FD)
require(ape)
require(phytools)
require(picante)
require(hypervolume)
require(dplyr)
require(usedist)

###################################################
########STEP 1: Functional & phylogenetic diversity calculation
###################################################

#### OCCURRENCE DATA PREP

#Read in data on species distributions
SpOcc <- read.csv("SpOcc_MtnElevBins.csv", header = T)

SpOcc$MtnNum <- as.factor(SpOcc$Mtn) #create new column to hold factors of mtn IDs
SpOcc$spId <- as.numeric(as.factor(SpOcc$Species)) #create new column to hold numeric species IDs
SpOcc$spphyId <- as.numeric(as.factor(SpOcc$Phylo_Name)) #create new column to hold numeric species IDs

colnames(SpOcc) <- c("poId", "MtnName", "elevation", "SR", "spp", "sppPhy", "moId", "spId", "spphyId")
head(SpOcc)
length(unique(SpOcc$poId))###poIds are a unique code for each mtn/elev bin combo (point location). moIDs are Mountain IDs, elevation is the minimum of each 100m bin, SR is species richness of a bin (from interpolated elevational ranges of species), spp and sppPhy are species names to match the trait and phylogenetic dataset, respectively.

###This is structured as a presence-only dataset.

###FOR SOME INDICES, HAVE TO REMOVE/IGNORE THE BINS WITH LESS THAN 3 SPECIES. - as such, need to make it uniform across analyses
SpOcc <- SpOcc[SpOcc$SR >= 3,] # >=3 removes 90 rows (=54 bins or poids = total of 968 poids); No species are removed by this action (ie. none are only found only in 1- or 2-sp comms)
#FEve cannot be calculated for communities with <3 functionally singular species.
#convex hulls require more than one continuous trait or dimension in 'x' (2 sp only allows for 1 dimension). If not, it's measured as the range, NOT as the convex hull volume; and m>n.
length(unique(SpOcc$poId))
length(unique(SpOcc$spp))

poIDs <- unique(SpOcc$poId) 
#################################################################

#### TRAIT DATA PREP

#read in trait matrix
sp.traits <- read.csv("SpTraits.csv", header=TRUE, row.names = 1)
#combine all Vertebrate diet subcategories
sp.traits <- mutate(sp.traits, Diet.Verts = Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk)

rownames(sp.traits) <- sp.traits[,1]
sp <- as.character(sp.traits[,1])
sp.ids <- sp.traits[,1:3]
head(sp.traits)


## Select traits of interest (vert diet combined & two very rare habitat types removed (caves... and marine.intertidal...))
funcdat <- select(sp.traits, CheekPouch, Diet.Inv, Diet.Verts, Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO, Activity.Nocturnal, Activity.Crepuscular, Activity.Diurnal, Body_mass_g, Forest, Savanna, Shrubland, Grassland, Wetland, Rocky.areas, Desert, HBL, RTL, RHFL, REL, Life.mode, Saltation) #Note that downstream methods adjust for range of individual traits and cross-correlation

#log transform Mass, and HBL to make more normal
##write over existing columns to keep same names
funcdat$Body_mass_g <- log10(funcdat$Body_mass_g)
funcdat$HBL <- log10(funcdat$HBL)

#set characters as factors
funcdat$Life.mode <- as.factor(funcdat$Life.mode)
funcdat$Saltation <- as.factor(funcdat$Saltation)
summary(funcdat)


###############################################
#use gawdis for distance matrix calculation to improve the equal weighting of traits/trait groupings

#designate trait groupings for gawdis (groups = different values)
Trait.groups = c(1, rep(2, 7), rep(3, 3), 4, rep(5, 7), 4, 6, 7, 8, 9, 10)
#check trait groups
funcdatcols <- colnames(funcdat)
cbind.data.frame(funcdatcols, Trait.groups)

#identify the trait groups that represent fuzzy coding or dummy variables - values matching the group values above.
fuzzy.groups = c(2, 3, 5)  #diet, activity time, habitats

#create functional dissimilarity matrix with equal contribution of traits and trait groups (fuzzy coded categories and log mass/head-body length (corr ~ 95%)) via gawdis
#must use optimized to avoid the negative weights of a few traits
funcdist <- gawdis(funcdat, w.type = "optimized", groups = Trait.groups, fuzzy = fuzzy.groups, opti.maxiter = 300) #300 is the default

#saveRDS(funcdist, file="gawdis_distmat.rds")  #must save out and use the same matrix because the result will vary slightly run to run (since the analytic solution won't complete) 
funcdist <- readRDS("gawdis_distmat.rds") 

#check out result
attr(funcdist, "correls")  #initial weights
attr(funcdist, "weights") #weights finally given to traits
attr(funcdist, "group.correls") #weights of each group

plot(attr(funcdist, "correls"))
plot(attr(funcdist, "weights"))
plot(attr(funcdist, "group.correls"))

############### 



#############################################
#### A. Functional richness - convex hull volume via Hypervolume package
poIDs <- unique(SpOcc$poId)

## Principal Coordinate Analysis
pcoa.calc <- pcoa(funcdist, correction="lingoes")
nbdim <- ncol(pcoa.calc$vectors)
eigen <- pcoa.calc$values
vect <- pcoa.calc$vectors[,1:nbdim]
rownames(vect) <- rownames(funcdat)
colnames(vect) <- paste("PC",1:nbdim,sep="")

FDH <- matrix(NA,length(poIDs),1)

for (i in 1:length(poIDs)){
  comm <- SpOcc[SpOcc$poId ==  poIDs[i],]
  sp.i <- unique(comm$spp)
  vect.i <- vect[rownames(vect) %in% sp.i,1:2]
  fdh.i <- expectation_convex(vect.i, check.memory=F)
  #saveRDS(fdh.i, file=paste0("HullVol_2D_",poIDs[i],".rds"))
 
  FDH[i,1] <- fdh.i@Volume # Volume of the convex hull
 }

FDHid <- as.data.frame(cbind(FDH,poIDs))
#FDHid$longitude <- 0
#FDHid$latitude <- 0
FDHid$elevation <- 0
FDHid$moIDs <- 0
FDHid$SR <- 0

for (k in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[k])
  FDHid[k,3:5] <- comm[1,2:4]
}

colnames(FDHid) <- c("FDH.vol","poIDs","moIDs","elevation", "SR")
head(FDHid)

saveRDS(FDHid, file="FDH_convexhull_2d.rds")



####################################################################
#### B. Indices of distance-based functional diversity

## Quantify non-abundance-weighted functional dispersion (FDis), functional evenness (FEve), and other indices via dbFD package for each point location (=elev bin)

DBFD <- matrix(NA,length(poIDs),5)

for (i in 1:length(poIDs)) {
  comm <- subset(SpOcc, poId == poIDs[i])
  nspec <- nrow(comm)
  comm <- setDF(comm)
  sp.i <- unique(comm$spp)
  funcdist.i <- dist_subset(funcdist, sp.i)  #package usedist extracts comm sp list to a subset dist matrix (doesn't work for N=1 bc doesn't have labels)  
  dbfd.i <- dbFD(x = funcdist.i, w.abun = F, corr = "lingoes", calc.FRic = F, calc.FGR = F, calc.CWM = F, calc.FDiv = F, print.pco = F) ##set corr = "lingoes" to make equivalent across all indices/packages used 
  
  #saveRDS(dbfd.i, file=paste0("dbfd_2D_",poIDs[i],".rds"))
  
  DBFD[i, 1] <- dbfd.i$FDis
  DBFD[i, 2] <- dbfd.i$FEve
  DBFD[i, 3] <- dbfd.i$RaoQ
  DBFD[i, 4] <- dbfd.i$nbsp
  DBFD[i, 5] <- dbfd.i$sing.sp
}

DBFDid <- as.data.frame(cbind(DBFD,poIDs))
DBFDid$moIDs <- 0
DBFDid$elevation <- 0
DBFDid$SR <- 0

for (i in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[i])
  DBFDid[i,7:9] <- comm[1,2:4]
}

colnames(DBFDid) <- c("FDis", "FEve", "RaoQ", "nbsp", "sing.sp", "poIDs","moIDs", "elevation", "SR")
head(DBFDid)
unique(DBFDid$nbsp == DBFDid$sing.sp) #checking SR calcs
unique(DBFDid$nbsp == DBFDid$SR)  #checking SR calcs
summary(DBFDid)

#write.csv(DBFDid, file = "DBFD_2D.csv") 
saveRDS(DBFDid, file="DBFD_2D.rds")


################
#############
###############################################################

#### C. Phylogenetic divesity (PD) based on branch lengths and topology directly
phylotree <- read.nexus("Phylogenies.nex") 
####need to match to the phylonames instead of trait names: changed column below
poIDs <- unique(SpOcc$poId)

PD <- matrix(NA,length(poIDs),100) # Phylogenetic diversity
PDSm <- matrix(NA,length(poIDs),100) # Mean phylogenetic distinctness
PDS <- list() # To hold species-level phylogenetic distinctness
PDSi <- list()

for (i in 1:100){	
  pdi <- phylotree[[i]]
  
  for (k in 1:length(poIDs)) {
    comm <- subset(SpOcc, poId == poIDs[k])
    nspec <- nrow(comm)
    comm <- setDF(comm)
    sp.id <- pdi$tip.label
    
    rmTip <- as.character(sp.id[! sp.id %in% comm$sppPhy]) # Remove missing tips 
    subTree <- drop.tip(pdi, tip=rmTip)  # remove those branches
    PD[k,i] <- sum(subTree$edge.length) # Sum the branch length
    di <- evol.distinct(subTree, type = "fair.proportion",scale = FALSE, use.branch.lengths = TRUE)
    PDSm[k,i] <- psych::geometric.mean(di[,2])
    PDSi[[k]] <- di[,2]
  }
  
  PDS[[i]] <- PDSi
}


PDid <- as.data.frame(cbind(PD,poIDs)) #PD proper (branch length sum)
PDid$moIDs <- 0
PDid$elevation <- 0
PDid$SR <- 0

PDSid <- as.data.frame(cbind(PDSm,poIDs))  #phylogenetic distinctness
PDSid$moIDs <- 0
PDSid$elevation <- 0
PDSid$SR <- 0

for (k in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[k])
  PDid[k,102:104] <- comm[1,2:4]
  PDSid[k,102:104] <- comm[1,2:4]
}

head(PDid[,99:104])
head(PDSid[,99:104])

#calculate mean values across the 100 trees randomly selected from the 10000 in Upham.
PDid$PD_mean <- rowMeans(PDid[,1:100])
PDSid$PDi_mean <- rowMeans(PDSid[,1:100])

head(PDid[,99:105])
head(PDSid[,99:105])

saveRDS(PDid, file="PD.rds")
saveRDS(PDSid, file="PDSm.rds")
saveRDS(PDS, file="PDS.rds")

#write.csv(PDid, file = "PD_GlobalResults.csv", row.names = F)

########################################################################
#####PEve using dbfd and MPD using picante (distance-based rather than direct branch lengths) for each point location (=elev bin)

PEve <- matrix(NA,length(poIDs),100) # Phylogenetic evenness
MPD <- matrix(NA,length(poIDs),100) # Mean pairwise phylogenetic distance

for (i in 1:100){	
  pdi <- phylotree[[i]]  #loop through each global tree
  phydist <- cophenetic.phylo(pdi) #global phylogenetic distance matrix
 
  for (k in 1:length(poIDs)) {
    comm <- subset(SpOcc, poId == poIDs[k])
    nspec <- nrow(comm)
    comm <- setDF(comm)
    sp.id <- unique(comm$sppPhy)
    #make a sample matrix that mpd can read (with species names and any numeric value to indicate presence)
    phylocomm.i <- select(comm, poId, SR, sppPhy) #comm name, sp abundance, sp name
    samp.i <- sample2matrix(phylocomm.i) #convert to matrix format
    
    phydist.i <- dist_subset(phydist, sp.id)  #package usedist extracts comm sp list to a subset dist matrix (doesn't work for N=1 bc doesn't have labels)  
   dbPd.i <- dbFD(x = phydist.i, w.abun = F, corr = "lingoes", calc.FRic = F, calc.FGR = F, calc.CWM = F, calc.FDiv = F, print.pco = F) 
   PEve[k,i] <- dbPd.i$FEve
  
   phydistmat.i <- as.matrix(phydist.i) #matrix format to match to sample format
   MPD[k,i] <- mpd(samp.i, phydistmat.i, abundance.weighted = F)
  }
}


PEveid <- as.data.frame(cbind(PEve,poIDs)) #PEve
PEveid$moIDs <- 0
PEveid$elevation <- 0
PEveid$SR <- 0

MPDid <- as.data.frame(cbind(MPD,poIDs))  #mean phylogenetic distance
MPDid$moIDs <- 0
MPDid$elevation <- 0
MPDid$SR <- 0

for (k in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[k])
  PEveid[k,102:104] <- comm[1,2:4]
  MPDid[k,102:104] <- comm[1,2:4]
}

head(PEveid[,99:104])
head(MPDid[,99:104])

#calculate mean values across the 100 trees randomly selected from the 10000 in Upham.
PEveid$PEve_mean <- rowMeans(PEveid[,1:100])
MPDid$MPD_mean <- rowMeans(MPDid[,1:100])

head(PEveid[,99:105])
head(MPDid[,99:105])

saveRDS(PEveid, file="PEve.rds")
saveRDS(MPDid, file="MPD.rds")

#write.csv(PEveid, file = "PEve_GlobalResults.csv", row.names = F)
#write.csv(MPDid, file = "MPD_GlobalResults.csv", row.names = F)
#################################################################################






###################################################
#######STEP 2: Null models
###################################################
#### Randomize species presences and calculate indices ---- see MainAnalysisNulls script



###################################################
#########STEP 3: INLA models 
###################################################
###see ModelTesting scripts


