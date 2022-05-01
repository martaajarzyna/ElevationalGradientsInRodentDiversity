###Rodent diversity patterns along elevational gradients  - phylogenetic signal tests - per mountain
##Kohli, Miyajima and Jarzyna, GEB

##last modified by BAK: 13 July 2021


require(data.table)
require(ape)
require(phytools)
require(picante)
require(dplyr)
require(usedist)
require(geomorph)



###############################################
#######Phylogenetic signal tests
#######################################




###
#PREP DATA FOR PHYLO SIGNAL TESTS
###

#Read back in the final tree
phylotree <- read.nexus("Phylogenies.nex") 


#read in all species trait data.
Alltraits <- read.csv("SpTraits.csv", header=TRUE, row.names = 1)
#combine all Vertebrate diet subcategories
Alltraits <- mutate(Alltraits, Diet.Verts = Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk)
#log transform Mass, and HBL to make more normal  ##write over existing columns to keep same names
Alltraits$Body_mass_g <- log10(Alltraits$Body_mass_g)
Alltraits$HBL <- log10(Alltraits$HBL)
#set characters as factors
Alltraits$Life.mode <- as.factor(Alltraits$Life.mode)
Alltraits$Saltation <- as.factor(Alltraits$Saltation)

## Select traits of interest (vert diet combined & two very rare habitat types removed (caves... and marine.intertidal...))
Alltraits <- select(Alltraits, Sp_Name, CheekPouch, Diet.Inv, Diet.Verts, Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO, Activity.Nocturnal, Activity.Crepuscular, Activity.Diurnal, Body_mass_g, Forest, Savanna, Shrubland, Grassland, Wetland, Rocky.areas, Desert, HBL, RTL, RHFL, REL, Life.mode, Saltation)
str(Alltraits)
head(Alltraits)

#read names file
names <- read.csv("RodentTaxonomy_AcrossDatabases.csv")
names <- select(names, Sp_Name, SpName, Phylo_Name, PhyloName)
head(names)

#join with phylo names in name file
Traits.PhyloNames <- left_join(names, Alltraits, by = "Sp_Name")
head(Traits.PhyloNames)
str(Traits.PhyloNames)
dim(Traits.PhyloNames)





#######
#Estimate phy sig

##################################################################

# source the evolutionary-model Mantel function - Debastiani and Duarte 2017
EM.mantel<-function(tree, traits, runs = 999, euclidean= TRUE, sqrtPhylo=FALSE, checkdata = TRUE, ...){
  phylo.dist<-cophenetic(tree)
  if(sqrtPhylo){
    phylo.dist<-sqrt(phylo.dist)
  }
  if(checkdata){
    if(is.null(tree$tip.label)){
      stop("\n Error in tip labels of tree\n")
    }
    if(is.null(rownames(traits))){
      stop("\n Error in row names of traits\n")
    }
    match.names <- match(rownames(traits),rownames(phylo.dist))
    if(sum(is.na(match.names)) > 0){
      stop("\n There are species from traits data that are not on phylogenetic tree\n")
    }
    phylo.dist <- phylo.dist[match.names, match.names]
  }
  if(length(tree$tip.label) > dim(traits)[1]){
    warning("Tree have more species that species in traits data")
  }
  if(dim(phylo.dist)[1] != dim(traits)[1] & checkdata == FALSE){
    stop("\n Different number of species in tree and in traits data, use checkdata = TRUE\n")
  }
  gow.dist<-gowdis(traits, ...)
  if(euclidean){
    gow.sim<-1-gow.dist
    gow.dist<-sqrt(1-gow.sim)
  }
  traits.attr<-attr(gow.dist, "Types", exact = TRUE)
  res.mantel<-mantel(phylo.dist,gow.dist,permutations=runs)
  res.BM<-matrix(NA,runs,1)
  for(k in 1:runs){
    traits_sim<-matrix(NA,length(tree$tip.label),dim(traits)[2])
    rownames(traits_sim)<-tree$tip.label
    for(i in 1:dim(traits)[2]){
      traits_sim[,i]<-rTraitCont(tree,model="BM")
    }
    traits_sim<-decostand(traits_sim,method="standardize",MARGIN=2)
    traits_sim<-as.data.frame(traits_sim)	
    for(i in 1:dim(traits)[2]){
      if(traits.attr[i] == "B" | traits.attr[i] == "A"){
        probs<-sum(traits[,i])/dim(traits)[1]
        threshold<-quantile(traits_sim[,i],probs=1-probs)
        traits_sim[,i]<-ifelse(traits_sim[,i]>=threshold,1,0)
      }
      if(traits.attr[i] == "N" | traits.attr[i] == "O"){
        n.levels<-length(levels(traits[,i]))
        traits.levels<-levels(traits[,i])
        probs<-cumsum(table(traits[,i]))/sum(table(traits[,i]))
        probs<-probs[1:(n.levels-1)]
        threshold<-quantile(traits_sim[,i],probs=probs)
        threshold<-c(min(traits_sim[,i]),threshold,max(traits_sim[,i]))
        temp<-matrix(NA,length(traits_sim[,i]),1)
        for(j in 1:n.levels){
          if(j < n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<threshold[j+1], traits.levels[j],temp)
          }
          if(j == n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<=threshold[j+1], traits.levels[j],temp)
          }
        }
        traits_sim[,i]<-as.factor(temp)
        if(traits.attr[i] == "O"){
          traits_sim[,i]<-ordered(temp,levels=levels(traits[,i]))
        }
      }
    }
    if(checkdata == TRUE){
      match.names <- match(rownames(traits),rownames(traits_sim))
      traits_sim<-traits_sim[match.names,,drop=FALSE]
    }
    gow.dist.BM<-gowdis(traits_sim, ...)
    if(euclidean){
      gow.sim.BM<-1-gow.dist.BM
      gow.dist.BM<-sqrt(1-gow.sim.BM)
    }
    res.mantel.BM<-mantel(phylo.dist,gow.dist.BM,permutations=0)
    res.BM[k,1]<-res.mantel.BM$statistic
  }
  p.BM<-(sum(ifelse(res.BM[,1]>=res.mantel$statistic,1,0))+1)/(runs+1)
  p.NULL<-res.mantel$signif
  r.Mantel<-res.mantel$statistic
  RES<-list(perm.NULL=res.mantel$perm,perm.BM=res.BM[,1],r.Mantel=r.Mantel,p.NULL=p.NULL,p.BM=p.BM)
  return(RES)
}

###################################################

#single, continuous traits
### 1.
Mass <- funcdat_phy$Body_mass_g
names(Mass)=rownames(funcdat_phy)
Mass

#geomorph K
physignal(Mass, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) 
#picante
#phylosignal(Mass1, phylotree[[1]], reps = 999, checkdata = T) 
#phytools
#phylosig(phylotree[[1]], Mass1, method="K", test=T, nsim=1000) 
#phylosig(phylotree[[1]], Mass1, method="lambda", test=T, nsim=1000)

### 2.
HBL1 <- funcdat_phy$HBL
names(HBL1)=rownames(funcdat_phy)
HBL1

#geomorph K
physignal(HBL1, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) 
#picante
#phylosignal(HBL1, phylotree[[1]], reps = 999, checkdata = T) 
#phytools
#phylosig(phylotree[[1]], HBL1, method="K", test=T, nsim=1000) 
#phylosig(phylotree[[1]], HBL1, method="lambda", test=T, nsim=1000)

### 3.
RTL1 <- funcdat_phy$RTL
names(RTL1)=rownames(funcdat_phy)
RTL1

#geomorph K
physignal(RTL1, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) 
#picante
#phylosignal(RTL1, phylotree[[1]], reps = 999, checkdata = T) 
#phytools
#phylosig(phylotree[[1]], RTL1, method="K", test=T, nsim=1000) 
#phylosig(phylotree[[1]], RTL1, method="lambda", test=T, nsim=1000)

### 4.
RHFL1 <- funcdat_phy$RHFL
names(RHFL1)=rownames(funcdat_phy)
RHFL1

#geomorph K
physignal(RHFL1, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) 
#picante
#phylosignal(RHFL1, phylotree[[1]], reps = 999, checkdata = T) 
#phytools
#phylosig(phylotree[[1]], RHFL1, method="K", test=T, nsim=1000) 
#phylosig(phylotree[[1]], RHFL1, method="lambda", test=T, nsim=1000)

### 5.
REL1 <- funcdat_phy$REL
names(REL1)=rownames(funcdat_phy)


#geomorph K
physignal(REL1, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) 
#picante
#phylosignal(REL1, phylotree[[1]], reps = 999, checkdata = T) 
#phytools
#phylosig(phylotree[[1]], REL1, method="K", test=T, nsim=1000) 
#phylosig(phylotree[[1]], REL1, method="lambda", test=T, nsim=1000)


####Multivariate
### 1.
DietMorph <- as.matrix(funcdat_phy[,2:8])
DietMorph
physignal(DietMorph, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) #multivariate Kmult

#run with em.mantel function code
diet <- select(funcdat_phy, Diet.Inv:Diet.PlantO)
diet
EM.mantel(phylotree[[1]], diet, runs = 999)
#EMout <- EM.mantel(phylotree[[1]], diet, runs = 9)
# $r.Mantel = The Mantel observed statistic.
# $p.NULL = The p value from no phylogenetic structure (standard p value in Mantel test).
# $p.BM = The p value under simulation of traits from Brownian phylogenetic structure.

### 2.
#HabMult <- as.matrix(funcdat_phy[,13:19])
#HabMult
#physignal(HabMult, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) #multivariate Kmult

#run with em.mantel function code - not appropriate with Kmult because binary
habitats <- select(funcdat_phy, Forest:Desert)
habitats
EM.mantel(phylotree[[1]], habitats, runs = 999)


### 3.
#ActivMult <- as.matrix(funcdat_phy[,9:11])
#ActivMult
#physignal(ActivMult, phylotree[[1]], iter = 999, seed = NULL, print.progress = TRUE) #multivariate Kmult

#run with em.mantel function code - not appropriate with Kmult because binary
activity <- select(funcdat_phy, Activity.Nocturnal:Activity.Diurnal)
activity
EM.mantel(phylotree[[1]], activity, runs = 999)


####Categorical traits (not binary)
## 1.
LifeMode1 <- select(funcdat_phy, Life.mode)
LifeMode1
EM.mantel(phylotree[[1]], LifeMode1, runs = 999)

## 2.
Salt1 <- select(funcdat_phy, Saltation)
Salt1
EM.mantel(phylotree[[1]], Salt1, runs = 999)





###############################################################################
#### OCCURRENCE DATA PREP

#Read in data on species distributions
SpOcc <- read.csv("SpOcc_MtnElevBins.csv", header = T)
head(SpOcc)
str(SpOcc)

SpOcc$moID <- as.factor(SpOcc$Mtn) #create new column to hold factors of mtn IDs
SpOcc$spId <- as.numeric(as.factor(SpOcc$Species)) #create new column to hold numeric species IDs
SpOcc$spphyId <- as.numeric(as.factor(SpOcc$Phylo_Name)) #create new column to hold numeric species IDs
colnames(SpOcc) <- c("poId", "MtnName", "elevation", "SR", "spp", "sppPhy", "moId", "spId", "spphyId")
head(SpOcc)
length(unique(SpOcc$moId))

SpOcc <- SpOcc[SpOcc$SR >= 3,] 
length(unique(SpOcc$moId))
length(unique(SpOcc$spp))

moIDs <- unique(SpOcc$moId) 

MassPS <- matrix(NA,length(moIDs),100)
MassPSp <- matrix(NA,length(moIDs),100)

HBLPS <- matrix(NA,length(moIDs),100)
HBLPSp <-matrix(NA,length(moIDs),100)

RTLPSp <- matrix(NA,length(moIDs),100)
RTLPS <- matrix(NA,length(moIDs),100)

RHFLPS <- matrix(NA,length(moIDs),100)
RHFLPSp <- matrix(NA,length(moIDs),100)

RELPS <- matrix(NA,length(moIDs),100)
RELPSp <- matrix(NA,length(moIDs),100)

DietPS <- matrix(NA,length(moIDs),100)
DietPSp <- matrix(NA,length(moIDs),100)

HabitatPS <- matrix(NA,length(moIDs),100)
HabitatPSp <- matrix(NA,length(moIDs),100)

ActivityPS <- matrix(NA,length(moIDs),100)
ActivityPSp <- matrix(NA,length(moIDs),100)

LifeModePS <- matrix(NA,length(moIDs),100)
LifeModePSp <- matrix(NA,length(moIDs),100)

SaltPS <- matrix(NA,length(moIDs),100)
SaltPSp <- matrix(NA,length(moIDs),100)

#cycle through each of 100 phylogenetic trees to calculate phy signal for each trait
for (i in 1:100){	
  pdi <- phylotree[[i]]
  
  for (k in 1:length(moIDs)) {
    comm <- subset(SpOcc, moId == moIDs[k])
    nspec <- nrow(comm)
    comm <- setDF(comm)
    sp.id <- pdi$tip.label
    
    rmTip <- as.character(sp.id[! sp.id %in% comm$sppPhy]) # Remove missing tips 
    subTree <- drop.tip(pdi, tip=rmTip)  # remove those branches
    
    #subset traits of just the mtn via phylonames
    subTraits <- left_join(comm, Traits.PhyloNames, by = c("sppPhy" = "Phylo_Name"))   # reorder the trait data matrix to match the tip labels by species name
    subTraits2 <- subTraits[match(subTree$tip.label, subTraits$sppPhy),]      
    row.names(subTraits2) <- subTraits2$sppPhy #set row names
    subFuncdat <- subTraits2[,13:ncol(subTraits2)]   #strip down to just traits
    
    ### 1. body mass
    #Mass1 <- subFuncdat$Body_mass_g
    #names(Mass1)=rownames(subFuncdat)
    #geomorph K
    #PSi <- physignal(Mass1, subTree, iter = 999, seed = NULL, print.progress = F) 
    #MassPS[k,i] <- PSi$phy.signal
    #MassPSp[k,i] <- PSi$pvalue
    
    ### 2. head body length
    HBL1 <- subFuncdat$HBL
    names(HBL1)=rownames(subFuncdat)
    #geomorph K
    PS2 <- physignal(HBL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    HBLPS[k,i] <- PS2$phy.signal
    HBLPSp[k,i] <- PS2$pvalue
    
    ### 3.relative tail length
    RTL1 <- subFuncdat$RTL
    names(RTL1)=rownames(subFuncdat)
    #geomorph K
    PS3 <- physignal(RTL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    RTLPS[k,i] <- PS3$phy.signal
    RTLPSp[k,i] <- PS3$pvalue
   
    ### 4. relative hind foot length
    RHFL1 <- subFuncdat$RHFL
    names(RHFL1)=rownames(subFuncdat)
    #geomorph K
    PS4 <- physignal(RHFL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    RHFLPS[k,i] <- PS4$phy.signal
    RHFLPSp[k,i] <- PS4$pvalue

    ### 5.relative ear length
    REL1 <- subFuncdat$REL
    names(REL1)=rownames(subFuncdat)
    #geomorph K
    PS5 <- physignal(REL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    RELPS[k,i] <- PS5$phy.signal
    RELPSp[k,i] <- PS5$pvalue
    
    
    ####Multivariate discrete traits
    # $r.Mantel = The Mantel observed statistic.
    # $p.NULL = The p value from no phylogenetic structure (standard p value in Mantel test).
    # $p.BM = The p value under simulation of traits from Brownian phylogenetic structure.
    
    ### 6.diet
    Diet <- select(subFuncdat, Diet.Inv:Diet.PlantO)
    PS6 <- EM.mantel(subTree, Diet, runs = 999)
    DietPS[k,i] <- PS6$r.Mantel
    DietPSp[k,i] <- PS6$p.BM
  
    ### 7.habitats
    Habitat <- select(subFuncdat, Forest:Desert)
    PS7 <- EM.mantel(subTree, Habitat, runs = 999)
    HabitatPS[k,i] <- PS7$r.Mantel
    HabitatPSp[k,i] <- PS7$p.BM
    
    ### 8. diel activty
    Activity <- select(subFuncdat, Activity.Nocturnal:Activity.Diurnal)
    PS8 <- EM.mantel(subTree, Activity, runs = 999)
    ActivityPS[k,i] <- PS8$r.Mantel
    ActivityPSp[k,i] <- PS8$p.BM
    
    ####single Categorical traits (not binary)
    ## 9. life mode (vertical stratum)
    LifeMode1 <- select(subFuncdat, Life.mode)
    PS9 <- EM.mantel(subTree, LifeMode1, runs = 999)
    LifeModePS[k,i] <- PS9$r.Mantel
    LifeModePSp[k,i] <- PS9$p.BM
    
    ## 10. saltatory locomotion
    Salt1 <- select(subFuncdat, Saltation)
    PS10 <- EM.mantel(subTree, Salt1, runs = 999)
    SaltPS[k,i] <- PS10$r.Mantel
    SaltPSp[k,i] <- PS10$p.BM
 
    ### 11. is binary cheek pouches - done separately.
    
    
  }
} #closing of tree bracket




#########################################################################

####
##summarize outputs
######

#1. Mass
MassPSid <- as.data.frame(cbind(MassPS,moIDs)) #phylosig values
MassPSid$mtn <- 0
MassPSpid <- as.data.frame(cbind(MassPSp,moIDs))  #phylosig signifs
MassPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  MassPSid[k,102] <- comm[2]
  MassPSpid[k,102] <- comm[2]
}
head(MassPSid[,99:102])
head(MassPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
MassPSid$MassPS_mean <- rowMeans(MassPSid[,1:100])
#tally how many trees had significant phylosig for the trait
MassPSpid$Massp_count <- rowSums(MassPSpid[,1:100] < 0.05)
head(MassPSid[,99:103])
head(MassPSpid[,99:103])
MassPS_summary <- as.data.frame(cbind(MassPSid[103], MassPSpid[103]))
write.csv(MassPSid, file = "MassPS.csv", row.names = F)
write.csv(MassPSpid, file = "MassPSp.csv", row.names = F)

#2. HBL
HBLPSid <- as.data.frame(cbind(HBLPS,moIDs)) #phylosig values
HBLPSid$mtn <- 0
HBLPSpid <- as.data.frame(cbind(HBLPSp,moIDs))  #phylosig signifs
HBLPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  HBLPSid[k,102] <- comm[2]
  HBLPSpid[k,102] <- comm[2]
}
head(HBLPSid[,99:102])
head(HBLPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
HBLPSid$HBLPS_mean <- rowMeans(HBLPSid[,1:100])
#tally how many trees had significant phylosig for the trait
HBLPSpid$HBLp_count <- rowSums(HBLPSpid[,1:100] < 0.05)
head(HBLPSid[,99:103])
head(HBLPSpid[,99:103])
HBLPS_summary <- as.data.frame(cbind(HBLPSid[103], HBLPSpid[103]))
write.csv(HBLPSid, file = "HBLPS.csv", row.names = F)
write.csv(HBLPSpid, file = "HBLPSp.csv", row.names = F)

#3. RTL
RTLPSid <- as.data.frame(cbind(RTLPS,moIDs)) #phylosig values
RTLPSid$mtn <- 0
RTLPSpid <- as.data.frame(cbind(RTLPSp,moIDs))  #phylosig signifs
RTLPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  RTLPSid[k,102] <- comm[2]
  RTLPSpid[k,102] <- comm[2]
}
head(RTLPSid[,99:102])
head(RTLPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
RTLPSid$RTLPS_mean <- rowMeans(RTLPSid[,1:100])
#tally how many trees had significant phylosig for the trait
RTLPSpid$RTLp_count <- rowSums(RTLPSpid[,1:100] < 0.05)
head(RTLPSid[,99:103])
head(RTLPSpid[,99:103])
RTLPS_summary <- as.data.frame(cbind(RTLPSid[103], RTLPSpid[103]))
write.csv(RTLPSid, file = "RTLPS.csv", row.names = F)
write.csv(RTLPSpid, file = "RTLPSp.csv", row.names = F)

#4. RHFL
RHFLPSid <- as.data.frame(cbind(RHFLPS,moIDs)) #phylosig values
RHFLPSid$mtn <- 0
RHFLPSpid <- as.data.frame(cbind(RHFLPSp,moIDs))  #phylosig signifs
RHFLPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  RHFLPSid[k,102] <- comm[2]
  RHFLPSpid[k,102] <- comm[2]
}
head(RHFLPSid[,99:102])
head(RHFLPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
RHFLPSid$RHFLPS_mean <- rowMeans(RHFLPSid[,1:100])
#tally how many trees had significant phylosig for the trait
RHFLPSpid$RHFLp_count <- rowSums(RHFLPSpid[,1:100] < 0.05)
head(RHFLPSid[,99:103])
head(RHFLPSpid[,99:103])
RHFLPS_summary <- as.data.frame(cbind(RHFLPSid[103], RHFLPSpid[103]))
write.csv(RHFLPSid, file = "RHFLPS.csv", row.names = F)
write.csv(RHFLPSpid, file = "RHFLPSp.csv", row.names = F)

#5. REL
RELPSid <- as.data.frame(cbind(RELPS,moIDs)) #phylosig values
RELPSid$mtn <- 0
RELPSpid <- as.data.frame(cbind(RELPSp,moIDs))  #phylosig signifs
RELPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  RELPSid[k,102] <- comm[2]
  RELPSpid[k,102] <- comm[2]
}
head(RELPSid[,99:102])
head(RELPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
RELPSid$RELPS_mean <- rowMeans(RELPSid[,1:100])
#tally how many trees had significant phylosig for the trait
RELPSpid$RELp_count <- rowSums(RELPSpid[,1:100] < 0.05)
head(RELPSid[,99:103])
head(RELPSpid[,99:103])
RELPS_summary <- as.data.frame(cbind(RELPSid[103], RELPSpid[103]))
write.csv(RELPSid, file = "RELPS.csv", row.names = F)
write.csv(RELPSpid, file = "RELPSp.csv", row.names = F)

#6. Diet
DietPSid <- as.data.frame(cbind(DietPS,moIDs)) #phylosig values
DietPSid$mtn <- 0
DietPSpid <- as.data.frame(cbind(DietPSp,moIDs))  #phylosig signifs
DietPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  DietPSid[k,102] <- comm[2]
  DietPSpid[k,102] <- comm[2]
}
head(DietPSid[,99:102])
head(DietPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
DietPSid$DietPS_mean <- rowMeans(DietPSid[,1:100])
#tally how many trees had significant phylosig for the trait
DietPSpid$Dietp_count <- rowSums(DietPSpid[,1:100] < 0.05)
head(DietPSid[,99:103])
head(DietPSpid[,99:103])
DietPS_summary <- as.data.frame(cbind(DietPSid[103], DietPSpid[103]))
write.csv(DietPSid, file = "DietPS.csv", row.names = F)
write.csv(DietPSpid, file = "DietPSp.csv", row.names = F)

#7. Habitat
HabitatPSid <- as.data.frame(cbind(HabitatPS,moIDs)) #phylosig values
HabitatPSid$mtn <- 0
HabitatPSpid <- as.data.frame(cbind(HabitatPSp,moIDs))  #phylosig signifs
HabitatPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  HabitatPSid[k,102] <- comm[2]
  HabitatPSpid[k,102] <- comm[2]
}
head(HabitatPSid[,99:102])
head(HabitatPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
HabitatPSid$HabitatPS_mean <- rowMeans(HabitatPSid[,1:100])
#tally how many trees had significant phylosig for the trait
HabitatPSpid$Habitatp_count <- rowSums(HabitatPSpid[,1:100] < 0.05)
head(HabitatPSid[,99:103])
head(HabitatPSpid[,99:103])
HabitatPS_summary <- as.data.frame(cbind(HabitatPSid[103], HabitatPSpid[103]))
write.csv(HabitatPSid, file = "HabitatPS.csv", row.names = F)
write.csv(HabitatPSpid, file = "HabitatPSp.csv", row.names = F)

#8. Activity
ActivityPSid <- as.data.frame(cbind(ActivityPS,moIDs)) #phylosig values
ActivityPSid$mtn <- 0
ActivityPSpid <- as.data.frame(cbind(ActivityPSp,moIDs))  #phylosig signifs
ActivityPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  ActivityPSid[k,102] <- comm[2]
  ActivityPSpid[k,102] <- comm[2]
}
head(ActivityPSid[,99:102])
head(ActivityPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
ActivityPSid$ActivityPS_mean <- rowMeans(ActivityPSid[,1:100])
#tally how many trees had significant phylosig for the trait
ActivityPSpid$Activityp_count <- rowSums(ActivityPSpid[,1:100] < 0.05)
head(ActivityPSid[,99:103])
head(ActivityPSpid[,99:103])
ActivityPS_summary <- as.data.frame(cbind(ActivityPSid[103], ActivityPSpid[103]))
write.csv(ActivityPSid, file = "ActivityPS.csv", row.names = F)
write.csv(ActivityPSpid, file = "ActivityPSp.csv", row.names = F)

#9. LifeMode
LifeModePSid <- as.data.frame(cbind(LifeModePS,moIDs)) #phylosig values
LifeModePSid$mtn <- 0
LifeModePSpid <- as.data.frame(cbind(LifeModePSp,moIDs))  #phylosig signifs
LifeModePSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  LifeModePSid[k,102] <- comm[2]
  LifeModePSpid[k,102] <- comm[2]
}
head(LifeModePSid[,99:102])
head(LifeModePSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
LifeModePSid$LifeModePS_mean <- rowMeans(LifeModePSid[,1:100])
#tally how many trees had significant phylosig for the trait
LifeModePSpid$LifeModep_count <- rowSums(LifeModePSpid[,1:100] < 0.05)
head(LifeModePSid[,99:103])
head(LifeModePSpid[,99:103])
LifeModePS_summary <- as.data.frame(cbind(LifeModePSid[103], LifeModePSpid[103]))
write.csv(LifeModePSid, file = "LifeModePS.csv", row.names = F)
write.csv(LifeModePSpid, file = "LifeModePSp.csv", row.names = F)

#10. Salt
SaltPSid <- as.data.frame(cbind(SaltPS,moIDs)) #phylosig values
SaltPSid$mtn <- 0
SaltPSpid <- as.data.frame(cbind(SaltPSp,moIDs))  #phylosig signifs
SaltPSpid$mtn <- 0
for (k in 1:length(moIDs)){
  comm <- subset(SpOcc, moId == moIDs[k])
  SaltPSid[k,102] <- comm[2]
  SaltPSpid[k,102] <- comm[2]
}
head(SaltPSid[,99:102])
head(SaltPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
SaltPSid$SaltPS_mean <- rowMeans(SaltPSid[,1:100])
#tally how many trees had significant phylosig for the trait
SaltPSpid$Saltp_count <- rowSums(SaltPSpid[,1:100] < 0.05)
head(SaltPSid[,99:103])
head(SaltPSpid[,99:103])
SaltPS_summary <- as.data.frame(cbind(SaltPSid[103], SaltPSpid[103]))
write.csv(SaltPSid, file = "SaltPS.csv", row.names = F)
write.csv(SaltPSpid, file = "SaltPSp.csv", row.names = F)





###########################################
######################################################################

###binary traits
require(caper)  #masks select from dplyr
#caper - phylo.d for binary traits (Fritz and Purvis) The D statistic is equal to 1 if the observed binary trait has a phylogenetically random distribution across the tips of the phylogeny and to 0 if the observed trait is as clumped as if it had evolved by Brownian motion under our threshold model.
####USE for CheekPouch


########
SpOcc2 <- SpOcc[SpOcc$moId != "Tucuman" & SpOcc$moId != "Itatiaia" & SpOcc$moId != "Manu"& SpOcc$moId != "Sneznik"& SpOcc$moId != "JiajinS"& SpOcc$moId != "Luoji"& SpOcc$moId != "Sejila"& SpOcc$moId != "Yulong"& SpOcc$moId != "AilaoE"& SpOcc$moId != "JiajinN"& SpOcc$moId != "KitangladRange"& SpOcc$moId != "Rwenzori"& SpOcc$moId != "KinabaluNP",]
moIDs2 <- unique(SpOcc2$moId) 


PouchPS <- matrix(NA,length(moIDs),100)
PouchPSp <- matrix(NA,length(moIDs),100)
#cycle through each of 100 phylogenetic trees to calculate phy signal for each trait
for (i in 1:100){	
  pdi <- phylotree[[i]]
    for (k in 1:length(moIDs2)) {
    comm <- subset(SpOcc2, moId == moIDs2[k])
    nspec <- nrow(comm)
    comm <- setDF(comm)
    sp.id <- pdi$tip.label
    
    rmTip <- as.character(sp.id[! sp.id %in% comm$sppPhy]) # Remove missing tips
    subTree <- drop.tip(pdi, tip=rmTip)  # remove those branches
        #subset traits of just the mtn via phylonames
    subTraits <- left_join(comm, Traits.PhyloNames, by = c("sppPhy" = "Phylo_Name"))   # reorder the trait data matrix to match the tip labels by species name
    subTraits2 <- subTraits[match(subTree$tip.label, subTraits$sppPhy),]      
    row.names(subTraits2) <- subTraits2$sppPhy #set row names
    subFuncdat <- subTraits2[,13:ncol(subTraits2)]   #strip down to just traits
        #DFs for binary trait - cheeck pouches
    Pouch1 <- as.data.frame(subFuncdat$CheekPouch)
    Pouch1$names <- rownames(subFuncdat)
    colnames(Pouch1) <- c("pouches", "names")
    #Fritz and Purvis D
    PS11 <- phylo.d(Pouch1, subTree, names.col = names, binvar = pouches, permut = 1000)
    PouchPS[k,i] <- PS11$DEstimate
    PouchPSp[k,i] <- PS11$Pval1
        #$DEstimate for Mantel test value
    #$Pval1 and Pval0 for probability of resulting from random str (ie. signif phylo signal) and brownian str.
      }
  }

####
##summarize output
######

#11. Pouch
PouchPSid <- as.data.frame(cbind(PouchPS,moIDs2)) #phylosig values
PouchPSid$mtn <- 0
PouchPSpid <- as.data.frame(cbind(PouchPSp,moIDs2))  #phylosig signifs
PouchPSpid$mtn <- 0
for (k in 1:length(moIDs2)){
  comm <- subset(SpOcc2, moId == moIDs2[k])
  PouchPSid[k,102] <- comm[2]
  PouchPSpid[k,102] <- comm[2]
}
head(PouchPSid[,99:102])
head(PouchPSpid[,99:102])
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
PouchPSid$PouchPS_mean <- rowMeans(PouchPSid[,1:100])
#tally how many trees had significant phylosig for the trait
PouchPSpid$Pouchp_count <- rowSums(PouchPSpid[,1:100] < 0.05)
head(PouchPSid[,99:103])
head(PouchPSpid[,99:103])
PouchPS_summary <- as.data.frame(cbind(PouchPSid[102:103], PouchPSpid[103]))
write.csv(PouchPSid, file = "PouchPS.csv", row.names = F)
write.csv(PouchPSpid, file = "PouchPSp.csv", row.names = F)






#######################################
###summarize ALL phylogenetic signal test outputs
MtnPhySig_moIDs1 <- as.data.frame(cbind(MassPS_summary, HBLPS_summary, RTLPS_summary, RHFLPS_summary, RELPS_summary, DietPS_summary, HabitatPS_summary, ActivityPS_summary, LifeModePS_summary, SaltPS_summary))
MtnPhySig_moIDs1
write.csv(MtnPhySig_moIDs1, file = "IndivMtnPhySigSummaryMinusCheeks.csv", row.names = F)

