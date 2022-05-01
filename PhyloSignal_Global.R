###Rodent diversity patterns along elevational gradients  - phylogenetic signal tests - global
##Kohli, Miyajima and Jarzyna, GEB
##last modified by BAK: 14 July 2021


require(data.table)
require(ape)
require(phytools)
require(picante)
require(dplyr)
require(usedist)
require(geomorph)
require(FD)



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
#####

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

#moIDs <- unique(SpOcc$moId) 

MassPS <- matrix(NA,1, 100)
MassPSp <- matrix(NA,1,100)

HBLPS <- matrix(NA,1,100)
HBLPSp <-matrix(NA,1,100)

RTLPSp <- matrix(NA,1,100)
RTLPS <- matrix(NA,1,100)

RHFLPS <- matrix(NA,1,100)
RHFLPSp <- matrix(NA,1,100)

RELPS <- matrix(NA,1,100)
RELPSp <- matrix(NA,1,100)

DietPS <- matrix(NA,1,100)
DietPSp <- matrix(NA,1,100)

HabitatPS <- matrix(NA,1,100)
HabitatPSp <- matrix(NA,1,100)

ActivityPS <- matrix(NA,1,100)
ActivityPSp <- matrix(NA,1,100)

LifeModePS <- matrix(NA,1,100)
LifeModePSp <- matrix(NA,1,100)

SaltPS <- matrix(NA,1,100)
SaltPSp <- matrix(NA,1,100)

#cycle through each of 100 phylogenetic trees to calculate phy signal for each trait
for (i in 1:100){	
  subTree <- phylotree[[i]]
  
    # reorder the trait data matrix to match the tip labels by species name
    subTraits2 <- Traits.PhyloNames[match(subTree$tip.label, Traits.PhyloNames$Phylo_Name),]      
    row.names(subTraits2) <- subTraits2$Phylo_Name #set row names
    subFuncdat <- subTraits2[,5:ncol(subTraits2)]   #strip down to just traits
    
    ### 1. body mass
    Mass1 <- subFuncdat$Body_mass_g
    names(Mass1)=rownames(subFuncdat)
    #geomorph K
    PSi <- physignal(Mass1, subTree, iter = 999, seed = NULL, print.progress = F) 
    MassPS[i] <- PSi$phy.signal
    MassPSp[i] <- PSi$pvalue

    ### 2. head body length
    HBL1 <- subFuncdat$HBL
    names(HBL1)=rownames(subFuncdat)
    #geomorph K
    PS2 <- physignal(HBL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    HBLPS[i] <- PS2$phy.signal
    HBLPSp[i] <- PS2$pvalue
    
    ### 3.relative tail length
    RTL1 <- subFuncdat$RTL
    names(RTL1)=rownames(subFuncdat)
    #geomorph K
    PS3 <- physignal(RTL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    RTLPS[i] <- PS3$phy.signal
    RTLPSp[i] <- PS3$pvalue
   
    ### 4. relative hind foot length
    RHFL1 <- subFuncdat$RHFL
    names(RHFL1)=rownames(subFuncdat)
    #geomorph K
    PS4 <- physignal(RHFL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    RHFLPS[i] <- PS4$phy.signal
    RHFLPSp[i] <- PS4$pvalue

    ### 5.relative ear length
    REL1 <- subFuncdat$REL
    names(REL1)=rownames(subFuncdat)
    #geomorph K
    PS5 <- physignal(REL1, subTree, iter = 999, seed = NULL, print.progress = F) 
    RELPS[i] <- PS5$phy.signal
    RELPSp[i] <- PS5$pvalue
    
    
    ####Multivariate discrete traits
    # $r.Mantel = The Mantel observed statistic.
    # $p.NULL = The p value from no phylogenetic structure (standard p value in Mantel test).
    # $p.BM = The p value under simulation of traits from Brownian phylogenetic structure.
    
    ### 6.diet
    Diet <- select(subFuncdat, Diet.Inv:Diet.PlantO)
    PS6 <- EM.mantel(subTree, Diet, runs = 999)
    DietPS[i] <- PS6$r.Mantel
    DietPSp[i] <- PS6$p.BM
  
    ### 7.habitats
    Habitat <- select(subFuncdat, Forest:Desert)
    PS7 <- EM.mantel(subTree, Habitat, runs = 999)
    HabitatPS[i] <- PS7$r.Mantel
    HabitatPSp[i] <- PS7$p.BM
    
    ### 8. diel activity
    Activity <- select(subFuncdat, Activity.Nocturnal:Activity.Diurnal)
    PS8 <- EM.mantel(subTree, Activity, runs = 999)
    ActivityPS[i] <- PS8$r.Mantel
    ActivityPSp[i] <- PS8$p.BM
    
    ####single Categorical traits (not binary)
    ## 9. life mode (vertical stratum)
    LifeMode1 <- select(subFuncdat, Life.mode)
    PS9 <- EM.mantel(subTree, LifeMode1, runs = 999)
    LifeModePS[i] <- PS9$r.Mantel
    LifeModePSp[i] <- PS9$p.BM
    
    ## 10. saltatory locomotion
    Salt1 <- select(subFuncdat, Saltation)
    PS10 <- EM.mantel(subTree, Salt1, runs = 999)
    SaltPS[i] <- PS10$r.Mantel
    SaltPSp[i] <- PS10$p.BM
 
    ### 11. is binary cheek pouches - done separately.
    
} #closing of tree bracket




#########################################################################

####
##summarize outputs
######

#1. Mass
MassPSid <- as.data.frame(MassPS) #phylosig values
MassPSpid <- as.data.frame(MassPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
MassPSid$MassPS_mean <- sum(MassPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
MassPSpid$Massp_count <- rowSums(MassPSpid[,1:100] < 0.05)
head(MassPSid[,99:101])
head(MassPSpid[,99:101])
MassPS_summary <- as.data.frame(cbind(MassPSid[101], MassPSpid[101]))
write.csv(MassPSid, file = "MassPS_global.csv", row.names = F)
write.csv(MassPSpid, file = "MassPSp_global.csv", row.names = F)

#2. HBL
HBLPSid <- as.data.frame(HBLPS) #phylosig values
HBLPSpid <- as.data.frame(HBLPSp)  #phylosig signifs
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
HBLPSid$HBLPS_mean <- sum(HBLPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
HBLPSpid$HBLp_count <- rowSums(HBLPSpid[,1:100] < 0.05)

head(HBLPSid[,99:101])
head(HBLPSpid[,99:101])
HBLPS_summary <- as.data.frame(cbind(HBLPSid[101], HBLPSpid[101]))
write.csv(HBLPSid, file = "HBLPS_global.csv", row.names = F)
write.csv(HBLPSpid, file = "HBLPSp)global.csv", row.names = F)

#3. RTL
RTLPSid <- as.data.frame(RTLPS) #phylosig values
RTLPSpid <- as.data.frame(RTLPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
RTLPSid$RTLPS_mean <- sum(RTLPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
RTLPSpid$RTLp_count <- rowSums(RTLPSpid[,1:100] < 0.05)
head(RTLPSid[,99:101])
head(RTLPSpid[,99:101])
RTLPS_summary <- as.data.frame(cbind(RTLPSid[101], RTLPSpid[101]))
write.csv(RTLPSid, file = "RTLPS_global.csv", row.names = F)
write.csv(RTLPSpid, file = "RTLPSp_global.csv", row.names = F)

#4. RHFL
RHFLPSid <- as.data.frame(RHFLPS) #phylosig values
RHFLPSpid <- as.data.frame(RHFLPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
RHFLPSid$RHFLPS_mean <- sum(RHFLPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
RHFLPSpid$RHFLp_count <- rowSums(RHFLPSpid[,1:100] < 0.05)
head(RHFLPSid[,99:101])
head(RHFLPSpid[,99:101])
RHFLPS_summary <- as.data.frame(cbind(RHFLPSid[101], RHFLPSpid[101]))
write.csv(RHFLPSid, file = "RHFLPS_global.csv", row.names = F)
write.csv(RHFLPSpid, file = "RHFLPSp_global.csv", row.names = F)

#5. REL
RELPSid <- as.data.frame(RELPS) #phylosig values
RELPSpid <- as.data.frame(RELPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
RELPSid$RELPS_mean <- sum(RELPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
RELPSpid$RELp_count <- rowSums(RELPSpid[,1:100] < 0.05)
head(RELPSid[,99:101])
head(RELPSpid[,99:101])
RELPS_summary <- as.data.frame(cbind(RELPSid[101], RELPSpid[101]))
write.csv(RELPSid, file = "RELPS_global.csv", row.names = F)
write.csv(RELPSpid, file = "RELPSp_global.csv", row.names = F)

#6. Diet
DietPSid <- as.data.frame(DietPS) #phylosig values
DietPSpid <- as.data.frame(DietPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
DietPSid$DietPS_mean <- sum(DietPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
DietPSpid$Dietp_count <- rowSums(DietPSpid[,1:100] < 0.05)
head(DietPSid[,99:101])
head(DietPSpid[,99:101])
DietPS_summary <- as.data.frame(cbind(DietPSid[101], DietPSpid[101]))
write.csv(DietPSid, file = "DietPS_global.csv", row.names = F)
write.csv(DietPSpid, file = "DietPSp_global.csv", row.names = F)

#7. Habitat
HabitatPSid <- as.data.frame(HabitatPS) #phylosig values
HabitatPSpid <- as.data.frame(HabitatPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
HabitatPSid$HabitatPS_mean <- sum(HabitatPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
HabitatPSpid$Habitatp_count <- rowSums(HabitatPSpid[,1:100] < 0.05)
head(HabitatPSid[,99:101])
head(HabitatPSpid[,99:101])
HabitatPS_summary <- as.data.frame(cbind(HabitatPSid[101], HabitatPSpid[101]))
write.csv(HabitatPSid, file = "HabitatPS_global.csv", row.names = F)
write.csv(HabitatPSpid, file = "HabitatPSp_global.csv", row.names = F)

#8. Activity
ActivityPSid <- as.data.frame(ActivityPS) #phylosig values
ActivityPSpid <- as.data.frame(ActivityPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
ActivityPSid$ActivityPS_mean <- sum(ActivityPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
ActivityPSpid$Activityp_count <- rowSums(ActivityPSpid[,1:100] < 0.05)
head(ActivityPSid[,99:101])
head(ActivityPSpid[,99:101])
ActivityPS_summary <- as.data.frame(cbind(ActivityPSid[101], ActivityPSpid[101]))
write.csv(ActivityPSid, file = "ActivityPS_global.csv", row.names = F)
write.csv(ActivityPSpid, file = "ActivityPSp_global.csv", row.names = F)

#9. LifeMode
LifeModePSid <- as.data.frame(LifeModePS) #phylosig values
LifeModePSpid <- as.data.frame(LifeModePSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
LifeModePSid$LifeModePS_mean <- sum(LifeModePSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
LifeModePSpid$LifeModep_count <- rowSums(LifeModePSpid[,1:100] < 0.05)
head(LifeModePSid[,99:101])
head(LifeModePSpid[,99:101])
LifeModePS_summary <- as.data.frame(cbind(LifeModePSid[101], LifeModePSpid[101]))
write.csv(LifeModePSid, file = "LifeModePS_global.csv", row.names = F)
write.csv(LifeModePSpid, file = "LifeModePSp_global.csv", row.names = F)

#10. Salt
SaltPSid <- as.data.frame(SaltPS) #phylosig values
SaltPSpid <- as.data.frame(SaltPSp)  #phylosig signifs

#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
SaltPSid$SaltPS_mean <- sum(SaltPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
SaltPSpid$Saltp_count <- rowSums(SaltPSpid[,1:100] < 0.05)
head(SaltPSid[,99:101])
head(SaltPSpid[,99:101])
SaltPS_summary <- as.data.frame(cbind(SaltPSid[101], SaltPSpid[101]))
write.csv(SaltPSid, file = "SaltPS_global.csv", row.names = F)
write.csv(SaltPSpid, file = "SaltPSp_global.csv", row.names = F)







#read back in mass summaries
MassPSid <- read.csv("MassPS_global.csv", header = F)
MassPSpid <- read.csv("MassPSp_global.csv", header = F)
MassPS_summary <- as.data.frame(cbind(MassPSid[,101], MassPSpid[,101]))





###########################################
######################################################################
# 
###binary traits
require(caper)  #masks select from dplyr
#caper - phylo.d for binary traits (Fritz and Purvis) The D statistic is equal to 1 if the observed binary trait has a phylogenetically random distribution across the tips of the phylogeny and to 0 if the observed trait is as clumped as if it had evolved by Brownian motion under our threshold model.
####USE for CheekPouch


########


PouchPS <- matrix(NA,1,100)
PouchPSp <- matrix(NA,1,100)
PouchPSpB <- matrix(NA,1,100)

#cycle through each of 100 phylogenetic trees to calculate phy signal for each trait
for (i in 1:100){	
  subTree <- phylotree[[i]]
  
  # reorder the trait data matrix to match the tip labels by species name
  subTraits2 <- Traits.PhyloNames[match(subTree$tip.label, Traits.PhyloNames$Phylo_Name),]      
  row.names(subTraits2) <- subTraits2$Phylo_Name #set row names
  subFuncdat <- subTraits2[,5:ncol(subTraits2)]   #strip down to just traits
  
    #DFs for binary trait - cheeck pouches
    Pouch1 <- as.data.frame(subFuncdat$CheekPouch)
    Pouch1$names <- rownames(subFuncdat)
    colnames(Pouch1) <- c("pouches", "names")
    #Fritz and Purvis D
    PS11 <- phylo.d(Pouch1, subTree, names.col = names, binvar = pouches, permut = 1000)
    PouchPS[i] <- PS11$DEstimate
    PouchPSp[i] <- PS11$Pval1
    PouchPSpB[i] <- PS11$Pval0
    #$DEstimate for Mantel test value
    #$Pval1 and Pval0 for probability of resulting from random str (ie. signif phylo signal) and brownian str.
    
  }

####
##summarize output
######

#11. Pouch
PouchPSid <- as.data.frame(PouchPS) #phylosig values
PouchPSpid <- as.data.frame(PouchPSp)  #phylosig signifs
PouchPSpBid <- as.data.frame(PouchPSpB)  #phylosig signifs
#calculate mean phylosig values across the 100 trees randomly selected from the 10000 in Upham.
PouchPSid$PouchPS_mean <- sum(PouchPSid[,1:100])/100
#tally how many trees had significant phylosig for the trait
PouchPSpid$Pouchp_count <- rowSums(PouchPSpid[,1:100] < 0.05)
PouchPSpBid$PouchpB_count <- rowSums(PouchPSpBid[,1:100] < 0.05)

head(PouchPSid[,99:101])
head(PouchPSpid[,99:101])
head(PouchPSpBid[,99:101])

PouchPS_summary <- as.data.frame(cbind(PouchPSid[101], PouchPSpid[101], PouchPSpBid[101]))
write.csv(PouchPSid, file = "PouchPS_global.csv", row.names = F)
write.csv(PouchPSpid, file = "PouchPSp_global.csv", row.names = F)
write.csv(PouchPSpBid, file = "PouchPSpBrownian_global.csv", row.names = F)






#######################################
###summarize ALL phylogenetic signal test outputs
MtnPhySig_glob1 <- as.data.frame(cbind(MassPS_summary, HBLPS_summary, RTLPS_summary, RHFLPS_summary, RELPS_summary, DietPS_summary, HabitatPS_summary, ActivityPS_summary, LifeModePS_summary, SaltPS_summary, PouchPS_summary))
MtnPhySig_glob1
write.csv(MtnPhySig_glob1, file = "IndivMtnPhySigSummary_global.csv", row.names = F)

