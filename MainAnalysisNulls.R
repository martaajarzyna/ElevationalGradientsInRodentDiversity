###Rodent diversity patterns along elevational gradients worldwide - Null models
##Kohli, Miyajima and Jarzyna, GEB
##last modified by BAK: 2 April 2021
#####

###Note: All null iterations were run via computing cluster using separate, single scripts due to their excessive run time, except for hypervolumes, but are included here together.


require(data.table)  
require(gawdis)      
require(FD)          
require(ape)         
require(phytools)    
require(picante)      
require(hypervolume) 
require(dplyr)        
require(usedist)      
require(vegan)        
require(moments)      

###################################################
########STEP 2: Functional & phylogenetic diversity null modelling
###################################################

#### OCCURRENCE DATA PREP

#Read in data on species distributions
SpOcc <- read.csv("SpOcc_MtnElevBins.csv", header = T)


SpOcc$MtnNum <- as.factor(SpOcc$Mtn) #create new column to hold factors of mtn IDs
SpOcc$spId <- as.numeric(as.factor(SpOcc$Species)) #create new column to hold numeric species IDs
SpOcc$spphyId <- as.numeric(as.factor(SpOcc$Phylo_Name)) #create new column to hold numeric species IDs

colnames(SpOcc) <- c("poId", "MtnName", "elevation", "SR", "spp", "sppPhy", "moId", "spId", "spphyId")
head(SpOcc)
length(unique(SpOcc$poId))###poIds are a unique code for each mtn/elev bin combo. moIDs are Mountain IDs, elevation is the minimum of each 100m bin, SR is species richness of a bin (from interpolated elevational ranges of species), spp and sppPhy are species names to match the trait and phylogenetic dataset, respectively.

###This is structured as a presence-only dataset.

###FOR SOME INDICES, WE REMOVE/IGNORE THE BINS WITH LESS THAN 3 SPECIES. - as such, need to make it uniform across analyses
SpOcc <- SpOcc[SpOcc$SR >= 3,] # Starting number of rows = 9488, starting number of poids = 1032. 
### SR >=3 gives 9400 rows and 979 poids. 

length(unique(SpOcc$poId))
length(unique(SpOcc$spp))  #retains all species 
head(SpOcc) #call all species occurrences

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
funcdat <- select(sp.traits, CheekPouch, Diet.Inv, Diet.Verts, Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO, Activity.Nocturnal, Activity.Crepuscular, Activity.Diurnal, Body_mass_g, Forest, Savanna, Shrubland, Grassland, Wetland, Rocky.areas, Desert, HBL, RTL, RHFL, REL, Life.mode, Saltation) #Note that gawdis methods adjust for range of individual traits and cross-correlation

#log transform Mass, and HBL to make more normal
##write over existing columns to keep same names
funcdat$Body_mass_g <- log10(funcdat$Body_mass_g)
funcdat$HBL <- log10(funcdat$HBL)

#set characters as factors
funcdat$Life.mode <- as.factor(funcdat$Life.mode)
funcdat$Saltation <- as.factor(funcdat$Saltation)

###############################################
#read in gawdis distance matrix
funcdist <- readRDS("gawdis_distmat.rds") 

############### 

####. Phylogenetic tree
phylotree <- read.nexus("Phylogenies.nex") #downloaded from vertlife.org
####when used, need to match to the phylonames instead of trait names: sppPhy

#####
#####Skew function
###functions and code for summarizing the skew of null randomizations to inform the reliability of SES values (see Botta-Dukat 2018 for reasoning and explanation)

#Created by: BA Kohli 
#on: 29 October 2020
#last modified: 3 March 2020

#calculates skew value of null distributions excluding first column, which is the observed values.
skew.func <- function(z){
  mean.skew <- mean(apply(z[,2:(ncol(z))], MARGIN = 1, skewness), na.rm = T)
  sd.skew <- sd(apply(z[,2:(ncol(z))], MARGIN = 1, skewness), na.rm = T)
  skew.out <- cbind.data.frame(mean.skew, sd.skew)
}
########################

###################################################
#######Constructing Null Assemblages 
###################################################


#### Randomize via Independent swap null (holds row and column totals constant; richness and species occurrence frequencies across the gradient)
######

head(SpOcc)
SpOcc.species <- select(SpOcc, spp, sppPhy, spId, spphyId)
SpOcc.speciesU <- distinct(SpOcc.species)
SpOcc.speciesU <- SpOcc.speciesU %>% mutate_if(is.numeric, as.factor)
moIDs <- unique(SpOcc$moId)
poIDs <- unique(SpOcc$poId)

for (i in 1:1000){
  comm.null.all <- list()
  
  for (k in 1:length(moIDs)) {
    comm.null.mo <- list()
    comm.mo <- subset(SpOcc, moId == moIDs[k]) #isolate all data for a given mountain
    #make a sample matrix that mpd can read (with species names and any numeric value to indicate presence)
    comm.mo.samp <- select(comm.mo, poId, SR, spId) #comm name, sp abundance, sp name
     #subset the data as needed to convert
    comm.mo.mat <- sample2matrix(comm.mo.samp) #convert to matrix format
    comm.mo.mat <- decostand(comm.mo.mat, method = "pa") #make binary (still the empirical matrix)
    comm.mo.mat.rdm <- randomizeMatrix(comm.mo.mat, null.model = "independentswap", iterations = 1000) #creates randomized occurrence matrix
     comm.mo.samp.rdm <- matrix2sample(comm.mo.mat.rdm)  #take the new randomized matrix to sample format
   colnames(comm.mo.samp.rdm) <- c("poId.null", "pres.null", "spId.null")
   
   #order the empirical and randomized data exactly the same and bind them together
 comm.mo.samp.rdm.order <-comm.mo.samp.rdm[order(comm.mo.samp.rdm$poId.null),]
   comm.mo.order <- comm.mo[order(comm.mo$poId),]
   comm.null.mo1 <- cbind(comm.mo.samp.rdm.order, comm.mo.order)
   #select and order column to match SpOcc
   comm.null.mo1 <- select(comm.null.mo1, poId, MtnName, elevation, SR, moId, spId.null) 
   #join with species info
   comm.null.mo.k <- inner_join(comm.null.mo1, SpOcc.speciesU, by = c("spId.null" = "spId"))
   comm.null.mo = comm.null.mo.k %>% select(1:4, 7:8, moId, spId.null, spphyId)

    comm.null.all <- rbind(comm.null.all, comm.null.mo)  #this step combines all the data for all the moids in the dataset = a complete random set of observations per bin per mountain.  It should have the same number of observations as the empiricial dataset (spocc).
  }
  saveRDS(comm.null.all, file=paste0("null_mountsp_IndSwap_",i,".rds"))
}

#look at an example
n1 <- readRDS("null_mountsp_IndSwap_1.rds")
#################################################
#################################################



##########################################################




#################################
###calculate nulls per index
############################

poIDs <- unique(SpOcc$poId)



######
#### A. #### Calculate null Functional richness - convex hull volume
######

## Principal Coordinate Analysis
pcoa.calc <- pcoa(funcdist, correction="lingoes")
nbdim <- ncol(pcoa.calc$vectors)
eigen <- pcoa.calc$values
vect <- pcoa.calc$vectors[,1:nbdim]
rownames(vect) <- rownames(funcdat)
colnames(vect) <- paste("PC",1:nbdim,sep="")

FDHnull <- matrix(NA,length(poIDs),1000)

#this takes a long time to complete#
for (i in 1:1000){
  n1 <- readRDS(file=paste0("null_mountsp_IndSwap_",i,".rds"))
  
  for (k in 1:length(poIDs)){
    comm <- n1[n1$poId ==  poIDs[k],]
    sp.i <- unique(comm$spp)
    vect.i <- vect[rownames(vect) %in% sp.i,1:2]
    fdh.i <- expectation_convex(vect.i, check.memory=F, verbose = F)
   
    FDHnull[k,i] <- fdh.i@Volume # Volume of the convex hull
   
  }
  saveRDS(FDHnull, file=("FDHnull_IndSwap1000_unproc.rds")) #saves the unprocessed matrix holding all null volume values
}



#
#
#
FDHnullid <- as.data.frame(cbind(FDHnull,poIDs))
FDHnullid$moIDs <- 0
FDHnullid$elevation <- 0
FDHnullid$SR <- 0

for (k in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[k])
  FDHnullid[k,1002:1004] <- comm[1,2:4]  
}

FDHnullid <- cbind(FDHnullid[1001:1004],FDHnullid[1:1000]) 
head(FDHnullid)
colSums(FDHnullid[5:1004])
colMeans(FDHnullid[5:1004])
rowSums(FDHnullid[5:1004])
range(rowMeans(FDHnullid[5:1004]))

saveRDS(FDHnullid, file="FDHnull_IndSwap1000.rds")

#### compare observed empirical value to null distribution of index 
FDH <- readRDS(file="FDH_convexhull_2d.rds")
FDHnull <- readRDS(file="FDHnull_IndSwap1000.rds")

#check order
FDH$poIDs == FDHnull$poIDs  #ensure they match

FDHnull.n <- cbind(FDHnull[,1:4],FDH[,1],FDHnull[,5:ncol(FDHnull)]) #the first column of values is the observed value, the other columns are random iterations' values, ##rows are for each site (bin).

head(FDHnull.n[1:10])
head(FDH)

## Rank observed value, calculate SES and p value
obs.null.output <- as.matrix(FDHnull.n[,5:(ncol(FDHnull.n))])
#quantify the rank or the observed value in the null distribution for each community
FDH$FDH.obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
#values range from 1 to the the number of iterations run +1

#quantify the Standardized Effect Size;  SES = (observed - null mean) / null sd
FDH$SES.FDH <- (obs.null.output[,1] - apply(obs.null.output[,2:(ncol(obs.null.output))], 1, mean)) / apply(obs.null.output[,2:(ncol(obs.null.output))], 1, sd)

#calculate p values 
FDH$FDH.p <- apply(cbind(obs.null.output[,1], obs.null.output[,2:(ncol(obs.null.output))]), MARGIN = 1, rank)[1,] / 1001 

head(FDH)
summary(FDH)
hist(FDH$FDH.vol)
hist(FDH$FDH.obs.rank)
hist(FDH$SES.FDH)

#########################investigating NA's
FDH[FDH$SES.FDH == "NaN",] #print any NAs:  Luoji 2600m, JiajinN 3600m, Oquirrh 2500m  #all the values are same across randomizations, thus divides by sd=0.  this is because the bin either contains EVERY species found on the mountain, or because the bin contains only species that are found IN EVERY bin, and thus composition is invariant across randomizations.
################################

#change SES NAs to 0
FDH$SES.FDH[is.na(FDH$SES.FDH)]<-0
summary(FDH)

#save out SES results bound to observed results
saveRDS(FDH, file="FDHnull_IndSwap_convexhull_2d.rds")
##################################################################################

############
####checking Skew

FDHnull_IndSwap_SKEW <- skew.func(obs.null.output)  #generates an overall skew summary for the whole index nulls across all poids

#or append it to the main summary file
FDH$FDH.skew <- apply(obs.null.output[,2:(ncol(obs.null.output))], MARGIN = 1, skewness)
head(FDH)
summary(FDH)
FDH[FDH$FDH.skew== "NaN",] 
FDH$FDH.skew[is.na(FDH$FDH.skew)]<-0
hist(FDH$FDH.skew)
quantile(FDH$FDH.skew, probs = c(.01, .025, .05, .95, .975, .99), na.rm = T)

#save out SES results bound to observed results
saveRDS(FDH, file="FDHnull_IndSwap1000_convexhull_2d.rds")
##########################################################################################






#########################################################################################
######
#### B. #### Calculate null Functional dispersion and evenness (dbFD)
######

# Quantify non-abundance-weighted functional dispersion (FDis), functional evenness (FEve), and other indices via dbFD package for each point location (=elev bin)
poIDs <- unique(SpOcc$poId)

DBFD_fdis <- matrix(NA,length(poIDs),1000)
DBFD_feve <- matrix(NA,length(poIDs),1000)

for (n in 1:1000){
  n1 <- readRDS(file=paste0("null_mountsp_IndSwap_",n,".rds"))
  
  for (i in 1:length(poIDs)) {
    comm <- subset(n1, poId == poIDs[i])
    nspec <- nrow(comm)
    #comm <- setDF(comm)
    sp.i <- unique(comm$spp)
    funcdist.i <- dist_subset(funcdist, sp.i)  #package usedist extracts comm sp list to a subset dist matrix (doesn't work for N=1 bc doesn't have labels)  
    dbfd.i <- dbFD(x = funcdist.i, w.abun = F, corr = "lingoes", calc.FRic = F, calc.FGR = F, calc.CWM = F, calc.FDiv = F, print.pco = F) ##set corr = "lingoes" to make equivalent across all indices/packages used 
    
    saveRDS(dbfd.i, file=paste0("dbfd_2D_",poIDs[i],".rds"))
    
    DBFD_fdis[i, n] <- dbfd.i$FDis
    DBFD_feve[i, n] <- dbfd.i$FEve
    #DBFD[i, 3] <- dbfd.i$RaoQ
    #DBFD[i, 4] <- dbfd.i$nbsp
    #DBFD[i, 5] <- dbfd.i$sing.sp
  }
}



##################################################

FDISnullid <- as.data.frame(cbind(DBFD_fdis,poIDs))
#FDISnullid$longitude <- 0
#FDISnullid$latitude <- 0
FDISnullid$moIDs <- 0
FDISnullid$elevation <- 0
FDISnullid$SR <- 0

for (i in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[i])
  FDISnullid[i,1002:1004] <- comm[1,2:4]  
}

FDISnullid <- cbind(FDISnullid[1001:1004], FDISnullid[1:1000]) 
head(FDISnullid)
saveRDS(FDISnullid, file="DBFD_2D_FDISnull_IndSwap.rds")

#   #   #   #

FEVEnullid <- as.data.frame(cbind(DBFD_feve,poIDs))
#FEVEnullid$longitude <- 0
#FEVEnullid$latitude <- 0
FEVEnullid$moIDs <- 0
FEVEnullid$elevation <- 0
FEVEnullid$SR <- 0

for (i in 1:length(poIDs)){
  comm <- subset(SpOcc, poId == poIDs[i])
  FEVEnullid[i,1002:1004] <- comm[1,2:4]
}

FEVEnullid <- cbind(FEVEnullid[1001:1004], FEVEnullid[1:1000]) 
head(FEVEnullid)
saveRDS(FEVEnullid, file="DBFD_2D_FEVEnull_IndSwap.rds")

################

#### compare observed empirical value to null distribution of index 
#FDis
dbfd.obs <- readRDS(file="DBFD_2D.rds")
FDISnull <- readRDS(file="DBFD_2D_FDISnull_IndSwap.rds")
###check order
table(dbfd.obs$poIDs == FDISnull$poIDs)  #

FDISnull.n <- cbind(FDISnull[,1:4],dbfd.obs[,1],FDISnull[,5:ncol(FDISnull)]) #the first column of values is the observed value, the other columns are random iterations' values, ##rows are for each site (bin).

head(FDISnull.n)
head(dbfd.obs)

## Rank observed value, calculate SES and p value
obs.null.output <- as.matrix(FDISnull.n[,5:(ncol(FDISnull.n))])
#quantify the rank or the observed value in the null distribution for each community
dbfd.obs$FDis.obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
#values range from 1 to the the number of iterations run +1

#quantify the Standardized Effect Size;  SES = (observed - null mean) / null sd
dbfd.obs$SES.FDis <- (obs.null.output[,1] - apply(obs.null.output[,2:(ncol(obs.null.output))], 1, mean)) / apply(obs.null.output[,2:(ncol(obs.null.output))], 1, sd)

#calculate p values 
dbfd.obs$FDis.p <- apply(cbind(obs.null.output[,1], obs.null.output[,2:(ncol(obs.null.output))]), MARGIN = 1, rank)[1,] / 1001 

head(dbfd.obs)
summary(dbfd.obs)
hist(dbfd.obs$FDis)
hist(dbfd.obs$FDis.obs.rank)
hist(dbfd.obs$SES.FDis)

#########################investigating NA's
dbfd.obs[dbfd.obs$SES.FDis == "NaN",] #print any NAs:  Luoji 2600m, JiajinN 3600m, Oquirrh 2500m  #all the values are same across randomizations, thus divides by sd=0.  this is because the bin either contains EVERY species found on the mountain, or because the bin contains only species that are found IN EVERY bin, and thus composition is invariant across randomizations.
################################

#change SES NAs to 0
dbfd.obs$SES.FDis[is.na(dbfd.obs$SES.FDis)]<-0
summary(dbfd.obs)
######

####checking Skew

null_IndSwap_SKEW <- skew.func(obs.null.output)  #generates an overall skew summary for the whole index nulls across all poids

#or append it to the main summary file
dbfd.obs$FDis.skew <- apply(obs.null.output[,2:(ncol(obs.null.output))], MARGIN = 1, skewness)
head(dbfd.obs)
summary(dbfd.obs)
dbfd.obs[dbfd.obs$FDis.skew== "NaN",] 
hist(dbfd.obs$FDis.skew)
quantile(dbfd.obs$FDis.skew, probs = c(.01, .025, .05, .95, .975, .99), na.rm = T)

###
####FEve
###

#observed are in dbfd.obs, read in above
FEVEnull <- readRDS(file="DBFD_2D_FEVEnull_IndSwap.rds")
FEVEnull.n <- cbind(FEVEnull[,1:4],dbfd.obs[,2],FEVEnull[,5:ncol(FEVEnull)]) #the first column of values is the observed value, the other columns are random iterations' values, ##rows are for each site (bin).
table(dbfd.obs$poIDs == FEVEnull$poIDs) 

head(FEVEnull.n)
head(dbfd.obs)

## Rank observed value, calculate SES and p value
obs.null.output <- as.matrix(FEVEnull.n[,5:(ncol(FEVEnull.n))])
#quantify the rank or the observed value in the null distribution for each community
dbfd.obs$FEve.obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
#values range from 1 to the the number of iterations run +1

#quantify the Standardized Effect Size;  SES = (observed - null mean) / null sd
dbfd.obs$SES.FEve <- (obs.null.output[,1] - apply(obs.null.output[,2:(ncol(obs.null.output))], 1, mean)) / apply(obs.null.output[,2:(ncol(obs.null.output))], 1, sd)

#calculate p values 
dbfd.obs$FEve.p <- apply(cbind(obs.null.output[,1], obs.null.output[,2:(ncol(obs.null.output))]), MARGIN = 1, rank)[1,] / 1001 

head(dbfd.obs)
summary(dbfd.obs)
hist(dbfd.obs$FEve)
hist(dbfd.obs$FEve.obs.rank)
hist(dbfd.obs$SES.FEve)

#########################investigating NA's
dbfd.obs[dbfd.obs$SES.FEve == "NaN",] #print any NAs:  Luoji 2600m, JiajinN 3600m, Oquirrh 2500m  #all the values are same across randomizations, thus divides by sd=0.  this is because the bin either contains EVERY species found on the mountain, or because the bin contains only species that are found IN EVERY bin, and thus composition is invariant across randomizations.
################################

#change SES NAs to 0
dbfd.obs$SES.FEve[is.na(dbfd.obs$SES.FEve)]<-0
summary(dbfd.obs)


null_IndSwap_SKEW <- skew.func(obs.null.output)  #generates an overall skew summary for the whole index nulls across all poids

#or append it to the main summary file
dbfd.obs$FEve.skew <- apply(obs.null.output[,2:(ncol(obs.null.output))], MARGIN = 1, skewness)
head(dbfd.obs)
summary(dbfd.obs)
dbfd.obs[dbfd.obs$FEve.skew== "NaN",] 
hist(dbfd.obs$FEve.skew)
quantile(dbfd.obs$FEve.skew, probs = c(.01, .025, .05, .95, .975, .99), na.rm = T)

#save out both FDis and FEve SES results bound to observed results
saveRDS(dbfd.obs, file="DBFD_2Dnull_IndSwap1000.rds")

##################################################################################





#################################################################################
######
#### c. #### Calculate null phylogenetic diversity (PD Faith)  - PD based on branch lengths and topology directly
#

#############################

for (j in 1:1000){ 
  
  cat(j)
  n1 <- readRDS(file=paste0("null_mountsp_IndSwap_",j,".rds")) #these are randomizations created before
  
  poIDs <- unique(n1$poId)
  PD <- matrix(NA,length(poIDs),100) #to hold PD from across 100 trees
  
  for (i in 1:100){	
    pdi <- phylotree[[i]]
    for (k in 1:nrow(PD)) {
      comm <- subset(n1, poId == poIDs[k])
      nspec <- nrow(comm)
      #comm <- setDF(comm)
      sp.id <- pdi$tip.label
      
      rmTip <- as.character(sp.id[! sp.id %in% comm$sppPhy])#removing missing tips 
      subTree <- drop.tip(pdi, tip=rmTip)  # remove those branches
      PD[k,i] <- sum(subTree$edge.length) #summing up the branch length
    }
  }
  
  #create ids and insert other info
  PDid <- cbind(PD,poIDs)
  PDid <- as.data.frame(PDid)
  PDid$moIDs <- 0
  PDid$elevation <- 0
  PDid$SR <- 0
  
  for (l in 1:nrow(PD)){
    comm <- subset(n1, poId == poIDs[l])
    PDid[l,102:104] <- comm[1,2:4] #dimensions of PDid to insert into will change depending on how many trees one uses
  }
  
  saveRDS(PDid, file=paste0("null_PD_IndSwap",j,".rds")) #this saves null values for each species randomization, you'll have as many files as randomizations
}
##################################################


############
##Calculate means across trees. read in each iteration (j) rds, and calculate the mean index values across the rows, bind all those together into a streamlined null values file (1000 cols)
PDmean <- matrix(NA,length(poIDs),1000)

for (z in 1:1000){ 
  
  z1 <- readRDS(file=paste0("null_PD_IndSwap",z,".rds"))
  
  PDmean[,z] <- rowMeans(z1[,1:100]) #calculate mean values across the N trees sampled
  PDmeanFinal <- cbind(PDid[101:104],PDmean)
  saveRDS(PDmeanFinal, file="PDmeanNull_IndSwap.rds")
}


###
####PD SES
###

#### compare observed empirical value to null distribution of index 
PD.obs <- readRDS(file="PD.rds")
PDnull <- readRDS(file="PDmeanNull_IndSwap.rds")
###check order
table(PD.obs$poIDs == PDnull$poIDs) #NOT THE SAME ORDER!
head(PD.obs[,101:105])
head(PDnull[,1:5])

#to unify the rows in the order of the observed data use left join, then select those columns desired, in desired order
PDObsNullJoin <- right_join(PD.obs[,101:105], PDnull2, by = "poIDs")
head(PDObsNullJoin)
#check what they match
unique(PD.obs$poIDs == PDObsNullJoin$poIDs)
unique(PDObsNullJoin$elevation.x == PDObsNullJoin$elevation.y)

PDnull.n <- cbind(PD.obs[,101:105],PDObsNullJoin[,9:ncol(PDObsNullJoin)]) #the first column of values is the observed value, the other columns are random iterations' values, ##rows are for each site (bin).

head(PDnull.n[,1:5])
head(PD.obs[,101:105])

## Rank observed value, calculate SES and p value
obs.null.output <- as.matrix(PDnull.n[,5:(ncol(PDnull.n))])
#quantify the rank or the observed value in the null distribution for each community
PD.obs$PD.obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
#values range from 1 to the the number of iterations run +1

#quantify the Standardized Effect Size;  SES = (observed - null mean) / null sd
PD.obs$SES.PD <- (obs.null.output[,1] - apply(obs.null.output[,2:(ncol(obs.null.output))], 1, mean)) / apply(obs.null.output[,2:(ncol(obs.null.output))], 1, sd)

#calculate p values 
PD.obs$PD.p <- apply(cbind(obs.null.output[,1], obs.null.output[,2:(ncol(obs.null.output))]), MARGIN = 1, rank)[1,] / 1001 

head(PD.obs)
summary(PD.obs)
hist(PD.obs$PD_mean)
hist(PD.obs$PD.obs.rank)
hist(PD.obs$SES.PD)

#########################investigating NA's
PD.obs[PD.obs$SES.PD == "NaN",] #print any NAs:  Luoji 2600m, JiajinN 3600m, Oquirrh 2500m  #all the values are same across randomizations, thus divides by sd=0.  this is because the bin either contains EVERY species found on the mountain, or because the bin contains only species that are found IN EVERY bin, and thus composition is invariant across randomizations.
################################

#change SES NAs to 0
PD.obs$SES.PD[is.na(PD.obs$SES.PD)]<-0
summary(PD.obs)


null_IndSwap_SKEW <- skew.func(obs.null.output)  #generates an overall skew summary for the whole index nulls across all poids

#or append it to the main summary file
PD.obs$PD.skew <- apply(obs.null.output[,2:(ncol(obs.null.output))], MARGIN = 1, skewness)
head(PD.obs)
summary(PD.obs)
PD.obs[PD.obs$PD.skew== "NaN",] 
PD.obs$PD.skew[is.na(PD.obs$PD.skew)]<-0
hist(PD.obs$PD.skew)
quantile(PD.obs$PD.skew, probs = c(.01, .025, .05, .95, .975, .99), na.rm = T)

#save out PD SES results bound to observed results
saveRDS(PD.obs, file="PDnull_IndSwap1000.rds")

##########################################################################################





#################################################################################
######
#### D. #### Calculate null phylogenetic evenness and dispersion via distances
######


for (j in 1:1000){ 
  cat(j)
  n1 <- readRDS(file=paste0("null_mountsp_IndSwap_",j,".rds")) #these are randomizations created before
  
  poIDs <- unique(n1$poId)
  PEve <- matrix(NA,length(poIDs),1000) # Phylogenetic evenness
  MPD <- matrix(NA,length(poIDs),1000) # Mean pairwise phylogenetic distance
  
  for (i in 1:100){	
    pdi <- phylotree[[i]]
    phydist <- cophenetic.phylo(pdi) #global phylogenetic distance matrix
    
    for (k in 1:nrow(PEve)) {
      comm <- subset(n1, poId == poIDs[k])
      nspec <- nrow(comm)
      #comm <- setDF(comm)
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
  
  #create ids and insert other info
  PEveid <- as.data.frame(cbind(PEve,poIDs)) #PEve
  PEveid$moIDs <- 0
  PEveid$elevation <- 0
  PEveid$SR <- 0
  
  MPDid <- as.data.frame(cbind(MPD,poIDs))  #mean phylogenetic distance
  MPDid$moIDs <- 0
  MPDid$elevation <- 0
  MPDid$SR <- 0
  
  for (l in 1:nrow(PEve)){
    comm <- subset(n1, poId == poIDs[l])
    PEveid[l,1002:1004] <- comm[1,2:4] #dimensions of PDid to insert into will change depending on how many trees one uses
    MPDid[l,1002:1004] <- comm[1,2:4]
  }
  
  saveRDS(PEveid, file=paste0("null_PEve_IndSwap",j,".rds")) #this saves 1000 null values for each species randomization, you'll have as many files as randomizations
  saveRDS(MPDid, file=paste0("null_MPD_IndSwap",j,".rds"))
}
##################################################




############
##Calculate means across trees. read in each iteration (j) rds, and calculate the mean index values across the rows, bind all those together into a streamlined null values file (1000 cols)
PEvemean <- matrix(NA,length(poIDs),1000) 
MPDmean <- matrix(NA,length(poIDs),1000) 

for (z in 1:1000){ 
  
  z1 <- readRDS(file=paste0("null_PEve_IndSwap",z,".rds"))
  
  PEvemean[,z] <- rowMeans(z1[,1:100]) #calculate mean values across the N trees sampled
  PEvemeanFinal <- cbind(PEveid[1001:1004],PEvemean)  #will need to change to 101:104 after doing manually/for other indices
  saveRDS(PEvemeanFinal, file="PEvemeanNull_IndSwap1000.rds")
  
}



for (y in 1:1000){ 
  
  y1 <- readRDS(file=paste0("null_MPD_IndSwap",y,".rds"))
  
  MPDmean[,y] <- rowMeans(y1[,1:100]) #calculate mean values across the N trees sampled
  MPDmeanFinal <- cbind(MPDid[1001:1004],MPDmean) #will need to change to 101:104 after doing manually/for other indices
  saveRDS(MPDmeanFinal, file="MPDmeanNull_IndSwap.rds")
}


###
####PEve SES
###

#### compare observed empirical value to null distribution of index 
PEve.obs <- readRDS(file="PEve.rds")
PEvenull <- readRDS(file= "PEvemeanNull_IndSwap1000.rds")
###check order
table(PEve.obs$poIDs == PEvenull$poIDs) #NOT THE SAME ORDER!
head(PEve.obs[,101:105])
head(PEvenull[,1:5])

#to unify the rows in the order of the observed data use left join, then select those columns desired, in desired order
PEveObsNullJoin <- right_join(PEve.obs[,101:105], PEvenull, by = "poIDs")
head(PEveObsNullJoin)
#check what they match
unique(PEve.obs$poIDs == PEveObsNullJoin$poIDs)
unique(PEveObsNullJoin$elevation.x == PEveObsNullJoin$elevation.y)

PEvenull.n <- cbind(PEve.obs[,101:105],PEveObsNullJoin[,9:ncol(PEveObsNullJoin)]) #the first column of values is the observed value, the other columns are random iterations' values, ##rows are for each site (bin).

head(PEvenull.n[,1:5])
head(PEve.obs[,101:105])

## Rank observed value, calculate SES and p value
obs.null.output <- as.matrix(PEvenull.n[,5:(ncol(PEvenull.n))])
#quantify the rank or the observed value in the null distribution for each community
PEve.obs$PEve.obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
#values range from 1 to the the number of iterations run +1

#quantify the Standardized Effect Size;  SES = (observed - null mean) / null sd
PEve.obs$SES.PEve <- (obs.null.output[,1] - apply(obs.null.output[,2:(ncol(obs.null.output))], 1, mean)) / apply(obs.null.output[,2:(ncol(obs.null.output))], 1, sd)

#calculate p values 
PEve.obs$PEve.p <- apply(cbind(obs.null.output[,1], obs.null.output[,2:(ncol(obs.null.output))]), MARGIN = 1, rank)[1,] / 1001 

head(PEve.obs)
summary(PEve.obs)
hist(PEve.obs$PEve_mean)
hist(PEve.obs$PEve.obs.rank)
hist(PEve.obs$SES.PEve)

#########################investigating NA's
PEve.obs[PEve.obs$SES.PEve == "NaN",] #print any NAs:  Luoji 2600m, JiajinN 3600m, Oquirrh 2500m  #all the values are same across randomizations, thus divides by sd=0.  this is because the bin either contains EVERY species found on the mountain, or because the bin contains only species that are found IN EVERY bin, and thus composition is invariant across randomizations.
################################

#change SES NAs to 0
PEve.obs$SES.PEve[is.na(PEve.obs$SES.PEve)]<-0
summary(PEve.obs)


null_IndSwap_SKEW <- skew.func(obs.null.output)  #generates an overall skew summary for the whole index nulls across all poids

#or append it to the main summary file
PEve.obs$PEve.skew <- apply(obs.null.output[,2:(ncol(obs.null.output))], MARGIN = 1, skewness)
head(PEve.obs)
summary(PEve.obs)
PEve.obs[PEve.obs$PEve.skew== "NaN",] 
PEve.obs$PEve.skew[is.na(PEve.obs$PEve.skew)]<-0
hist(PEve.obs$PEve.skew)
quantile(PEve.obs$PEve.skew, probs = c(.01, .025, .05, .95, .975, .99), na.rm = T)


#save out PEve SES results bound to observed results
saveRDS(PEve.obs, file="PEvenull_IndSwap1000.rds")
##########################################################################################



###
####MPD SES
###

#### compare observed empirical value to null distribution of index 
MPD.obs <- readRDS(file="MPD.rds")
MPDnull <- readRDS(file= "MPDmeanNull_IndSwap1000.rds")
###check order
MPD.obs$poIDs == MPDnull$poIDs #NOT THE SAME ORDER!
head(MPD.obs[,101:105])
head(MPDnull[,1:5])

#to unify the rows in the order of the observed data use left join, then select those columns desired, in desired order
MPDObsNullJoin <- right_join(MPD.obs[,101:105], MPDnull, by = "poIDs")
head(MPDObsNullJoin)
#check what they match
unique(MPD.obs$poIDs == MPDObsNullJoin$poIDs)
unique(MPDObsNullJoin$elevation.x == MPDObsNullJoin$elevation.y)

MPDnull.n <- cbind(MPD.obs[,101:105],MPDObsNullJoin[,9:ncol(MPDObsNullJoin)]) #the first column of values is the observed value, the other columns are random iterations' values, ##rows are for each site (bin).

head(MPDnull.n[,1:5])
head(MPD.obs[,101:105])

## Rank observed value, calculate SES and p value
obs.null.output <- as.matrix(MPDnull.n[,5:(ncol(MPDnull.n))])
#quantify the rank or the observed value in the null distribution for each community
MPD.obs$MPD.obs.rank <- apply(obs.null.output, MARGIN = 1, rank)[1,]
#values range from 1 to the the number of iterations run +1

#quantify the Standardized Effect Size;  SES = (observed - null mean) / null sd
MPD.obs$SES.MPD <- (obs.null.output[,1] - apply(obs.null.output[,2:(ncol(obs.null.output))], 1, mean)) / apply(obs.null.output[,2:(ncol(obs.null.output))], 1, sd)

#calculate p values 
MPD.obs$MPD.p <- apply(cbind(obs.null.output[,1], obs.null.output[,2:(ncol(obs.null.output))]), MARGIN = 1, rank)[1,] / 1001 

head(MPD.obs)
summary(MPD.obs)
hist(MPD.obs$MPD_mean)
hist(MPD.obs$MPD.obs.rank)
hist(MPD.obs$SES.MPD)

#########################investigating NA's
MPD.obs[MPD.obs$SES.MPD == "NaN",] #print any NAs:  Luoji 2600m, JiajinN 3600m, Oquirrh 2500m  #all the values are same across randomizations, thus divides by sd=0.  this is because the bin either contains EVERY species found on the mountain, or because the bin contains only species that are found IN EVERY bin, and thus composition is invariant across randomizations.
################################

#change SES NAs to 0
MPD.obs$SES.MPD[is.na(MPD.obs$SES.MPD)]<-0
summary(MPD.obs)


null_IndSwap_SKEW <- skew.func(obs.null.output)  #generates an overall skew summary for the whole index nulls across all poids

#or append it to the main summary file
MPD.obs$MPD.skew <- apply(obs.null.output[,2:(ncol(obs.null.output))], MARGIN = 1, skewness)
head(MPD.obs)
summary(MPD.obs)
MPD.obs[MPD.obs$MPD.skew== "NaN",] 
MPD.obs$MPD.skew[is.na(MPD.obs$MPD.skew)]<-0
hist(MPD.obs$MPD.skew)
quantile(MPD.obs$MPD.skew, probs = c(.01, .025, .05, .95, .975, .99), na.rm = T)


#save out MPD SES results bound to observed results
saveRDS(MPD.obs, file="MPDnull_IndSwap1000.rds")
##########################################################################################