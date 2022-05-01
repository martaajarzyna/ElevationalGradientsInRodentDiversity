###Rodent diversity patterns along elevational gradients  - INLA Models per each individual mountain
##Kohli, Miyajima and Jarzyna, GEB
##last modified by BAK: 20 April 2021


###################################################
#########STEP 3: INLA models per mountain
###################################################
require(INLA)
require(dplyr)
#### INLA global models

#read in the master summary of all indices and mtn chars created in the Figures script.
FDall <- read.csv(file="GlobalDivMtn_IndexSummary.csv")
str(FDall)

#subset to what we need in the model for simplicity
FD <- select(FDall, 1:19, AI.0.5deg, 32:35,40:43,49,56:59)
str(FD)

FD$moIDs <- as.factor(FD$Mtn) #create new column to hold factors of mtn IDs

moIDs <- unique(FD$moIDs)
n.loc <- length(unique(FD$moIDs))
FD$moID1 <- as.numeric(FD$moIDs)
FD$moID2 <- FD$moID1
FD$moID3 <- FD$moID1
##Log transform raw FD indices that are not beta distributions
FD$PDlog <- log(FD$PD)
FD$MPDlog <- log(FD$MPD)
FD$Vollog <- log(FD$Vol+1)
FD$FDislog <- log(FD$FDis+1)



#####################################################
######
#per mountain models of just elevation
#######
MtnWAIC <- matrix(NA,length(moIDs),2)
###########SR 
#linear
formIndiv1 <- SR ~ elev_norm

for (i in 1:length(moIDs)){
  resultIndiv1 <- inla(formIndiv1, data=FD[FD$Mtn == moIDs[i],], family="poisson", control.compute = list(dic = TRUE, waic = TRUE))
  #summary(resultIndiv) 
 # saveRDS(resultIndiv1, file=paste0("SR_elevation_",moIDs[i],"_1.rds"))
  MtnWAIC[i,1] <- resultIndiv1$waic$waic
}

#Mid-elevation peak
formIndiv2 <- SR ~ elev_norm + I(elev_norm^2) 

for (i in 1:length(moIDs)){
  resultIndiv2 <- inla(formIndiv2, data=FD[FD$Mtn == moIDs[i],], family="poisson", control.compute = list(dic = TRUE, waic = TRUE))
  #summary(resultIndiv) 
 # saveRDS(resultIndiv2, file=paste0("SR_elevation_",moIDs[i],"_2.rds"))
  MtnWAIC[i,2] <- resultIndiv2$waic$waic
}

MtnWAIC_moIDs <- as.data.frame(cbind(as.character(moIDs),MtnWAIC))
MtnWAIC_moIDs
write.csv(MtnWAIC_moIDs, file = "IndivMtnWAIC_SR.csv", row.names = F)

###############evenness 
MtnWAIC <- matrix(NA,length(moIDs),2)
#linear
formIndiv1 <- FEve ~ elev_norm

for (i in 1:length(moIDs)){
  resultIndiv1 <- inla(formIndiv1, data=FD[FD$Mtn == moIDs[i],], family="beta", control.compute = list(dic = TRUE, waic = TRUE))
  #summary(resultIndiv)  
  #saveRDS(resultIndiv1, file=paste0("PEve_elevation_",moIDs[i],"_1.rds"))
  MtnWAIC[i,1] <- resultIndiv1$waic$waic
}

#Mid-elevation peak
formIndiv2 <- FEve ~ elev_norm + I(elev_norm^2) 

for (i in 1:length(moIDs)){
  resultIndiv2 <- inla(formIndiv2, data=FD[FD$Mtn == moIDs[i],], family="beta", control.compute = list(dic = TRUE, waic = TRUE))
  #summary(resultIndiv) 
  #saveRDS(resultIndiv2, file=paste0("PEve_elevation_",moIDs[i],"_2.rds"))
  MtnWAIC[i,2] <- resultIndiv2$waic$waic
}

MtnWAIC_moIDs <- as.data.frame(cbind(as.character(moIDs),MtnWAIC))
MtnWAIC_moIDs
write.csv(MtnWAIC_moIDs, file = "IndivMtnWAIC_FEve.csv", row.names = F)

################raw FD/PD richness and dispersion and ALL SES indices
MtnWAIC <- matrix(NA,length(moIDs),2)
#linear
formIndiv1 <- SES.MPD ~ elev_norm

  for (i in 1:length(moIDs)){
    resultIndiv1 <- inla(formIndiv1, data=FD[FD$Mtn == moIDs[i],], family="gaussian", verbose = F, control.compute = list(dic = TRUE, waic = TRUE))
  #summary(resultIndiv)  
  #saveRDS(resultIndiv1, file=paste0("MPDlog_elevation_",moIDs[i],"_1.rds"))
  MtnWAIC[i,1] <- resultIndiv1$waic$waic
}

#Mid-elevation peak
formIndiv2 <- SES.MPD ~ elev_norm + I(elev_norm^2) 

for (i in 1:length(moIDs)){
  resultIndiv2 <- inla(formIndiv2, data=FD[FD$Mtn == moIDs[i],], family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))
  #summary(resultIndiv) 
 # saveRDS(resultIndiv2, file=paste0("MPDlog_elevation_",moIDs[i],"_2.rds"))
  MtnWAIC[i,2] <- resultIndiv2$waic$waic
}

MtnWAIC_moIDs <- as.data.frame(cbind(as.character(moIDs),MtnWAIC))
MtnWAIC_moIDs
write.csv(MtnWAIC_moIDs, file = "IndivMtnWAIC_SES.MPD.csv", row.names = F)
#################################################################