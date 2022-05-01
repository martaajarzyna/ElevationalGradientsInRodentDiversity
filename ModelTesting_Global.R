###Rodent diversity patterns along elevational gradients  - INLA Models
##Kohli, Miyajima and Jarzyna, GEB

##last modified by BAK: 19 May 2021



###################################################
#########STEP 3: INLA models 
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

##Log transform raw FD indices (can't do for SES) to make gaussian appropriate - for Vol, FDis, all values are below 1 (well below), so used a constant of 1 to avoid negative values. Not a problem with PD or MPD.


FD$PDlog <- log(FD$PD)
FD$MPDlog <- log(FD$MPD)
FD$Vollog <- log(FD$Vol+1)
FD$FDislog <- log(FD$FDis+1)

summary(FD)

#######
#only indices
FDsimp <- select(FD, 4:16, 38:41)
#correlation between indices
IndexCorrs <- cor(FD[4:16])


#######

##############################################################
#FULL MODEL generic structure
    #RESPONSE variables: Raw and SES indices (SR + 12 FD) 
#3 main effect PREDICTOR variables (=ecological variables): per-mtn-normalized (then scaled) elevation (and squared), scaled aridity (from 0.5 degree grid cells), scaled absolute latitude  
    #simple covariate PREDICTORS (=methodological variables): DataSource, EffortGrade [no interactions]
    #random effects: mtn id

####Error distribution structures
names(inla.models()$likelihood)
#SR =  poisson
#FEve/PEve = beta 
#others = gaussian
#######################################################################





##################################
#GLOBAL - Full model
#####

#other indices (have to change the response variable name for each but model structure the same)
#Complete global model - linear
form <- SES.MPD ~ elev_normScaledGlob + aridity_scaled + LatitudeAbs_scaled + elev_normScaledGlob:aridity_scaled + elev_normScaledGlob:LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE)) 

summary(result) 
saveRDS(result, file="DivGlobal_SES.MPD_FullModel_1.rds")

#Complete global model - mid-elev peak
form <- SES.MPD ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + aridity_scaled + LatitudeAbs_scaled + elev_normScaledGlob:aridity_scaled + I(elev_normScaledGlob^2):aridity_scaled + elev_normScaledGlob:LatitudeAbs_scaled + I(elev_normScaledGlob^2):LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE)) 

summary(result)  
saveRDS(result, file="DivGlobal_SES.MPD_FullModel_2.rds")
#########

#raw evenness indices

#Complete global model - linear
form <- PEve ~ elev_normScaledGlob + aridity_scaled + LatitudeAbs_scaled + elev_normScaledGlob:aridity_scaled + elev_normScaledGlob:LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE)) 

summary(result) 
saveRDS(result, file="DivGlobal_PEve_FullModel_1.rds")

#Complete global model - mid-elev peak
form <- PEve ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + aridity_scaled + LatitudeAbs_scaled + elev_normScaledGlob:aridity_scaled + I(elev_normScaledGlob^2):aridity_scaled + elev_normScaledGlob:LatitudeAbs_scaled + I(elev_normScaledGlob^2):LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE)) 

summary(result)  
saveRDS(result, file="DivGlobal_PEve_FullModel_2.rds")
########################

####richness
#Complete global model - linear
form <- SR ~ elev_normScaledGlob + aridity_scaled + LatitudeAbs_scaled + elev_normScaledGlob:aridity_scaled + elev_normScaledGlob:LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")

result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE, link = 1)) 
summary(result)  
saveRDS(result, file="DivGlobal_SRpois_FullModel_1.rds")


#Complete global model - mid-elev peak
form <- SR ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + aridity_scaled + LatitudeAbs_scaled + elev_normScaledGlob:aridity_scaled + I(elev_normScaledGlob^2):aridity_scaled + elev_normScaledGlob:LatitudeAbs_scaled + I(elev_normScaledGlob^2):LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")

result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE), control.predictor = list(compute = TRUE, link = 1)) 
summary(result)  
saveRDS(result, file="DivGlobal_SRpois_FullModel_2.rds")

#############

##
#Global - elevation only
###


####
#SR
####

# Linear relationship
form <- SR ~ elev_normScaledGlob + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE))

summary(result)  
saveRDS(result, file="DivGlobal_SRpois_ElevOnly_1.rds")

##NOTE: Files ending in _1 are with a linear elevation component only
#Files ending in _2 include a quadratic component
#Reponse variable is in the file name.

#if need to read model results back in
#x1 <- readRDS(file="DivGlobal_SRpois_ElevOnly_1.rds")
#summary(x1)

# Mid-elevation peak
form <- SR ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")  #random effects are f()
result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE))

summary(result) 
saveRDS(result, file="DivGlobal_SRpois_ElevOnly_2.rds")

#######################################################################

####
###evenness (raw FEve, PEve)
####

#linear
form <- PEve ~ elev_normScaledGlob + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE))  

summary(result) 
saveRDS(result, file="DivGlobal_PEve_ElevOnly_1.rds")

# Mid-elevation peak
form <- PEve ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")  #random effects are f()
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE))

summary(result)  
saveRDS(result, file="DivGlobal_PEve_ElevOnly_2.rds")
###########################################################################

####
###raw FD/PD richness and dispersion and ALL SES indices
####

# Linear relationship
form <- SES.MPD ~ elev_normScaledGlob + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))  

summary(result)  
saveRDS(result, file="DivGlobal_SES.MPD_ElevOnly_1.rds")

# Mid-elevation peak
form <- SES.MPD ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")  #random effects are f()
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))

summary(result)  
saveRDS(result, file="DivGlobal_SES.MPD_ElevOnly_2.rds")
#######################################################################



#####
#Latitude 

### INLA models with latitude as an interaction term
####

####
###SR 
####
form <- SR ~ elev_normScaledGlob + LatitudeAbs_scaled + elev_normScaledGlob:LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE))
summary(result)
saveRDS(result, file="DivGlobal_SRpois_ElevLat_1.rds")

# Mid-elevation peak
form <- SR ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + LatitudeAbs_scaled + elev_normScaledGlob:LatitudeAbs_scaled + I(elev_normScaledGlob^2):LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE))
summary(result) #for scaled elev: WAIC 2891.00 (not an improvement)

saveRDS(result, file="DivGlobal_SRpois_ElevLat_2.rds")
#####################################################################

####
###evenness (raw FEve, PEve)
####
form <- PEve ~ elev_normScaledGlob + LatitudeAbs_scaled + elev_normScaledGlob:LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE))

summary(result)
saveRDS(result, file="DivGlobal_PEve_ElevLat_1.rds")

# Mid-elevation peak
form <- PEve ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + LatitudeAbs_scaled + elev_normScaledGlob:LatitudeAbs_scaled + I(elev_normScaledGlob^2):LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE))

summary(result)
saveRDS(result, file="DivGlobal_PEve_ElevLat_2.rds")
#####################################################################

#all others - gaussian
# Linear relationship
form <- SES.MPD ~ elev_normScaledGlob + LatitudeAbs_scaled + elev_normScaledGlob:LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))

summary(result) 
saveRDS(result, file="DivGlobal_SES.MPD_ElevLat_1.rds")

# Mid-elevation peak
form <- SES.MPD ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + LatitudeAbs_scaled + elev_normScaledGlob:LatitudeAbs_scaled + I(elev_normScaledGlob^2):LatitudeAbs_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))

summary(result) 
saveRDS(result, file="DivGlobal_SES.MPD_ElevLat_2.rds")

####################################################################

#####
#aridity
####


### INLA models with aridity as an interaction term

####
###SR 
####
#linear
form <- SR ~ elev_normScaledGlob + aridity_scaled + elev_normScaledGlob:aridity_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE))
summary(result)
saveRDS(result, file="DivGlobal_SRpois_ElevAI_1.rds")

# Mid-elevation peak
form <- SR ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + aridity_scaled + elev_normScaledGlob:aridity_scaled + I(elev_normScaledGlob^2):aridity_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="poisson", control.compute = list(dic = TRUE, waic = TRUE))
summary(result) 
saveRDS(result, file="DivGlobal_SRpois_ElevAI_2.rds")
#####################################################################

####
###evenness (raw FEve, PEve)
####

#linear
form <- PEve ~ elev_normScaledGlob + aridity_scaled + elev_normScaledGlob:aridity_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE))

summary(result)
saveRDS(result, file="DivGlobal_PEve_ElevAI_1.rds")

# Mid-elevation peak
form <- PEve ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + aridity_scaled + elev_normScaledGlob:aridity_scaled + I(elev_normScaledGlob^2):aridity_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="beta", control.compute = list(dic = TRUE, waic = TRUE))

summary(result) 
saveRDS(result, file="DivGlobal_PEve_ElevAI_2.rds")
#####################################################################

#all others

# Linear relationship
form <- SES.MPD ~ elev_normScaledGlob + aridity_scaled + elev_normScaledGlob:aridity_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))

summary(result) 
saveRDS(result, file="DivGlobal_SES.MPD_ElevAI_1.rds")

# Mid-elevation peak
form <- SES.MPD ~ elev_normScaledGlob + I(elev_normScaledGlob^2) + aridity_scaled + elev_normScaledGlob:aridity_scaled + I(elev_normScaledGlob^2):aridity_scaled + DataSource + EffortGrade + f(moID1, model = "iid") + f(moID2, elev_normScaledGlob, copy="moID1") + f(moID3, I(elev_normScaledGlob^2), copy="moID1")
result <- inla(form, data=FD, family="gaussian", control.compute = list(dic = TRUE, waic = TRUE))

summary(result) 
saveRDS(result, file="DivGlobal_SES.MPD_ElevAI_2.rds")
#
##################################################################
