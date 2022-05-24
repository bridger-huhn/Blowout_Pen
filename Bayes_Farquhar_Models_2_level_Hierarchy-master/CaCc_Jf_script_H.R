## Bayesian Estimate of PS paramters from Brassica ACi curves ###
setwd("/Users/bridgerhuhn/Documents/Research/Blowout_Pen/Bayes_Farquhar_Models_2_level_Hierarchy-master")   ### Change 4 server
require(rjags)### Change 4 server
#source("farQ_TPU.R")  ### Use output pars to predict net assimilation
source("CaCc_Jf_model_H.R") ### Bayesian Modelw/ Temp Dependency & farQ limitation

#### Set up for rjags #########
parameters = c("Vcmax", "Rd","gammaS" ,"Kc", "Ko",
"mu.Vcmax", "mu.Rd","mu.gammaS",
"tau.Vcmax", "tau.Rd","tau.gammaS", "tau") ### pars to be monitored

adaptSteps = 100000             # Number of steps to "tune" the samplers.
burnInSteps = 200000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps=50000        # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per

######################
######################
#####DATA SET UP######
######################
ACdat<-read.delim("~/Documents/ACiP4/Data/ACdata_reduced1_clean.txt")  ### Change 4 server
#useID<-read.delim("~/Documents/ACiP4/Data/IDs_to_use-9-9-17.txt")
#gID<-useID$ID
#IDgeno<-as.character(useID$geno)
group <- length(unique(d$comment))
ACi<-d
####  DATA SET  UP  ############
#### PULL OUT Data for each genotype   ###
ID<-sort(unique(ACi$comment))

#### ID Data
IDdat<-vector("list", group)
for(i in 1:group){
    IDdat[[i]]<-as.factor(ACi[ACi$geno==ID[i],]$AC.ID)
}


#### Photo Data
A = vector("list", group)
for(i in 1:group){
    A[[i]]<-ACi[ACi$geno==ID[i],]$Photo
}
## CP = CO2 concentraion intercelluar space (in partial pressure)
ACi$CP<-ACi$Ci*ACi$Press/1000   ## calcs CO2 concentraion intercelluar space (in partial pressure)
CP = vector("list", group)   ### pulls for each ID
for(i in 1:group){
    CP[[i]]<-ACi[ACi$geno==ID[i],]$CP
}

### Leaf Temp (C) ###
T = vector("list", group) ### pulls for each ID
for(i in 1:group){
    T[[i]]<-ACi[ACi$geno==ID[i],]$Tleaf
}


### Leaf Temp (K) ###
K = vector("list", group)  ### pulls for each ID
for(i in 1:group){
    K[[i]]<-ACi[ACi$geno==ID[i],]$Tleaf+273.15
}


###  O2 in pp from atmospheric pressure (kPa)
### PP(Pa) =  MolFact * pressure of the gas mixture
ACi$O<-(.21)*(ACi$Press*1000)
O = vector("list", group)   ### pulls for each ID
for(i in 1:group){
    O[[i]]<-ACi[ACi$geno==ID[i],]$O
}

### PARi
Jf = vector("list", group)   ### pulls for each ID
for(i in 1:group){
    Jf[[i]]<-ACi[ACi$geno==ID[i],]$ETR
}

## CP = CO2 concentraion intercelluar space (in partial pressure)
ACi$CP<-(ACi$CO2S)/ 1.0e6 *(ACi$Press*1000.0)   ## calcs CO2 concentraion sample chamber in partial pressure)
CP = vector("list", group)   ### pulls for each ID
for(i in 1:group){
    CP[[i]]<-ACi[ACi$geno==ID[i],]$CP
}

### Cond converted to cond to Co2 umol m-s s-1 Pa-1 ###
#   ACi$CndCO2   --- Conductane to CO2 in (mol/m2/s)
# convert mol/m2/s to umol/m2/s/Pa   = CndCO2 * 10e6 / Pressure(Pa)
ACi$Cond_Conv<-ACi$CndCO2* 1.0e6 /(ACi$Press*1000.0)
g = vector("list", group)   ### pulls for each ID
for(i in 1:group){
    g[[i]]<-ACi[ACi$geno==ID[i],]$Cond_Conv
}


#### Sample size
#N<-36
## Contants ##
Kref=298.15; R=0.008314; Tref=25
for (i in 1:group) {
    assign(paste0("datalist", i), list(N=length(A[[1]]), NID=length(unique(IDdat[[1]])), An=A[[1]], Ca=CP[[1]], g=g[[1]],  Jf=Jf[[1]], O=O[[1]], ID=IDdat[[1]], Kref=Kref, R=R, Tref=Tref))
}


#################################
##################################
##### Impliment model in JAGS ####
#################################
### running each curve
source("CaCc_Jf_model_H.R") ### Bayesian Modelw/ no Temp Dependency, No Gm and chl flo to describe farQ light limitation

print("initialize models")
model1 <- jags.model(textConnection(CaCc_Jf),
data = datalist1, n.chains=nChains , n.adapt=adaptSteps)

model2 <- jags.model(textConnection(CaCc_Jf),
data = datalist2, n.chains=nChains , n.adapt=adaptSteps)
model3 <- jags.model(textConnection(CaCc_Jf),
data = datalist3, n.chains=nChains , n.adapt=adaptSteps)
model4 <- jags.model(textConnection(CaCc_Jf),
data = datalist4, n.chains=nChains , n.adapt=adaptSteps)
model5 <- jags.model(textConnection(CaCc_Jf),
data = datalist5, n.chains=nChains , n.adapt=adaptSteps)
model6 <- jags.model(textConnection(CaCc_Jf),
data = datalist6, n.chains=nChains , n.adapt=adaptSteps)


#################################
print("updating")
update(model1, burnInSteps) # Burnin for burnInSteps samples
update(model2, burnInSteps)
update(model3, burnInSteps)
update(model4, burnInSteps)
update(model5, burnInSteps)
update(model6, burnInSteps)
##########################################

###### SAMPLE all 18 individual curves ###

##########################################


##########################################
print("sampling chains")
##### mcmc_samples  model 1 IDs  #####
#####   ID order   1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760
##########################################

mcmc_samples1<- coda.samples(model1,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
####### Plot results #####
#plot(mcmc_samples1)
mcmcChain_r301_3 = as.matrix( mcmc_samples1)
chainLength = NROW(mcmcChain_r301_3)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain_r301_3[, "tau" ] )
mcmcChain_r301_3 = as.data.frame(cbind( mcmcChain_r301_3, sigma ))
g1<-gelman.diag(mcmc_samples1)
####
mcmc_samples2<- coda.samples(model2,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_r46_3 = as.matrix( mcmc_samples2)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain_r301_3[, "tau" ] )
mcmcChain_r46_3 = as.data.frame(cbind( mcmcChain_r46_3, sigma ))
g2<-gelman.diag(mcmc_samples2)
######
mcmc_samples3 <- coda.samples(model3,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_bro_3 = as.matrix( mcmc_samples3)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_bro_3[, "tau" ] )
mcmcChain_bro_3 = as.data.frame(cbind( mcmcChain_bro_3, sigma ))
g3<-gelman.diag(mcmc_samples3)
######
mcmc_samples4 <- coda.samples(model4,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_cab_3 = as.matrix( mcmc_samples4)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_cab_3[, "tau" ] )
mcmcChain_cab_3 = as.data.frame(cbind( mcmcChain_cab_3, sigma ))
g4<-gelman.diag(mcmc_samples4)
######
mcmc_samples5 <- coda.samples(model5,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_oil_3 = as.matrix( mcmc_samples5)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_oil_3[, "tau" ] )
mcmcChain_oil_3 = as.data.frame(cbind( mcmcChain_oil_3, sigma ))
g5<-gelman.diag(mcmc_samples5)
######
mcmc_samples6 <- coda.samples(model6,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
mcmcChain_tur_3 = as.matrix( mcmc_samples6)
# Convert precision (tau) to SD###
sigma = 1  / sqrt( mcmcChain_tur_3[, "tau" ] )
mcmcChain_tur_3 = as.data.frame(cbind( mcmcChain_tur_3, sigma ))
g6<-gelman.diag(mcmc_samples6)

##########################################

##########################################
######### Gelman potential scale reduction factor summary for all par and multivariate #############
gelmin<-c(min(g1$psrf[,1]), min(g2$psrf[,1]),min(g3$psrf[,1]), min(g4$psrf[,1]),
min(g5$psrf[,1]), min(g6$psrf[,1]))


gelmax<-c(max(g1$psrf[,1]), max(g2$psrf[,1]),max(g3$psrf[,1]), max(g4$psrf[,1]),
max(g5$psrf[,1]), max(g6$psrf[,1]))

gelmulit<-c(g1$mpsrf, g2$mpsrf, g3$mpsrf, g4$mpsrf,g5$mpsrf, g6$mpsrf)

gelsum<-as.data.frame(cbind(gelmin, gelmax,gelmulit))
colnames(gelsum)<-c("gel_min", "gel_max","gel_multi")

##########################################
#########  SAVE Samples ####################
# ID 1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760
###
print("writing samples")


setwd("~/Documents/ACiP4/Post_Data_H")
write.table(mcmcChain_r301_3,"mcmcChain_r301_3", sep="\t", col.name=TRUE)
write.table(mcmcChain_r46_3,"mcmcChain_r46_3", sep="\t", col.name=TRUE)
write.table(mcmcChain_bro_3,"mcmcChain_bro_3", sep="\t", col.name=TRUE)
write.table(mcmcChain_cab_3,"mcmcChain_cab_3", sep="\t", col.name=TRUE)
write.table(mcmcChain_oil_3,"mcmcChain_oil_3", sep="\t", col.name=TRUE)
write.table(mcmcChain_tur_3,"mcmcChain_tur_3", sep="\t", col.name=TRUE)
#########################################################################

#########################################################################

#################    PARAMETER ESTIMATES      ###########################

#########################################################################

#########################################################################
###  for revision 3 (2_24_18) unbalanced data results in a different number of ind par est ###
###  for simplicity select only 3 indivudals for saving in data structure
common_cols <- colnames(mcmcChain_oil_3)
mcmcChain_r301_3<-subset(mcmcChain_r301_3, select = common_cols)
mcmcChain_r46_3<-subset(mcmcChain_r46_3, select = common_cols)
mcmcChain_bro_3<-subset(mcmcChain_bro_3, select = common_cols)
mcmcChain_cab_3<-subset(mcmcChain_cab_3, select = common_cols)
mcmcChain_tur_3<-subset(mcmcChain_tur_3, select = common_cols)
#############################################


#############################################

EstPars1<-apply(mcmcChain_r301_3, 2, median)
EstPars2<-apply(mcmcChain_r46_3, 2, median)
EstPars3<-apply(mcmcChain_bro_3, 2, median)
EstPars4<-apply(mcmcChain_cab_3, 2, median)
EstPars5<-apply(mcmcChain_oil_3, 2, median)
EstPars6<-apply(mcmcChain_tur_3, 2, median)

ParEst<-rbind(EstPars1,EstPars2,EstPars3,EstPars4,EstPars5,EstPars6)




#########################################################################

#########################################################################

########################################################################

#########################################################################

#################   Calc DIC with pD criteria  ###########################

#########################################################################

#########################################################################
# ID 1845 1894 1902 1937 2157 2208 2304 2319 2320 2349 2603 2629 2655 2696 2712 2723 2737 2760

print("sampling for DIC")
dic1 <- dic.samples(model1, DICsteps, "pD")
dic1a <- dic.samples(model1, DICsteps, "popt")
dic_r301_3<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model2, DICsteps, "pD")
dic1a <- dic.samples(model2, DICsteps, "popt")
dic_r46_3<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model3, DICsteps, "pD")
dic1a <- dic.samples(model3, DICsteps, "popt")
dic_bro_3<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model4, DICsteps, "pD")
dic1a <- dic.samples(model4, DICsteps, "popt")
dic_cab_3<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model5, DICsteps, "pD")
dic1a <- dic.samples(model5, DICsteps, "popt")
dic_oil_3<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))

dic1 <- dic.samples(model6, DICsteps, "pD")
dic1a <- dic.samples(model6, DICsteps, "popt")
dic_tur_3<-c(sum(dic1$deviance), sum(dic1$penalty),sum(dic1$deviance, dic1$penalty),
sum(dic1a$deviance), sum(dic1a$penalty), sum(dic1a$deviance, dic1a$penalty))


dics<-rbind(dic_r301_3, dic_r46_3, dic_bro_3, dic_cab_3, dic_oil_3,dic_tur_3)
colnames(dics)<-c("pD_dev", "pD_pen", "pD_DIC", "popt_dev", "popt_pen", "popt_DIC")



modelname<-rep("CaCc_Jf",6)

###################################################

###################################################

#   FINAL TABLE W/ median estimates and DICS ########

###################################################
print("writing PAR Estimates")
ParEstFull<-as.data.frame(cbind(modelname, ID, ParEst,dics, gelsum))
write.table(ParEstFull,"ParEstFull_model_3", sep="\t", col.name=TRUE,row.names = F)

print("Finished")



require(ggmcmc)

pdf("~/Documents/ACiP4/Trace_plots/density_r301_3.pdf",
colormodel='cmyk' ,width=6.25, height=(80))
ggs_density(ggs(mcmc_samples1))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_46_3.pdf",
colormodel='cmyk',width=6.25, height=(80))
ggs_density(ggs(mcmc_samples2))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_bro_3.pdf", colormodel='cmyk',width=6.25, height=(80))
ggs_density(ggs(mcmc_samples3))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_cab_3.pdf", colormodel='cmyk',width=6.25, height=(80))
ggs_density(ggs(mcmc_samples4))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_oil_3.pdf", colormodel='cmyk',width=6.25, height=(80))
ggs_density(ggs(mcmc_samples5))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/density_tur_3.pdf", colormodel='cmyk',width=6.25, height=(80))
ggs_density(ggs(mcmc_samples6))
dev.off()




pdf("~/Documents/ACiP4/Trace_plots/cor_r301_3.pdf",
colormodel='cmyk',width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples1))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_46_3.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples2))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_bro_3.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples3))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_cab_3.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples4))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_oil_3.pdf", colormodel='cmyk',,width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples5))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/cor_tur_3.pdf", colormodel='cmyk',width=10, height=(10))
ggs_crosscorrelation(ggs(mcmc_samples6))
dev.off()


pdf("~/Documents/ACiP4/Trace_plots/Rhat_r301_3.pdf",
colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples1))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_46_3.pdf",
colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples2))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_bro_3.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples3))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_cab_3.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples4))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_oil_3.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples5))
dev.off()

pdf("~/Documents/ACiP4/Trace_plots/Rhat_tur_3.pdf", colormodel='cmyk',width=6.25, height=(20))
ggs_Rhat(ggs(mcmc_samples6))
dev.off()

































