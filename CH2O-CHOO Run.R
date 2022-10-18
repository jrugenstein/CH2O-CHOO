######## CH2O-CHOO Model #############
# This model solves the carbon cycle through time and as a function of different
# runoff sensitivities. It is based on the CLiBeSO-BFF model, originally presented in 
# Rugenstein et al. (2019)--Nature.

# For more detailed explanation, please review Rugenstein and Winkler (in review). 

# NOTE: This script uses parallel programming to speed up model simulations and test a wide range of model parameters. However, the parallelization code provided here only works on Mac computers. Please contact Jeremy Rugenstein for a version of the script that functions on Windows computers.

# Created on: 16-11-2017
# Last Modified: 26-09-2022
# Last Modified by: Jeremy K. C. Rugenstein (CSU; jeremy.rugenstein@colostate.edu)

library(deSolve)
library(seacarb) 
library(plotrix)
library(nls2)
library(np)
library(parallel)
library(parallelMap)
library(tidyverse)
library(reshape2)
library(colorspace)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source(file="CH2O-CHOO functions.R")

#### Initial Carbonate System Values ####
oceanV <- 1.4e21 # L of ocean
d13C.i <- 0 # initial DIC permil
depth <- 1 # meters
temp.o.i <- 15 # initial ocean surface temp in °C

salinity <- 35 # salinity in salinity units
sal <- (salinity/1000)+1 # salinity in kg/L
pCO2.i <- 400 # ppm
pH.i <- 8.15
Mg.Ca.i <- 5 # seawater Mg/Ca ratio
Ca.i <- 0.01028 # mol/kg
Mg2.i <- 0.05654 #Mg.Ca.i*Ca.i # mol/kg
carb.parms.i <- co2sys.init(pH = pH.i,pCO2 = pCO2.i,temp=temp.o.i)
carb.parms.i
omega.i <- carb.parms.i$OmegaCalcite
RCO2.i <- pCO2.i/pCO2.i
DIC.i <- carb.parms.i$DIC*oceanV
DIC.i
Alk.i <- carb.parms.i$ALK*oceanV

#### Model parameters ####
p <- paramfunc()

# Initial Reservoirs
S.i <- 10e6*p['Fwsulf.i'] # mol sulfur; assumes a 10 Ma residence time of S in the ocean; modern-day estimates are 5e19 (Broecker and Peng 1982), but this includes a substantial gypsum component
# Note that, while CH2O-CHOO has a sulfur sub-cycle, it is not activated in this version of the code. 

y0  <- c(DIC=DIC.i, Alk=Alk.i, d13C=d13C.i, S=S.i)
names(y0) <- c(names(y0)[1:3],"S")

#### Time set up and Perturbation ####
#### Initial Set-up 
# This is used to set-up an initial perturbation and is called when the land surface reactivity changes prior to a model simulation
duration.init <- 6e6 # Duration of model run
dt.init <- 100 # time steps (in years)
t.init <- seq(0, duration.init, dt.init)

pert.dur.init <- 1e2 # yrs (must be a factor of 100)
C.extra.init <- 0 # Pg converted to mols C
pert.time.init <- 1e4 # time of perturbation

C.increase.init <- c(rep(0,pert.time.init/dt.init),rep(C.extra.init,pert.dur.init/dt.init),
                rep(0,(length(t.init)-pert.dur.init/dt.init-pert.time.init/dt.init)))[1:length(t.init)]
C.perturbation.init <- approxfun(x=t.init,y=C.increase.init,method="linear",rule=2)
length(C.increase.init)
sum(C.increase.init)
results.exclude <- length(c(rep(0,pert.time.init/dt.init),rep(C.extra.init,pert.dur.init/dt.init)))
results.exclude

# Perturbation set-up
# This is used during a normal perturbation to estimate the pCO2 recovery time 
duration <- 4e6 # Duration of model run
dt <- 100 # time steps (in years)
t <- seq(0, duration, dt)

pert.dur <- 1e2 # yrs (must be a factor of 100)
C.extra <- 1000e15/pert.dur/12 # Pg converted to mols C
pert.time <- 1e4 # time of perturbation

C.increase <- c(rep(0,pert.time/dt),rep(C.extra,pert.dur/dt),
                rep(0,(length(t)-pert.dur/dt-pert.time/dt)))[1:length(t)]
C.perturbation <- approxfun(x=t,y=C.increase,method="linear",rule=2)
d13C.extra <- -25 # ‰
d13C.increase <- c(rep(0,pert.time/dt),d13C.extra,
                   rep(0,(length(t)-pert.time/dt)))[1:length(t)]
d13C.perturbation <- approxfun(x=t,y=d13C.increase,method="linear",rule=2)

#### Test Run of CH2O-CHOO ####
p <- paramfunc(remove = "init")
Conc0.in <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'])$Conc0
p <- c(p,q.sens=0.04,Conc0=Conc0.in,init=FALSE)
p

output <- as.data.frame(lsodes(y=y0, times=t, func=CLiBeSO, parms=p,rtol=1e-6))
colnames(output) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb',
                      'Fborg','Fwsil','RCO2','Temp','Fworg',
                      "Fwsulf","Fbsulf","q.1","Conc.1","Dw.1")

y.time <- output$time

par(mfrow=c(4,1),mar=c(1,5,1,1))
plot(x=y.time/1e6,y=output$pCO2,bty="n",type="l", # CO2
     ylab=expression(paste("Atmospheric CO"[2])),
     xaxt="n")
plot(x=y.time,y=output$d13C,bty='n',type='l',xaxt="n",ylab=expression(paste(delta^13,"C"))) # d13C
plot(x=y.time,y=output$q,bty="n",type="l",xaxt="n",ylab="Runoff (m/yr)") # Runoff
par(mar=c(5,5,1,1))
plot(x=y.time,y=output$Conc,bty='n',type='l',ylab=expression(paste("[C]"[sil])),
     xlab="Time (yr)") # Concentration

#### Test of Different Parameterizations and Weathering Functions ####
# Overall q Sensitivity Vector
q.sens.vec <- c(seq(-0.1,-0.04,0.01),
                seq(-0.03,0.03,0.005),
                seq(0.04,0.1,0.01))

# Standard parameters (ie, MC standard parameters)
y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

p <- paramfunc()
Conc0.in <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'])$Conc0
p <- c(p,Conc0=Conc0.in)
p.add <- list()
  for (i in 1:length(q.sens.vec)){
    p.add[[i]] <- c(p,q.sens=q.sens.vec[i])
  }

detectCores()
parallelStart("multicore",cpus=detectCores()-2)
results1 <- list()
  results1 <- parallelLapply(p.add,CLiBeSO.uber)
parallelStop()

# Ts Sensitivity
# First, achieve a new steady-state given the change in Ts. Second, analyze a perturbation
Ts.vec <- seq(1e4,10e4,2e4)
results.Ts <- list()

for (j in 1:length(Ts.vec)) {
  p <- paramfunc(remove = c("Ts","init"))
  p <- c(p,Conc0=Conc0.in,q.sens=0.04,Ts=Ts.vec[j],init=TRUE)
  p
  
  init2.results <- as.data.frame(lsodes(y=y0, times=t.init, func=CLiBeSO, parms=p,rtol=1e-6))
  colnames(init2.results) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb',
                               'Fborg','Fwsil','RCO2','Temp','Fworg',
                               "Fwsulf","Fbsulf","q.1","Conc.1","Dw.1")
  
  temp.o.i <- 15
  pCO2.i <- tail(init2.results$pCO2)[5]
  y0  <- c(DIC=tail(init2.results$DIC)[5], Alk=tail(init2.results$Alk)[5], 
           d13C=tail(init2.results$d13C)[5], S=tail(init2.results$S)[5])
  names(y0) <- c(names(y0)[1:3],"S")
  new.q0 <- tail(init2.results$q.1)[5]
  
  p <- paramfunc(remove=c("Ts","temp.a.i","q0"))
  p.add <- list()
  for (i in 1:length(q.sens.vec)){
    p.add[[i]] <- c(p,q.sens=q.sens.vec[i],Ts=Ts.vec[j],
                    Conc0=tail(init2.results$Conc.1)[5],q0=new.q0,
                    temp.a.i=tail(init2.results$Temp)[5])
  }

  #### Parallel Execution ####
  detectCores()
  parallelStart("multicore",cpus=detectCores()-3)
    results.Ts[[j]] <- parallelLapply(p.add,CLiBeSO.uber)
  parallelStop()
}

results.Ts <- flatten(results.Ts)

# West 2012 Parameterization
p <- paramfunc(remove = c("Ts","init"))
p <- c(p,Conc0=Conc0.in,q.sens=0.04,Ts=6e4,init=TRUE)
p

y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

init3.results <- as.data.frame(lsodes(y=y0, times=t.init, func=CLiBeSO.west, parms=p,rtol=1e-6))
colnames(init3.results) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb','Fborg','Fwsil','Fextra','Conc','Conc0','RCO2','Temp','Fworg',"Fwsulf","Fbsulf","q.1","Conc.1","Dw.1")

pCO2.i <- tail(init3.results$pCO2)[5]
y0  <- c(DIC=tail(init3.results$DIC)[5], Alk=tail(init3.results$Alk)[5], 
         d13C=tail(init3.results$d13C)[5], S=tail(init3.results$S)[5])
names(y0) <- c(names(y0)[1:3],"S")
new.q0 <- tail(init3.results$q.1)[5]
new.temp.a.i <- tail(init3.results$Temp)[5]

p <- paramfunc(remove = c("Ts","q0","temp.a.i","Fwcarb.i","Fbcarb.i","Fborg.i","Fwsil.i"))
p.add <- list()
for (i in 1:length(q.sens.vec)){
  p.add[[i]] <- c(p,q.sens=q.sens.vec[i],Ts=6e4,q0=new.q0,temp.a.i=new.temp.a.i,
                  Fwcarb.i=tail(init3.results$Fwcarb)[5],
                  Fbcarb.i=tail(init3.results$Fbcarb)[5],
                  Fborg.i=tail(init3.results$Fborg)[5],
                  Fwsil.i=tail(init3.results$Fwsil)[5])
}
length(p.add)
p.add[[1]]

detectCores()
parallelStart("multicore",cpus=detectCores()-2)
results4 <- list()
  results4 <- parallelLapply(p.add,CLiBeSO.west.uber)
parallelStop()

# Constant Concentrations
p <- paramfunc()
p.add <- list()
for (i in 1:length(q.sens.vec)){
  p.add[[i]] <- c(p,Conc0=Conc0.in,q.sens=q.sens.vec[i])
}
length(p.add)
p.add[[1]]

y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

detectCores()
parallelStart("multicore",cpus=detectCores()-2)
results5 <- list()
  results5 <- parallelLapply(p.add,CLiBeSO.chem.uber)
parallelStop()

# Seafloor weathering
p <- paramfunc()
p.add <- list()
for (i in 1:length(q.sens.vec)){
  p.add[[i]] <- c(p,q.sens=q.sens.vec[i])
}
length(p.add)
p.add[[1]]

y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

detectCores()
parallelStart("multicore",cpus=detectCores()-2)
results10 <- list()
  results10 <- parallelLapply(p.add,CLiBeSO.sfw.uber)
parallelStop()

# Two box sensitivity experiment
# NOTE: For the two experiments below, there are actually 3 boxes in case these are needed. However, one of the boxes is effectively "turned off". 

# One Area (ie, Box) of High Sensitivity
C.perturbation.init <- approxfun(x=t.init,y=C.increase.init,method="linear",rule=2)
y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

p <- paramfunc()

Conc0.in.1 <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'])$Conc0
Conc0.in.2 <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'])$Conc0
Conc0.in.3 <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'],Ts=1e3)$Conc0

p <- paramfunc(remove=c("Trop.area","Ts","init"))
p <- c(p,Conc0=c(Conc0.in.1,Conc0.in.2,Conc0.in.3),q.sens=c(0.04,0.04,0.04),Trop.area=0.5,Ts=c(3.9e4,3.9e4,1e3),
       q0=c(0.3,0.3,0.3),init=TRUE)

init4.results <- as.data.frame(lsodes(y=y0, times=t.init, func=CLiBeSO.3, parms=p,rtol=1e-6))
colnames(init4.results) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb', 'Fborg','Fwsil','Fwsil.1','Fwsil.3','RCO2','Temp','Fworg',"Fwsulf","Fbsulf","q.1","q.2","q.3","q.avg","Conc.1","Conc.2","Conc.3","Conc.avg","Dw.1","Dw.2","Dw.3")

pCO2.i <- tail(init4.results$pCO2)[5]
y0  <- c(DIC=tail(init4.results$DIC)[5], Alk=tail(init4.results$Alk)[5], 
         d13C=tail(init4.results$d13C)[5], S=tail(init4.results$S)[5])
names(y0) <- c(names(y0)[1:3],"S")
new.q0 <- c(tail(init4.results$q.1)[5],tail(init4.results$q.2)[5],tail(init4.results$q.3)[5])
new.Conc0 <- c(tail(init4.results$Conc.1)[5],tail(init4.results$Conc.2)[5],tail(init4.results$Conc.3)[5])

p <- paramfunc(remove=c("Trop.area","Ts","temp.a.i","q0"))
p.add <- list()
for (i in 1:length(q.sens.vec)){
  p.add[[i]] <- c(p,Trop.area=0.5,Ts=c(3.9e4,3.9e4,1e3),Conc0=new.Conc0,
                  q.sens=c(q.sens.vec[i],q.sens.vec[i],0.03),
                  q0=new.q0,temp.a.i=tail(init4.results$Temp)[5])
}
length(p.add)
p.add[[1]]

detectCores()
parallelStart("multicore",cpus=detectCores()-2)
results7 <- list()
  results7 <- parallelLapply(p.add,CLiBeSO.3.uber)
parallelStop()

# One Area of High Sensitivity with high runoff sensitivity
C.perturbation.init <- approxfun(x=t.init,y=C.increase.init,method="linear",rule=2)
y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

p <- paramfunc()

Conc0.in.1 <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'],Ts=3.9e4)$Conc0
Conc0.in.2 <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'],Ts=3.9e4)$Conc0
Conc0.in.3 <- MCfunc(pCO2.i,0.3,0.3,Temperature = p['temp.a.i'],Ts=1e3)$Conc0

p <- paramfunc(remove=c("Trop.area","Ts","init"))
p <- c(p,Conc0=c(Conc0.in.1,Conc0.in.2,Conc0.in.3),q.sens=c(0.04,0.04,0.04),Trop.area=0.5,Ts=c(3.9e4,3.9e4,1e3),
       q0=c(0.3,0.3,0.3),init=TRUE)

init6.results <- as.data.frame(lsodes(y=y0, times=t.init, func=CLiBeSO.3, parms=p,rtol=1e-6))
colnames(init6.results) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb', 'Fborg','Fwsil','Fwsil.1','Fwsil.3','RCO2','Temp','Fworg',"Fwsulf","Fbsulf","q.1","q.2","q.3","q.avg","Conc.1","Conc.2","Conc.3","Conc.avg","Dw.1","Dw.2","Dw.3")

pCO2.i <- tail(init6.results$pCO2)[5]
y0  <- c(DIC=tail(init6.results$DIC)[5], Alk=tail(init6.results$Alk)[5], 
         d13C=tail(init6.results$d13C)[5], S=tail(init6.results$S)[5])
names(y0) <- c(names(y0)[1:3],"S")
new.q0 <- c(tail(init6.results$q.1)[5],tail(init6.results$q.2)[5],tail(init6.results$q.3)[5])
new.Conc0 <- c(tail(init6.results$Conc.1)[5],tail(init6.results$Conc.2)[5],tail(init6.results$Conc.3)[5])

p <- paramfunc(remove=c("Trop.area","Ts"))
p.add <- list()
for (i in 1:length(q.sens.vec)){
  p.add[[i]] <- c(p,Trop.area=0.5,Ts=c(3.9e4,3.9e4,1e3),Conc0=new.Conc0,
                  q.sens=c(q.sens.vec[i],q.sens.vec[i],0.06),
                  q0=new.q0,temp.a.i=tail(init4.results$Temp)[5])
}

detectCores()
parallelStart("multicore",cpus=detectCores()-2)
results9 <- list()
  results9 <- parallelLapply(p.add,CLiBeSO.3.uber)
parallelStop()

#### Calculate return times ####
y0 <- init.carb()$y0
pCO2.i <- init.carb()$pCO2

pCO2.recov <- 1 # ie, this sets the "tolerance" for recovery for pCO2

# Standard Parameters
return.time1 <- array(dim=c(length(results1),3))
for (i in 1:nrow(return.time1)){
  time.result <- results1[[i]][[2]]$time[-c(1:150)]
  pCO2.result <- results1[[i]][[2]]$pCO2[-c(1:150)]
  d13C.result <- results1[[i]][[2]]$d13C[-c(1:150)]
  max.q.result <- results1[[i]][[2]]$q[-c(1:150)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result,max.q.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C","q.max")#,"Conc")
  return.time1[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<results1[[i]][[2]]$pCO2[1]+pCO2.recov)]) 
  return.time1[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C,1)==0)]) 
  return.time1[i,3] <- max(pCO2.result$q.max)
}
return.time1 <- as.data.frame(cbind(return.time1,q.sens.vec))
colnames(return.time1) <- c("pCO2","d13C","q","Q")

# Ts Senstivity Experiment
return.time2 <- array(dim=c(length(results.Ts),2))
for (i in 1:nrow(return.time2)){
  time.result <- results.Ts[[i]][[2]]$time[-c(1:250)]
  pCO2.result <- results.Ts[[i]][[2]]$pCO2[-c(1:250)]
  d13C.result <- results.Ts[[i]][[2]]$d13C[-c(1:250)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C")
  return.time2[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<results.Ts[[i]][[2]]$pCO2[1]+(pCO2.recov))]) 
  return.time2[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C[500:length(pCO2.result$d13C)],1)==round(results.Ts[[i]][[2]]$d13C[1],2))]) 
}
return.time2 <- as.data.frame(return.time2)
colnames(return.time2) <- c("pCO2","d13C")
head(return.time2)

q.sens.vec.2 <- numeric(length=length(results.Ts))
Ts.vec.2 <- numeric(length=length(results.Ts))
Dw.vec.2 <- numeric(length=length(results.Ts))
for (i in 1:length(q.sens.vec.2)){
  q.sens.vec.2[i] <- results.Ts[[i]][[1]]['q.sens']
  Ts.vec.2[i] <- results.Ts[[i]][[1]]['Ts']
  Dw.vec.2[i] <- results.Ts[[i]][[2]]$Dw[1]
}

return.time2.mat <- as_tibble(cbind(return.time2,q.sens.vec.2*100,Ts.vec.2,Dw.vec.2))
colnames(return.time2.mat) <- c("pCO2","d13C","Q","Ts","Dw")
return.time2.mat %>% filter(Ts<8e4) -> return.time2.mat
return.time2.mat

# West (2012) Parameterization
return.time4 <- array(dim=c(length(results4),3))
for (i in 1:nrow(return.time4)){
  time.result <- results4[[i]][[2]]$time[-c(1:100)]
  pCO2.result <- results4[[i]][[2]]$pCO2[-c(1:100)]
  d13C.result <- results4[[i]][[2]]$d13C[-c(1:100)]
  max.q.result <- results4[[i]][[2]]$q[-c(1:100)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result,max.q.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C","q.max")
  return.time4[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<pCO2.i+pCO2.recov)]) 
  return.time4[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C,1)==0)]) 
  return.time4[i,3] <- max(pCO2.result$q.max)
}
return.time4 <- as_tibble(cbind(as.data.frame(return.time4),q.sens.vec))
colnames(return.time4) <- c("pCO2","d13C","q","Q")

# Constant Concentrations (ie, chemostatic)
return.time5 <- array(dim=c(length(results5),3))
for (i in 1:nrow(return.time5)){
  time.result <- results5[[i]][[2]]$time[-c(1:100)]
  pCO2.result <- results5[[i]][[2]]$pCO2[-c(1:100)]
  d13C.result <- results5[[i]][[2]]$d13C[-c(1:100)]
  max.q.result <- results5[[i]][[2]]$q[-c(1:100)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result,max.q.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C","q.max")
  return.time5[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<pCO2.i+pCO2.recov)]) 
  return.time5[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C,1)==0)]) 
  return.time5[i,3] <- max(pCO2.result$q.max)
}
return.time5 <- as_tibble(cbind(as.data.frame(return.time5),q.sens.vec))
colnames(return.time5) <- c("pCO2","d13C","q","Q")

# Seafloor Weathering
return.time10 <- array(dim=c(length(results10),3))
for (i in 1:nrow(return.time10)){
  time.result <- results10[[i]][[2]]$time[-c(1:100)]
  pCO2.result <- results10[[i]][[2]]$pCO2[-c(1:100)]
  d13C.result <- results10[[i]][[2]]$d13C[-c(1:100)]
  max.q.result <- results10[[i]][[2]]$q[-c(1:100)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result,max.q.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C","q.max")
  return.time10[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<pCO2.i+pCO2.recov)]) 
  return.time10[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C,1)==0)]) 
  return.time10[i,3] <- max(pCO2.result$q.max)
}
return.time10 <- as_tibble(cbind(as.data.frame(return.time10),q.sens.vec))
colnames(return.time10) <- c("pCO2","d13C","q","Q")

# One Area (ie, Box) of High Sensitivity
return.time7 <- array(dim=c(length(results7),4))
for (i in 1:nrow(return.time7)){
  time.result <- results7[[i]][[2]]$time[-c(1:100)]
  pCO2.result <- results7[[i]][[2]]$pCO2[-c(1:100)]
  d13C.result <- results7[[i]][[2]]$d13C[-c(1:100)]
  max.q.result <- results7[[i]][[2]]$q.avg[-c(1:100)]
  temp.result <- results7[[i]][[2]]$Temp[-c(1:100)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C")
  q.sens.time <- ((max.q.result[2]-results7[[i]][[2]]$q.avg[1])/(temp.result[2]-results7[[i]][[1]]['temp.a.i']))*100
  return.time7[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<results7[[i]][[2]]$pCO2[1]+pCO2.recov)]) 
  return.time7[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C,1)==0)]) 
  return.time7[i,3] <- max(max.q.result)
  return.time7[i,4] <- q.sens.time #mean(q.sens.time[c(101:2000)])
}
return.time7 <- as_tibble(as.data.frame(return.time7))
colnames(return.time7) <- c("pCO2","d13C","qmax","Q")

# One Area of High Sensitivity with high runoff sensitivity
return.time9 <- array(dim=c(length(results9),4))
for (i in 1:nrow(return.time9)){
  time.result <- results9[[i]][[2]]$time[-c(1:100)]
  pCO2.result <- results9[[i]][[2]]$pCO2[-c(1:100)]
  d13C.result <- results9[[i]][[2]]$d13C[-c(1:100)]
  max.q.result <- results9[[i]][[2]]$q.avg[-c(1:100)]
  temp.result <- results9[[i]][[2]]$Temp[-c(1:100)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C")
  q.sens.time <- ((max.q.result[2]-results9[[i]][[2]]$q.avg[1])/(temp.result[2]-p['temp.a.i']))*100
  return.time9[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<results9[[i]][[2]]$pCO2[1]+pCO2.recov)]) 
  return.time9[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C,1)==0)]) 
  return.time9[i,3] <- max(max.q.result)
  return.time9[i,4] <- q.sens.time #mean(q.sens.time[c(101:2000)])
}
return.time9 <- as_tibble(as.data.frame(return.time9))
colnames(return.time9) <- c("pCO2","d13C","qmax","Q")

return.time.p <- ggplot() +
  geom_vline(xintercept = 0,size=1,linetype=3) +
  geom_hline(yintercept = 10,size=1,linetype=3) +
  geom_vline(xintercept = -10,size=1,linetype=2) +
  
  geom_polygon(aes(x=c(2.2,4,4,2.2),y=c(-3,-3,43,43)),fill="gray90") +
  
  geom_line(data=filter(return.time2.mat,pCO2 < 5e6),aes(x=Q,y=pCO2/1e5,colour=Ts/1e4,group=Ts/1e4),size=3,linetype=2) +
  
  scale_color_distiller(type = "seq",direction = -1,palette = "Greys",values = scales::rescale(return.time2.mat$Ts),
                        name="Soil Age (10^4 yr)") +
  
  # geom_line(data=filter(return.time1,pCO2 < 5e6),aes(x=Q*100,y=pCO2/1e5),size=3,linetype=1,colour="firebrick3") + # Standard
  geom_line(data=filter(return.time4,pCO2 < 5e6),aes(x=Q*100,y=pCO2/1e5),size=3,linetype=1,colour="firebrick3") + # West
  geom_line(data=filter(return.time5,pCO2 < 5e6),aes(x=Q*100,y=pCO2/1e5),size=3,linetype=1,colour="dodgerblue3") + # Chemostatic
  geom_line(data=filter(return.time10,pCO2 < 5e6),aes(x=Q*100,y=pCO2/1e5),size=3,linetype=1,colour="darkgreen") + # Seafloor
  
  # geom_line(data=return.time7,aes(x=Q,y=pCO2/1e5),size=3,linetype=1,colour="dodgerblue3") +
  # geom_line(data=return.time9,aes(x=Q,y=pCO2/1e5),size=3,linetype=1,colour="firebrick3") +
  
  scale_x_continuous(expand=expand_scale()) +
  scale_y_continuous(expand=expand_scale()) +
  ylab(expression(paste("pCO"[2]," Recovery Time (10"^5," yr)"))) +
  xlab(expression(paste(Delta,"Q (%/K)"))) +
  theme_linedraw(base_size = 20)
  
return.time.p

#### Calculate for PETM ####
gj <- as_tibble(read.csv("gutjahr.csv",header=T))
gj

gj.org <- as_tibble(read.csv("gutjahr-orgC.csv",header=T))
gj.org

# PETM Initial Conditions
salinity <- 35 # salinity in salinity units
sal <- (salinity/1000)+1 # salinity in kg/L
pCO2.i <- 1000 # ppm
pH.i <- 7.8
Mg.Ca.i <- 2 # seawater Mg/Ca ratio
Ca.i <- 20e-3 # mol/kg
Mg2.i <- Mg.Ca.i*Ca.i # mol/kg
# carb.parms.i <- carb(flag=21, var1=pCO2.i, var2=pH.i,
#                      # Mg2=Mg2.i*2, # For use with seacarbMg.R
#                      S=salinity, T=temp.o.i, P=0)
carb.parms.i <- co2sys.init(pH = pH.i, pCO2 = pCO2.i,
                       temp=temp.o.i,sal=salinity,p=1,Ca2=Ca.i,Mg=Mg2.i)
carb.parms.i
omega.i <- carb.parms.i$OmegaCalcite[1]
RCO2.i <- pCO2.i/pCO2.i
DIC.i <- carb.parms.i$DIC[1]*oceanV
Alk.i <- carb.parms.i$ALK[1]*oceanV

# Set-up PETM perturbation
pert.time <- 1e4 # time of perturbation

# No org burial
t.gj <- seq(1,round(max(gj$Time))*10,1)
C.emission.gj <- approx(x=gj$Time*10,y=gj$Emissions*1e15/12,xout=t.gj,rule=2)$y
C.d13C.gj <- approx(x=gj$Time*10,y=gj$d13C,xout=t.gj,rule=2)$y
C.orgB <- 0
C.orgB.d13C <- 0

# With org burial
t.gj.org <- seq(1,round(max(gj.org$Time))*10,1)
C.emission.gj. <- approx(x=gj.org$Time*10,y=gj.org$Emissions*1e15/12,xout=t.gj.org,rule=2)$y
C.d13C.gj <- approx(x=gj.org$Time*10,y=gj.org$d13C,xout=t.gj.org,rule=2)$y
C.orgB <- approx(x=gj.org$Time*10,y=gj.org$OrgC*1e15/12*-1,xout=t.gj.org,rule=2)$y
C.orgB.d13C <- approx(x=gj.org$Time*10,y=gj.org$OrgCd13C,xout=t.gj.org,rule=2)$y

C.increase <- c(rep(0,pert.time/dt),C.emission.gj,
                rep(0,(length(t)-pert.time/dt-length(C.emission.gj))))[1:length(t)]
C.perturbation <- approxfun(x=t,y=C.increase,method="linear",rule=2)

d13C.increase <- c(rep(0,pert.time/dt),C.d13C.gj,
                   rep(0,(length(t)-pert.time/dt-length(C.d13C.gj))))[1:length(t)]
d13C.perturbation <- approxfun(x=t,y=d13C.increase,method="linear",rule=2)

Corg.increase <- c(rep(0,pert.time/dt),C.orgB,
                   rep(0,(length(t)-pert.time/dt-length(C.orgB))))[1:length(t)]
Corg.perturbation <- approxfun(x=t,y=Corg.increase,method="linear",rule=2)

d13C.Corg.increase <- c(rep(0,pert.time/dt),C.orgB.d13C,
                   rep(0,(length(t)-pert.time/dt-length(C.orgB.d13C))))[1:length(t)]
d13C.Corg.perturbation <- approxfun(x=t,y=d13C.Corg.increase,method="linear",rule=2)

y0  <- c(DIC=DIC.i, Alk=Alk.i, d13C=d13C.i, S=S.i)
names(y0) <- c(names(y0)[1:3],"S")

# Test PETM version of the model
p <- paramfunc()
p <- c(p,q.sens=0.04)

output <- as.data.frame(lsodes(y=y0, times=t, func=CLiBeSO.PETM, parms=p,rtol=1e-6))
colnames(output) <- c('time','DIC', 'Alk','d13C',"S",
                      'pCO2','omega', 'Fwcarb',
                      'Fbcarb', 'Fborg','Fwsil',
                      'RCO2','Temp','Fworg',"Fwsulf","Fbsulf","q","Conc","Dw")

y.time <- output$time

par(mfrow=c(4,1),mar=c(1,5,1,1))
plot(x=y.time/1e6,y=output$pCO2,bty="n",type="l", # CO2
     ylab=expression(paste("Atmospheric CO"[2])),
     xaxt="n")
plot(x=y.time,y=output$d13C,bty='n',type='l',xaxt="n",ylab=expression(paste(delta^13,"C"))) # d13C
plot(x=y.time,y=output$q,bty="n",type="l",xaxt="n",ylab="Runoff (m/yr)") # Runoff
par(mar=c(5,5,1,1))
plot(x=y.time,y=output$Conc,bty='n',type='l',ylab=expression(paste("[C]"[sil])),
     xlab="Time (yr)") # Concentration

#### Model Run ####
Ts.vec <- seq(1e4,5e4,1e4)
PETM.results <- list()

for (j in 1:length(Ts.vec)) {
  p <- paramfunc(remove = c("Ts","init"))
  p <- c(p,Conc0=Conc0.in,q.sens=0.04,Ts=Ts.vec[j],init=TRUE)
  p

  init2.results <- as.data.frame(lsodes(y=y0, times=t.init, func=CLiBeSO, parms=p,rtol=1e-6))
  colnames(init2.results) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb','Fborg','Fwsil','RCO2','Temp','Fworg',"Fwsulf","Fbsulf","q.1","Conc.1","Dw.1")

  temp.o.i <- 15#tail(init2.results$Temp)[5]
  pCO2.i <- tail(init2.results$pCO2)[5]
  y0  <- c(DIC=tail(init2.results$DIC)[5], Alk=tail(init2.results$Alk)[5], 
           d13C=tail(init2.results$d13C)[5], S=tail(init2.results$S)[5])
  names(y0) <- c(names(y0)[1:3],"S")
  new.q0 <- tail(init2.results$q.1)[5]

  p <- paramfunc(remove=c("Ts","temp.a.i","q0"))
  q.sens.vec <- seq(-0.1,0.2,0.04)
  p.add <- list()
    for (i in 1:length(q.sens.vec)){
      p.add[[i]] <- c(p,q.sens=q.sens.vec[i],Ts=Ts.vec[j],
                      Conc0=tail(init2.results$Conc.1)[5],q0=new.q0,
                      temp.a.i=tail(init2.results$Temp)[5])
    }

#### Parallel Execution ####
  detectCores()
  parallelStart("multicore",cpus=detectCores()-3)
    PETM.results[[j]] <- parallelLapply(p.add,CLiBeSO.uber)
  parallelStop()
}

PETM.results <- flatten(PETM.results)

#### Calculate return times ####
pCO2.recov <- 1
return.time <- array(dim=c(length(PETM.results),2))
for (i in 1:nrow(return.time)){
  time.result <- PETM.results[[i]][[2]]$time[-c(1:250)]
  pCO2.result <- PETM.results[[i]][[2]]$pCO2[-c(1:250)]
  d13C.result <- PETM.results[[i]][[2]]$d13C[-c(1:250)]
  pCO2.result <- as.data.frame(cbind(time.result,pCO2.result,d13C.result))
  colnames(pCO2.result) <- c("time","pCO2","d13C")
  return.time[i,1] <- min(pCO2.result$time[which(pCO2.result$pCO2<PETM.results[[i]][[2]]$pCO2[1]+(pCO2.recov))]) 
  return.time[i,2] <- min(pCO2.result$time[which(round(pCO2.result$d13C[500:length(pCO2.result$d13C)],1)==round(PETM.results[[i]][[2]]$d13C[1],2))]) 
}
return.time <- as.data.frame(return.time)
colnames(return.time) <- c("pCO2","d13C")

q.sens.vec.2 <- numeric(length=length(PETM.results))
Ts.vec.2 <- numeric(length=length(PETM.results))
Dw.vec.2 <- numeric(length=length(PETM.results))
for (i in 1:length(q.sens.vec.2)){
  q.sens.vec.2[i] <- PETM.results[[i]][[1]]['q.sens']
  Ts.vec.2[i] <- PETM.results[[i]][[1]]['Ts']
  Dw.vec.2[i] <- PETM.results[[i]][[2]]$Dw[1]
}

return.time.mat <- as_tibble(cbind(return.time,q.sens.vec.2*100,Ts.vec.2,Dw.vec.2))
colnames(return.time.mat) <- c("pCO2","d13C","Q","Ts","Dw")

d13C.breaks <- seq(0,1e6,5e4)
d13C.labels <- c(0,rep("",3),2e5,rep("",3),4e5,rep("",3),6e5,rep("",3),
                   8e5,rep("",3),1e6)

return.time.mat %>% filter(d13C!=Inf) -> blank

RTplot <- ggplot(data=return.time.mat,aes(x=Q,y=Ts,fill=d13C)) +
  geom_raster() +
  # geom_hline(aes(yintercept=6e4),size=3) +
  # geom_hline(aes(yintercept=2e4),size=2,linetype=2) +
  # geom_vline(aes(xintercept=0),size=3) +
  geom_contour(aes(x=Q,y=Ts,z=d13C),breaks=c(2e5),colour="black",size=5) +
  # geom_contour(aes(x=Q,y=Ts,z=pCO2),breaks=c(1e6),colour="black",size=5) +
  
  scale_x_continuous(expand=expand_scale()) +
  scale_y_continuous(expand=expand_scale()) +
  # scale_fill_continuous(name=expression(paste("Recovery\nTime (yr)"))) +
  scale_fill_viridis_c(name=expression(paste("Recovery\nTime (yr)")),
                       trans="log") + # ,breaks=d13C.breaks,labels=d13C.labels) +
  # labs(title=expression(paste(delta^13,"C"))) +
  # labs(title="pCO2") +
  xlab("Runoff Sensitivity (%/K)") + 
  # ylab("Concentration Sensitivity (%/K)") +
  ylab("Global Average Weathering Zone Age (yrs)") +
  # ylab("Concentration Sensitivity (%/CO2 doubling)") +
  theme_linedraw(base_size = 20)
RTplot

return.time.mat %>% filter(d13C!=Inf) -> return.time.mat.2

RTplot.2 <- ggplot(data=return.time.mat,aes(x=Q,y=d13C/1e5)) +
  geom_hline(yintercept = 1.4,linetype="dotted") +
  geom_hline(yintercept = 2,linetype="dotted") +
  geom_vline(xintercept = 0,linetype="dashed") +
  
  geom_line(aes(color=Ts/1e4,group=Ts/1e4),size=4) +
  
  scale_color_distiller(type = "seq",direction = -1,palette = "Greys",values = scales::rescale(return.time.mat.2$Ts),
                        name="Soil Age (10^4 yr)") +
  
  scale_y_continuous(trans="log10") +
  xlab("Runoff Sensitivity (%/K)") + ylab(expression(paste(delta^13,"C Recovery Time (10"^5," yr)"))) +
  theme_linedraw(base_size = 20) +
  theme(panel.grid = element_blank())
RTplot.2
