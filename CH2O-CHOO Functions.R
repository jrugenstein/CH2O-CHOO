# Functions used in the CH2O-CHOO Model. Please load these functions first. 

# For more detailed explanation, please review Rugenstein and Winkler (in review). 

# Jeremy Rugenstein acknowledges Kimberly Lau and Adam Jost for help in 
# using package('deSolve').

# Created on: 04-03-2020
# Last Modified: 14-06-2022
# Last Modified by: Jeremy K. C. Rugenstein (CSU)

#### Model Parameter Function ####
paramfunc <- function(remove=NA) {
  p <- c( #parameter vector
    
    # Initialization Run?
    init = FALSE,
    
    # Physical Constants
    Clim.Sens = 4, # Earth System Sensitivity (°C/CO2 doubling) (Knuti et al., 2017--Nature Geoscience)
    temp.a.i = temp.o.i, # Temperature of the atmosphere
    
    # Initial Sulfur fluxes
    Fwsulf.i = 0.5e12, # from Calmels et al. 2007 (Geology)
    Fbsulf.i = 0.5e12, # To balance the initial sulfur cycle
    
    # Initial Carbon and Alkalinity fluxes (mol/yr)
    Fvolc.i = 6e12, # solid Earth degassing CO2
    Fworg.i = 6e12, # organic C weathering
    Fwcarb.i = 12e12, # carbonate weathering
    Fborg.i = 6e12, # burial of organic carbon
    Fbcarb.i = 18e12, # carbonate burial
    Fwsil.i = 6e12, # silicate weathering flux
    
    # Initial d13C Values
    Fwcarb.d13C = 0,
    Fworg.d13C = -22,
    Fvolc.d13C = -5, 
    DB = 27, # Cap-delta between inorganic carbon and organic carbon
    
    # Phosphorus Parameters (default parameters from Shield and Mills 2017)
    Fwp.i=4.7e10,
    Fbp.i=4.7e10,
    sil.P = 0.58, # P derived from silicate weathering
    carb.P = 0.21, # P derived from carbonate weathering
    org.P = 0.21, # P derived from organic carbon weathering
    
    # Other scaling exponents and factors
    borg.carb = 1, # Scaling between organic carbon burial and rain rate
    carb.clim.exp = 0.05, # Exponent carbonate weathering Arrhenius reaction (not implemented as a changeable parameter)
    Rk.i = 1/(log2(pCO2.i/280)+1), # Silicate reactivity of land surface (k in Eq. 3) (Caves et al., 2016--EPSL)
    Trop.area = 1/3, # Effective "area" for tropical weathering
    Conc.sens = 0.01, # (% change per K)--only used in CLiBeSO.qC
    q0 = 0.3, # Initial runoff
    
    # Parameters used in MC 2014 Function
    keff.ref = 8.7e-6, # mol/m2/yr (reference reaction rate) (from Maher and Chamberlain SI Table S1)
    m = 270, # g/mol (molar mass) (from Maher and Chamberlain SI Table S1)
    A = 0.1, # m2/g (specific surface area) (from Maher and Chamberlain SI Table S1)
    Rn.max.ref = 1085, # umol SiO2/L/yr # Maher and Chamberlain 2014
    L.phi = 0.1, # Reactive length scale times effective porosity (Maher and Chamberlain 2014)
    Ts = 2e4, # Timescale for Solids (Maher and Chamberlain 2014)
    Ceq = 374, # Equilibrium Concentration for SiO2 (Maher and Chamberlain 2014)
    soil.CO2.flag = 1, # Flag to activate soil CO2 module
      # 1: Soil CO2 calculated using Volk 1987
      # 2: Soil CO2 equals atmospheric CO2
      # 3: Soil CO2 equals atmospheric CO2 multiplied by magnification factor (soil.CO2.mag)
    soil.CO2.mag = 10, # Soil CO2 magnification factor
    Ea = 38, # (kJ/mol) Activation Energy (Maher and Chamberlain 2014)
    
    # Parameters used in West 2012 weathering function
    K.west = 2.6e-4, # Range: 7.6e-6 to 1.2e-3
    Kw.west = 7.6e-5, # Range: 1.5e-6 to 3.0e-4
    z.west = 8.9, # Range: 0.26 to 41 t/m2
    Xm.west = 0.09, # Range: 0.04 to 0.13
    sigma.west = 0.89 # Range: 0.66 to 1.13
  )
  
  if (is.na(remove)==TRUE) {return.p <- p}
  else {for (i in 1:length(remove)){
    p <- p[-(which(names(p)==c(remove[i])))]}
    return.p <- p
  }
  return.p
}

#### Carbonate Solvers ####
# From Adam Jost

# Currently not used
co2sys.old <- function(DIC, TA, temp=25, sal=35, p=1, Ca2=0.01028, Mg=0.05654){
  
  Mg.Ca <- Mg/Ca2
  
  t=273.15+temp #convert to K
  k0=0.02839188 #from seacarb, assumes sal=35, temp=25, p=0
  
  #calculate k1, k2, kw, and kb based on the method by DOE 1994 (as shown in Zeebe and Wolf-Gladrow)
  k1.1 = exp(2.83655 - 2307.1266/t - 1.5529413*log(t) - (0.207608410 + 4.0484/t)*sqrt(sal) + 0.0846834*sal - 0.00654208*sal^(3/2) + log(1-0.001005*sal))
  k1 <- k1.1*((0.15505*(Mg-0.05654)))+k1.1
  
  k2.1 = exp(-9.226508 - 3351.6106/t - 0.2005743*log(t) - (0.106901773 + 23.9722/t)*sqrt(sal) + 0.1130822*sal - 0.00846934*sal^(3/2) + log(1 - 0.001005*sal))
  k2 <- k2.1*((0.44224*(Mg-0.05654)))+k2.1
  
  kw = exp(148.96502 - 13847.26/t - 23.6521*log(t) + (118.67/t - 5.977 + 1.0495*log(t))*sal^(1/2) - 0.01615*sal)
  kb = exp((-8966.90 - 2890.53*sal^(1/2) - 77.942*sal + 1.728*sal^(3/2) - 0.0996*sal^2)/t + 148.0248 +137.1942*sal^(1/2) + 1.62142*sal - (24.4344 + 25.085*sal^(1/2) + 0.2474*sal)*log(t) + 0.053105*sal^(1/2)*T)
  
  bt = 4.16e-4*sal/35 #eq. A.7.15 in Zeebe and Wolf-Gladrow 2001
  
  #calculate ksp for calcite and aragonite from Muci (1983) as described in Zeebe and Wolf-Gladrow, 2001
  ksp.c.1 = 10^(-171.9065 - 0.077993*t + 2839.319/t + 71.595*log10(t) + (-0.77712 + 0.0028426*t + 178.34/t)*sal^(1/2) - 0.07711*sal + 0.0041249*sal^1.5)
  ksp.c <- ksp.c.1 - 3.655e-8*(5.14-Mg.Ca)
  #ksp.a
  
  #creating a function to solve for [H+] (referred to as 'h' here) given DIC and TA. Refer to Zeebe and Wolf-Gladrow 2001 for instructions
  h.solve = function(DIC, TA){
    #polynomials
    a = -1
    b = -TA - kb - k1
    c = DIC*k1 - TA*kb - TA*k1 + kb*bt + kw - kb*k1 - k1*k2
    d = DIC*kb*k1 + 2*DIC*k1*k2 - TA*kb*k1 - TA*k1*k2 + kb*bt*k1 + kw*kb + kw*k1 - kb*k1*k2 
    e = 2*DIC*kb*k1*k2 - TA*kb*k1*k2 + kb*bt*k1*k2 + kw*kb*k1 + kw*k1*k2
    f = kw*kb*k1*k2
    p = c(f,e,d,c,b,a)
    r = polyroot(p)
    h=max(Re(r))
    return(h)
  }
  
  h = h.solve(DIC, TA)
  pH = -log10(h)
  CO2 = DIC/(1 + k1/h + k1*k2/(h^2))
  HCO3 = DIC/(1 + h/k1 + k2/h)
  CO32 = DIC/(1 + h/k2 + (h^2)/k1*k2)
  OmegaCalcite = Ca2*CO32/ksp.c
  #omega.a = Ca2*CO32/ksp.a
  pCO2=CO2/k0*1e6
  list = as.data.frame(pH)
  list$DIC=DIC; list$ALK=TA; list$CO2=CO2; list$HCO3=HCO3; list$CO32=CO32; list$OmegaCalcite=OmegaCalcite; list$pCO2=pCO2#; list$omega.a=omega.a
  list = as.data.frame(list)
  list
}

# Used only to set initial conditions
co2sys.init <- function(pH, pCO2, temp=25, sal=35, p=1, Ca2=0.01028, Mg=0.05654, SO42=0.0282){
  
  # Inputs should be in mol/kg
  
  Mg.Ca <- Mg/Ca2
  
  # From Wikipedia
  Ca2.m <- 0.01028
  Mg.m <- 0.05654
  SO42.m <- 0.0282
  
  # From Zeebe and Tyrrell (2019; GCA)--Table 2
  Ca.K1 <- 5e-3
  Mg.K1 <- 17e-3
  SO42.K1 <- 208e-3
  
  Ca.K2 <- 157e-3
  Mg.K2 <- 420e-3
  SO42.K2 <- 176e-3
  
  Ca.Ksp <- 185e-3
  Mg.Ksp <- 518e-3
  SO42.Ksp <- 106e-3
  
  t=273.15+temp #convert to K
  k0=0.02839188 #from seacarb, assumes sal=35, temp=25, p=0
  
  #calculate k1, k2, kw, and kb based on the method by DOE 1994 (as shown in Zeebe and Wolf-Gladrow)
  k1.m <- exp(2.83655 - 2307.1266/t - 1.5529413*log(t) - (0.207608410 + 4.0484/t)*sqrt(sal) + 0.0846834*sal - 0.00654208*sal^(3/2) + log(1-0.001005*sal))
  k1 <- k1.m*(1+((Ca.K1*(Ca2/Ca2.m-1))+(Mg.K1*(Mg/Mg.m-1))+(SO42.K1*(SO42/SO42.m-1))))
  
  k2.m = exp(-9.226508 - 3351.6106/t - 0.2005743*log(t) - (0.106901773 + 23.9722/t)*sqrt(sal) + 0.1130822*sal - 0.00846934*sal^(3/2) + log(1 - 0.001005*sal))
  k2 <- k2.m*(1+((Ca.K2*(Ca2/Ca2.m-1))+(Mg.K2*(Mg/Mg.m-1))+(SO42.K2*(SO42/SO42.m-1))))
  
  kw = exp(148.96502 - 13847.26/t - 23.6521*log(t) + (118.67/t - 5.977 + 1.0495*log(t))*sal^(1/2) - 0.01615*sal)
  kb = exp((-8966.90 - 2890.53*sal^(1/2) - 77.942*sal + 1.728*sal^(3/2) - 0.0996*sal^2)/t + 148.0248 +137.1942*sal^(1/2) + 1.62142*sal - (24.4344 + 25.085*sal^(1/2) + 0.2474*sal)*log(t) + 0.053105*sal^(1/2)*T)
  
  bt = 4.16e-4*sal/35 #eq. A.7.15 in Zeebe and Wolf-Gladrow 2001
  
  #calculate ksp for calcite and aragonite from Muci (1983) as described in Zeebe and Wolf-Gladrow, 2001
  ksp.c.m = 10^(-171.9065 - 0.077993*t + 2839.319/t + 71.595*log10(t) + (-0.77712 + 0.0028426*t + 178.34/t)*sal^(1/2) - 0.07711*sal + 0.0041249*sal^1.5)
  ksp.c <- ksp.c.m*(1+((Ca.Ksp*(Ca2/Ca2.m-1))+(Mg.Ksp*(Mg/Mg.m-1))+(SO42.Ksp*(SO42/SO42.m-1))))
  
  h <- 10^(-pH) # convert to [H+]
  CO2 <- pCO2*k0/1e6
  DIC <- CO2*(1+k1/h+k1*k2/h^2)
  HCO3 <- DIC/(1+h/k1+k2/h)
  CO32 <- DIC/(1+h/k2+h^2/k1*k2)
  OmegaCalcite = Ca2*CO32/ksp.c
  TA <- CO2*(k1/h+2*k1*k2/h^2)+bt*kb/(kb+h)+kw/h-h
  
  # list <- c(DIC=DIC,ALK=TA,CO2=CO2,HCO3=HCO3,CO32=CO32,OmegaCalcite=OmegaCalcite)
  results <- as.data.frame(array(data=c(DIC,TA,CO2,HCO3,CO32,OmegaCalcite),dim=c(1,6)))
  colnames(results) <- c("DIC","ALK","CO2","HCO3","CO32","OmegaCalcite")
  results
}

# Used within the CLiBeSO functions
co2sys <- function(DIC, TA, temp=25, sal=35, p=1, Ca.in=0.01028, Mg.in=0.05654, 
                   SO42.in=0.0282){
  
  # Inputs should be in mol/kg
  
  Mg.Ca <- Mg.in/Ca.in
  
  # From Wikipedia
  Ca.m <- 0.01028
  Mg.m <- 0.05654
  SO42.m <- 0.0282
  
  # From Zeebe and Tyrrell (2019; GCA)--Table 2
  Ca.K1 <- 5e-3
  Mg.K1 <- 17e-3
  SO42.K1 <- 208e-3
  
  Ca.K2 <- 157e-3
  Mg.K2 <- 420e-3
  SO42.K2 <- 176e-3
  
  Ca.Ksp <- 185e-3
  Mg.Ksp <- 518e-3
  SO42.Ksp <- 106e-3
  
  t=273.15+temp #convert to K
  k0=0.02839188 #from seacarb, assumes sal=35, temp=25, p=0
  
  #calculate k1, k2, kw, and kb based on the method by DOE 1994 (as shown in Zeebe and Wolf-Gladrow)
  k1.m <- exp(2.83655 - 2307.1266/t - 1.5529413*log(t) - (0.207608410 + 4.0484/t)*sqrt(sal) + 0.0846834*sal - 0.00654208*sal^(3/2) + log(1-0.001005*sal))
  k1 <- k1.m*(1+((Ca.K1*(Ca.in/Ca.m-1))+(Mg.K1*(Mg.in/Mg.m-1))+(SO42.K1*(SO42.in/SO42.m-1))))
  
  k2.m = exp(-9.226508 - 3351.6106/t - 0.2005743*log(t) - (0.106901773 + 23.9722/t)*sqrt(sal) + 0.1130822*sal - 0.00846934*sal^(3/2) + log(1 - 0.001005*sal))
  k2 <- k2.m*(1+((Ca.K2*(Ca.in/Ca.m-1))+(Mg.K2*(Mg.in/Mg.m-1))+(SO42.K2*(SO42.in/SO42.m-1))))
  
  kw = exp(148.96502 - 13847.26/t - 23.6521*log(t) + (118.67/t - 5.977 + 1.0495*log(t))*sal^(1/2) - 0.01615*sal)
  kb = exp((-8966.90 - 2890.53*sal^(1/2) - 77.942*sal + 1.728*sal^(3/2) - 0.0996*sal^2)/t + 148.0248 +137.1942*sal^(1/2) + 1.62142*sal - (24.4344 + 25.085*sal^(1/2) + 0.2474*sal)*log(t) + 0.053105*sal^(1/2)*T)
  
  bt = 4.16e-4*sal/35 #eq. A.7.15 in Zeebe and Wolf-Gladrow 2001
  
  #calculate ksp for calcite and aragonite from Muci (1983) as described in Zeebe and Wolf-Gladrow, 2001
  ksp.c.m = 10^(-171.9065 - 0.077993*t + 2839.319/t + 71.595*log10(t) + (-0.77712 + 0.0028426*t + 178.34/t)*sal^(1/2) - 0.07711*sal + 0.0041249*sal^1.5)
  ksp.c <- ksp.c.m*(1+((Ca.Ksp*(Ca.in/Ca.m-1))+(Mg.Ksp*(Mg.in/Mg.m-1))+(SO42.Ksp*(SO42.in/SO42.m-1))))
  
  #creating a function to solve for [H+] (referred to as 'h' here) given DIC and TA. Refer to Zeebe and Wolf-Gladrow 2001 for instructions
  h.solve = function(DIC, TA){
    #polynomials
    a = -1
    b = -TA - kb - k1
    c = DIC*k1 - TA*kb - TA*k1 + kb*bt + kw - kb*k1 - k1*k2
    d = DIC*kb*k1 + 2*DIC*k1*k2 - TA*kb*k1 - TA*k1*k2 + kb*bt*k1 + kw*kb + kw*k1 - kb*k1*k2 
    e = 2*DIC*kb*k1*k2 - TA*kb*k1*k2 + kb*bt*k1*k2 + kw*kb*k1 + kw*k1*k2
    f = kw*kb*k1*k2
    p = c(f,e,d,c,b,a)
    r = polyroot(p)
    h=max(Re(r))
    return(h)
  }
  
  h = h.solve(DIC, TA)
  pH = -log10(h)
  CO2 = DIC/(1 + k1/h + k1*k2/(h^2))
  HCO3 = DIC/(1 + h/k1 + k2/h)
  CO32 = DIC/(1 + h/k2 + (h^2)/k1*k2)
  OmegaCalcite = Ca.in*CO32/ksp.c
  #omega.a = Ca2*CO32/ksp.a
  pCO2=CO2/k0*1e6
  
  results <- as.data.frame(array(data=c(DIC,TA,CO2,HCO3,CO32,OmegaCalcite,pH,pCO2),dim=c(1,8)))
  colnames(results) <- c("DIC","ALK","CO2","HCO3","CO32","OmegaCalcite","pH","pCO2")
  results
}

# Initial Carbonate System Conditions Function

init.carb <- function(){
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
  omega.i <- carb.parms.i$OmegaCalcite
  RCO2.i <- pCO2.i/pCO2.i
  DIC.i <- carb.parms.i$DIC*oceanV
  Alk.i <- carb.parms.i$ALK*oceanV

# Initial Reservoirs
  S.i <- 10e6*p['Fwsulf.i'] # mol sulfur; assumes a 10 Ma residence time of S in the ocean; modern-day estimates are 5e19 (Broecker and Peng 1982), but this includes a substantial gypsum component
  PO4.i <- 3.1e15 # mol P

  y0  <- c(DIC=DIC.i, Alk=Alk.i, d13C=d13C.i, S=S.i
           # ,P=P.i
  )
  names(y0) <- c(names(y0)[1:3],"S"
                 # ,"P"
  )

  carb.init.output <- list(pCO2.i,y0)
  names(carb.init.output) <- c("pCO2.i","y0")
  carb.init.output
}


#### Weathering Functions ####
soilCO2func <- function(pCO2,prod.max=2,pCO2.min=100) {
  
  # From Tyler Volk (1989)
  
  prod0 <- 1
  pCO2.soil.i <- pCO2.i*10
  pCO2.50 <- (prod.max/prod0-1)*(pCO2.i-pCO2.min)
  prod <- prod.max*(pCO2-pCO2.min)/(pCO2.50+(pCO2-pCO2.min))
  pCO2.soil <- (prod/prod0*(1-pCO2.i/pCO2.soil.i)+pCO2/pCO2.soil.i)*pCO2.soil.i+(pCO2-pCO2.i)
  RCO2.soil <- pCO2.soil/pCO2.soil.i
  
  RCO2.soil
}

# Arrhenius Function
arrhenius <- function(Temp.C,Ea.in,keff.ref.in,Ts.in){
  m <- p['m'] # g/mol (molar mass) (from Maher and Chamberlain SI Table S1)
  A <- p['A'] # m2/g (specific surface area) (from Maher and Chamberlain SI Table S1)
  R <- 8.314 # J/K/mol
  Ea.used <- Ea.in*1000
  
  # mkeffa <- m*A*keff.ref.in*((Ea.used/R)*((1/(Temp.C+273.15))-(1/(p['temp.a.i']+273.15))))
  # mkeffa <- m*A*keff.ref.in*exp((Ea.used/R)*((1/(p['temp.a.i']+273.15))-(1/(Temp.C+273.15))))
  mkeffa <- exp((Ea.used/R)*((1/(15+273.15))-(1/(Temp.C+273.15))))
  # mkeffa <- exp((Ea.used/R)*exp(-keff.ref.in*Ts.in)*((1/(p['temp.a.i']+273.15))-(1/(Temp.C+273.15))))
  mkeffa
}

arrhenius.grp <- function(Temp.C,Ea.in,keff.ref.in){ # Arrhenius from Graham and Pierrehumbert (2020)
  R <- 8.314 # J/K/mol
  Ea.used <- Ea.in*1000 # J/mol
  Temp.K <- Temp.C+273.15
  Temp.ref <- p['temp.a.i']+273.15
  
  Te <- (Temp.ref^2)*R/Ea.used
  
  keff <- keff.ref.in*exp((Temp.K-Temp.ref)/Te)
  keff
}

# Maher and Chamberlain (2014) equation
MCfunc <- function(pCO2,q0,q,Ea=p['Ea'],Ts=p['Ts'],keff.ref=p['keff.ref'],
                   Temperature=Temp,soilCO2.use=p['soil.CO2.flag'],
                   mag.factor=p['soil.CO2.mag']){
  
  # Ea is input in kJ/mol
  
  RCO2 <- pCO2/pCO2.i
  
  if (soilCO2.use==1){
    RCO2.soil <- soilCO2func(pCO2=pCO2)
  }
  if (soilCO2.use==2){
    RCO2.soil <- RCO2
  }
  if (soilCO2.use==3){
    RCO2.soil <- RCO2*mag.factor
  }
  
  m.init <- p['m'] # g/mol (molar mass) (from Maher and Chamberlain SI Table S1)
  A.init <- p['A'] # m2/g (specific surface area) (from Maher and Chamberlain SI Table S1)
  
  Rn.max <- p['Rn.max.ref']*(keff.ref/8.7e-6) # Maximum net reaction rate scaled to keff (umol/L/yr); 8.7e-6 is the reference reaction rate in MC14

  # Calculation Initial Concentration
  Dw0 <- ((p['L.phi']*Rn.max*(1/(1+m.init*A.init*keff.ref*Ts)))/p['Ceq'])
  Conc0 <- p['Ceq']*(((exp(1)^2)*Dw0/q0)/(1+(exp(1)^2)*Dw0/q0))
  
  # mkeffA <- arrhenius(Temp.C=Temperature,Ea.in=Ea,keff.ref.in=keff.ref)
  Ceq <- ((RCO2.soil)^(0.316))*p['Ceq'] #1000
  # fw <- 1/(1+mkeffA*Ts)
  fw <- 1/(1+p['m']*p['A']*keff.ref*Ts)
  Dw <- (p['L.phi']*Rn.max*fw/Ceq)*arrhenius(Temp.C=Temperature,Ea.in=Ea,keff.ref.in=keff.ref,Ts.in = Ts)
  Conc <- Ceq*(((exp(1)^2)*Dw/q)/(1+(exp(1)^2)*Dw/q))
  
  results <- c(q0,q,Conc0,Conc,Ceq,Dw,RCO2.soil)
  names(results) <- c("q0","q","Conc0","Conc","Ceq","Dw","RCO2.soil")
  results <- as.data.frame(t(as.matrix(results,nrow=1,ncol=length(results))))
  results
}

# West (2012) equation
west.W <- function(q,Temperature,Temperature.0=15,#E.rate=100){
                   E.rate = 1/p['Ts']){
  # q in m/yr; Temperature in °C; E.rate in m/yr
  
  q.used <- q*1000 # Convert to mm/yr
  Temp.K <- Temperature+273.15 # K
  Temp.K0 <- Temperature.0+273.15
  erosion <- E.rate*1e2*2.8*1e10/1e6 # tons/km2/yr
  Ea.used <- p['Ea']*1000 # J/mol
  R <- 8.314 # J/K/mol (Universal gas constant)
  
  # Convert z from t/m2 to t/km2
  z <- p['z.west']*1e6
  
  runoff.term <- 1-exp(-p['Kw.west']*q.used)
  Arrhenius.term <- exp((Ea.used/(R*Temp.K0))-(Ea.used/(R*Temp.K)))
  z.term <- ((z/erosion)^p['sigma.west'])/p['sigma.west']
  exp.term <- 1-exp(-p['K.west']*runoff.term*Arrhenius.term*z.term)
  Weathering.flux <- erosion*p['Xm.west']*exp.term
  Weathering.flux # tons/km2/yr
}

#### Base CLiBeSO Functions ####
CLiBeSO = function(t,y,p){
  with(as.list(c(p,y)),{
    
    if (p['init']==FALSE) {
      Fextra <- C.perturbation(t)
      d13C.Fextra <- d13C.perturbation(t)
    }

    if (p['init']==TRUE) {
      Fextra <- C.perturbation.init(t)
      d13C.Fextra <- d13C.perturbation(t)
    }
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)

    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- p['q0']
    q <- q0*(1+p['q.sens'])^(Temp-p['temp.a.i'])
    
    # Concentrations
    # Conc <- Conc0*(1+p['Conc.sens'])^(Temp-p['temp.a.i'])
    # Conc <- (log2(RCO2)*p['Conc.sens']*Conc0)+Conc0 #(RCO2*(1+p['Conc.sens']))/RCO2
    # Conc <- exp(p['Ea']/R*(1/(p['temp.a.i']+273.15)-1/(Temp+273.15))) # Activation Energy
    
    # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
    # open-system calculation for Ceq
    Conc.results <- MCfunc(pCO2 = pCO2.i,q0 = q0, q = q, Ea = p['Ea'], Ts=p['Ts'],
                           keff.ref=p['keff.ref'],Temperature=Temp,
                           soilCO2.use=p['soil.CO2.flag'])
    Conc0 <- p['Conc0']
    # Conc0 <- Conc.results$Conc0
    Conc <- Conc.results$Conc
    Dw <- Conc.results$Dw

    # Silicate Weathering
    Fwsil <- p['Fwsil.i']*q/q0*Conc/Conc0

    # Carbonate Fluxes
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

CLiBeSO.3 = function(t,y,p){ # This function contains three land "boxes"
  with(as.list(c(p,y)),{

    if (p['init']==FALSE) {
      Fextra <- C.perturbation(t)
      d13C.Fextra <- d13C.perturbation(t)
    }
    
    if (p['init']==TRUE) {
      Fextra <- C.perturbation.init(t)
      d13C.Fextra <- d13C.perturbation(t)
    }
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)
    
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # BOX 1--Polar/Low weatherability (Ts)
      area.1 <- (1-p['Trop.area'])/2
      
      Temp.1 <- Temp
      
      # Runoff    
      q0.1 <- p['q01']
      q.1 <- q0.1*(1+p['q.sens1'])^(Temp-p['temp.a.i'])
      
      # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
      # open-system calculation for Ceq
      Conc.results.1 <- MCfunc(pCO2 = pCO2.i,q0 = q0.1, q = q.1, Ea = p['Ea'], Ts=p['Ts1'],
                               keff.ref=p['keff.ref'],Temperature = Temp.1,
                               soilCO2.use=p['soil.CO2.flag'])
      Conc0.1 <- p['Conc01'] #Conc.results.1$Conc0
      Conc.1 <- Conc.results.1$Conc
      Dw.1 <- Conc.results.1$Dw
      
      # Silicate and Carbonate Weathering
      Fwsil.1 <- p['Fwsil.i']*area.1*q.1/q0.1*Conc.1/Conc0.1
      Fwcarb.1 <- p['Fwcarb.i']*area.1*q.1/q0.1*Conc.1/Conc0.1

    # BOX 2--Mid-latitudes/Medium weatherability (Ts)
      area.2 <- (1-p['Trop.area'])/2
      
      Temp.2 <- Temp
      
      # Runoff    
      q0.2 <- p['q02']
      q.2 <- q0.2*(1+p['q.sens2'])^(Temp-p['temp.a.i'])
      
      # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
      # open-system calculation for Ceq
      Conc.results.2 <- MCfunc(pCO2 = pCO2.i,q0 = q0.2, q = q.2, Ea = p['Ea'], Ts=p['Ts2'], Temperature = Temp.2,
                               soilCO2.use=p['soil.CO2.flag'])
      Conc0.2 <- p['Conc02'] #Conc.results.2$Conc0
      Conc.2 <- Conc.results.2$Conc
      Dw.2 <- Conc.results.2$Dw
      
      # Silicate and Carbonate Weathering
      Fwsil.2 <- p['Fwsil.i']*area.2*q.2/q0.2*Conc.2/Conc0.2
      Fwcarb.2 <- p['Fwcarb.i']*area.2*q.2/q0.2*Conc.2/Conc0.2
      
    # BOX 3--Tropics/High weatherability (Ts)
      area.3 <- p['Trop.area']
      
      Temp.3 <- Temp
      
      # Runoff    
      q0.3 <- p['q03']
      q.3 <- q0.3*(1+p['q.sens3'])^(Temp-p['temp.a.i'])
      
      # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
      # open-system calculation for Ceq
      Conc.results.3 <- MCfunc(pCO2 = pCO2.i,q0 = q0.3, q = q.3, Ea = p['Ea'], Ts=p['Ts3'], Temperature = Temp.3,
                               soilCO2.use=p['soil.CO2.flag'])
      Conc0.3 <- p['Conc03'] #Conc.results.3$Conc0
      Conc.3 <- Conc.results.3$Conc
      Dw.3 <- Conc.results.3$Dw
      
      # Silicate and Carbonate Weathering
      Fwsil.3 <- p['Fwsil.i']*area.3*q.3/q0.3*Conc.3/Conc0.3
      Fwcarb.3 <- p['Fwcarb.i']*area.3*q.3/q0.3*Conc.3/Conc0.3

    # Sum weathering fluxes
    Fwsil <- Fwsil.1 + Fwsil.2 + Fwsil.3
    Fwcarb <- Fwcarb.1 + Fwcarb.2 + Fwcarb.3
    
    # Average global q (assumes each "box" is 1/3 of total area)
    q.avg <- q.1*area.1 + q.2*area.2 + q.3*area.3
    
    # Weighted average global [C] (weighted by 'Trop.area')
    Conc.avg <- Conc.1*area.1 + Conc.2*area.2 + Conc.3*area.3

    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil, Fwsil.1,Fwsil.3,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q.1, q.2, q.3, q.avg, Conc.1, 
         Conc.2, Conc.3, Conc.avg, Dw.1, Dw.2, Dw.3)
  })
}

CLiBeSO.chem = function(t,y,p){
  with(as.list(c(p,y)),{
    Fextra <- C.perturbation(t)
    d13C.Fextra <- d13C.perturbation(t)
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)
    
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- 0.3
    # q <- (log2(RCO2)*p['q.sens']*q0)+q0
    q <- q0*(1+p['q.sens'])^(Temp-p['temp.a.i'])
    
    # Concentrations
    Conc0 <- 1
    Conc <- 1
    Dw <- 1

    # Silicate Weathering
    Fwsil <- p['Fwsil.i']*q/q0*Conc/Conc0

    # Carbonate Fluxes
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

CLiBeSO.west = function(t,y,p){ # Uses the West 2012 equation to solve for Fwsil
  with(as.list(c(p,y)),{

    Temp.0 <- p['temp.a.i']

      if (p['init']==FALSE) {
      Fextra <- C.perturbation(t)
      d13C.Fextra <- d13C.perturbation(t)
      Conc0 <- west.W(q=q0,Temperature=Temp.0,Temperature.0 = Temp.0,E.rate = 1/p['Ts'])
    }
    
    if (p['init']==TRUE) {
      Fextra <- C.perturbation.init(t)
      d13C.Fextra <- d13C.perturbation(t)
      Conc0 <- west.W(q=q0,Temperature=temp.o.i,Temperature.0=temp.o.i,E.rate = 1/2e4)
    }
        
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- Temp.0+p['Clim.Sens']*log2(RCO2.temp)

    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- Temp.0+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- p['q0']
    q <- q0*(1+p['q.sens'])^(Temp-Temp.0)
    
    # Weathering Fluxes from West Formulation
    # Conc0 <- west.W(q=q0,Temperature=15,E.rate = 1/2e4)
    Conc <- west.W(q=q,Temperature=Temp,Temperature.0=Temp.0,E.rate = 1/p['Ts'])
    Dw <- 1
    
    # Silicate Weathering
    Fwsil <- p['Fwsil.i']*Conc/Conc0

    # Carbonate Fluxes
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,Fextra,Conc,Conc0,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

CLiBeSO.P = function(t,y,p){ # Calculate P weathering and burial fluxes and Forgb is dependent upon P
  with(as.list(c(p,y)),{
    Fextra <- C.perturbation(t)
    d13C.Fextra <- d13C.perturbation(t)
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)
    
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- 0.3
    # q <- (log2(RCO2)*p['q.sens']*q0)+q0
    q <- q0*(1+p['q.sens'])^(Temp-p['temp.a.i'])
    
    # Concentrations
    # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
    # open-system calculation for Ceq
    Conc.results <- MCfunc(pCO2 = pCO2.i,q0 = q0, q = q, Ea = p['Ea'], Ts=p['Ts'],
                           keff.ref=p['keff.ref'],Temperature=Temp,
                           soilCO2.use=p['soil.CO2.flag'])
    Conc0 <- Conc.results$Conc0
    Conc <- Conc.results$Conc
    Dw <- Conc.results$Dw

    # Silicate Weathering
    Fwsil <- p['Fwsil.i']*q/q0*Conc/Conc0
    
    # Carbonate Fluxes
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    Fwp <- p['Fwp.i']*(((p['sil.P'])*(Fwsil/p['Fwsil.i']))+((p['carb.P'])*(Fwcarb/p['Fwcarb.i']))+((p['org.P'])*(Fworg/p['Fworg.i'])))
    
    # Organic and P burial fluxes
    # Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    Fborg <- p['Fborg.i']*(y['P']/PO4.i)^2
    Fbp <- p['Fbp.i']*(Fborg/p['Fborg.i'])
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    dP <- Fwp - Fbp
    
    results <- c(dDIC, dAlk, dd13C, dS, dP)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

CLiBeSO.sfw = function(t,y,p){
  with(as.list(c(p,y)),{
    Fextra <- C.perturbation(t)
    d13C.Fextra <- d13C.perturbation(t)
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)
    
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- 0.3
    q <- q0*(1+p['q.sens'])^(Temp-p['temp.a.i'])
    
    # Concentrations
    # Conc <- Conc0*(1+p['Conc.sens'])^(Temp-p['temp.a.i'])
    # Conc <- (log2(RCO2)*p['Conc.sens']*Conc0)+Conc0 #(RCO2*(1+p['Conc.sens']))/RCO2
    # Conc <- exp(p['Ea']/R*(1/(p['temp.a.i']+273.15)-1/(Temp+273.15))) # Activation Energy
    
    # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
    # open-system calculation for Ceq
    Conc.results <- MCfunc(pCO2 = pCO2.i,q0 = q0, q = q, Ea = p['Ea'], Ts=p['Ts'],
                           keff.ref=p['keff.ref'],Temperature=Temp,
                           soilCO2.use=p['soil.CO2.flag'])
    Conc0 <- Conc.results$Conc0
    Conc <- Conc.results$Conc
    Dw <- Conc.results$Dw
    
    # Seafloor Weathering
    Fsfw <- p['Fwsil.i']/5*(RCO2)^0.32 # Equation and exponent from Brady and Gislason 1997
    
    # Total Silicate Weathering
    Fwsil <- (p['Fwsil.i']*0.8)*q/q0*Conc/Conc0 + Fsfw
    
    # Carbonate Fluxes--equations are from Shields and Mills 2017 (PNAS)
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

CLiBeSO.qC = function(t,y,p){
  with(as.list(c(p,y)),{
    Fextra <- C.perturbation(t)
    d13C.Fextra <- d13C.perturbation(t)
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)
    
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- 0.3
    q <- q0*(1+p['q.sens'])^(Temp-p['temp.a.i'])
    
    # Concentrations
    Conc0 <- 1
    Conc <- Conc0*(1+p['Conc.sens'])^(Temp-p['temp.a.i'])
    # Conc <- (log2(RCO2)*p['Conc.sens']*Conc0)+Conc0 #(RCO2*(1+p['Conc.sens']))/RCO2
    # Conc <- exp(p['Ea']/R*(1/(p['temp.a.i']+273.15)-1/(Temp+273.15))) # Activation Energy
    
    Dw <- 1
    
    # Silicate Weathering
    Fwsil <- p['Fwsil.i']*q/q0*Conc/Conc0
    
    # Carbonate Fluxes--equations are from Shields and Mills 2017 (PNAS)
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

CLiBeSO.PETM = function(t,y,p){
  with(as.list(c(p,y)),{
    Fextra <- C.perturbation(t)
    d13C.Fextra <- d13C.perturbation(t)
    ForgB.extra <- Corg.perturbation(t)
    d13C.ForgB.extra <- d13C.Corg.perturbation(t)
    
    #Calculate carbon system parameters (two-step iterative procedure to get effect of temperature change correct)
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=temp.o.i,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    pCO2.temp <- carb.parms$pCO2
    RCO2.temp <- pCO2.temp/pCO2.i
    Temp.temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2.temp)
    
    carb.parms <- co2sys(DIC=y['DIC']/oceanV,TA=y['Alk']/oceanV,temp=Temp.temp,
                         sal=salinity,p=1,Ca.in=Ca.i,Mg.in=Mg2.i)
    
    pCO2 <- carb.parms$pCO2
    omega <- carb.parms$OmegaCalcite
    RCO2.sil = pCO2/280
    RCO2 <- pCO2/pCO2.i
    
    # Temperature
    Temp <- p['temp.a.i']+p['Clim.Sens']*log2(RCO2)
    
    # Volcanism
    Fvolc <- p['Fvolc.i']
    
    # Runoff    
    q0 <- 0.3
    q <- q0*(1+p['q.sens'])^(Temp-p['temp.a.i'])
    
    # Concentrations
    # Conc <- Conc0*(1+p['Conc.sens'])^(Temp-p['temp.a.i'])
    # Conc <- (log2(RCO2)*p['Conc.sens']*Conc0)+Conc0 #(RCO2*(1+p['Conc.sens']))/RCO2
    # Conc <- exp(p['Ea']/R*(1/(p['temp.a.i']+273.15)-1/(Temp+273.15))) # Activation Energy
    
    # Using a Maher and Chamberlain-style concentration feedback with a Winnick and Maher 
    # open-system calculation for Ceq
    Conc.results <- MCfunc(pCO2 = pCO2.i,q0 = q0, q = q, Ea = p['Ea'], Ts=p['Ts'],
                           keff.ref=p['keff.ref'],Temperature=Temp,
                           soilCO2.use=p['soil.CO2.flag'])
    Conc0 <- Conc.results$Conc0
    Conc <- Conc.results$Conc
    Dw <- Conc.results$Dw
    
    # Silicate Weathering
    Fwsil <- p['Fwsil.i']*q/q0*Conc/Conc0
    
    # Carbonate Fluxes
    Fwcarb <- p['Fwcarb.i']*q/q0*Conc/Conc0
    
    Fbcarb <- p['Fbcarb.i']*(omega/omega.i) #carbonate burial, mol/yr
    
    # Alkalinity Fluxes    
    Falkin = 2*Fwcarb+2*Fwsil
    Falkout = 2*Fbcarb
    
    # Organic and P weathering fluxes
    Fworg <- p['Fworg.i']
    
    # Organic and P burial fluxes
    Fborg <- p['Fborg.i']*(Fbcarb/p['Fbcarb.i'])^p['borg.carb']
    
    # Sulfur Fluxes and associated Carbon fluxes
    Fwsulf <- p['Fwsulf.i']
    Fbsulf <- p['Fbsulf.i']*y['S']/S.i
    
    Fwsulf.C <- 2*Fwsulf # Assumes 1:2 molar ratio (S -> C)
    Fbsulf.C <- 2*Fbsulf # Assumes 1:2 molar ratio (S -> C)
    
    dDIC <- Fvolc + Fworg + Fwcarb - Fborg - Fbcarb + Fwsulf.C - Fbsulf.C + Fextra - ForgB.extra
    dAlk <- Falkin - Falkout
    dd13C <- (Fworg*(p['Fworg.d13C']-y['d13C']) + Fwcarb*(p['Fwcarb.d13C']-y['d13C']) + Fvolc*(p['Fvolc.d13C']-y['d13C']) + p['DB']*Fborg + Fextra*(d13C.Fextra-y['d13C']) - ForgB.extra*(d13C.ForgB.extra-y['d13C']))/y['DIC']
    dS <- Fwsulf - Fbsulf
    
    results <- c(dDIC, dAlk, dd13C, dS)
    
    list(results, pCO2, omega, Fwcarb, Fbcarb, Fborg, Fwsil,
         RCO2, Temp, Fworg, Fwsulf, Fbsulf, q, Conc, Dw)
  })
}

#### Uber CLiBeSO Functions ####
# These functions are necessary for parallelizing the above CLiBeSO functions
CLiBeSO.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",
                            'pCO2','omega', 'Fwcarb',
                            'Fbcarb', 'Fborg','Fwsil',
                            'RCO2','Temp','Fworg',"Fwsulf",
                            "Fbsulf","q","Conc","Dw")
  list(parameters,final.results)
}

CLiBeSO.3.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO.3,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",'pCO2','omega', 'Fwcarb','Fbcarb', 
                            'Fborg','Fwsil','Fwsil.1','Fwsil.3','RCO2','Temp','Fworg',"Fwsulf","Fbsulf",
                            "q.1","q.2","q.3",
                            "q.avg","Conc.1","Conc.2","Conc.3","Conc.avg","Dw.1","Dw.2","Dw.3")
  list(parameters,final.results)
}

CLiBeSO.chem.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO.chem,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",
                            'pCO2','omega', 'Fwcarb',
                            'Fbcarb', 'Fborg','Fwsil',
                            'RCO2','Temp','Fworg',"Fwsulf",
                            "Fbsulf","q","Conc","Dw")
  list(parameters,final.results)
}

CLiBeSO.west.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO.west,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",
                            'pCO2','omega', 'Fwcarb',
                            'Fbcarb', 'Fborg','Fwsil',"Fextra","Conc","Conc0",
                            'RCO2','Temp','Fworg',"Fwsulf",
                            "Fbsulf","q","Conc","Dw")
  list(parameters,final.results)
}

CLiBeSO.sfw.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO.sfw,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",
                            'pCO2','omega', 'Fwcarb',
                            'Fbcarb', 'Fborg','Fwsil',
                            'RCO2','Temp','Fworg',"Fwsulf",
                            "Fbsulf","q","Conc","Dw")
  list(parameters,final.results)
}

CLiBeSO.qC.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO.qC,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",
                            'pCO2','omega', 'Fwcarb',
                            'Fbcarb', 'Fborg','Fwsil',
                            'RCO2','Temp','Fworg',"Fwsulf",
                            "Fbsulf","q","Conc","Dw")
  list(parameters,final.results)
}

CLiBeSO.PETM.uber <- function(parameters){
  final.results <- as.data.frame(lsodes(y=y0,times=t,func=CLiBeSO.PETM,parms=parameters,rtol=1e-6))
  names(final.results) <- c('time','DIC', 'Alk','d13C',"S",
                            'pCO2','omega', 'Fwcarb',
                            'Fbcarb', 'Fborg','Fwsil',
                            'RCO2','Temp','Fworg',"Fwsulf",
                            "Fbsulf","q","Conc","Dw")
  list(parameters,final.results)
}

print("Functions loaded. This train is bound for glory!")