#
#More sophisticated C3 model with nitrogen assimilation
# (c) 2019 D.S. Mackay
#
#This model will compute C3 photosynthesis (A), stomatal conductance (gc), 
#  and intercellular (Ci) and chloroplast (Cc) CO2 partial pressures. It also
#  computes nitrogen assimilation (glycine and serine). These are solved
#  simultaneously using nested binary search algorithms for efficiency.
#
#It is run by providing a proposed maximum gc, atmospheric CO2 (Ca). 
#  leaf temperature (T), and incoming PAR (I)
#
#C3 parameters are entered as constants within this code; change at will
#
#The model will find the lowest gs that gives the maximum A if run using coupledA_gc()
#
#This is implemented using von Caemmerer's (2013) equations that approximate both enzyme
#  and transport limited C3 photosynthesis. 
#
#References: von Caemmerer 2013 Plant, Cell & Environment
#            Busch et al 2017 Nature Plants
#
#
calcC3photo = function(Ca, I, Tl, gc)
{
  #This is a set of parameters (at 25 degrees C) and drivers
  Vcmax25 = 100 #Maximum rubisco activity, umol m-2 s-1
  Jmax25 = 1.7 * Vcmax25 #Maximum electron transport rate, umol electrons m-2 s-1
  theta = 0.85 #unitless curvature parameter
  Kc25 = 259 # , Michaelis-Menten constant for carboxylase, ubar
  Ko25 = 179 # , M-M constant for oxygenase, mbar
  R = 8.31 #universal gas constant, J K-1 mol-1
  E_Vcmax = 58.52 #activation energy maximum carboxylation rate, kJ mol-1
  E_Kc = 59.36 #activation energy for the M-M carboxylation constant, kJ mol-1
  E_Ko = 35.94 #activation energy for the M-M oxygenation constant, kJ mol-1
  E_Jmax = 37.0 #activation energy for electron transport, kJ mol-1
  E_gammaStar = 23.4 #activation energy for CO2 compensation point, kJ mol-1
  E_Rd = 66.3 #activation energy for leaf mitochondrial respiration, kJ mol-1
  Pa = 1013 #atmospheric pressure, mbar
  O = 0.21*Pa #O2 partial pressure, mbar
  gammaStar25 = 38.6 # , ubar
  f = 0.15 #correction for spectral quality of light
  x = 0.4 #partitioning factor of electron transport rate
  absorptance = 0.85 #fraction of irradiance, I, absorbed
  I2 = I * absorptance * (1-f)/2
  Rd25 = Vcmax25/100 #leaf mitochondrial respiration
  gm = 0.35 #mesophyll conductance, mol m-2 s-1
  gt = 1/ (1/gc+1/gm) #diffusive conductance from mesophyll to leaf surface, mol m-2 s-1
  Tlk = Tl + 273 #absolute temperature
  Vcmax = Vcmax25 * exp((Tlk-298)*E_Vcmax*1000/(298*R*Tlk))
  Kc = Kc25 * exp((Tlk-298)*E_Kc*1000/(298*R*Tlk))
  Ko = Ko25 * exp((Tlk-298)*E_Ko*1000/(298*R*Tlk))
  gammaStar = gammaStar25 * exp((Tlk-298)*E_gammaStar*1000/(298*R*Tlk))
  Rd = Rd25 * exp((Tlk-298)*E_Rd*1000/(298*R*Tlk))
  Jmax = Jmax25 * exp((Tlk-298)*E_Jmax*1000/(298*R*Tlk))
  Tp = 0.167*Vcmax #Triosphoshate limitation constant, umol m-2 s-1
  
#These constants are for nitrogen assimilation, Busch et al 2017 Nature Plants
  alphaGmax = 0.0 #fraction of glycolate carbon diverted from photorespiration to glycine
  alphaSmax = 0.0 #fraction of glycolate carbon diverted from photorespiration to serine
  Nmax = 1.5   #maximum rate of de novo nitrogen supply to the chloroplast, umol N m-2 s-1
  
  stop = 0
  Cc_max = Ca
  Cc_min = 0
  
  #Iterative solution of coupled Cc, Ci, and A
  #Correct solution is when gas exchange Aref equals the C3 A
  #Aref declines as Cc increases, while C3 A increases with Cm, and so there
  #is an equilibrium point, which is found efficiently with a binary search 
  while (stop==0)
  {
    Cc = 0.5 * (Cc_max + Cc_min) #select mid-point between current limits
    #Note: From leaf surface to inter-cell, A = gc(Ca - Ci)
    #      and from inter-cell to chloroplast, A = gm(Ci = Cc)
    #      so from leaf surface to chloroplast, A = gt(Ca - Cc)
    #      where gt = 1/(1/gc + 1/gm)
    
    #(1) gas exchange solution for photosynthesis
    Aref = gt * (Ca - Cc) 
    
    #
    #(2) Biochemical C and N assimilation
    #
    #(2a) Enzyme limited with N assimilation
    #Note: enzyme limited carboxylation is independent of N assimilation
    #Using von Caemmerer approximation
    
    beta = 0
    if (alphaGmax > 0)
    {
      beta = 3*alphaGmax/(3*alphaGmax + 2*alphaSmax)
    }
      
    PHI = 2.0*gammaStar/Cc
    Vcc = Cc*Vcmax/(Cc+Kc*(1+O/Ko)) 
    Voc = PHI * Vcc
    alphaGc = min(alphaGmax, Nmax*beta/Voc)
    alphaSc = min(alphaSmax, 3*Nmax*(1-beta)/(2*Voc))
    gammaStarG = gammaStar * (1 - alphaGc)
    Ac = Vcmax*(Cc - gammaStar)/(Cc + Kc*(1+O/Ko)) - Rd
    
    #(2b) Transport limited with N assimilation
    J = (I2 + Jmax - sqrt((I2+Jmax)^2- 4 * theta * I2 * Jmax))/(2*theta)
    if (J > Nmax * (2*beta+6))
    {
      alphaGj = min(alphaGmax,4*Nmax*beta*(1/PHI+1)/(J-Nmax*(2*beta+6)))
      alphaSj = min(alphaSmax,6*Nmax*(1-beta)*(1/PHI+1)/(J-Nmax*(2*beta+6)))
    }
    else
    {
      alphaGj = alphaGmax
      alphaSj = alphaSmax
    }
    gammaStarG = gammaStar * (1 - alphaGj)
    Aj = J * (Cc - gammaStarG)/(4*Cc + 8*gammaStarG * (1 + 2*alphaGj + alphaSj)) - Rd
    
    #(2c) Trios phosphate limited with N assimilation
    alphaGp = min(alphaGmax, Nmax*beta*(2/PHI-1)/(6*Tp+3*Nmax*(2-beta)))
    alphaSp = min(alphaSmax, (3/2)*Nmax*(1-beta)*(2/PHI-1)/(6*Tp+3*Nmax*(2-beta)))
    gammaStarG = gammaStar * (1 - alphaGp)
    Ap = 3*Tp*(Cc - gammaStarG)/(Cc - gammaStarG*(1 + 3*alphaGp + 4*alphaSp)) - Rd
    
    #Net photosynthesis is the minimum of three pathways
    A = min(Ac, Aj, Ap)
    #print(c(Ac,Aj,Ap))
    
    #Next compute the export of glycine and serine, umol N m-2 s-1
    if (A == Ac)
    {
      NG = alphaGc * Voc
      NS = 2/3 * alphaSc * Voc
    }
    else if (A == Aj)
    {
      Voj = J / (4/PHI + (4 + 8 * alphaGj + 4 * alphaSj))
      NG = alphaGj * Voj
      NS = 2/3 * alphaSj * Voj
    }
    else
    {
      Vop = 3 * Tp / (1/PHI - 0.5 * (1 + 3 * alphaGp + 4 * alphaSp))
      NG = alphaGp * Vop
      NS = 2/3 * alphaSp * Vop
    }
    
    #Back out Ci
    Ci = Cc + A/gm
    
    if (Aref > A)
    {
      Cc_min = Cc #solution is not lower than current Cm
    }
    else
    {
      Cc_max = Cc #solution is no higher than current Cm
    }
    if ((A/Aref) > 0.9999 && (A/Aref) < 1.0001 || A < 0)
    {
      stop = 1 #solution found
    }
  }
  return(c(Ci,Cc,A,NG,NS))
}

#
#Solves the coupled C3 photosynthesis - stomatal conductance
#
#The algorithm works by reducing gc to its lowest value without reducing A. It starts
# with a high gc (e.g., maximum supported by soil-xylem hydraulics) and finds the 
# smallest gc for the reference A, which is either enzyme or transport limited. This
# function calls calcCphoto() to compute A, Ci, and Cm. The algorithm uses a binary
# search for efficiency.
#
#Why does this matter? At low light, A is transport limited and so gc can be reduced
# from its maximum without affecting A. At high atmospheric CO2 concentration, A is
# may already be enzyme limited and so gs can be reduced without affecting A.
#
coupledA_gc = function(Ca, I, Tl, gc0)
{
  stop = 0
  gc_min = min(gc0,0.01)
  gc_max = gc0*2.0
  result = calcC3photo(Ca, I, Tl, gc0)
  A_ref = result[3]
  
  if (A_ref <= 0.0)
  {
    gc = gc_min
    stop = 1
  }
  
  while (stop == 0)
  {
    gc = 0.5 * (gc_min + gc_max) #use the mid-point between the two limits
    result = calcC3photo(Ca, I, Tl, gc)
    A = result[3]
    
    if (A < 0.998*A_ref)
    {
      gc_min = gc #solution is no lower than current gc
    }
    else
    {
      gc_max = gc #solution is no higher than the current gc
    }
    if (gc_min/gc_max > 0.9999 && gc_min/gc_max < 1.0001)
    {
      stop = 1 #solution found
    }
  }
  Ci = result[1]
  Cc = result[2]
  A = result[3]
  NG = result[4]
  NS = result[5]
  return(c(Ci,Cc,A,gc,NG,NS))
}


#
#This code computes A, Ci, Cm, and gc for a set of atmospheric CO2 partial pressures
#
#Use this to find out what gc is needed at low input CO2 partial pressure
# or how much A is limited by gas exchange at high PAR
#
gc_vec = c(0.001,0.01,0.02,0.04,0.08,0.16,0.32,0.64,1.28) 
Cc_vec = array(data=0,dim=c(1,length(gc_vec)))
Ci_vec = array(data=0,dim=c(1,length(gc_vec)))
A_vec = array(data=0,dim=c(1,length(gc_vec)))
NG_vec = array(data=0,dim=c(1,length(gc_vec)))
NS_vec = array(data=0,dim=c(1,length(gc_vec)))
#gc_vec = array(data=0,dim=c(1,length(gc_vec)))
I = 1000 #PAR irradiance, umol m-2 s-1
Tl = 25 #leaf temperature, degrees C
Ca = 2000 #start very high

#loop through all input CO2 partial pressures
for (i in 1:length(gc_vec))
{
  gc = gc_vec[i]
  result = calcC3photo(Ca, I, Tl, gc)
  Ci_vec[i] = result[1]
  Cc_vec[i] = result[2]
  A_vec[i] = result[3]
  NG_vec[i] = result[4]
  NS_vec[i] = result[5]
  #gc_vec[i] = result[4]
}

#
#at different light levels, predict stomatal conductane and photosynthesis
#
I_vec = c(15,27,30,40,80,120,180,240,320,400,480,600,720,750,780,820,860,900,940,1200,1600,2000)
Cc_vec2 = array(data=0,dim=c(1,length(I_vec)))
Ci_vec2 = array(data=0,dim=c(1,length(I_vec)))
A_vec2 = array(data=0,dim=c(1,length(I_vec)))
gc_vec2 = array(data=0,dim=c(1,length(I_vec)))
NG_vec2 = array(data=0,dim=c(1,length(I_vec)))
NS_vec2 = array(data=0,dim=c(1,length(I_vec)))
N_vec2 = array(data=0,dim=c(1,length(I_vec)))
Ca = 280
Tl = 25 #leaf temperature, degrees C
gc = 1.0
#loop through PAR levels
for (i in 1:length(I_vec))
{
  I = I_vec[i]
  result = coupledA_gc(Ca, I, Tl, gc)
  Ci_vec2[i] = result[1]
  Cc_vec2[i] = result[2]
  A_vec2[i] = result[3]
  gc_vec2[i] = result[4]
  NG_vec2[i] = result[5]
  NS_vec2[i] = result[6]
  N_vec2[i] = result[5]+result[6]
}

#
#again with elevated atmospheric CO2
#
I_vec = c(15,27,30,40,80,120,180,240,320,400,480,600,720,750,780,820,860,900,940,1200,1600,2000)
Cc_vec3 = array(data=0,dim=c(1,length(I_vec)))
Ci_vec3 = array(data=0,dim=c(1,length(I_vec)))
A_vec3 = array(data=0,dim=c(1,length(I_vec)))
gc_vec3 = array(data=0,dim=c(1,length(I_vec)))
NG_vec3 = array(data=0,dim=c(1,length(I_vec)))
NS_vec3 = array(data=0,dim=c(1,length(I_vec)))
N_vec3 = array(data=0,dim=c(1,length(I_vec)))
Ca = 3*280
Tl = 25 #leaf temperature, degrees C
gc = 1.0
#loop through gc levels
for (i in 1:length(I_vec))
{
  I = I_vec[i]
  result = coupledA_gc(Ca, I, Tl, gc)
  Ci_vec3[i] = result[1]
  Cc_vec3[i] = result[2]
  A_vec3[i] = result[3]
  gc_vec3[i] = result[4]
  NG_vec3[i] = result[5]
  NS_vec3[i] = result[6]
  N_vec3[i] = result[5]+result[6]
}

#
#Plots A-Ci and light response output
#
exp_list <- c(as.expression(bquote(italic(A)[N]~"("~mu ~"mol" ~m^-2~s^-1~")" )),
              as.expression(bquote(italic(C)[i] ~"("~ mu~"bar)")),
              as.expression(bquote(italic(I)[In]~ "("~mu ~"mol"~ m^2~s^1~")")),
              as.expression(bquote(italic(g)[C]~"(mmol"~m^2~s^1~")")),
              as.expression(bquote(italic(N)[A]~"(umol N"~m^2~s^1~")")))

par(ps=c(24), cex=3.0)
par(bg="white")
nf <- layout(matrix(seq(1,4,1), 2, 2, byrow = TRUE))
layout.show(nf)
par(oma=c(1.0,1.0,1.0,0.1))
par(mar=c(8,10,1,1))
par(mgp=c(5.0,1.0,0))

#plot A-Ci curve
ylow0 <- min(c(Ci_vec2,Ci_vec3))
yhigh0 <- max(c(Ci_vec3,Ci_vec3))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(c(I_vec))
xhigh0 <- max(c(I_vec))
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(I_vec,Ci_vec2,type="l", lwd=4, col="black", xlab=exp_list[3], ylab=exp_list[2],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
lines(I_vec, Ci_vec3, lty=2, lwd=4, col="black")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,400,800,1200,1600,2000),labels=c("0","400","800","1200","1600","2000"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,100,200,300,400,500,600,700,800,900,1000),
     labels=c("0","100","200","300","400","500","600","700","800","900","1000"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
#legend(130,35, lty=c(1,1), col=c("black","blue"),lwd=c(4,4),
 #      c("A-Ci","A-Cm"), bty="n")
box(lwd=4)

#plot light response curve
ylow0 <- min(c(A_vec2,A_vec3))
yhigh0 <- max(c(A_vec2,A_vec3))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(I_vec)
xhigh0 <- max(I_vec)
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(I_vec,A_vec2,type="l", lwd=4, col="black", xlab=exp_list[3], ylab=exp_list[1],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
lines(I_vec, A_vec3, lty=2, lwd=4, col="black")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,400,800,1200,1600,2000),labels=c("0","400","800","1200","1600","2000"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,5,10,15,20,25,30,35),labels=c("0","5","10","15","20","25","30","35"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
box(lwd=4)

#plot gc response to light level
ylow0 <- min(c(gc_vec2,gc_vec3))
yhigh0 <- max(c(gc_vec2,gc_vec3))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(I_vec)
xhigh0 <- max(I_vec)
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(I_vec,gc_vec2,type="l", lwd=4, col="black", xlab=exp_list[3], ylab=exp_list[4],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
lines(I_vec, gc_vec3, lty=2, lwd=4, col="black")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,400,800,1200,1600,2000),labels=c("0","400","800","1200","1600","2000"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5),
     labels=c("0","150","300","450","600","750","900","1050","1200","1350","1500"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
box(lwd=4)

#plot N versus I
ylow0 <- min(c(0))
yhigh0 <- max(c(N_vec2,N_vec3))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(c(I_vec))
xhigh0 <- max(c(I_vec))
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(I_vec,NG_vec2,type="l", lwd=4, col="black", xlab=exp_list[3], ylab=exp_list[5],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
lines(I_vec, NS_vec2, lty=1, lwd=4, col="green")
lines(I_vec, N_vec2, lty=1, lwd=4, col="darkgreen")
lines(I_vec, NG_vec3, lty=2, lwd=4, col="black")
lines(I_vec, NS_vec3, lty=2, lwd=4, col="green")
lines(I_vec, N_vec3, lty=2, lwd=4, col="darkgreen")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,400,800,1200,1600,2000),labels=c("0","400","800","1200","1600","2000"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,0.25,0.5,0.75,1.0,1.25),labels=c("0","0.25","0.50","0.75","1.0","1.25"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
box(lwd=4)


#
#This code computes A, Ci, Cm, and gc for a set of atmospheric CO2 partial pressures
#
#Use this to find out what gc is needed at low input CO2 partial pressure
# or how much A is limited by gas exchange at high PAR
#
Ca_vec = c(50,75,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000) 
Cc_vec = array(data=0,dim=c(1,length(Ca_vec)))
Ci_vec = array(data=0,dim=c(1,length(Ca_vec)))
A_vec = array(data=0,dim=c(1,length(Ca_vec)))
gc_vec = array(data=0,dim=c(1,length(Ca_vec)))
NG_vec = array(data=0,dim=c(1,length(Ca_vec)))
NS_vec = array(data=0,dim=c(1,length(Ca_vec)))
I = 2000 #PAR irradiance, umol m-2 s-1
Tl = 25 #leaf temperature, degrees C
gc = 0.8 #start very high

#loop through all input CO2 partial pressures
for (i in 1:length(Ca_vec))
{
  Ca = Ca_vec[i]
  result = coupledA_gc(Ca, I, Tl, gc)
  Ci_vec[i] = result[1]
  Cc_vec[i] = result[2]
  A_vec[i] = result[3]
  gc_vec[i] = result[4]
  NG_vec[i] = result[5]
  NS_vec[i] = result[6]
}

#
#Plots A-Ci and light response output
#
exp_list <- c(as.expression(bquote("C3 "~italic(A)[N]~"("~mu ~"mol" ~m^-2~s^-1~")" )),
              as.expression(bquote(italic(C)[i] ~"("~ mu~"bar)")),
              as.expression(bquote(italic(C)[a]~ "("~mu ~"bar")),
              as.expression(bquote(italic(g)[C]~"(mmol"~m^2~s^1~")")),
              as.expression(bquote(italic(N)[A]~"(umol N"~m^2~s^1~")")))

par(ps=c(24), cex=3.0)
par(bg="white")
nf <- layout(matrix(seq(1,4,1), 2, 2, byrow = TRUE))
layout.show(nf)
par(oma=c(1.0,1.0,1.0,0.1))
par(mar=c(8,10,1,1))
par(mgp=c(5.0,1.0,0))

#plot A-Ci curve
ylow0 <- min(c(A_vec,A_vec))
yhigh0 <- max(c(A_vec,A_vec))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(c(Ci_vec,Ci_vec))
xhigh0 <- max(c(Ci_vec,Ci_vec))
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(Ci_vec,A_vec,type="l", lwd=4, col="black", xlab=exp_list[2], ylab=exp_list[1],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
#lines(Ci_vec3, A_vec3, lty=2, lwd=4, col="red")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,300,600,900,1200,1500,1800),
     labels=c("0","300","600","900","1200","1500","1800"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,5,10,15,20,25),labels=c("0","5","10","15","20","25"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
#legend(130,35, lty=c(1,1), col=c("black","blue"),lwd=c(4,4),
#      c("A-Ci","A-Cm"), bty="n")
box(lwd=4)

#plot light response curve
ylow0 <- min(c(NG_vec,NS_vec))
yhigh0 <- max(c(NG_vec,NS_vec))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(Ci_vec)
xhigh0 <- max(Ci_vec)
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(Ci_vec,NG_vec,type="l", lwd=4, col="black", xlab=exp_list[2], ylab=exp_list[5],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
lines(Ci_vec, NS_vec, lty=1, lwd=4, col="green")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,300,600,900,1200,1500,1800),
     labels=c("0","300","600","900","1200","1500","1800"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,0.25,0.5,0.75,1.0,1.25),labels=c("0","0.25","0.50","0.75","1.0","1.25"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
box(lwd=4)

#plot gc response to light level
ylow0 <- min(c(gc_vec,gc_vec))
yhigh0 <- max(c(gc_vec,gc_vec))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(Ci_vec)
xhigh0 <- max(Ci_vec)
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(Ci_vec,gc_vec,type="l", lwd=4, col="black", xlab=exp_list[2], ylab=exp_list[4],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
#lines(I_vec, gc_vec3, lty=2, lwd=4, col="red")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,300,600,900,1200,1500,1800),
     labels=c("0","300","600","900","1200","1500","1800"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5),
     labels=c("0","150","300","450","600","750","900","1050","1200","1350","1500"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
box(lwd=4)

#plot A versus gc
ylow0 <- min(c(A_vec,A_vec))
yhigh0 <- max(c(A_vec,A_vec))
ylow <- ylow0-0.05*(max(abs(ylow0),abs(yhigh0)))
yhigh <- yhigh0+0.12*(max(abs(ylow0),abs(yhigh0)))
xlow0 <- min(c(gc_vec,gc_vec))
xhigh0 <- max(c(gc_vec,gc_vec))
xlow <- xlow0-0.05*(max(abs(xlow0),abs(xhigh0)))
xhigh <- xhigh0+0.12*(max(abs(xlow0),abs(xhigh0)))

par(xaxt="n")
par(yaxt="n")

plot(gc_vec,A_vec,type="l", lwd=4, col="black", xlab=exp_list[4], ylab=exp_list[1],
     cex.axis=1.5,cex.lab=1.5,xlim=c(xlow,xhigh), ylim=c(ylow,yhigh))
#lines(gc_vec3, A_vec3, lty=2, lwd=4, col="red")
par(xaxt="s")
par(yaxt="s")
axis(1, at=c(0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5),
     labels=c("0","150","300","450","600","750","900","1050","1200","1350","1500"),
     cex.axis=1.25,lwd.ticks=4,padj=0.35, tck=0.025)
axis(2, at=c(0,5,10,15,20,25),labels=c("0","5","10","15","20","25"),
     cex.axis=1.25,lwd.ticks=4, tck=0.025)
box(lwd=4)

