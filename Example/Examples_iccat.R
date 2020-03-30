library(JABBApkg)
setwd("C:/Work/Research/GitHub/JABBApkg_testing/example")

#><>><>><>><>><>><>><>><>><>><>><>
# Bigeye Tuna ICCAT
#><>><>><>><>><>><>><>><>><>><>><>

assessment = "BETiccat"
scenario = "S1"
dir.create(assessment,showWarnings = F)
output.dir = file.path(getwd(),assessment)
setwd(output.dir)
#------------------------------------------------------
# Simple example fit JABBA to Catch and CPUE with SEs
#-------------------------------------------------------
data(iccat)
bet = iccat$bet
jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment=assessment,scenario = "FitCPUE",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01)
jabba = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)

# Make individual plots
jbplot_catch(jabba)
jbplot_catcherror(jabba)
jbplot_ppdist(jabba)
jbplot_mcmc(jabba)
jbplot_residuals(jabba)
jbplot_cpuefits(jabba)
jbplot_runstest(jabba)
jbplot_logfits(jabba)
jbplot_procdev(jabba)
jbplot_bprior(jabba)

# Trajectories
jbplot_trj(jabba,type="B")
jbplot_trj(jabba,type="F")
jbplot_trj(jabba,type="BBmsy")
jbplot_trj(jabba,type="FFmsy")
jbplot_trj(jabba,type="BB0")

jbplot_spphase(jabba)
jbplot_kobe(jabba)

# Write all as png
jabba_plots(jabba,output.dir = output.dir)

#------------------------------------------------------
# Estimate shape m as function of Bmsy/K
#-------------------------------------------------------

jbinput = build_jabba(catch=bet$catch,cpue=bet$cpue,se=bet$se,assessment=assessment,scenario = "Fit_shape",model.type = "Pella_m",BmsyK=0.37,shape.CV = 0.3,sigma.est = FALSE,fixed.obsE = 0.01)
jabba = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)



#----------------
# Catch-Only
#----------------
jbinput = build_jabba(catch=bet$catch,model.type = "Fox",assessment=assessment,scenario =  "CatchOnly" ,b.prior=c(0.7,0.2,2010,"bbmsy"))
jabba = fit_jabba(jbinput,save.jabba=TRUE,output.dir=output.dir)
jbplot_catcherror(jabba)
jbplot_trj(jabba,type="B")
jbplot_trj(jabba,type="BBmsy")
jbplot_bprior(jabba)

jabba_plots(jabba,output.dir = output.dir)

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# White Marlin (Maurato et al. 2019)
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>




Scenario = Scenarios[s] 

# Add CV on Catch of 0.2 (lognormal error)
add.catch.CV = TRUE
catch.cv = 0.2

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Suplus Production model specifications
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

# Choose model type: 
# 1: Schaefer
# 2: Fox  
# 3: Pella-Tomlinsson  

Model = 3

Mod.names = c("Schaefer","Fox","Pella")[Model]

# Depensation opiton:
# Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25) 
# Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models 
Plim = 0

# Required specification for Pella-Tomlinson (Model = 3)
BmsyK = 0.39 # Set Surplus Production curve inflection point
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>  

#--------------------------------------------------
# Read csv files
#--------------------------------------------------

# Use SEs from csv file for abudance indices (TRUE/FALSE)
SE.I = TRUE

# Load assessment data
catch = read.csv(paste0(File,"/",assessment,"/catch",assessment,".csv"),sep=";")
cpue = read.csv(paste0(File,"/",assessment,"/cpue",assessment,".csv"),sep=";")#

if(SE.I ==TRUE){
  se =  read.csv(paste0(File,"/",assessment,"/se",assessment,".csv"),sep=";")
}

#><> NEW Scale all indices to a mean of 0.3
for(i in 2:ncol(se)){
  se[,i] = se[,i]*0.3/mean(se[,i], na.rm=T)
}
# Check
apply(se[,-1],2,mean, na.rm=T)

names(cpue)
names(catch)

#--------------------------------------------------
# option to exclude CPUE time series or catch year
#--------------------------------------------------
# Exclude SP LL
if(s==1){  
  cpue = cpue[,-c(13)] 
  se = se[,-c(13)] 
}
# Remove first 3 years for JPN
cpue[1:3,2] = NA

#------------------------------------------------------
# Option use mean CPUE from state-space cpue averaging
#-----------------------------------------------------
meanCPUE = FALSE
#------------------------------------------------
# Prior for unfished biomass K
#------------------------------------------------
# The option are: 
# a) Specify as a lognormal prior with mean and CV 
# b) Specify as range to be converted into lognormal prior

K.dist = c("lnorm","range")[1]

# if lnorm use mean and CV; if range use lower,upper bound
K.prior = c(25000,2)  

#-----------------------------------------------------------
# mean and CV and sd for Initial depletion level P1= SB/SB0
#-----------------------------------------------------------
# Set the initial depletion prior B1/K 
# To be converted into a lognormal prior (with upper bound at 1.1)

psi.dist= c("lnorm","beta")[1]
# specify as mean and CV 
psi.prior = c(1,0.25) 

#--------------------------------------------------------------
# Determine estimation for catchability q and observation error 
#--------------------------------------------------------------
# Assign q to CPUE
sets.q = 1:(ncol(cpue)-1) 


#----------------------------------------------------
# Determine r prior
#----------------------------------------------------
# The option are: 
# a) Specifying a lognormal prior 
# b) Specifying a resiliance category after Froese et al. (2017; CMSY)
# Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)

# use [1] lognormal(mean,stdev) or [2] range (min,max) or
r.dist = c("lnorm","range")[1] 
r.prior = c(0.181,0.180) #h=0.6




#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
# Observation Error
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>

#To Estimate additional observation variance set sigma.add = TRUE
sigma.est = TRUE

# Series
sets.var = 1:(ncol(cpue)-1) # estimate individual additional variace

# As option for data-weighing
# minimum fixed observation error for each variance set (optional choose 1 value for both)
fixed.obsE = c(0.01) # Important if SE.I is not availble

# Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
# Process Error
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
#Estimate set sigma.proc == True
sigma.proc =  TRUE #ifelse(s<4,FALSE,TRUE)
# Determines if process error deviation are estimated for all years (TRUE)  
# or only from the point the first abundance index becomes available (FALSE)
proc.dev.all = FALSE 

#------------------------------------------
if(sigma.proc == TRUE){
  #><>igamma = c(4,0.01) #specify inv-gamma parameters
  igamma = c(0.001,0.001) # use unimformative prior
  # Process error check
  gamma.check = 1/rgamma(1000,igamma[1],igamma[2])
  # check mean process error + CV
  mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))
  
  # check CV
  round(c(mu.proc,CV.proc),3)
  quantile(sqrt(gamma.check),c(0.1,0.9))
}else{
  sigma.proc = 0.07 #IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
}




