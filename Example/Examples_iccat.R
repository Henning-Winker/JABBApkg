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





