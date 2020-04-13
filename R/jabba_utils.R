#' gets beta prior paramenters
#'
#' gets beta prior paramenters and plots distribution
#' @param mu mean of beta distribution
#' @param CV cv of beta distribution
#' @param Min min of x-axis range (default = 0)
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return a and b parameter (shape and scale)
#' @export
get_beta <- function(mu,CV,Min=0,Prior="x",Plot=FALSE){
  a = seq(0.0001,1000,0.001)
  b= (a-mu*a)/mu
  s2 = a*b/((a+b)^2*(a+b+1))
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = (a-mu*a)/mu
  x = seq(Min,1,0.001)
  pdf = dbeta(x,a,b)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(a,b))
}

#' gets gamma prior paramenters
#'
#' gets gamma prior paramenters and plots distribution
#' @param mu mean of gamma distribution
#' @param CV cv of gamma distribution
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return a and b parameter (shape and scale)
#' @export
get_gamma <- function(mu,CV,Prior="x", Plot=FALSE){
  a = seq(0.00001,10000,0.0001)
  b = a/mu
  s2 = (a/b^2)
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = a/mu
  x = sort(rgamma(1000,a,b))
  pdf = dgamma(x,a,b)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(a,b))
}


#' gets lognormal prior paramenters
#'
#' gets lognormal prior paramenters mean and log.sd with bias correction and plots distribution
#' @param mu mean of lognormal distribution
#' @param CV cv of lognoram distribution
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return mean and lod.sd
#' @export
plot_lnorm <- function(mu,CV,Prior="x",Plot=FALSE){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu)-0.5*sdev^2,sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))
  pdf = dlnorm(x,log(mu),sdev)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(mu,sdev))
}

#' Gets lognormal prior paramenters
#'
#' Gets lognormal prior paramenters mean and log.sd with bias correction and plots distribution
#' @param mu mean of lognormal distribution
#' @param CV cv of lognoram distribution
#' @param Prior X-axis label
#' @param Plot c(TRUE,FALSE)
#' @return mean and lod.sd
#' @export
plot_lnorm <- function(mu,CV,Prior="x",Plot=FALSE){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu)-0.5*sdev^2,sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))
  pdf = dlnorm(x,log(mu),sdev)
  if(Plot==TRUE){
    plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  }
  return(c(mu,sdev))
}

#' Function kobeJabba for FLR
#'
#' Function to convert kobe posteriors into KOBE FLR input object
#' @param x posterior array dims(iter,year,stock,harvest)
#' @export
kobeJabba<-function(x){

  out=cbind(reshape::melt(x[,,2]),c(x[,,3]))
  names(out)=c("iter","year","stock","harvest")
  out$year=out$year
  out}

#' Function to convert JABBA projections into  FLR object
#'
#' Function to convert kobe projection matrix posteriors into Kobe FLR input object
#' @param x posterior array dims(iter,year,tac,stock,harvest)
#' @export
kobeJabbaProj<-function(x){

  out=cbind(reshape::melt(x[,,,2]),c(x[,,,3]))
  names(out)=c("iter","year","tac","stock","harvest")
  out$year=out$year

  out}

#' Function to do runs.test and 3 x sigma limits
#'
#' runs test is conducted with library(snpar)
#' @param x residuals from CPUE fits
#' @param type only c("resid","observations")
#' @return runs p value and 3 x sigma limits
#' @export
runs_sig3 <- function(x,type=NULL) {
  if(is.null(type)) type="resid"
  if(type=="resid"){mu = 0}else{mu = mean(x, na.rm = TRUE)}
  # Average moving range
  mr  <- abs(diff(x - mu))
  amr <- mean(mr, na.rm = TRUE)
  # Upper limit for moving ranges
  ulmr <- 3.267 * amr
  # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
  mr  <- mr[mr < ulmr]
  amr <- mean(mr, na.rm = TRUE)
  # Calculate standard deviation, Montgomery, 6.33
  stdev <- amr / 1.128
  # Calculate control limits
  lcl <- mu - 3 * stdev
  ucl <- mu + 3 * stdev
  if(nlevels(factor(sign(x)))>1){
    runstest = snpar::runs.test(x)
    pvalue = round(runstest$p.value,3)} else {
      pvalue = 0.001
    }

  return(list(sig3lim=c(lcl,ucl),p.runs= pvalue))
}

#' jabba_hindcast()
#'
#' Wrapper to coduct histcasts for retrospective analysis and cross-validation
#' @param jbinput List of input variables as output by build_jabba()
#' MCMC settings
#' @param ni number of iterations
#' @param nt thinning interval of saved iterations
#' @param nb burn-in
#' @param nc number of mcmc chains
#' Initial values
#' @param init.values = FALSE,
#' @param init.K = NULL,
#' @param init.r = NULL,
#' @param init.q = NULL,# vector
#' @param peels sequence of retro spective peels default 0:5
#' @param save.jabba Save individual JABBA run outputs
#' @param output.dir path to save plot. default is getwd()
#' @param save.hc Save hindcast list output as .rdata
#' @param plotall if TRUE makes jabba_plots() for each run    
#' @param speedup Reduces MCMC after setting runs 2+ inits to first "full" reference run
#' @return hc = list() containing estimates of key joint results from all hindcast run 
#' @export
jabba_hindcast = function(jbinput,
                          # MCMC settings
                          ni = 30000, # Number of iterations
                          nt = 5, # Steps saved
                          nb = 5000, # Burn-in
                          nc = 2, # number of chains
                          # init values
                          init.values = FALSE,
                          init.K = NULL,
                          init.r = NULL,
                          init.q = NULL,# vector
                          peels = 0:5, # retro peel option
                          save.jabba = FALSE,
                          output.dir = getwd(),
                          save.hc = FALSE,
                          plotall = FALSE,
                          speedup = TRUE 
){
  # hindcast object define object  
  hc = list(yr=jbinput$data$yr,catch=jbinput$jagsdata$TC,peels=peels,timeseries = NULL,refpts=NULL,pfunc=NULL,diags=NULL,settings=NULL)
  hc$settings$cols = jbinput$settings$cols
  hc$settings$harvest = jbinput$settings$harvest.label
  hc$settings$catch.metric = jbinput$settings$catch.metric
  
  Scenario = jbinput$settings$scenario
  for(i in 1:length(peels)){
    jbinput$settings$scenario = peels[i]
    if(i == 1 | speedup==FALSE){
      mci = ni
      mct = nt
      mcb = nb
      mcc = nc
      Kin = init.K
      rin = init.r
      qin = init.q
    } else {
      mci = 5500
      mct = 1
      mcb = 500
      mcc = 2
      init.values=TRUE
    } 
    if(peels[i]%in%peels[2]){
      
      Kin = fithc$pars[1,1]
      rin = fithc$pars[2,1]
      qin = fithc$pars[(3:(2+fithc$settings$nq)),1]
    }
    
    fithc = fit_jabba(jbinput,save.jabba=save.jabba,output.dir=output.dir,
                      ni = mci, # Number of iterations
                      nt = mct, # Steps saved
                      nb = mcb, # Burn-in
                      nc = mcc, # number of chains
                      init.values = init.values,
                      init.K = Kin,
                      init.r = rin,
                      init.q = qin,# vector
                      peels = peels[i]) # retro peel option
    
    hc$timeseries$mu = rbind(hc$timeseries$mu,data.frame(factor=fithc$diags[1,1],level=fithc$diags[1,2],fithc$timeseries[,"mu",])) 
    hc$timeseries$lci = rbind(hc$timeseries$lci,data.frame(factor=fithc$diags[1,1],level=fithc$diags[1,2],fithc$timeseries[,"lci",])) 
    hc$timeseries$uci = rbind(hc$timeseries$uci,data.frame(factor=fithc$diags[1,1],level=fithc$diags[1,2],fithc$timeseries[,"uci",])) 
    hc$diags = rbind(hc$diags,fithc$diags)
    hc$refpts= rbind(hc$refpts,fithc$refpts[1,])
    hc$pfunc= rbind(hc$pfunc ,fithc$pfunc)
    
    if(plotall==TRUE){
      jabba_plots(fithc,output.dir = output.dir)
    }
    
  } # end of loop
  if(save.hc==TRUE){
    save(hc,file=paste0(output.dir,"/hc_",Scenario,".rdata"))  
  }
  return(hc)
} #end of hindcast function 
#------------------------------------------------------------------------------------------------------------------------------------


