#' Install missing R packages
#' 
#' checks if all required R packages are installed
#' 
#' @export

jabba_libs <- function(){
  # Install required packages if missing
  list.of.packages <- c("gplots", "coda","rjags","R2jags","fitdistrplus","reshape","mvtnorm","scales",'snpar')
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  }  