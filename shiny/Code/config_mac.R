#!/usr/bin/env Rscript 
args <- commandArgs(TRUE)

tryCatch({
  if (!"shiny" %in% installed.packages()[,"Package"]) install.packages("shiny", repos="http://cran.rstudio.com/")
}, error= function(e) {
  print("Cannot install 'Shiny'. Possible cause - No internet connection or CRAN repository inaccessible. Please check your internet connection or try installing Shiny manually.")
})

if (!"shiny" %in% installed.packages()[,"Package"]) stop("'Shiny' package missing. Execution aborted.")

list.of.packages <- c("ggplot2", "reshape2", "tools", "corrplot", "foreach", "parallel", "doParallel", "plyr", "seqinr", "rmngb", "ineq", "ttutils")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
tryCatch({
  if(length(new.packages)) install.packages(new.packages, repos="http://cran.rstudio.com/")
}, error= function(e) {
  print(paste0("Cannot install required packages. Possible cause - No internet connection or CRAN repository inaccessible. Please check your internet connection or try installing following package manually: ", new.packages))
})

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) stop("Required packages missing. Execution aborted.")

require(shiny)
library(methods)
SEQualyzer_path <- args[1]
runApp(SEQualyzer_path, launch.browser=TRUE)