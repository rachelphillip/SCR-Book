## this to create:
## load("/Users/dlb/Dropbox/kruger (1)/popn.RData")
## kruger.popn=popn
## save(kruger.popn,file="./scrmlebook/data/kruger-popn.RData")
##
##
#' @name kruger-popn
#' @title Simulated leopard population in part of Kruger National Park.
#' @docType data
#' @description Simulated leopard population in part of Kruger National Park, created from fits to
#' real survey data. Used in Chapter 'Spatial Encounter Rates and Detection Probabilities'. An object
#' of class 'popn' (see package \code{secr}).
#' @usage data('kruger-popn')
#' @references 
#' Maputla, N.W. 2014. Drivers of leopard population dynamics in the Kruger National Park, 
#' South Africa. PhD Thesis, University of Pretoria, Pretoria, RSA.
#' @examples
#'  library(scrmlebook)
#'  data('kruger-popn')
#'  plot(kruger.popn)
#'  data('ERdetfun-data') # to get camera locations
#'  plot(kruger.cams,add=TRUE)
NULL
