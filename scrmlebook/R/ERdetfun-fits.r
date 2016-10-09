#' @name ERdetfun-fits
#' @title Fitted models for Chapter 'Spatial Encounter Rates and Detection Probabilities'.
#' @docType data
#' @description Contains these objects of class \code{secr}, created using package \code{secr} (see that package 
#' for object formats and content):
#'  \describe{
#'    \item{\code{er.hn.fit}:}{Encounter rate model with 'half-normal' (Gaussian) form, 
#'    fitted to \code{kruger.capt} (\code{detector(traps(er.hn.fit$capthist))=="count"}).}
#'    \item{\code{er.hr.fit}:}{Encounter rate model with 'hazard-rate' form, 
#'    fitted to \code{kruger.capt} (\code{detector(traps(er.hn.fit$capthist))=="count"}).}
#'    \item{\code{er.bin.hn.fit}:}{Encounter rate model with 'half-normal' (Gaussian) form, 
#'    fitted to \code{kruger.capt.bin} (\code{detector(traps(er.hn.fit$capthist))=="proximity"}).}
#'    \item{\code{er.bin.hr.fit}:}{Encounter rate model with 'hazard-rate' form, 
#'    fitted to \code{kruger.capt.bin} (\code{detector(traps(er.hn.fit$capthist))=="proximity"}).}
#'    \item{\code{er.hn.fit}:}{Detection function model with 'half-normal' (Gaussian) form, 
#'    fitted to \code{kruger.capt} (\code{detector(traps(er.hn.fit$capthist))=="count"}).}
#'    \item{\code{er.hr.fit}:}{Detection function model with 'hazard-rate' form, 
#'    fitted to \code{kruger.capt} (\code{detector(traps(er.hn.fit$capthist))=="count"}).}
#'    \item{\code{er.bin.hn.fit}:}{Detection function model with 'half-normal' (Gaussian) form, 
#'    fitted to \code{kruger.capt.bin} (\code{detector(traps(er.hn.fit$capthist))=="proximity"}).}
#'    \item{\code{er.bin.hr.fit}:}{Detection function model with 'hazard-rate' form, 
#'    fitted to \code{kruger.capt.bin} (\code{detector(traps(er.hn.fit$capthist))=="proximity"}).}
#'  }
#'  
#' @usage data('ERdetfun-fits')
#' @examples
#'  library(secr)
#'  library(scrmlebook)
#'  data('ERdetfun-fits')
#'  er.hr.fit
#'  plot(er.hr.fit)
NULL