## this to create the above:
#load("./analysis/data/ERdetfun-data.RData")
## then fit stuff, then save as below:
#savelist=c("kruger.capt","kruger.capt.bin","kruger.cams","kruger.cams.bin",
#           "kruger.mask",
#           "er.hn.fit","df.hn.fit","er.bin.hn.fit","df.bin.hn.fit",
#           "er.hr.fit","df.hr.fit","er.bin.hr.fit","df.bin.hr.fit")
#
#' @name ERdetfun-data
#' @title Datasets for Chapter 'Spatial Encounter Rates and Detection Probabilities'.
#' @docType data
#' @description Contains these objects created using package \code{secr} (see that package 
#' for data formats):
#'  \describe{
#'    \item{\code{kruger.capt}:}{A capture history object (class \code{capthist}) 
#'    containing the numbers of captures of each individual at each of 62 detectors.}
#'    \item{\code{kruger.capt.bin}:}{A capture history object (class \code{capthist}) 
#'    containing binary indicators of capture/not of each individual at each of 62 
#'    detectors.}
#'    \item{\code{kruger.cams}:}{A traps object (class \code{traps}) containing the 
#'    numbers of captures for each individual at each of 62 detectors of type 'count' 
#'    (\code{detector(traps(capthist))=="count"}). This object is contained in 
#'    \code{traps(kruger.capt)}.}
#'    \item{\code{kruger.cams.bin}:}{A traps object (class \code{traps}) containing 
#'    binary indicators of capture/not for each individual at each of 62 detectors of 
#'    type 'proximity' (\code{detector(traps(capthist))=="proximity"}). This object 
#'    is contained in \code{traps(kruger.capt.bin)}.}
#'    \item{\code{kruger.mask}:}{A mask object (class \code{mask}) 
#'    containing a mask suitable for use with the above objects.}
#'  }
#'  All of the above have these covariates attached:
#'  \describe{
#'    \item{\code{sname}:}{Name of the array. A factor with levels 'Malelane' or 'Pretoriuskop'.}
#'    \item{\code{habitat.cov}:}{Covariate indexing habitat suitability (a continuous 
#'    variable).}
#'    \item{\code{landscape}:}{Habitat class. A factor with levels 'Pretoriuskop Sourveld',
#'    'Mixed Bushwillow Woodlands', 'Malelane Mountain Bushveld', 'Thorn Veld'.}
#'    \item{\code{water}:}{A measure of annual water. A continuous variable.}
#'    \item{\code{high.water}:}{A binary variable which is 1 where \code{water} is above a 
#'    threshold value.}
#'    \item{\code{dist.to.water}:}{Distance to closest \code{high.water} value of 1.}
#'  }
#'  
#' @usage data('ERdetfun-data')
#' @references 
#' Maputla, N.W. 2014. Drivers of leopard population dynamics in the Kruger National Park, 
#' South Africa. PhD Thesis, University of Pretoria, Pretoria, RSA.
#' @examples
#'  library(scrmlebook)
#'  data('ERdetfun-data')
#'  plotcovariate(kruger.mask,covariate="habitat.cov",asp=1,bty="n"))
#'  plot(kruger.cams,add=TRUE)
NULL