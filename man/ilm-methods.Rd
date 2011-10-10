\name{ilm-methods}
\docType{methods}
\alias{getIntens}
\alias{getIntens-methods}
\alias{getIntens,ILM-method}
\alias{getIntens,ILM,missing-method}
\alias{getIntens,ILM,character-method}
\alias{getProbeConcs}
\alias{getProbeConcs-methods}
\alias{getProbeConcs,ILM-method}
\alias{getProbeConcs,ILM,missing-method}
\alias{getProbeConcs,ILM,character-method}
\alias{getExprSummary}
\alias{getExprSummary-methods}
\alias{getExprSummary,ILM-method}
\alias{getExprSummary,ILM,missing-method}
\alias{getExprSummary,ILM,character-method}
\alias{getSDSummary}
\alias{getSDSummary-methods}
\alias{getSDSummary,ILM-method}
\alias{getSDSummary,ILM,missing-method}
\alias{getSDSummary,ILM,character-method}
\alias{[,ILM-method}
\alias{show}
\alias{show,ILM-method}

\title{
Methods to access the results of ilm.
}
\description{
These methods allow to access the results of the function \code{\link[affyILM:ilm]{ilm}} stored in an object of type \code{\linkS4class{ILM}}.
}

\usage{
getIntens(object,y)
getProbeConcs(object,y)
getExprSummary(object,y,z)
getSDSummary(object,y,z)
}
\arguments{
  \item{object}{An object of type \code{\linkS4class{ILM}} }

  \item{y}{A character string or vector of probe set name(s)}
  
  \item{z}{A character string specifying the type of set that is requested}
}

\details{

\code{"getIntens()"} is used to access the intensity values.
 
The probe concentrations are calculated on behalf of the Langmuir model. For each probe 
(of a probeset) the concentration is estimated in picoMolar and can be accessed via \code{"getProbeConcs()"}. Medianpolish, Transposed Medianpolish or Median can then be used to compute probeset summarized expression values. 

The results can be accessed via \code{"getExprSummary()"} and the associated standard deviation is accessed via \code{"getSDSummary()"}. z, if specified, can take two values : Probe.Set and Cluster.Set (the last one has been created for the definition of cluster sets in the last generation of arrays, whole-gene).

If y=NULL, the results for all probe sets are shown. 

}

\value{
A \code{"matrix"}

An object of the class \code{\linkS4class{ILM}} when subsetting \code{"["}
}

\author{Myriam Kroll, Fabrice Berger and Enrico Carlon}

\seealso{
\code{\linkS4class{ILM}}}

\examples{
## Locate and read in CEL-file
path <- system.file("rawData/FusionSDK_HG-Focus/HG-Focus/2.Calvin/",package="AffymetrixDataTestFiles");
file1 <- paste(path,"HG-Focus-1-121502.CEL", sep="")
## Calculation of background estimates and expression values (concentrations)
result <- ilm(file1)
## Background intensities of all probes
getIntens(result)
## Background intensities for one or more probesets
getIntens(result,"203561_at") 
getIntens(result,c("203561_at","40359_at"))

## Get concentrations of all probesets
getProbeConcs(result)
## Get concentrations (in picoMolar)
getProbeConcs(result,"203561_at")
getProbeConcs(result,c("203561_at","40359_at"))

## Subsetting
result["203561_at"] 

}

\keyword{ methods }
