\name{ilm-methods}
\docType{methods}
\alias{getIntens}
\alias{getIntens-methods}
\alias{getIntens,ILM-method}
\alias{getIntens,ILM,missing-method}
\alias{getIntens,ILM,character-method}
\alias{getConcs}
\alias{getConcs-methods}
\alias{getConcs,ILM-method}
\alias{getConcs,ILM,missing-method}
\alias{getConcs,ILM,character-method}
\alias{getParams}
\alias{getParams-methods}
\alias{getParams,ILM-method}
\alias{[,ILM-method}
\alias{show}
\alias{show,ILM-method}

\title{
Methods to access the results of ilm
}
\description{
These methods allow to access the results of the function \code{\link[affyILM:ilm]{ilm}} stored in an object of type \code{\linkS4class{ILM}}.
}

\usage{
getIntens(object,y)
getConcs(object,y)
getParams(object)
}

\arguments{
  \item{object}{An object of type \code{\linkS4class{ILM}} }

  \item{y}{A character string or vector of probe set name(s). If y is missing, values are retrieved for all probesets.}
}

\details{A linear function is used to approximate the shape of the background intensities. By  
means of the linear least square method the 50 parameters characterizing this function are
determined. The optimized set of parameters can be accessed using \code{getParams()}.

\code{getIntens()} is used to access the calculated background intensities.
 
The concentrations are calculated on behalf of the Langmuir model. For each probe 
(of a probeset) the concentration is estimated in picoMolar. Its median value is then considered
as the concentration/probeset. The results can be accessed via \code{getConcs()}.

If y is missing, the results for all probe sets are shown. 

}

\value{
A \code{matrix} of dimension (number of probes) x (2x(number of CEL-files)) for \code{getIntens()} 
A \code{matrix} of dimension 50 x (number of CEL-files) and for \code{getParams()}
A \code{matrix} of dimension (number of probesets) x (number of CEL-files) for \code{getConcs()}

An object of the class \code{\linkS4class{ILM}} when subsetting \code{[}
}

\author{Myriam Kroll, Fabrice Berger and Enrico Carlon}

\seealso{
\code{\linkS4class{ILM}}}

\examples{
## Locate and read in CEL-file
file1 <- system.file("rawData","FusionSDK_HG-Focus","HG-Focus","2.Calvin","HG-Focus-1-121502.CEL",package="AffymetrixDataTestFiles")
## Calculation of background estimates and expression values (concentrations)
result <- ilm(file1)
## Background intensities of all probes
getIntens(result)
## Background intensities for one or more probesets
getIntens(result,"203561_at") 
getIntens(result,c("203561_at","40359_at"))

## Get concentrations of all probesets
getConcs(result)
## Get concentrations (in picoMolar)
getConcs(result,"203561_at")
getConcs(result,c("203561_at","40359_at"))

## Parameter
getParams(result)

## Subsetting
result["203561_at"] 

}

\keyword{ methods }
