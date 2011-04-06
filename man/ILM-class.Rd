\name{ILM-class}
\Rdversion{1.1}
\docType{class}
\alias{ILM-class}
\alias{ILM}

\title{Class to contain the results of an ilm calculation}
\description{ S4 object to contain background intensities, results 
from parameter fit and concentrations}
\section{Slots}{
	 \describe{
    \item{\code{intens}:}{A \code{"matrix"} of size (number of probes) x (2x number of CEL-files). A set of two columns corresponds to one CEL-file where the first column contains the raw PM and the second one the background intensities per feature.}
    \item{\code{concs}:}{A \code{"matrix"} of size (number of probe sets) x (number of CEL-files) holding the concentration per probe set (picoMolar) calculated according to the Langmuir model.} 
    \item{\code{params}:}{A \code{"matrix"} holding the 50 optimized parameter values/fitted parameters for each CEL-file (number of columns)}
    \item{\code{satLim}:}{Is the \code{"numeric"} saturation limit of the intensities of the Langmuir Isotherm, i.e. where the concentration is high or the probe-target binding free energies are large. The default value of \code{satLim} is 10000. }

  }
}

\section{Methods}{
  \describe{
    \item{\code{[}}{Subset }
    \item{getIntens}{Get background intensity (of a particular probe set) }
    \item{getConcs}{Get concentration in pM}
    \item{getParams}{Optimized parameters}
    \item{plotIntens}{Plot PM and calculated background intensities per probe set}
	 }
}

\author{Myriam Kroll, Fabrice Berger and Enrico Carlon}


\examples{
showClass("ILM")
}
\keyword{classes}
