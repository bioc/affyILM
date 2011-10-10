\name{ILM-class}
\Rdversion{1.1}
\docType{class}
\alias{ILM-class}
\alias{ILM}

\title{Class to contain the results of an ilm calculation}
\description{ S4 object to contain intensities, probe concentrations, annotations, and summarized expression values.}
\section{Slots}{
	 \describe{
    \item{\code{Ipm}:}{A \code{"matrix"} of size (number of probes) x (number of CEL-files). Each columns corresponds to one CEL-file with the raw PM.}
    \item{\code{I0}:}{A \code{"matrix"} of size 1x1 or (number of probes) x (number of CEL-files). In the current release, I0 is set to 0 (default value). In the future releases, I0 will include estimates of background intensities.}
    \item{\code{probe.concs}:}{A \code{"matrix"} of size (number of probes) x (number of CEL-files) holding the concentration for each probe (picoMolar) computed according to the Langmuir model.} 
    \item{\code{exprSummary}:}{A \code{"list"} with 1 slots, \code{Probe.Set} and \code{Cluster.Set}, each containing either NA of a \code{matrix} of size (number of sets) x (number of CEL-files). Each slot has been included to provide summarized expression values (using medianpolish or transposed medianpolish), computed at the level of Probe Sets (old and new generation arrays), or Cluster Sets (only for whole gene arrays, if probe table contains \code{"Cluster.Set.Name"} definition, added manually)}
    \item{\code{se.exprSummary}:}{A \code{"list"} with 1 slots, \code{Probe.Set} and \code{Cluster.Set}, each containing either NA of a \code{matrix} of size (number of sets) x (number of CEL-files). Each slot has been included to provide standard deviation estimates for summarized expression values, computed for each Probe Set (old and new generation arrays), or Cluster Set (only for whole gene arrays, if probe table contains \code{"Cluster.Set.Name"} definition, added manually)}     
    \item{\code{satLim}:}{Is the \code{"numeric"} saturation limit of the intensities of the Langmuir Isotherm, i.e. where the concentration is high or the probe-target binding free energies are large. The default value of \code{A} is 10000. }
    \item{\code{deltaG.pm}:}{A \code{"matrix"} providing the DNA-RNA hybridization free energies (deltaG) computed for each PM probe sequence. }
    \item{\code{deltaGp.pm}:}{A \code{"matrix"} providing the RNA-RNA hybridization free energies (deltaG) for each PM probe sequence. }
    \item{\code{info}:}{A \code{"list"} containing several annotations : \code{"ncol"} and \code{"nrow"} provide the chip dimensions, \code{"cdfName"} provides information on the array model,\code{"pmindex"} contains the indices of the probes together with the name of the Probe Set, \code{"alpha"} is a matrix with intermediate value used by the Langmuir Model (used by the function \code{plotILM}), and \code{"probe.table"} contains the probe annotation table.}    
  }
}

\section{Methods}{
  \describe{
    \item{\code{[}}{Subset ILM objects}
    \item{getIntens}{Get background intensity (of a particular probe set) }
    \item{getProbeConcs}{Get probes concentrations in pM}
    \item{getExprSummary}{Get summarized expression values}
    \item{getSDSummary}{Get estimates for summarized expression values}
    \item{plotIntens}{Plot probe intensities of selected probe set}
    \item{plotILM}{Plot Langmuir Isotherm for selected probe set}
	 }
}

\author{Myriam Kroll, Fabrice Berger and Enrico Carlon}


\examples{
showClass("ILM")
}
\keyword{classes}
