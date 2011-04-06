\name{ilm}
\Rdversion{1.1}
\alias{ilm}
\title{affyILM}

\description{
This function is the working horse of the package and is used as an overall
function to calculate the background intensities as well as the concentrations.
}
\usage{
ilm(celfiles, threshold = 350, satLim = 10000)
}
\arguments{
  \item{celfiles}{
filenames of CEL-files
}
  \item{threshold}{
Criterium to filter probes used for background calculation. By default set to
350,~i.e.~use only probes with intensity values below threshold value.
}
  \item{satLim}{
Saturation Limit of the Langmuir isotherm (determined by scanner). 
Default value 10000.
}
}

\details{
A linear function with 50 parameters is used to describe the background intensity.
On behalf of the linear-least square method these 50 parameters are determined of which
18 reflect the influence of the neighboring spots on the background intensity of a particular
feature. The next 16 parameters incorporate the nearest-neighbor free energies (sequence dependence) 
with and the last 16 parameters modify the strength of the sequence dependence based upon
the position of a nucleotide along the sequence.
The parameters are determined through the linear least square method making use of a standard algorithm 
implementing singular value decomposition.
(SVD).
}

\value{An object of class \code{\linkS4class{ILM}}}

\note{The AffymetrixDataTestFiles-package must be installed to run examples.}

\references{
KM Kroll, E Carlon and GT Barkema (2009), Linear method for fast background subtraction in oligonucleotide microarrays, Algorithms for Molecular Biology 2009, 4:15

G Mulders, GT Barkema and E Carlon, Inverse Langmuir method for oligonucleotide microarray analysis, BMC Bioinformatics (2009) 10, 64

E. Carlon and T. Heim, Thermodynamics of RNA/DNA hybridization in high-density oligonucleotide microarrays, Physica A (2006), 362, 433
}

\author{Myriam Kroll, Fabrice Berger, Gerard Barkema and Enrico Carlon}


\seealso{\code{\link{getIntens}}, \code{\link{getConcs}}}

\examples{
## Locate and read in CEL-file
file1 <- system.file("rawData","FusionSDK_HG-Focus","HG-Focus","2.Calvin","HG-Focus-1-121502.CEL",package="AffymetrixDataTestFiles")
## Calculation of background estimates and expression values (concentrations)
result <- ilm(file1)
## show all
show(result) 

## per probeset (example probeset randomly chosen)
result["AFFX-r2-P1-cre-5_at"]

## Analysis of two files
file2 <- system.file("rawData","FusionSDK_HG-Focus","HG-Focus","2.Calvin","HG-Focus-2-121502.CEL",package="AffymetrixDataTestFiles")
result2 <- ilm(c(file1,file2))
}

\keyword{ manip }
