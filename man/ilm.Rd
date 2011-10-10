\name{ilm}
\Rdversion{1.1}
\alias{ilm}
\title{affyILM}

\description{
This function is the working horse of the package and is used as an overall
function to calculate the background intensities as well as the concentrations.
}
\usage{
result <-ilm(celfiles,satLim=10000,scale.method="linear",scale.target="concs",cdf.name=NULL,probe.table=NULL,probe.name=NULL,na.replace=list(median,max),bgcorrect=F,summarize.level="none",summary.method="none",summary.na.rm=TRUE,dgDRpairs=NULL,dgRRpairs=NULL,beta=NULL)
}
\arguments{
  \item{celfiles}{
filenames of CEL-files
}
  \item{satLim}{
Saturation Limit of the Langmuir isotherm (determined by scanner). 
Default value 10000.
}
  \item{scale.method}{
Method used to scale the values between arrays. "linear" (default) uses a coefficient to adjust slopes of all pairwise comparisons to unity, and "linear.stat" makes use of both median and median absolute deviation to scale the values.
}
  \item{scale.target}{
Specify on which value the scaling step has to be performed. Authorized values are "intens" and "concs", respectively refering to probe intensities and probe concentrations.
}
  \item{cdf.name}{
Name of a package providing the chip definition file. By default, this information is extracted from the chip-model provided by the CEL files. This parameter can be used to specify alternative definition files.
}
  \item{probe.table}{
A matrix providing probe annotation table, with the same structure than the probe package provided on Bioconductor. Manual creation of probe tables, or subsets, can be used in place of alternative definition files and probe packages. An additional column, named "Cluster.Set.Name", can be added for new whole-gene arrays.
}
  \item{probe.name}{
Name of a package providing the probe annotation table. By default, this information is extracted from the chip-model provided by the CEL files. This parameter can be used to specify alternative probe table packages.
}
  \item{na.replace}{
Experimental. The methods that can be used to replace missing values, for each probe set. na.replace should be a list with 2 slots, each containing a function. The first slot should contain the function used to missing probe concentrations due to missing probe intensities (i.e. median), and the second slot should contain the function used to compute probe concentrations when the missing value is due to intensity values that are higher than the saturation limit defined by satLim (i.e. max). Those function are then used to compute replacement values based on concentration values of the same probe set.
}
  \item{bgcorrect}{
Not used in the current release. Default value is FALSE. This parameter will be use in the next releases to compute the background correction.
}
  \item{summarize.level}{
Parameter for summarization (character). Default value is "none", to avoid this step. Other accepted values are "probeset" for traditionale definition of probe sets, "cluster" if a Cluster.Set.Name column has been manually added to the probe table and provided with parameter \code{probe.table}, or "both" to run summarization using both definitions.
}
  \item{summary.method}{
Parameter for summarization (character). Default value is "none", to avoid this step. Other accepted values are "median" to use the median of probe concentrations, "medpolish" to use the medianpolish on the log2 of concentrations, and "tmedpolish" for the transposed medianpolish procedure.
}
  \item{summary.na.rm}{
Parameter for summarization (logical). Default value is TRUE, sothat to be able to compute summarized values when data contains missing values.
}
  \item{dgDRpairs}{
Experts only. Allow to tune the computation of DNA/RNA hybridization free energies. dgDRpairs is used to provide list of deltaG values for di-nucleotides.
}
  \item{dgRRpairs}{
Experts only. Allow to tune the computation of RNA/RNA hybridization free energies. dgRRpairs is used to provide list of deltaG values for di-nucleotides.
}
  \item{beta}{
Experts only. Allow to tune the computation of concentrations using the Langmuir Isotherm. beta = 1/RT allow to change the effective hybridization temperature.
}

}

\details{
The Langmuir Isotherm is used as a model to estimate probe concentrations from measured PM intensities, thanks to the computation of sequence-specific DNA/RNA hybridization free energies.
}

\value{An object of class \code{\linkS4class{ILM}}}

\note{The AffymetrixDataTestFiles-package must be installed to run examples.}

\references{
KM Kroll, E Carlon and GT Barkema (2009), Linear method for fast background subtraction in oligonucleotide microarrays Algorithms for Molecular Biology 2009, 4:15
G Mulders, GT Barkema and E Carlon, Inverse Langmuir method for oligonucleotide microarray analysis, BMC Bioinformatics (2009) 10, 64

}

\author{Myriam Kroll, Fabrice Berger, Gerard Barkema and Enrico Carlon}


\seealso{\code{\link{getIntens}}, \code{\link{getProbeConcs}}, \code{\link{getExprSummary}}, \code{\link{getSDSummary}}}

\examples{
## Locate and read in CEL-file
path <- system.file("rawData/FusionSDK_HG-Focus/HG-Focus/2.Calvin/",package="AffymetrixDataTestFiles");
file1 <- paste(path,"HG-Focus-1-121502.CEL", sep="")
## Calculation of background estimates and expression values (concentrations)
result <- ilm(file1)
## show all
show(result) 

## per probeset (example probeset randomly chosen)
result["AFFX-r2-P1-cre-5_at"]

## Analysis of two files
file2 <- paste(path,"HG-Focus-2-121502.CEL", sep="")
result2 <- ilm(c(file1,file2))
}

\keyword{ manip }
