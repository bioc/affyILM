\name{plotILM}
\docType{methods}
\alias{plotILM}
\alias{plotILM,ILM-method}
\alias{plotILM,ILM,character-method}
\alias{plotILM,ILM,missing-method}
\title{
plotILM
}
\description{
Illustrate the Langmuir Isotherm for selected probe set.
}

\usage{
plotILM(object,y,z,...)
}
\arguments{
  \item{object}{An object of type \code{\linkS4class{ILM}}}
  \item{y}{A probe set}
  \item{z}{The name of a sample}
  \item{...}{Graphical parameters can be given as arguments to \code{\link[graphics]{par}}}
}

\details{\code{plotIntens} plots one graph for each CEL-file. Note that it is only possible to plot one probe set at a time (if y (or z) is a vector, only the first value is used). If y=NULL or z=NULL, there is no output.  Optional plot.error argument (numeric) can be used to define the illustration of the variability. plot.error=1 computes the error on the log scale. plot.error=2 (default value) computes the error on the concentrations. plot.error=3 computes the error on the concentrations, for 2 intervals. Note that error curves can only be plotted if concentration is higher than error (negative concentration does not exist!).}

\author{Myriam Kroll, Fabrice Berger and Enrico Carlon}

\seealso{
\code{\linkS4class{ILM}},\code{\link{plotIntens}}}

\examples{
path <- system.file("rawData/FusionSDK_HG-Focus/HG-Focus/2.Calvin/",package="AffymetrixDataTestFiles");
file1 <- paste(path,"HG-Focus-1-121502.CEL", sep="")
file2 <- paste(path,"HG-Focus-2-121502.CEL", sep="")
result2 <- ilm(c(file1,file2))
## plot output
plotILM(result2,y="203561_at",z="HG-Focus-2-121502.CEL")
}

\keyword{ methods }

