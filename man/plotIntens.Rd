\name{plotIntens}
\docType{methods}
\alias{plotIntens}
\alias{plotIntens,ILM-method}
\alias{plotIntens,ILM,character-method}
\alias{plotIntens,ILM,missing-method}
\title{
plotIntens
}
\description{
Plot the PM intensity and the calculated background of one probe set.
}

\usage{
plotIntens(object,y,z,...)
}
\arguments{
  \item{object}{An object of type \code{\linkS4class{ILM}}}
  \item{y}{A probe set}
  \item{z}{The name of a sample}
  \item{...}{Graphical parameters can be given as arguments to \code{\link[graphics]{par}}}
}

\details{\code{plotIntens} plots one graph for each CEL-file. Note that it is only possible to plot one probe set at a time (if y (or z) is a vector, only the first value is used). If y=NULL or z=NULL, there is no output.}

\author{Myriam Kroll, Fabrice Berger and Enrico Carlon}

\seealso{
\code{\linkS4class{ILM}},\code{\link{plotILM}}}

\examples{
path <- system.file("rawData/FusionSDK_HG-Focus/HG-Focus/2.Calvin/",package="AffymetrixDataTestFiles");
file1 <- paste(path,"HG-Focus-1-121502.CEL", sep="")
file2 <- paste(path,"HG-Focus-2-121502.CEL", sep="")
result2 <- ilm(c(file1,file2))
## plot output
plotIntens(result2,y="203561_at",z="HG-Focus-2-121502.CEL")
}

\keyword{ methods }

