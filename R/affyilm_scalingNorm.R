# Normalization by rescaling: rescale intensities of chips to median intensity value of reference chip (by default 1st chip read in)

.scalingNormalization <- function(allIntensities,NIntens){
# Median values per column ie CEL-file
	medallIntens <- apply(allIntensities,2,median,na.rm=TRUE)
# Scaling
	sf1 <- .scaleIntens(medallIntens)
	Iscaled1 <- allIntensities*rep(sf1,rep(NIntens[[1]],NIntens[[2]]))

	return(Iscaled1)
}



.scaleIntens <- function(scaleVal){
	refval <- scaleVal[1]
	sf1 <- mapply(function(A,ref){ref/A},A=scaleVal,ref=refval)

	return(sf1)
}




