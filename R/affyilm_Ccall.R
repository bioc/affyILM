# Call to C function to do ILM


.cGetConcsILM <- function(Intensities,chipdim,seq_tot,zaehler,Nfiles,threshold,v,satLim){
	NIntens <- length(Intensities)
	getConcs <- .C("GetConcsILM",
			as.integer(chipdim),
			as.integer(NIntens),
			as.integer(zaehler),
			as.integer(Nfiles),
			as.double(threshold),
			Intens=as.double(Intensities),
			Seq=as.character(seq_tot),
			Params=double(length=50),
			BackI=double(length=NIntens),
			DGDR=double(length=NIntens),
			DGRR=double(length=NIntens),
			concs=double(length=NIntens),
			as.integer(v),
			aLim=as.double(satLim),
			PACKAGE="affyILM"
		      )

	return(list(getConcs$BackI,getConcs$Params,getConcs$concs))
}





















