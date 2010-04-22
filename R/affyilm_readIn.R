.getChipInfo <- function(chiptype) {

## Following few lines adapted from gcrma-package
# Load library, download automatically if not installed 
	cdf.tag <- cleancdfname(chiptype,addcdf=FALSE)
	probepackagename <- paste(cdf.tag,"probe",sep="")
	getProbePackage(probepackagename)
	cat(paste("Probepackage", probepackagename,"loaded \n"),sep="")
	probetype <- get(probepackagename)

### In case cdf-package is needed... ###
#	cdfpackagename <- paste(cdf.tag,"cdf",sep="") 
#	getCDF(cdfpackagename)

	return(list(probetype,probepackagename))

}


.identifyChip <- function(celfiles) {
	headers <- lapply(celfiles, readCelHeader)
	chiptype <- unique(sapply(headers, "[[", i="chiptype"))

	if(length(chiptype)!=1){
		stop("Different chiptypes!")
	} else {
		ncols <- headers[[1]]$cols
		nrows <- headers[[1]]$rows
		if(all.equal(ncols,nrows)){
			chipdim <- ncols
			cat(paste("Chip dimension",chipdim," x ",chipdim, "\n"),sep="")
		} else {
			stop("Sorry, in the current version of 
affyILM, only chips with ncol=nrow can be analyzed...")
		}
	}

	return(list(chiptype,chipdim,ncols,nrows))
}






















