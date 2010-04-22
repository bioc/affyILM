### Calculation of background intensities for each feature and estimation of
### concentration per gene (in picoMolar) ###

## The definition of the local parameters is not uniform, ie
## the global variables ix[9] and iy[9] define the positioning as follows

##            -------------
##            | 8 | 9 | 7 |  PM
##            -------------
##            | 3 | 1 | 2 |  MM
##            -------------
##            | 5 | 6 | 4 |  PM corresponding to MM above
##            -------------

###########################################################

ilm <- function(celfiles,threshold=350,satLim=10000){

	Ncel <- length(celfiles)
	filenames <- sub("^/?([^/]*/)*", "", celfiles) ## R-versions <=2.11.0 additional argument extended=TRUE
	missing <- !file.exists(celfiles)
	if(any(missing)) {
		missing <- paste(celfiles[missing], collapse=", ")
		stop("Cannot read CEL files. Some files not found: ", missing)
	}
	threshold <- as.double(threshold)

# Get chip information
	idChip <- .identifyChip(celfiles)
	chiptype <- idChip[[1]]
	chipdim <- idChip[[2]]
	ncols <- idChip[[3]]
	nrows <- idChip[[4]]

	moreChipInfo <- .getChipInfo(chiptype)
	probetype <- moreChipInfo[[1]]


# Select the correct probe file and sort it according to 
# the global index. We notice that bioconductor does not always
# provide "probe variables" which are sorted according to the
# same criterion. Try for instance library(hgu133aprobe) or
# library(hgu95aprobe)	
	index <- order(probetype$y*ncols+probetype$x+1)
	index_sort <- sort(probetype$y*ncols+probetype$x+1)
	index_all <- c(index_sort,index_sort+ncols)

	v <- c()
	v <- rep(0,times=ncols*nrows)
	
	yt <- probetype$y[index] #all y-values 
	yt[yt%%2!=0] <- 0 # assign 0 to all ODD (normal) y-values of PM's 
	yt[yt!=0] <- 1 # assign 1 to all EVEN y-values (not normal)

	yt -> v[index_sort]
	yt -> v[index_sort+ncols]

# Sequence information 
	seq_null <- paste(rep(LETTERS[24],25),collapse="")
	seq_tot <- rep(seq_null,times=ncols*nrows)
	seqsPM <- probetype$sequence[index]
	seqsPM -> seq_tot[index_sort]
	seqsPM -> seq_tot[index_sort+ncols]

# Probesets	
	setsPM <- probetype$Probe.Set.Name[index]
	genes <- levels(factor(setsPM))

# Read Intensities
# value: intensity matrix of dimension (#features)x(#celfiles)
	cat(paste("Reading intensities..."))
	allIntensities <- readCelIntensities(celfiles)
	NIntens <- dim(allIntensities)
	dimnames(allIntensities)[[2]] <- filenames
	cat(paste("done\n"))
	intensColNames<-c()
# Only if more than one CEL-file is analyzed: between-chip-normalization (scaling)
	if(Ncel!=1){
		scaledIntens <- .scalingNormalization(allIntensities,NIntens)

## Allocate matrices
		Params <- matrix(NA,nrow=50,ncol=Ncel,dimnames=list(c(), filenames ))
		intensMat <- matrix(NA,nrow=length(setsPM),ncol=2*Ncel,dimnames=list(Probeset=setsPM,rep(c("IPM","I0PM"),times=Ncel) ))
		concMat <- matrix(NA,nrow=length(levels(factor(setsPM))),ncol=Ncel,dimnames=list(genes=genes,filenames) )
		for(i in (1:Ncel)) {
			cat(paste("File analyzed: ", filenames[i],"\n"))
			Intensities <- scaledIntens[,i]

			ilm_alles <- .cGetConcsILM(Intensities,chipdim,seq_tot,i,Ncel,threshold,v,satLim)
			I0 <- ilm_alles[[1]]
			Params[,i] <- as.matrix(ilm_alles[[2]])
			Concs <- ilm_alles[[3]]
		
			Concs <- Concs[index_sort]
			Concs[Concs<=0.0] <- NA
			
			Concs <- tapply(Concs,factor(setsPM),median,na.rm=TRUE)
			concMat[,i] <- Concs
	
			intensMat[,2*i-1] <- round(Intensities[index_sort],2)
			intensMat[,2*i] <- round(I0[index_sort],2)
			intensColNames<-c(intensColNames,paste(c("IPM","I0PM"),filenames[i],sep="-"))
		}
	} else {
		i <- 1
		intensMat <- matrix(NA,nrow=length(setsPM),ncol=2,dimnames=list(Probeset=setsPM,Intensities=c("IPM","I0PM")))
		concMat <- matrix(NA,nrow=length(levels(factor(setsPM))),ncol=1,dimnames= list(genes=genes,filenames))

		Intensities <- allIntensities

		ilm_alles <- .cGetConcsILM(Intensities,chipdim,seq_tot,i,Ncel,threshold,v,satLim)
		I0 <- ilm_alles[[1]]
		Params <- as.matrix(ilm_alles[[2]])
		Concs <- ilm_alles[[3]]

		IPM <- allIntensities[index_sort]
		I0PM <- I0[index_sort]
		Concs <- Concs[index_sort]
		Concs[Concs<=0.0] <- NA

		intensMat[,1] <- IPM
		intensMat[,2] <- round(I0PM,2)
		Concs <- tapply(Concs,factor(setsPM),median,na.rm=TRUE)
		concMat[,1] <- Concs
		intensColNames<-c(intensColNames,paste(c("IPM","I0PM"),filenames,sep="-"))
	}
		colnames(intensMat)<-intensColNames
		colnames(Params)<-filenames
		new(Class="ILM",intens=intensMat,params=Params,concs=concMat,satLim=satLim)
}









