.computeDG<-function(sequence,dgDRpairs,dgRRpairs,dgInit=-3.1)
{
	message("Computing hybridization free energies...")
	nnuc<-nchar(sequence[1])
	deltaGpairsDR<-rep(dgInit,length(sequence))
	deltaGpairsRR<-rep(0,length(sequence))
	for(i in 1:(nnuc-1))
		{
		pairs<-substr(sequence,start=i,stop=i+1)
		deltaGpairsDR<-deltaGpairsDR+dgDRpairs[pairs]
		deltaGpairsRR<-deltaGpairsRR+dgRRpairs[pairs]
		}
	message("...done")
	return(cbind(deltaG=deltaGpairsDR,deltaGp=deltaGpairsRR))
}

.computeConcs<-function(Ipm,I0,dgDR,dgRR,sat,beta)
{
	message("Computing concentrations...")
	alpha<-1/(1+exp(-0.686*(46.5-dgRR)))
	ImI0<-Ipm-I0
	temp<-log(alpha)/beta
	temp<-(dgDR+temp)*-beta
	probe.concs<-(10^12)*exp(temp)*ImI0/(sat-ImI0)
	mask<-ImI0<0
	probe.concs[mask]<-(-1)
	mask<-ImI0>sat
	probe.concs[mask]-(-2)
	message("...done")
	return(list(alpha=alpha,probe.concs=probe.concs))
} 

.listCbind<-function(myList,item)
{
	values<-cbind(myList[[item]])
	return(values)
}

#corriger pour les cas Inf et -Inf
#faire 2 fonctions et virer les tests if()
.summarizeMedpolish<-function(pset,data,na.rm)
{
	try(temp<-medpolish(data[pset,],na.rm=na.rm,trace.iter=FALSE),silent=TRUE)
	temp2<-rep(NA,ncol(data))
	try(temp2<-temp$overall+temp$col,silent=TRUE)
	se.temp2<-rep(NA,ncol(data))
	se.temp2<-try(sqrt(rowSums(apply(data[pset,],1,'-',temp2)^2)/(length(temp2)-1)))
	return(c(concMat=temp2,se.concMat=se.temp2))
}

.summarizeTrMedpolish<-function(pset,data,na.rm)
{
	try(temp<-medpolish(t(data[pset,]),na.rm=na.rm,trace.iter=FALSE),silent=TRUE)
	temp2<-rep(NA,ncol(data))
	try(temp2<-temp$overall+temp$row)
	se.temp2<-rep(NA,ncol(data))
	se.temp2<-try(sqrt(rowSums(apply(data[pset,],1,'-',temp2)^2)/(length(temp2)-1)))
	return(c(concMat=temp2,se.concMat=se.temp2))
}

.getPsetIndex<-function(pset.id,psetnames)
{
	return(which(psetnames==pset.id))
}

.summarizeConcs<-function(i,concMatProbe.scaled,setsPM,genes,filenames)
{
	concMat.scaled <- matrix(NA,nrow=length(unique(setsPM)),ncol=1,dimnames= list(genes=genes,filenames[i]))
	se.concMat.scaled <- matrix(NA,nrow=length(unique(setsPM)),ncol=1,dimnames=list(genes=genes,filenames[i]))
	se.Concs.temp <- tapply(concMatProbe.scaled[,i],factor(setsPM),mad,na.rm=TRUE)*1.4826
	Concs.temp <- tapply(concMatProbe.scaled[,i],factor(setsPM),median,na.rm=TRUE)
	concMat.scaled[,1] <- Concs.temp
	se.concMat.scaled[,1] <- se.Concs.temp
	return(list(concMat=concMat.scaled,se.concMat=se.concMat.scaled))
}

.summarizeAffy<-function(concMatProbe.scaled,cdfName,summary.method,pmindex,nrow,ncol)
{
#	methods<-c("medianpolish","avgdiff","mas","liwong","playerout","median")
#	if(!sum(methods==method))
#	{
#		method<-"median"
#	}
#	if(method=="median")
#	{
#
#	}
#	else
#	{
	nfeatures<-nrow*ncol
	probeData<-array(NA,c(nfeatures,ncol(concMatProbe.scaled)))
	colnames(probeData)<-colnames(concMatProbe.scaled)
	rownames(probeData)<-c(1:nfeatures)
	probeData[pmindex,]<-concMatProbe.scaled
	temp.affybatch<-.createAffyBatch(probeData,cdfName=cdfName)
	results<-computeExprSet(temp.affybatch,pmcorrect.method="pmonly",summary.method=summary.method)
	concMat<-exprs(results)
	se.concMat<-NA
#	}
	return(list(ConcMat=concMat,se.concMat=se.concMat))
}

.createAffyBatch<-function (probeData,cdfName) 
{
    nsamples <- ncol(probeData)
    nprobes<-nrow(probeData)
    samplenames<-colnames(probeData)
    #phenoData
	phenoIndex <- data.frame(sample = c(1:nsamples), row.names = samplenames)
    phenoData=new("AnnotatedDataFrame")
    phenoData <- new("AnnotatedDataFrame", data = phenoIndex,varMetadata = data.frame(labelDescription = "index",row.names = "sample"))
    #description
    description <- new("MIAME")
    preproc(description)$filenames <- samplenames
    preproc(description)$affyversion <- NA
    notes(description) <- ""
#    headdetails <- read.celfile.header(as.character(filenames[[1]]))
#    dim.intensity <- headdetails[[2]]
    protocol <- new("AnnotatedDataFrame", data = data.frame(ScanDate = rep(NA,nsamples),row.names = sampleNames(phenoData), stringsAsFactors = FALSE),dimLabels = c("sampleNames", "sampleColumns"))
    exprs<-probeData
    rownames(exprs)<-c(1:(nprobes))
    colnames(exprs)<-samplenames
    data.affy<-new("AffyBatch", exprs = exprs, cdfName = cdfName,phenoData = phenoData, nrow = nprobes, ncol = nsamples,annotation = cleancdfname(cdfName, addcdf = FALSE), protocolData = protocol, description = description,notes = notes)
    return(data.affy)
}

.ilmNAreplace<-function(concMatProbeVec,na.replace,setsPM)
{
	concMatProbeNew<-tapply(concMatProbeVec,factor(setsPM),.ilmNAreplacePset,na.replace,simplify=T)
	concMatProbeVecNew<-unlist(concMatProbeNew,use.names=F)
	index<-tapply(c(1:length(concMatProbeVec)),factor(setsPM),"*",1,simplify=T)
	indexVec<-unlist(index,use.names=F)
	concMatProbe<-concMatProbeVec
	concMatProbe[indexVec]<-concMatProbeVecNew
	return(concMatProbe)
}

.ilmNAreplacePset<-function(values,na.replace)
{
	temp<-values
	mask1<-values==-1
	mask2<-values==-2
	maskall<-as.logical(mask1+mask2)
	if(is.function(na.replace[[1]]))
	{
		na.replace.function<-na.replace[[1]]
		temp[mask1]<-na.replace.function(values[!maskall])
	}else{
		temp[mask1]<-na.replace[[1]]
	}
	if(is.function(na.replace[[2]]))
	{
		na.replace.function<-na.replace[[2]]		
		temp[mask2]<-na.replace.function(values[!maskall])
	}else{
		temp[mask2]<-na.replace[[2]]
	}
	return(temp)
}

.sgrep<-function(pattern,vect)
{
	result<-grep(vect,pattern=pattern,fixed=T)
	return(result)
}
