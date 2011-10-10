### Main function ###

ilm<-function(celfiles,satLim=10000,scale.method="linear",scale.target="concs",cdf.name=NULL,probe.table=NULL,probe.name=NULL,na.replace=NULL,bgcorrect=F,summarize.level="none",summary.method="none",summary.na.rm=TRUE,dgDRpairs=NULL,dgRRpairs=NULL,beta=NULL){

#check1<-proc.time()
	Ncel <- length(celfiles)
	filenames <- sub("^/?([^/]*/)*", "", celfiles) ## R-versions <=2.11.0 additional argument extended=TRUE
	missing <- !file.exists(celfiles)
	if(any(missing)) {
		missing <- paste(celfiles[missing], collapse=", ")
		stop("Cannot read CEL files. Some files not found: ", missing)
	}
# Get chip information
	idChip <- .identifyChip(celfiles)
	#should we use alternatives cdf and probe packages?
	if(is.null(cdf.name))
	{
		chiptype <- idChip[[1]]
	}else{
		chiptype <- cdf.name
	}
	chipdim <- idChip[[2]]
	ncols <- idChip[[3]]
	nrows <- idChip[[4]]
	#is the probe table provided by a data frame? an alternative package? or the package associated to the cdf?
	if(is.null(probe.table))
	{
		if(!is.null(probe.name))
		{
		getProbePackage(probe.name)
		probe.table <- get(probe.name)
		}else{
		moreChipInfo <- .getChipInfo(chiptype)
		probe.table <- moreChipInfo[[1]]
		}
	}

# Probe index
	orig<-probe.table$y*ncols+probe.table$x+1	

# Sequence information (put X everywhere then replace with probe sequences)
	seqsPM <- probe.table$sequence
	
# Probesets	
	setsPM <- probe.table$Probe.Set.Name
	genes <- levels(factor(setsPM))

# Read Intensities
# value: intensity matrix of dimension (#features)x(#celfiles)
	cat(paste("Reading intensities..."))
	rawIntensities <- cbind(readCelIntensities(celfiles))
#check2<-proc.time()
#print("1")
#print(check2-check1)
	Ipm<-cbind(rawIntensities[orig,])
	print("...done",quote=F)
	if(scale.target=="intens" && Ncel>1)
	{
		message("Scaling intensities between arrays...")
		Ipm<-cbind(.affyILM_scale(rawIntensities,method=scale.method))
				message("...done")
	}
	colnames(Ipm)<-filenames
	rownames(Ipm)<-probe.table$Probe.Set.Name
	I0<-0
	rm(rawIntensities)
#check3<-proc.time()
#print("2")
#print(check3-check2)
	#computing DeltaG and DeltaGp
	if(is.null(dgDRpairs))
	{
		#sugimoto 37degres
		dgDRpairs<-c(2.906853,2.696777,1.101192,1.313245,1.692877,2.104984,0.890439,0.905607,1.592308,1.498023,0.211268,0.604984,1.811945,2.085401,0.887893,1.008153)
		names(dgDRpairs)<-c("CC","GC","AC","TC","CG","GG","AG","TG","CA","GA","AA","TA","CT","GT","AT","TT")
	}
	if(is.null(dgRRpairs))
	{
		#xia 37degres
		dgRRpairs<-c(3.248749,3.436203,2.251165,2.360775,2.359529,3.248749,2.075477,2.097503,2.097503,2.360775,0.927530,1.332335,2.075477,2.251165,1.099529,0.927530)
		names(dgRRpairs)<-c("CC","GC","AC","TC","CG","GG","AG","TG","CA","GA","AA","TA","CT","GT","AT","TT")
	}
	if(is.null(beta))
	{
		beta<-0.74
	}
	DeltaG<-.computeDG(seqsPM,dgDRpairs,dgRRpairs)
	colnames(DeltaG)<-c("dgDR","dgRR")
	result<-.computeConcs(Ipm,I0,dgDR=DeltaG[,"dgDR"],dgRR=DeltaG[,"dgRR"],sat=satLim,beta)
	probe.concs<-cbind(result$probe.concs)
	alpha<-result$alpha
	names(alpha)<-probe.table$Probe.Set.Name
	rownames(probe.concs)<-probe.table$Probe.Set.Name
#check4<-proc.time()
#print("3")
#print(check4-check3)
	if((is.null(na.replace)|sum(is.na(na.replace))))
	{
		probe.concs[probe.concs<=0.0] <- NA
	}else{
		probe.concs<-cbind(apply(probe.concs,2,.ilmNAreplace,na.replace,setsPM))
		rownames(probe.concs)<-probe.table$Probe.Set.Name
		colnames(probe.concs)<-filenames
	}
if(scale.target=="concs" && Ncel>1)
			{
				message("Scaling concentrations between arrays...")
				probe.concs.scaled<-.affyILM_scale(probe.concs,method=scale.method)
				rownames(probe.concs.scaled)<-rownames(probe.concs)
				colnames(probe.concs.scaled)<-colnames(probe.concs)
			}else{
				probe.concs.scaled<-probe.concs				
			}
#check6<-proc.time()
#print("5")
#print(check6-check5a)
	exprSummary<-list(Probe.Set=NA,Cluster.Set=NA)
	se.exprSummary<-list(Probe.Set=NA,Cluster.Set=NA)
	if(summarize.level!="none")
	{
	message("Summarizing ...")
	}
	if((summarize.level=="probeset" | summarize.level=="both") & summary.method=="median")
	{
	message("...Probe.Set level")
	summary.list<-lapply(c(1:Ncel),.summarizeConcs,probe.concs.scaled,setsPM,genes,filenames)
	concs.scaled.list<-lapply(summary.list,.listCbind,"concMat")
	concs.scaled<-do.call(cbind,concs.scaled.list)
	se.concs.scaled.list<-lapply(summary.list,.listCbind,"se.concMat")
	se.concs.scaled<-do.call(cbind,se.concs.scaled.list)
	rownames(concs.scaled)<-paste("ps",genes,"id",sep=".")
	colnames(concs.scaled)<-filenames
	rownames(se.concs.scaled)<-rownames(concs.scaled)
	colnames(se.concs.scaled)<-colnames(concs.scaled)
	probe.concs<-probe.concs.scaled
	probenames.mat<-cbind(rownames(probe.concs.scaled),as.character(setsPM))
	probenames<-apply(probenames.mat,1,paste,collapse=".")
	rownames(probe.concs.scaled)<-probenames
	exprSummary[["Probe.Set"]]<-concs.scaled
	rm(concs.scaled)
	se.exprSummary[["Probe.Set"]]<-se.concs.scaled
	rm(se.concs.scaled)
	message("... done")
	}
	if((summarize.level=="probeset" | summarize.level=="both") & sum(summary.method==c("medpolish","tmedpolish")))
	{
	message("...Probe.Set level")
	index.psets<-lapply(genes,.getPsetIndex,setsPM)
	names(index.psets)<-genes
	if(summary.method=="medpolish")
	{
	message("......MedianPolish")
	data.exprs<-lapply(index.psets,.summarizeMedpolish,log2(probe.concs.scaled),na.rm=summary.na.rm)
	}else{
	message("......Transposed MedianPolish")
	data.exprs<-lapply(index.psets,.summarizeTrMedpolish,log2(probe.concs.scaled),na.rm=summary.na.rm)
	}
	data.exprs<-do.call(rbind,data.exprs)
	probenames.mat<-cbind(rownames(probe.concs.scaled),as.character(setsPM))
	probenames<-apply(probenames.mat,1,paste,collapse=".")
	rownames(probe.concs.scaled)<-probenames
	rownames(data.exprs)<-paste("ps",genes,"id",sep=".")
	colnames(data.exprs)<-c(filenames,filenames)
	exprSummary[["Probe.Set"]]<-data.exprs[,1:Ncel]
	se.exprSummary[["Probe.Set"]]<-data.exprs[,c((Ncel+1):ncol(data.exprs))]	
	}
	if((summarize.level=="cluster" | summarize.level=="both") & sum(colnames(probe.table)=="Cluster.Set.Name") & summary.method=="median")
	{
	message("...Cluster.Set level")
	clustersetsPM<-probe.table$Cluster.Set.Name
	probenames.mat<-cbind(rownames(probe.concs.scaled),as.character(clustersetsPM))
	probenames<-apply(probenames.mat,1,paste,collapse=".")
	rownames(probe.concs.scaled)<-probenames
	clustergenes<-levels(factor(clustersetsPM))
	summary.list<-lapply(c(1:Ncel),.summarizeConcs,probe.concs.scaled,setsPM=clustersetsPM,genes=clustergenes,filenames)
	concs.scaled.list<-lapply(summary.list,.listCbind,"concMat")
	concs.scaled<-do.call(cbind,concs.scaled.list)
	se.concs.scaled.list<-lapply(summary.list,.listCbind,"se.concMat")
	se.concs.scaled<-do.call(cbind,se.concs.scaled.list)				
	rownames(concs.scaled)<-paste("ps",clustergenes,"id",sep=".")
	colnames(concs.scaled)<-filenames
	rownames(se.concs.scaled)<-rownames(concs.scaled)
	colnames(se.concs.scaled)<-colnames(concs.scaled)
	probe.concs<-probe.concs.scaled
	exprSummary[["Cluster.Set"]]<-concs.scaled
	se.exprSummary[["Cluster.Set"]]<-se.concs.scaled
	message("... done")
	}
	if((summarize.level=="cluster" | summarize.level=="both") & sum(colnames(probe.table)=="Cluster.Set.Name") & sum(summary.method==c("medpolish","tmedpolish")))
	{
	message("...Cluster.Set level")
	index.psets<-lapply(clustergenes,.getPsetIndex,clustersetsPM)
	names(index.psets)<-clustergenes
	if(summary.method=="medpolish")
	{
	message("......MedianPolish")
	data.exprs<-lapply(index.psets,.summarizeMedpolish,log2(probe.concs.scaled),na.rm=summary.na.rm)
	}else{
	message("......Transposed MedianPolish")
	data.exprs<-lapply(index.psets,.summarizeTrMedpolish,log2(probe.concs.scaled),na.rm=summary.na.rm)	
	}
	data.exprs<-do.call(rbind,data.exprs)
	rownames(data.exprs)<-paste("ps",clustergenes,"id",sep=".")
	colnames(data.exprs)<-c(filenames,filenames)
	probenames.mat<-cbind(rownames(probe.concs.scaled),as.character(clustersetsPM))
	probenames<-apply(probenames.mat,1,paste,collapse=".")
	rownames(probe.concs.scaled)<-probenames
	exprSummary[["Cluster.Set"]]<-data.exprs[,1:Ncel]
	se.exprSummary[["Cluster.Set"]]<-data.exprs[,c((Ncel+1):ncol(data.exprs))]	
	}
					
#check7<-proc.time()
#print("6")
#print(check7-check6)

pmindex<-orig
rownames(probe.concs)<-paste("ps",rownames(probe.concs),"id",sep=".")
colnames(probe.concs)<-filenames
names(pmindex)<-probe.table$Probe.Set.Name
rownames(Ipm)<-rownames(probe.concs)
rownames(DeltaG)<-rownames(probe.concs)

#check8<-proc.time()
#print("7")
#print(check8-check7)

		new(Class="ILM",Ipm=Ipm,I0=as.matrix(I0),exprSummary=exprSummary,se.exprSummary=se.exprSummary,probe.concs=probe.concs,satLim=satLim,deltaG.pm=cbind(DeltaG[,"dgDR"]),deltaGp.pm=cbind(DeltaG[,"dgRR"]),info=list(ncol=ncols,nrow=nrows,cdfName=chiptype,pmindex=pmindex,alpha=alpha,probe.table=probe.table))

}
