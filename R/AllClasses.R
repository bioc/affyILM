## Definition of class ILM and validation function		

setClass("ILM",
		representation(
		Ipm="matrix", 
		I0="matrix", 
		exprSummary="list",
		se.exprSummary="list",
		probe.concs="matrix",
		satLim="numeric",
		deltaG.pm="matrix",
		deltaGp.pm="matrix",
		info="list"),
		validity=function(object){
			msg<-NULL
			#checking that all slots are presents
			valid.slots<-c("Ipm","I0","exprSummary","se.exprSummary","probe.concs","satLim","deltaG.pm","deltaGp.pm","info")
			object.slots<-slotNames(object)
			flag<-setequal(valid.slots,object.slots)
			if(!flag){
				msg<-c(msg,"The list of slots does not match ILM object definition")
			}
			#checking that the probeset names are the same between the slots
			pmnames.list<-sapply(rownames(object@Ipm),strsplit,split=".",fixed=T)
			pmnames.mat<-do.call(rbind,pmnames.list)
			if(names(object@exprSummary)[1]=="Cluster.Set")
			{
				exprSummary.test<-object@exprSummary[["Cluster.Set"]]
			}else{
				exprSummary.test<-object@exprSummary[["Probe.Set"]]				
			}
			if(!is.na(exprSummary.test))
			{
			psetnames.list<-sapply(rownames(exprSummary.test),strsplit,split=".",fixed=T)
			psetnames.mat<-do.call(rbind,psetnames.list)
			flag<-setequal(unique(pmnames.mat[,2]),psetnames.mat[,2])
			if(!flag){
				msg<-c(msg,"Probeset names (row names) differ between slots @Ipm and @exprSummary")
				}
			}
			probenames.list<-sapply(rownames(object@probe.concs),strsplit,split=".",fixed=T)
			probenames.mat<-do.call(rbind,probenames.list)
			flag<-setequal(pmnames.mat[,2],probenames.mat[,2])
			if(!flag){
				msg<-c(msg,"Probeset names (row names) differ between slots @Ipm and @probe.concs")
			}
			#checking that the column names and number are correct between slots
			pm.colnames<-colnames(object@Ipm)
			if(!is.na(exprSummary.test))
			{
			exprSummary.names<-colnames(exprSummary.test)
			flag<-(length(pm.colnames)==length(exprSummary.names))
			if(!flag){
				msg<-c(msg,"Incorrect number of columns in slot @Ipm : should equal the number of samples in slot @exprSummary")
				}
			flag<-sum(pm.colnames==exprSummary.names)
			if(!flag){
				msg<-c(msg,"Filenames (col names) differ between slots @Ipm and @exprSummary")
				}
			}			
			if(is.null(msg)){
				return(TRUE)
			}else{
				return(msg)			
			}
		}
	)
