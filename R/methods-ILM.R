## Methods of class ILM

setMethod("getIntens",
		signature="ILM",
		function(object,y){
			if(missing(y)){
				if(length(object@I0)==1)
				{
					I0<-cbind(array(rep(object@I0,nrow(object@Ipm)*ncol(object@Ipm)),c(nrow(object@Ipm),ncol(object@Ipm))))
					colnames(I0)<-colnames(object@Ipm)
					rownames(I0)<-rownames(object@Ipm)
				}else{
					I0<-cbind(object@I0)
					rownames(I0)<-rownames(object@Ipm)
				}
				return(cbind(Ipm=object@Ipm,I0=I0))	
			}else{
				y<-paste(".",y,".",sep="")
				if(length(object@I0)==1)
				{
					I0<-cbind(array(rep(object@I0,nrow(object@Ipm)*ncol(object@Ipm)),c(nrow(object@Ipm),ncol(object@Ipm))))
					colnames(I0)<-colnames(object@Ipm)
					rownames(I0)<-rownames(object@Ipm)
				}else{
					I0<-cbind(object@I0)
					rownames(I0)<-rownames(object@Ipm)
				}				
				idx<- unique(c(unlist(sapply(y,.sgrep,rownames(object@Ipm)))))
				I0<-cbind(I0[idx,,drop=FALSE])
			   .object <- cbind(Ipm=object@Ipm[idx,,drop=FALSE],I0=I0)
		       return(.object)
			}
		}
	)


setMethod("plotIntens",
		signature="ILM",
		function(object,y,z,...,drop=FALSE){
			if(!missing(y)){
				n.psets <- length(y) 
				if(n.psets > 1){
					y <- y[1] 
					message("Warning: Several probesets have been provided!") 
					message("plotIntens() only uses the first one...") 
				}
				yback<-y
				y<-paste(".",y,".",sep="")
				intensColNames <- colnames(getIntens(object)) 
				n <- length(intensColNames)
				if(!missing(z)){
					nf <- length(z) 
					files <- z 
				}else{
					nf <- as.integer(n/2)
					files <- colnames(getExprSummary(object))
				}
				if(nf>1){
					files <- files[1] 
					message("Too many cel-files have been provided!") 
					message("plotIntens() only uses the first array...") 
					nf <- 1 
				}
				ipm <- getIntens(object)[grep(row.names(getIntens(object)),pattern=y,fixed=T),grep(intensColNames,pattern=files[1])[1]]
				i0pm <- getIntens(object)[grep(row.names(getIntens(object)),pattern=y,fixed=T),grep(intensColNames,pattern=files[1])[2]]
				plot(ipm,type="b",col="red",ylim=c(0,max(ipm)),main=paste(yback),sub=paste(files[1]),ylab="Intensity",...)
				if(unique(i0pm)!=0)
				{
				lines(i0pm,type="o",col="blue",ylim=c(0,max(ipm)),...)
				legend("topright", c("IPM","I0"),col=c("red","blue"), pch=c(1,2))
				}else{
				legend("topright", c("IPM"),col=c("red"), pch=c(1,2))					
				}
		
			} else {
				message("No probeset entered!")
			}
		}
)

setMethod("plotILM",
		signature="ILM",
		function(object,y,z,plot.error,...,drop=FALSE){
			if(!missing(y)){
				n.psets <- length(y) 
				if(n.psets > 1){
					y <- y[1] 
					message("Warning: Several probesets have been provided!") 
					message("plotILM() only uses the first one...") 
				}
				yback<-y
				y<-paste(".",y,".",sep="")
				intensColNames <- colnames(getIntens(object))
				deltaGColNames <- colnames(object@deltaG.pm)
				deltaGpColNames <- colnames(object@deltaGp.pm)
				concsColNames <- colnames(object@Ipm)
				n <- length(concsColNames)
				if(!missing(z)){
					nf <- length(z) 
					files <- z 
				}else{
					nf <- n
					files <- colnames(getExprSummary(object))
				}
#				if(nf>1){
#					files <- files[1] 
#					message("Too many cel-files have been provided!") 
#					message("plotILM() only uses the first array...") 
#					nf <- 1 
#				}
				satLim<-object@satLim
				beta<-0.74
				ipm <- getIntens(object)[grep(row.names(getIntens(object)),pattern=y,fixed=T),grep(intensColNames,pattern=files[1])[1]]
				i0pm <- getIntens(object)[grep(row.names(getIntens(object)),pattern=y,fixed=T),grep(intensColNames,pattern=files[1])[2]]
				deltaG <- object@deltaG.pm[grep(row.names(object@deltaG.pm),pattern=y,fixed=T),]
				deltaGp <- object@deltaGp.pm[grep(row.names(object@deltaGp.pm),pattern=y,fixed=T),]
				alpha=1/(1+exp((46.5-deltaGp)*(-0.686)))
				xvals<-(deltaG+(log(alpha)/beta))
				plot(c(20,40),c(50,15000),log="y",type='n',main=paste(yback,files[1],sep=" - "),xlab="deltaG + RT*log(alpha)",ylab="Intensities",...)
				text(x=xvals,y=ipm-i0pm,col="blue",labels=c(1:length(xvals)),cex=0.6)
				flag=unique(i0pm);
				if(length(flag)!=1 & nf==1)
				{
				text(x=xvals,y=ipm,col="grey",labels=c(1:length(xvals)),cex=0.6)
				legend(x="bottomright",legend=c("Ipm : Raw PM Intensities","Ipm-I0 : PM Intens. with background correction"),text.col=c("grey","blue"),col=c("grey","blue"),pch="*")
				}else{
					legend(x="bottomright",legend=c("Ipm : Raw PM Intensities"),text.col=c("blue"),col=c("blue"),pch="*")
				}
				if(nf>1)
				{
					for(i in c(2:length(files)))
					{
					ipm.tmp <- getIntens(object)[grep(row.names(getIntens(object)),pattern=y,fixed=T),grep(intensColNames,pattern=files[i])[1]]
					i0pm.tmp <- getIntens(object)[grep(row.names(getIntens(object)),pattern=y,fixed=T),grep(intensColNames,pattern=files[i])[2]]
					deltaG.tmp <- object@deltaG.pm[grep(row.names(object@deltaG.pm),pattern=y,fixed=T),]
					deltaGp.tmp <- object@deltaGp.pm[grep(row.names(object@deltaGp.pm),pattern=y,fixed=T),]
					alpha.tmp=1/(1+exp((46.5-deltaGp.tmp)*(-0.686)))
					xvals.tmp<-(deltaG.tmp+(log(alpha.tmp)/beta))
					text(x=xvals.tmp,y=ipm.tmp-i0pm.tmp,labels=c(1:length(xvals.tmp)),col="blue",cex=0.6)
					ipm<-cbind(ipm,ipm.tmp);
					i0pm<-cbind(i0pm,i0pm.tmp);
					deltaG<-cbind(deltaG,deltaG.tmp);
					deltaGp<-cbind(deltaGp,deltaGp.tmp);
					alpha<-cbind(alpha,alpha.tmp);
					xvals<-cbind(xvals,xvals.tmp);
					}
				}
				if(nf==1)
				{
				concs <- object@probe.concs[grep(row.names(object@probe.concs),pattern=y,fixed=T),grep(concsColNames,pattern=files[1])]
				}else{
				concs <- object@probe.concs[grep(row.names(object@probe.concs),pattern=y,fixed=T),files]					
				}
				concs.med<-median(concs,na.rm=T)
				concs.mad<-mad(concs,na.rm=T)
				cat(paste("Median = ",round(concs.med,0),"\n"))
				cat(paste("M.a.d. = ",round(concs.mad,2),"\n"))
				#theoretical curve
				if(missing(plot.error))
				{
					plot.error<-2
				}
				dg<-seq(15,40,by=0.2)
				expdg<-(concs.med/1.e12)*exp(dg*beta)
				lang1<-10000*expdg/(1+expdg)
				lines(dg,lang1,lwd=2.0)
				#plot.error =1 is computed on the log, =2 is computed on the concs, =3 is computed on the concs for 2 intervals
				if(plot.error==1)
				{
					errcm<-mad(log(concs),na.rm=T)
					if(!(is.na(errcm) | is.null(errcm)))
					{
						cm<-concs.med
						cmp<-cm*exp(errcm)/1.e12;expdg<-cmp*exp(dg*beta)
						lang1<-10000*expdg/(1+expdg)
						lines(dg,lang1,lty=2,lwd=1.2)
						cmm<-cm*exp(-errcm)/1.e12;expdg<-cmm*exp(dg*beta)
						if(cmm>0)
						{
							lang1<-10000*expdg/(1+expdg)
							lines(dg,lang1,lty=2,lwd=1.2)
						}
					}
				}else{
				#new test
					errcm<-log(mad(concs,na.rm=T))
					if(!(is.na(errcm) | is.null(errcm)))
					{
						cm<-concs.med
						cmp<-(cm+exp(errcm))/1.e12;expdg<-cmp*exp(dg*beta)
						lang1<-10000*expdg/(1+expdg)
						lines(dg,lang1,lty=2,lwd=1.2)
						cmm<-(cm-exp(errcm))/1.e12;expdg<-cmm*exp(dg*beta)
						if(cmm>0)
						{
							lang1<-10000*expdg/(1+expdg)
							lines(dg,lang1,lty=2,lwd=1.2)
						}
					}
					if(plot.error==3)
					{
						errcm<-log(2*mad(concs,na.rm=T))
						if(!(is.na(errcm) | is.null(errcm)))
						{	
							cmp<-(cm+exp(errcm))/1.e12;expdg<-cmp*exp(dg*beta)
							lang1<-10000*expdg/(1+expdg)
							lines(dg,lang1,lty=3,lwd=1.2)
							cmm<-(cm-exp(errcm))/1.e12;expdg<-cmm*exp(dg*beta)
							if(cmm>0)
							{
								lang1<-10000*expdg/(1+expdg)
								lines(dg,lang1,lty=3,lwd=1.2)
							}
						}
					}
				}
				text(35,1000,paste("c = ",floor(cm),"pM"),cex=0.8)
				return(list(Ipm=ipm,I0pm=i0pm,ImI0=ipm-i0pm,deltaG=deltaG,deltaGp=deltaGp,alpha=alpha,deltaGpRTlogA=xvals,Concs=concs))
			} else {
				message("No probeset entered!")
				return()
			}
		}
)

setMethod("getExprSummary",
		signature="ILM",
		function(object,y,z){
			types<-c("Probe.Set","Cluster.Set")
			mask<-0
			if(!missing(z))
				{
					mask<-(types==z)
					idx<-which(mask)
					subset<-types[idx]
				}
			if(missing(z) | sum(mask)==0)
				{
				subset="Probe.Set"
				}
			if(missing(y)){
				return(object@exprSummary[[subset]])					
			}else{
				y<-paste(".",y,".",sep="")
				if(!is.na(object@exprSummary[[subset]]))
				{
#				idx<-grep(rownames(object@exprSummary[[subset]]),pattern=y)
				idx<- unique(c(unlist(sapply(y,.sgrep,rownames(object@exprSummary[[subset]])))))
			    .object <- object@exprSummary[[subset]][idx,,drop=FALSE]
				}else{
					.object<-NA
				}
				return(.object)
			}
		}
)

setMethod("getSDSummary",
		signature="ILM",
		function(object,y,z){
			types<-c("Probe.Set","Cluster.Set")
			mask<-0
			if(!missing(z))
				{
					mask<-(types==z)
					idx<-which(mask)
					subset<-types[idx]
				}
			if(missing(z) | sum(mask)==0)
				{
				subset="Probe.Set"
				}
			if(missing(y)){
				return(object@se.exprSummary[[subset]])					
			}else{
				y<-paste(".",y,".",sep="")
				if(!is.na(object@exprSummary[[subset]]))
				{				
#				idx<-grep(rownames(object@se.exprSummary[[subset]]),pattern=y)
				idx<- unique(c(unlist(sapply(y,.sgrep,rownames(object@se.exprSummary[[subset]])))))
			    .object <- object@se.exprSummary[[subset]][idx,,drop=FALSE]
				}else{
					.object<-NA
				}
				return(.object)
			}
		}
)


setMethod("getProbeConcs",
		signature="ILM",
		function(object,y){
			if(missing(y)){
				return(object@probe.concs)
			}else{
				y<-paste(".",y,".",sep="")
				idx<- unique(c(unlist(sapply(y,.sgrep,rownames(object@probe.concs)))))
			    .object <- object@probe.concs[idx,,drop=FALSE]
				return(.object)
			}
		}
)

setMethod("show",
		signature="ILM",
	        function(object){
        		message("\n","Object of class 'ILM'")
		        message("containing the following slots:","\n")
		        nIntens<-nrow(getIntens(object))
		        message("@Ipm : Probe Intensities","\n")
		        if(nIntens!=0){
            			if(nIntens>40){
			                message("**** printout limited to 40 entries ****","\n")
			                nIntens <- 40
            			}
            			Ipm<-cbind(object@Ipm[1:nIntens,])
            			colnames(Ipm)<-colnames(object@Ipm)
			            print(Ipm,digits=4,quote=F)
		        }else{
			            message("\n","Empty slot")
        		}
		        message("@I0 : Probe Intensities","\n")
		        if(nIntens!=0){
		        	if(nrow(object@I0)!=1)
		        	{
            			if(nIntens>40){
			                message("**** printout limited to 40 entries ****","\n")
			                nIntens <- 40
            			}
            			I0<-cbind(object@I0[1:nIntens,])
						colnames(I0)<-colnames(object@I0)
			            print(I0,digits=4,quote=F)
		        	}else{
		        		print(object@I0,digits=4,quote=F)
		        	}
		        }else{
			            message("\n","Empty slot")
        		}
		        message("\n","@exprSummary : expression values","\n")
				Exprs.sets<-names(object@exprSummary)
				if(sum(Exprs.sets=="Probe.Set"))
					{
						message("\n","...$Probe.Set","\n")
						if(is.matrix(object@exprSummary$Probe.Set))
						{
			        	nExprs<-nrow(object@exprSummary$Probe.Set)
			        	if(nExprs!=0)
			        		{
       						if(nExprs>40)
       							{
		    	            	message("**** printout limited to 40 entries ****","\n")
                				nExprs <- 40
		         				}
				        	print(object@exprSummary$Probe.Set[1:nExprs,],digits=4,quote=F)					
							}
						}else{
							message("\n","Empty slot")
						}
					}
				if(sum(Exprs.sets=="Cluster.Set"))
					{
						message("\n","...$Cluster.Set","\n")
						if(is.matrix(object@exprSummary$Cluster.Set))
						{
			        	nExprs<-nrow(object@exprSummary$Cluster.Set)
			        	if(nExprs!=0)
			        		{
       						if(nExprs>40)
       							{
		    	            	message("**** printout limited to 40 entries ****","\n")
                				nExprs <- 40
		         				}
				        	print(object@exprSummary$Cluster.Set[1:nExprs,],digits=4,quote=F)					
							}
						}else{
							message("\n","Empty slot")
						}
					}
				nsets<-length(Exprs.sets)
				if(!nsets)
				{
			        message("\n","Empty slot")
        		}
		        message("\n","@se.exprSummary : Variability estimates for expression values","\n")
				se.Exprs.sets<-names(object@se.exprSummary)
				if(sum(se.Exprs.sets=="Probe.Set"))
					{
						message("\n","...$Probe.Set","\n")
						if(is.matrix(object@exprSummary$Probe.Set))
						{
			        	nExprs<-nrow(object@exprSummary$Probe.Set)
			        	if(nExprs!=0)
			        		{
       						if(nExprs>40)
       							{
		    	            	message("**** printout limited to 40 entries ****","\n")
                				nExprs <- 40
		         				}
				        	print(object@se.exprSummary$Probe.Set[1:nExprs,],digits=4,quote=F)					
							}
						}else{
							message("\n","Empty slot")
						}
					}
				if(sum(Exprs.sets=="Cluster.Set"))
					{
						message("\n","...$Cluster.Set","\n")
						if(is.matrix(object@exprSummary$Cluster.Set))
						{
			        	nExprs<-nrow(object@exprSummary$Cluster.Set)
			        	if(nExprs!=0)
			        		{
       						if(nExprs>40)
       							{
		    	            	message("**** printout limited to 40 entries ****","\n")
                				nExprs <- 40
		         				}
				        	print(object@se.exprSummary$Cluster.Set[1:nExprs,],digits=4,quote=F)					
							}
						}else{
							message("\n","Empty slot")
						}
					}
				nsets<-length(se.Exprs.sets)
				if(!nsets)
				{
			        message("\n","Empty slot")
        		}
			        message("\n","@satLim : Saturation value","\n")
		        print(object@satLim,quote=F)
		        message("\n","**** saturation limit ****","\n")
		        message("**** End Show (ILM) ****","\n")
        	}
	)


setMethod(f="[",
		signature="ILM",
		function(x,i=NULL,j,...,drop){ 
			tmp.intens <- getIntens(x,i)
			nintens<-ncol(tmp.intens)/2
			tmp.exprSummary<-list(Probe.Set=getExprSummary(x,i,"Probe.Set"),Cluster.Set=getExprSummary(x,i,"Cluster.Set"))
			tmp.se.exprSummary<-list(Probe.Set=getSDSummary(x,i,"Probe.Set"),Cluster.Set=getSDSummary(x,i,"Cluster.Set"))
			tmp.probeconcs<-getProbeConcs(x,i)
			Ipm<-cbind(tmp.intens[,c(1:nintens)])
			I0<-cbind(tmp.intens[,c((nintens+1):(2*nintens))])
			colnames(Ipm)<-colnames(tmp.intens)[c(1:nintens)]
			colnames(I0)<-colnames(tmp.intens)[c((nintens+1):(2*nintens))]			
			y <- new("ILM",Ipm=Ipm,I0=I0,exprSummary=tmp.exprSummary,se.exprSummary=tmp.se.exprSummary,probe.concs=tmp.probeconcs,deltaG.pm=x@deltaG.pm,deltaGp.pm=x@deltaGp.pm,info=x@info,satLim=x@satLim)
			return(y)
		}
	)







