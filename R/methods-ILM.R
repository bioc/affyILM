## Methods of class ILM



setMethod("getParams",
		signature="ILM",
		function(object){
			return(object@params)
		}
	 )


setMethod("getIntens",
		signature="ILM",
		function(object,y){
			if(missing(y)){
				return(object@intens)	
			}else{
		        	m <- match(row.names(object@intens), y)
			        hasmatch <- !is.na(m)
			        idx <- which(hasmatch)[order(m[hasmatch], na.last=NA)]
			        return(object@intens[idx, ])
			}
		}
	)


setMethod(f="plotIntens",
		signature=c(object="ILM"),
		function(object,y,z,...,drop=FALSE){
			if(!missing(y)){
				n.psets <- length(y) 
				if(n.psets > 1){
					y <- y[1] 
					message("Warning: Several probesets have been provided!") 
					message("plotIntens() only uses the first one...") 
				}
				intensColNames <- colnames(getIntens(object)) 
				n <- length(intensColNames)
				if(!missing(z)){
					nf <- length(z) 
					files <- z 
				}else{
					nf <- as.integer(n/2)
					files <- colnames(getParams(object))
				}
				if(nf>1){
					files <- files[1] 
					message("Too many cel-files have been provided!") 
					message("plotIntens() only uses the first array...") 
					nf <- 1 
				}
				ipm <- object@intens[row.names(object@intens)==y,grep(intensColNames,pattern=files[1])[1]]
				i0pm <- object@intens[row.names(object@intens)==y,grep(intensColNames,pattern=files[1])[2]]
				plot(ipm,type="b",col="red",ylim=c(0,max(ipm)),main=paste(y),sub=paste(files[1]),ylab="Intensity",...)
				lines(i0pm,type="o",col="blue",ylim=c(0,max(ipm)),...)
				legend("topright", c("IPM","I0"),col=c("red","blue"), pch=c(1,2))
		
			} else {
				message("No probeset entered!")
			}
		}
)


setMethod("getConcs",
		signature="ILM",
		function(object,y){
			if(missing(y)){
				return(object@concs)
			}else{
				m <- match(row.names(object@concs), y)
			        hasmatch <- !is.na(m)
			        idx <- which(hasmatch)[order(m[hasmatch], na.last=NA)]
			        .object <- object@concs[idx,,drop=FALSE]
				return(.object)
			}
		}
)

setMethod("show",
		signature="ILM",
	        function(object){
        		message("\n","Object of class 'ILM'")
		        message("containing the following slots:","\n")
		        nIntens<-nrow(object@intens)
		        nConcs<-nrow(object@concs)
		        message("@intens : Probe Intensities","\n")
		        if(nIntens!=0){
            			if(nIntens>40){
			                message("**** printout limited to 40 entries ****","\n")
			                nIntens <- 40
            			}
			            print(object@intens[1:nIntens,],digits=4,quote=F)
		        }else{
			            message("\n","Empty slot")
        		}

		        message("\n","@concs : Probeset Concentrations","\n")
		        if(nConcs!=0){
       				if(nConcs>40){
		                message("**** printout limited to 40 entries ****","\n")
                		nConcs <- 40
		         	}
			        print(object@concs[1:nConcs,],digits=4,quote=F)
        		}else{
			        message("\n","Empty slot")
        		}

		        message("\n","@params : Optimized parameters used to compute background 
intensities","\n")
		        print(object@params,digits=4,quote=F)
		        message("\n","**** in total 50 optimized values per CEL-file 
**** ","\n")
		        message("@satLim : Saturation value","\n")
		        print(object@satLim,quote=F)
		        message("\n","**** saturation limit ****","\n")
		        message("**** End Show (ILM) ****","\n")
        	}
	)


setMethod(f="[",
		signature="ILM",
		function(x,i=NULL,j,...,drop){ 
			tmp.intens <- getIntens(x,i)
			tmp.concs <- getConcs(x,i)
			y <- new("ILM",intens=tmp.intens,concs=tmp.concs,params=x@params,satLim=x@satLim)
			return(y)
		}
	)















