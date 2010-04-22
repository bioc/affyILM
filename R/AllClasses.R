## Definition of class ILM and validation function

.validILMobject = function(object){
		if(!setequal(unique(rownames(object@intens)),rownames(object@concs))){
			stop("Probeset names of slots intens and concs must be identical!")
		}
		return(TRUE)
		}
		



setClass("ILM",
		representation(
		intens="matrix", 
		concs="matrix",
		params="matrix",
		satLim="numeric"),
		validity=.validILMobject
	)




		






