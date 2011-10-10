## Definition of all generic functions

setGeneric("getIntens",signature="object",
		function(object,y)standardGeneric("getIntens")
)

setGeneric("getExprSummary",signature="object",
		function(object,y,z)standardGeneric("getExprSummary")
)

setGeneric("getSDSummary",signature="object",
		function(object,y,z)standardGeneric("getSDSummary")
)

setGeneric("getProbeConcs",signature="object",
		function(object,y)standardGeneric("getProbeConcs")
)

setGeneric("plotIntens", signature="object",
		function(object,y,z,...)standardGeneric("plotIntens")
)

setGeneric("plotILM", signature="object",
		function(object,y,z,...)standardGeneric("plotILM")
)



