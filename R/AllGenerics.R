## Definition of all generic functions

setGeneric("getParams",signature="object",
		function(object)standardGeneric("getParams")
)

setGeneric("getIntens",signature="object",
		function(object,y)standardGeneric("getIntens")
)

setGeneric("getConcs",signature="object",
		function(object,y)standardGeneric("getConcs")
)

setGeneric("plotIntens", signature="object",
		function(object,y,z,...)standardGeneric("plotIntens")
)






