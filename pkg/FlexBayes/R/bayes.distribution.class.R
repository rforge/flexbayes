#########################################################
#  Definitions and Method Assignments for Class "bhnm"  #
#  Insightful: NIH Bayes II                             #
#  Kjell Konis 7/29/2                                   #
#########################################################


setClass("bayes.distribution", representation(
	name = "character",
	parameters = "list")
)

setMethod("show", "bayes.distribution", function(object)
	print.bayes.distribution(object))


