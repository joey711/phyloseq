################################################################################
# Validity methods:
#
# These are delicate, because they are effectively at the S4 infrastructure
# level, in between "new" and the constructor. Some of the issues that might
# otherwise go here for a check are handled by the constructors. In many
# cases it desirable to let the constructor handle this, because it allows
# greater flexibility and transparency. These tests should be limited to 
# conditions that are not fixed automatically by the constructors, and/or
# could not be because the deficiency/error is too fundamental. By design,
# we expect the validity errors to cause a fault before (nearly) any action
# by the constructor. 
#
# This is a special case where the accessors are not-used, in favor of the
# S4 @tags. E.g. object@otu_table instead of otu_table(object). This is to avoid
# any complications with the accessors interacting with objects early on. 
# Perhaps this is a mistake, but its a very limited case and won't be difficult
# to change.
#
# Also, for now these are not documented at all at the user-level, 
# and are not expected to ever
# be at the "user-level", so formal documentation probably unnecessary. Lots
# of comments throughout this code will need to compensate.
################################################################################
########################################
# otu_table:
# # # * all values must be numeric (otu_table()-constructor should probably round values by default))
# # # * all values must be >= 0 (no negative abundances)
########################################
validotu_table <- function(object){
	# Both dimensions must have non-zero length.
	if( any(dim(object)==0) ){
		return("\n OTU abundance data must have non-zero dimensions.")
	}
	# Verify that it is numeric matrix
	if( !is.numeric(object@.Data[, 1]) ){
		return("\n Non-numeric matrix provided as OTU abundance table.\n Abundance must be numeric (usually integer).")
	}
	# Further verify all values are positive
	if( any(object@.Data < 0) ){
		return("\n At least one abundance value in otu_table less than zero.\n Abundance must be greater than zero.")
	}
	return(TRUE)
}
## assign the function as the validity method for the otu_table class
setValidity("otu_table", validotu_table)
########################################
########################################
# sample_data:
########################################
validsample_data <- function(object){
	if( any(dim(object)==0) ){
		return("Sample Data must have non-zero dimensions.")
	}
	return(TRUE)
}
## assign the function as the validity method for the sample_data class
setValidity("sample_data", validsample_data)
########################################
########################################
# taxonomyTable:
########################################
# # # * all values must be a character
# # # * at least some non-NULL (or equiv) values
# taxonomyTable validity function
########################################
validTaxonomyTable <- function(object){
	# Both dimensions must have non-zero length.
	if( any(dim(object)==0) ){
		return("\n Taxonomy Table must have non-zero dimensions.")
	}
	# Verify that it is character matrix
	if( !is.character(object@.Data[, 1]) ){
		return("\n Non-character matrix provided as Taxonomy Table.\n Taxonomy must be characters.")
	}
	# Further verify at least one value is non-NA
	if( all(is.na(object@.Data)) ){
		return("\n All values in Taxonomy Table are NA.\n Must have at least one informative value.")
	}	
	return(TRUE)
}
## assign the function as the validity method for the sample_data class
setValidity("taxonomyTable", validTaxonomyTable)
########################################
########################################
# tree:
########################################
# # (Any rules about trees appropriate in this context?)

########################################
########################################
# phyloseq-class:
########################################
# Because data-index complete-matching is checked/enforced by the phyloseq() constructor,
# it should not be checked here, or the constructor will fail validity tests before
# it gets the chance to groom the objects.
# Instead, the validity test can check if there is any intersection of the species names
# and/or sample names, prior to any attempt by the constructor to prune (which would end)
# in a mysterious index error, anyway
########################################
validphyloseq <- function(object){
	# There must be an otu_table
	if( is.null(object@otu_table) ){
		return("\n An otu_table is required for most analysis / graphics in the phyloseq-package")
	}
	# intersection of species-names must have non-zero length
	if( length(intersect_taxa(object)) <= 0 ){
		return(paste("\n Component taxa/OTU names do not match.\n",
			" Taxa indices are critical to analysis.\n Try taxa_names()", sep=""))
	}
	# If there is sample data, check that sample-names overlap
	if( !is.null(object@sam_data) ){
		if( length(intersect(sample.names(object@sam_data), sample.names(object@otu_table))) <= 0 ){
			return("\n Component sample names do not match.\n Try sample_names()")
		}
	}
	return(TRUE)
}
## assign the function as the validity method for the otu_table class
setValidity("phyloseq", validphyloseq)
########################################
