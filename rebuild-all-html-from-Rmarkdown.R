################################################################################
# Script to re-build ALL .rmd files in the present directory. 
# This is intended to help with maintaining/updating a 
# directory that has many separate .rmd-originating 
# web content that needs to be updated on a frequent basis...
################################################################################
# Remove the cache and figures directory.
# figure/ is superfluous when creating HTML5 files with everything embedded.
# cache/ is only needed if you are developing a page and will need
#  to re-build several times as you go. For one-off builds don't need it,
#  and for reproducibility, you don't want it hanging around in commmits.
unlink("cache", TRUE)
unlink("figure", TRUE)

# Create vector of .rmd files in this directory
files = list.files(pattern=".rmd")

# Load and run knit2html on all of the .rmd files
library("knitr")
for( i in files ){
	knit2html(i)
}

# Remove the cache and figures directory again.
unlink("cache", TRUE)
unlink("figure", TRUE)
################################################################################