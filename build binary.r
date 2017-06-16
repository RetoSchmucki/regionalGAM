## helpful instruction: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

library(devtools)
library(roxygen2)
setwd("parent_directory")

## example
create("package_name")
## Edit Description
## add R function in the R folder with extention ".R"
## add documentation for each function 

	#' the Function
	#'
	#' This function allows you to express.
	#' @param love Do you love cats? Defaults to TRUE.
	#' @keywords cats
	#' @export
	#' @examples
	#' the_function()

setwd("./package_name")
document()
devtools::use_package("ncdf4")
devtools::use_package("chron")
devtools::use_package("zoo")
devtools::use_package("FNN")
devtools::use_package("tcltk")


setwd("..")

install("package_name")


##ALSO create BINARY for DIFFERENT OS

## WINDOWS
dir.create("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/regionalGAM/binary/WIN",recursive=TRUE)
devtools::build("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/regionalGAM",path="C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/regionalGAM/binary/WIN",vignette=FALSE,manual=TRUE)

## MacOS
dir.create("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/regionalGAM/binary/MacOS",recursive=TRUE)
devtools::build("/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/regionalGAM",path="/Users/retoschmucki/CEH-OneDriveBusiness/OneDrive - Natural Environment Research Council/regionalGAM/binary/MacOS",vignette=FALSE,manual=TRUE)

