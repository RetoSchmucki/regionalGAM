## regionalGAM

With a rapid expansion of monitoring efforts and the usefulness of conducting integrative analyses to inform conservation initiatives, the choice of a robust abundance index is crucial to adequately assess the species status. Butterfly Monitoring Schemes (BMS)operate in increasing number of countries with broadly the same methodology, yet they differ in their observation frequencies and often in the method used to compute annual abundance indices.

Here we implemented the method for computing an abundance index with the *regional GAM* approach, an extension of the two-stages model introduced by Dennis et al. (2013). This index offers the best performance for a wide range of sampling frequency, providing greater robustness and unbiased estimates then the popular linear interpolation approach (Schmucki et al. 2015).

#### Installation

To install this package from GitHub, you will fist need to install the package `devtools` that is available from CRAN. From there, simply use the the function `install_github()` to install the `RegionalGAM` pacakge on your system. Note that this package was build with R 3.2, so you might you might have to update your R installation. If you are unable to install this package, you might consider sourcing the R script that can be found 

```
	install.packages("devtools")

	library(devtools)
	install_github("RetoSchmucki/regionalGAM")
```


	library(RegionalGAM)


	# load count data for the gatekeeper in the Cold temperate and Moist region
	data("gatekeeper_CM")
	head(gatekeeper_CM)
	
	# format data to compute the flight curve
	dataset1 <- gatekeeper_CM[,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")]

	# compute the annual flight curve for the regional dataset
	# WARNING, this may take some time.
	pheno <- flight_curve(dataset1)
	
	# plot pheno for year 2005
	plot(pheno$DAYNO[pheno$year==2005],pheno$nm[pheno$year==2005],pch=19,cex=0.7,type='o',col='red',xlab="day",ylab="relative abundance")

	# format data to compute abundance indices
		# Note: here we restrict index computation to the sites used in the trend analysis in the paper entitled,
		# A Regionally informed abundance index for supporting integrative analyses across butterfly monitoring schemes
		# But this can be extended to all sites
	dataset2 <- gatekeeper_CM[gatekeeper_CM$TREND==1,c("SPECIES","SITE","YEAR","MONTH","DAY","COUNT")]
	
	# compute the annual abundance indices
	data.index <- abundance_index(dataset2, pheno)



