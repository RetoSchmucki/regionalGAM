
flight_curve <- function(your_dataset) {

if("mgcv" %in% installed.packages() == "FALSE") {
   print("mgcv package is not installed.")
   x <- readline("Do you want to install it? Y/N")
   if (x == 'Y') { 
        install.packages("mgcv")}
   if (x == 'N') {
        stop("flight curve can not be computed without the mgcv package, sorry")
    }
}
your_dataset$DAYNO <- strptime(paste(your_dataset$DAY, your_dataset$MONTH,
    your_dataset$YEAR, sep = "/"), "%d/%m/%Y")$yday + 1
dataset <- your_dataset[, c("SPECIES", "SITE", "YEAR", "MONTH",
    "DAY", "DAYNO", "COUNT")]
sample_year <- unique(dataset$YEAR)
sample_year <- sample_year[order(sample_year)]
if (length(sample_year) >1 ) {
    for (y in sample_year) {
        dataset_y <- dataset[dataset$YEAR == y, ]
        nsite <- length(unique(dataset_y$SITE))
        # Determine missing days and add to dataset
        sp_data_all <- year_day_func(dataset_y)
        if (nsite > 200) {
            sp_data_all <- sp_data_all[as.character(sp_data_all$SITE) %in% as.character(unique(dataset_y$SITE)[sample(1:nsite,
                200, replace = F)]), ]
            sp_data_all <- sp_data_all
        }
        sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        print(paste("Fitting the GAM for",as.character(sp_data_all$SPECIES[1]),"and year",y,"with",length(unique(sp_data_all$SITE)),"sites :",Sys.time()))
        if(length(unique(sp_data_all$SITE))>1){
           gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE) -1,
              data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        } 
        else {
            gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1,
                data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        }
        # Give a second try if the GAM does not converge the first time
        if (class(gam_obj_site)[1] == "try-error") {
        # Determine missing days and add to dataset
            sp_data_all <- year_day_func(dataset_y)
            if (nsite > 200) {
            	sp_data_all <- sp_data_all[as.character(sp_data_all$SITE) %in% as.character(unique(dataset_y$SITE)[sample(1:nsite,
            	200, replace = F)]), ]        
            } 
            else {
                sp_data_all <- sp_data_all
            }
            sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        print(paste("Fitting the GAM for",sp_data_all$SPECIES[1],"at year", y,"with",length(unique(sp_data_all$SITE)),"sites :",Sys.time(),"second try"))
        if(length(unique(sp_data_all$SITE))>1){
            gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE) -1,
            	data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        }
        else {
            gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1,
                data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        }
        if (class(gam_obj_site)[1] == "try-error") {
            print(paste("OUPS! Flight period for",sp_data_all$SPECIES[1],"at year", y, " was not computed, the GAM did not converge after two trials"))
            sp_data_all[, "FITTED"] <- NA
            sp_data_all[, "COUNT_IMPUTED"] <- NA
            sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- NA
            sp_data_all[, "NM"] <- NA
        } 
        else {
            # Generate a list of values for all days from the additive model and use
            # these value to fill the missing observations
            sp_data_all[, "FITTED"] <- mgcv::predict.gam(gam_obj_site, newdata = sp_data_all[,
                c("trimDAYNO", "SITE")], type = "response")
            # force zeros at the beginning end end of the year
            sp_data_all[sp_data_all$trimDAYNO < 60, "FITTED"] <- 0
            sp_data_all[sp_data_all$trimDAYNO > 305, "FITTED"] <- 0
            sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
            sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
            # Define the flight curve from the fitted values and append them over
            # years (this is one flight curve per year for all site)
            site_sums <- aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE),
                FUN = sum)
            # Rename sum column
            names(site_sums)[names(site_sums) == "x"] <- "SITE_YR_FSUM"
            # Add data to sp_data data.frame (ensure merge does not sort the data!)
            sp_data_all = merge(sp_data_all, site_sums, by <- c("SITE"),
                all = TRUE, sort = FALSE)
            # Calculate normalized values
            sp_data_all[, "NM"] <- sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        }
        } 
        else {
            # Generate a list of values for all days from the additive model and use
            # these value to fill the missing observations
            sp_data_all[, "FITTED"] <- mgcv::predict.gam(gam_obj_site, newdata = sp_data_all[,
                c("trimDAYNO", "SITE")], type = "response")
            # force zeros at the beginning end end of the year
            sp_data_all[sp_data_all$trimDAYNO < 60, "FITTED"] <- 0
            sp_data_all[sp_data_all$trimDAYNO > 305, "FITTED"] <- 0
            sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
            sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
            # Define the flight curve from the fitted values and append them over
            # years (this is one flight curve per year for all site)
            site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE),
                FUN = sum)
            # Rename sum column
            names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"
            # Add data to sp_data data.frame (ensure merge does not sort the data!)
            sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE,
                sort = FALSE)
            # Calculate normalized values
            sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        }
        sp_data_filled <- sp_data_all
        flight_curve <- data.frame(species = sp_data_filled$SPECIES, year = sp_data_filled$YEAR,
            week = sp_data_filled$WEEK, DAYNO = sp_data_filled$DAYNO, DAYNO_adj = sp_data_filled$trimDAYNO,
            nm = sp_data_filled$NM)[!duplicated(paste(sp_data_filled$YEAR,
            sp_data_filled$DAYNO, sep = "_")), ]
        flight_curve <- flight_curve[order(flight_curve$DAYNO), ]
        # bind if exist else create
        if ("flight_pheno" %in% ls()) {
            flight_pheno <- rbind(flight_pheno, flight_curve)
        } 
        else {
        flight_pheno <- flight_curve
        }
    }  # end of year loop
} 
else {
    y <- unique(dataset$YEAR)
    dataset_y <- dataset[dataset$YEAR == y, ]
    nsite <- length(unique(dataset_y$SITE))
    # Determine missing days and add to dataset
    sp_data_all <- year_day_func(dataset_y)
    if (nsite > 200) {
        sp_data_all <- sp_data_all[as.character(sp_data_all$SITE) %in% as.character(unique(dataset_y$SITE)[sample(1:nsite,
        200, replace = F)]), ]
    } 
    else {
        sp_data_all <- sp_data_all
    }
    sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
    print(paste("Fitting the GAM for",sp_data_all$SPECIES[1],"at year", y,":",Sys.time()))
    if(length(unique(sp_data_all$SITE))>1){
        gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE) -1,
    	      data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
    } 
    else {
        gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1,
        data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
    }
    # Give a second try if the GAM does not converge the first time
    if (class(gam_obj_site)[1] == "try-error") {
        # Determine missing days and add to dataset
        sp_data_all <- year_day_func(dataset_y)

        if (nsite > 200) {
            sp_data_all <- sp_data_all[as.character(sp_data_all$SITE) %in% as.character(unique(dataset_y$SITE)[sample(1:nsite,
            	200, replace = F)]), ]
        } else {
            sp_data_all <- sp_data_all
        }
        sp_data_all$trimDAYNO <- sp_data_all$DAYNO - min(sp_data_all$DAYNO) + 1
        print(paste("Fitting the GAM for",sp_data_all$SPECIES[1],"at year", y,"with",length(unique(sp_data_all$SITE)),"sites :",Sys.time(),"second try"))
        if(length(unique(sp_data_all$SITE))>1){
            gam_obj_site <- try(mgcv::bam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE) - 1,
            data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        } 
        else {
               gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1,
               data = sp_data_all, family = poisson(link = "log")), silent = TRUE)
        }
        if (class(gam_obj_site)[1] == "try-error") {
            print(paste0("OUPS! Flight period for",sp_data_all$SPECIES[1],"at year", y, " was not computed, the GAM did not converge after two trial"))
            sp_data_all[, "FITTED"] <- NA
            sp_data_all[, "COUNT_IMPUTED"] <- NA
            sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- NA
            sp_data_all[, "NM"] <- NA
        } 
        else {
            # Generate a list of values for all days from the additive model and use
            # these value to fill the missing observations
            sp_data_all[, "FITTED"] <- mgcv::predict.gam(gam_obj_site, newdata = sp_data_all[,
            c("trimDAYNO", "SITE")], type = "response")
            # force zeros at the beginning end end of the year
            sp_data_all[sp_data_all$trimDAYNO < 60, "FITTED"] <- 0
            sp_data_all[sp_data_all$trimDAYNO > 305, "FITTED"] <- 0
            sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
            sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
            # Define the flight curve from the fitted values and append them over
            # years (this is one flight curve per year for all site)
            site_sums <- aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE),
                FUN = sum)
            # Rename sum column
            names(site_sums)[names(site_sums) == "x"] <- "SITE_YR_FSUM"
            # Add data to sp_data data.frame (ensure merge does not sort the data!)
            sp_data_all = merge(sp_data_all, site_sums, by <- c("SITE"),
                all = TRUE, sort = FALSE)
            # Calculate normalized values
            sp_data_all[, "NM"] <- sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
            }

        } 
        else {
            # Generate a list of values for all days from the additive model and use
            # these value to fill the missing observations
            sp_data_all[, "FITTED"] <- mgcv::predict.gam(gam_obj_site, newdata = sp_data_all[,
            c("trimDAYNO", "SITE")], type = "response")
            # force zeros at the beginning end end of the year
            sp_data_all[sp_data_all$trimDAYNO < 60, "FITTED"] <- 0
            sp_data_all[sp_data_all$trimDAYNO > 305, "FITTED"] <- 0
            sp_data_all[, "COUNT_IMPUTED"] <- sp_data_all$COUNT
            sp_data_all[is.na(sp_data_all$COUNT), "COUNT_IMPUTED"] <- sp_data_all$FITTED[is.na(sp_data_all$COUNT)]
            # Define the flight curve from the fitted values and append them over
            # years (this is one flight curve per year for all site)
            site_sums = aggregate(sp_data_all$FITTED, by = list(SITE = sp_data_all$SITE),
            FUN = sum)
            # Rename sum column
            names(site_sums)[names(site_sums) == "x"] = "SITE_YR_FSUM"
            # Add data to sp_data data.frame (ensure merge does not sort the data!)
            sp_data_all = merge(sp_data_all, site_sums, by = c("SITE"), all = TRUE,
            sort = FALSE)
            # Calculate normalized values
            sp_data_all[, "NM"] = sp_data_all$FITTED/sp_data_all$SITE_YR_FSUM
        }
    sp_data_filled <- sp_data_all
    flight_curve <- data.frame(species = sp_data_filled$SPECIES, year = sp_data_filled$YEAR,
    week = sp_data_filled$WEEK, DAYNO = sp_data_filled$DAYNO, DAYNO_adj = sp_data_filled$trimDAYNO,
    nm = sp_data_filled$NM)[!duplicated(paste(sp_data_filled$YEAR,
    sp_data_filled$DAYNO, sep = "_")), ]
    flight_curve <- flight_curve[order(flight_curve$DAYNO), ]
        # bind if exist else create
        if ("flight_pheno" %in% ls()) {
            flight_pheno <- rbind(flight_pheno, flight_curve)
        } 
        else {
            flight_pheno <- flight_curve
        }
    }
    return(flight_pheno)
}
