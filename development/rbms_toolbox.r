##==========================================
## Name Convention in the rbms package
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   02.01.2018
##
##      rbms_toolbox: useful functions for the rbms package
##
##==========================================

#' check_package
#' Internal function to verified the required package is installed
#' @param pkgName A string with the package name
#' @param message1 A string to inform about the dependency
#' @param message2 A string to inform what happen if not installed
#' @return If package is not installed, the function ask to install the package.
#' @author Reto Schmucki - retoschm[at]ceh.ac.uk
#' check_package()

check_package <- function(pkgName=NULL, message1='you need to install the package ',message2='This version requires '){
        if (!requireNamespace(pkgName)) {
            print(paste(message1, pkgName))
            x <- readline(paste("Do you want to install",pkgName,"? Y/N"))            
            if (toupper(x) == 'Y') { 
                    install.packages(pkgName)
            }
            if (toupper(x) == 'N') {
                print(paste(message2,pkgName))
            }
        }
    }


#' check_names
#' Verify for the required column names in the data
#' @param x a data.table object with column names
#' @param y vector with the required variable names
#' @return Verify if column names listed in \code{y} vector are found in the data set \code{x}, if not, a message identifies the 
#' missing column name and stops.
#' @details This function is not case sensitive, but it does not accept different names or spelling. 
#' @author Reto Schmucki - retoschm[at]ceh.ac.uk
#' @examples
#' DF <- data.frame(DAY=c(1:5),MONTH=rep(month.name[2],5),YEAR=rep(format(Sys.Date(),'%Y'),5))
#' check_names(DT,c('DAY','month','Years'))
#' check_names()
#'

check_names <- function(x, y){
        dt_names <- y %in% names(x)
        if(sum(dt_names)!=length(y)) {
            stop(paste('You need to have a variable named -', paste(y[!dt_names],collapse=' & '), '- in table',deparse(substitute(x)),'\n'))
        }
    }


#' ts_date_seq
#' Generate a time-series of dates (per day) from the beginning of a starting year to the end of an ending years.
#' @param InitYear start year of the time-series, 4 numbers format (e.g 1987)
#' @param LastYear end year of the time-series, if not provided, current year is used instead 
#' @return return a POSIXct vector with the format 'YYYY-MM-DD HH:MM:SS'
#' @keywords time series
#' @author Reto Schmucki - retoschm[at]ceh.ac.uk
#' @examples
#' ts_date_seq()
#'

ts_date_seq <- function(InitYear=1970,LastYear=format(Sys.Date(),"%Y")) {

        init_date <- as.Date(paste((InitYear-1), "01-01", sep = "-"))
        last_date <- as.Date(paste((as.numeric(LastYear)+1), "12-31", sep = "-"))

        date_series <- as.POSIXct(seq(from=init_date, to= last_date, by = "day"), format = "%Y-%m-%d")
        date_series <- date_series[!format(date_series,'%Y') %in% c((InitYear-1),as.numeric(LastYear)+1)]

        return(date_series)

    }


### check_pheno function check for the flight curve of a specific year and if missing impute the nearest available 
###             within a span of 5 years 

check_pheno <- function(sp_count_flight_y,sp_count_flight){

        if(sp_count_flight_y[is.na(NM),.N]>0){
            tr<-1
            z <- rep(1:5,rep(2,5))*c(-1,1)
            search_op <-sp_count_flight[,unique(as.integer(M_YEAR))]
            valid_y <- c(y+z)[c(y+z)>min(search_op) & c(y+z)<max(search_op)]
            alt_flight <- unique(sp_count_flight[as.integer(M_YEAR)==y,.(M_YEAR,trimDAYNO,NM)])

            while(alt_flight[is.na(NM),.N]>0 & tr<=length(valid_y)){
                alt_flight <- unique(sp_count_flight[as.integer(M_YEAR)==valid_y[tr],.(M_YEAR,trimDAYNO,NM)])
                tr<-tr+1
            }
            
            if(alt_flight[is.na(NM),.N]>0){
                next(paste("No reliable flight curve available within a 5 year horizon of",sp_count_flight_y[1,M_YEAR,]))
            }else{
                warning(paste("We used the flight curve of",alt_flight[1,M_YEAR],"to compute abundance indices for year",sp_count_flight_y[1,M_YEAR,]))
                sp_count_flight_y[,trimDAYNO:=DAY_SINCE-min(DAY_SINCE)+1]
                data.table::setnames(alt_flight,'NM','NMnew')
                alt_flight[,M_YEAR:=NULL]
                data.table::setkey(sp_count_flight_y,trimDAYNO)
                data.table::setkey(alt_flight,trimDAYNO)
                sp_count_flight_y <- merge(sp_count_flight_y,alt_flight,by='trimDAYNO',all.x=TRUE)
                sp_count_flight_y[,NM:=NMnew][,NMnew:=NULL]
            }
        }
        return(sp_count_flight_y)
    }


### fit_glm function to fit and predict daily butterfly counts using the flight curve and the glm method provided in stats package 
### 

fit_glm <- function(sp_count_flight_y,non_zero,FamilyGlm){

            if(sp_count_flight_y[unique(SITE_ID),.N]>1){
                glm_obj_site <- try(glm(COUNT ~ factor(SITE_ID) + offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,],
                family=FamilyGlm, control=list(maxit=100)),silent=TRUE)
            } else {
                glm_obj_site <- try(glm(COUNT ~ offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,],
                family=FamilyGlm, control=list(maxit=100)),silent=TRUE)
            }
             
            if (class(glm_obj_site)[1] == "try-error") {
                sp_count_flight_y[SITE_ID %in% non_zero,c("FITTED","COUNT_IMPUTED"):=.(NA,NA)]
                print(paste("Computation of abundance indices for year",sp_count_flight_y[1,M_YEAR,],"failed with the RegionalGAM, verify the data you provided for that year"))
                next()
            }else{
                sp_count_flight_y[SITE_ID %in% non_zero,FITTED:= predict.glm(glm_obj_site,newdata=sp_count_flight_y[SITE_ID %in% non_zero,],type = "response")]
            }

            sp_count_flight_mod_y <- list(sp_count_flight_y=sp_count_flight_y,glm_obj_site=glm_obj_site)

        return(sp_count_flight_mod_y)
    }

fit_glm.nb <- function(sp_count_flight_y,non_zero){

            if(sp_count_flight_y[unique(SITE_ID),.N]>1){
                glm_obj_site <- try(MASS::glm.nb(COUNT ~ factor(SITE_ID)+offset(NM),data=sp_count_flight_y[SITE_ID %in% non_zero,]),silent=TRUE)
            } else {
                glm_obj_site <- try(MASS::glm.nb(COUNT ~ offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,]),silent=TRUE)
            }

            if (class(glm_obj_site)[1] == "try-error") {
                sp_count_flight_y[SITE_ID %in% non_zero,c("FITTED","COUNT_IMPUTED"):=.(NA,NA)]
                print(paste("Computation of abundance indices for year",sp_count_flight_y[1,M_YEAR,],"failed with the RegionalGAM, verify the data you provided for that year"))
                next()
            }else{
                sp_count_flight_y[SITE_ID %in% non_zero,FITTED:= predict.glm(glm_obj_site,newdata=sp_count_flight_y[SITE_ID %in% non_zero,],type = "response")]
            }

            sp_count_flight_mod_y <- list(sp_count_flight_y=sp_count_flight_y,glm_obj_site=glm_obj_site)

        return(sp_count_flight_mod_y)
    }

### fit_speedglm function to fit and predict daily butterfly counts using the flight curve and the speedglm method provided in the speedglm package
### 

fit_speedglm <- function(sp_count_flight_y,non_zero,FamilyGlm){

            if(sp_count_flight_y[unique(SITE_ID),.N]>1){
                glm_obj_site <- try(speedglm::speedglm(COUNT ~ factor(SITE_ID) + offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,],
                family=FamilyGlm, control=list(maxit=100)),silent=TRUE)
            } else {
                glm_obj_site <- try(speedglm::speedglm(COUNT ~ offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,],
                family=FamilyGlm, control=list(maxit=100)),silent=TRUE)
            }
             
            if (class(glm_obj_site)[1] == "try-error") {
                sp_count_flight_y[SITE_ID %in% non_zero,c("FITTED","COUNT_IMPUTED"):=.(NA,NA)]
                print(paste("Computation of abundance indices for year",sp_count_flight_y[1,M_YEAR,],"failed with the RegionalGAM, verify the data you provided for that year"))
            }else{
                sp_count_flight_y[SITE_ID %in% non_zero,FITTED:= predict(glm_obj_site,newdata=sp_count_flight_y[SITE_ID %in% non_zero,],type = "response")]
            }

            sp_count_flight_mod_y <- list(sp_count_flight_y=sp_count_flight_y,glm_obj_site=glm_obj_site)

        return(sp_count_flight_mod_y)
    }


### impute_count function to compute the Abundance Index across sites and years from 
###                 your count dataset and the regional flight curve

impute_count <- function(ts_season_count,ts_flight_curve,FamilyGlm=quasipoisson(),CompltSeason=TRUE,
                                    SelectYear=NULL,SpeedGlm=FALSE) {
        
        ts_flight_curve <- ts_flight_curve$f_pheno

        check_package('data.table')
        if(isTRUE(SpeedGlm)){
            check_package('speedglm', message2='glm() will be used instead')
        }
        
        if(isTRUE(CompltSeason)){
            ts_season_count <- ts_season_count[COMPLT_SEASON==1]
        }
            
        sp_ts_season_count <- data.table::copy(ts_season_count)
        sp_ts_season_count[,SPECIES:=ts_flight_curve$SPECIES[1]]
        data.table::setkey(sp_ts_season_count,DATE)
        data.table::setkey(ts_flight_curve,DATE)
        sp_count_flight <- merge(sp_ts_season_count,ts_flight_curve[,.(DATE,trimDAYNO,NM)],all.x=TRUE)
        data.table::setkey(sp_count_flight,M_YEAR,DATE,SITE_ID)

        glmMet <- "glm()"
        if(isTRUE(SpeedGlm)){
            glmMet <- "speedglm()"
        }

        if( FamilyGlm[1]=='nb' & isTRUE(SpeedGlm)){
            glmMet <- "glm()"
            SpeedGlm <- FALSE
            cat('SpeedGlm is not implemented with Negative Binomial, we will use glm.nb() from the MASS package instead /n')
        }

        if(is.null(SelectYear)){
            year_series <- ts_season_count[,unique(as.integer(M_YEAR))]
        } else {
            year_series <- ts_season_count[M_YEAR %in% SelectYear,unique(as.integer(M_YEAR))]
        }

        for(y in year_series){
            
            sp_count_flight_y <-  data.table::copy(sp_count_flight[as.integer(M_YEAR)==y,])
            sp_count_flight_y <- check_pheno(sp_count_flight_y,sp_count_flight)

            print(paste("Computing abundance indices for species",sp_count_flight_y[1,SPECIES],"monitored in year", sp_count_flight_y[1,M_YEAR],"across",sp_count_flight_y[unique(SITE_ID),.N],"sites, using",glmMet,":",Sys.time()))

            sp_count_flight_y[M_SEASON==0L,COUNT:=NA]
            sp_count_flight_y[M_SEASON!=0L & NM==0,NM:=0.000001]
            non_zero <- sp_count_flight_y[,sum(COUNT,na.rm=TRUE),by=(SITE_ID)][V1>0,SITE_ID]
            zero <- sp_count_flight_y[,sum(COUNT,na.rm=TRUE),by=(SITE_ID)][V1==0,SITE_ID]

            if(length(non_zero)>=1){
                if(isTRUE(SpeedGlm)){
                    sp_count_flight_l <- fit_speedglm(sp_count_flight_y,non_zero,FamilyGlm)             
                    sp_count_flight_y <- sp_count_flight_l$sp_count_flight_y
                    sp_count_flight_mod <- sp_count_flight_l$glm_obj_site 
                }else{
                    if(FamilyGlm[1]=='nb'){
                    sp_count_flight_l <- fit_glm.nb(sp_count_flight_y,non_zero)    
                    sp_count_flight_y <- sp_count_flight_l$sp_count_flight_y
                    sp_count_flight_mod <- sp_count_flight_l$glm_obj_site
                    }else{
                    sp_count_flight_l <- fit_glm(sp_count_flight_y,non_zero,FamilyGlm)    
                    sp_count_flight_y <- sp_count_flight_l$sp_count_flight_y
                    sp_count_flight_mod <- sp_count_flight_l$glm_obj_site
                    }  
                }
            }

            sp_count_flight_y[SITE_ID %in% zero,FITTED:=0]
            sp_count_flight_y[is.na(COUNT),COUNT_IMPUTED:=FITTED][!is.na(COUNT),COUNT_IMPUTED:=as.numeric(COUNT)][M_SEASON==0L,COUNT_IMPUTED:=0] 

            data.table::setkey(sp_ts_season_count,SITE_ID,DAY_SINCE)
            data.table::setkey(sp_count_flight_y,SITE_ID,DAY_SINCE)

            if("FITTED" %in% names(sp_ts_season_count)){
                sp_ts_season_count[sp_count_flight_y,':='(trimDAYNO=i.trimDAYNO,NM=i.NM,FITTED=i.FITTED,COUNT_IMPUTED=i.COUNT_IMPUTED)]
            }else{
                sp_ts_season_count <- merge(sp_ts_season_count, sp_count_flight_y[,.(DAY_SINCE,SITE_ID,trimDAYNO,NM,FITTED,COUNT_IMPUTED)], all.x=TRUE) 
            }

           if ("imp_glm_model" %in% ls()) {
            glm_model <- list(sp_count_flight_mod)
            names(glm_model) <- paste0('imput_glm_mod_',sp_count_flight_y[1,M_YEAR])
            imp_glm_model <- c(imp_glm_model,glm_model)
           } else { 
            imp_glm_model <- list(sp_count_flight_mod)
            names(imp_glm_model) <- paste0('imput_glm_mod_',sp_count_flight_y[1,M_YEAR])
           }
        }

    if(!is.null(SelectYear)){
        return(list(sp_ts_season_count=sp_ts_season_count[M_YEAR %in% SelectYear,],glm_model=imp_glm_model))
    } else {
        return(list(sp_ts_season_count=sp_ts_season_count,glm_model=imp_glm_model))
    }
} 


### butterfly_day function to count cumulative butterfly count observed over one monitoring season.

butterfly_day <- function(sp_ts_season_count){

            b_day <- sp_ts_season_count[,sum(COUNT_IMPUTED),by=.(SPECIES,M_YEAR,SITE_ID)]
            data.table::setnames(b_day,"V1","BUTTERFLY_DAY")
        
        return(b_day)
    }


### SIMULATE DATA

### sim_emerg_curve() estimates an emergence curve shape following a logistic distribution along a time series [t_series], with peak position [peak_pos] relative along a  , 
### vector using the percentile and a standard deviation around the peak [sd_peak] in days, with two shape parameters
### sigma [sig] for left or right skewness (left when > 1, right when < 1, logistic when  = 1) and bet [bet] for a scale parameter. 

### return a vector of relative emergence along a vector of length t.

### Calabrese, J.M. (2012) How emergence and death assumptions affect count-based estimates of butterfly abundance and lifespan. Population Ecology, 54, 431–442.


sim_emerg_curve <- function (t_series, PeakPos=50, sdPeak=1, sigE=0.15, betE=3) {
            
            u1 <- round(PeakPos * length(t_series)/100) + rnorm(1, 0, sdPeak)
            fE1 <- (sigE * (exp((t_series - u1)/betE)))/(betE * (1 + exp((t_series - u1)/betE)))^(sigE + 1)
            sdfe <- fE1/sum(fE1)  

        return(sdfe)       
    }

### sim_emerg_count simulates emergence of n adults [TotalEmerg] according an emergence curve [sdfe] using a Poisson process

sim_emerg_nbr <- function(sdfe, TotalEmerg=100) {

            n_emerg <- unlist(lapply(TotalEmerg*sdfe,function(x) {rpois(1,x)}))

        return(n_emerg)
    }

### sim_adult_count simulates the number of adults in a population, according to individual maximum life span [max_life] and daily
### mortality risk based on a beta distribution with ShapeA and ShapeB parameters

sim_adult_nbr <- function(n_emerg, MaxLife=15, ShapeA=0.5, ShapeB=0.2){

           c_mat <- matrix(0,nrow=length(n_emerg),ncol=length(n_emerg)+MaxLife+1)
           m_hazard <- c(0,diff(pbeta(seq(0,1,(1/MaxLife)),ShapeA,ShapeB)),1)
        
           for (i in seq_along(n_emerg)){
                s <- n_emerg[i]
                y=2
                l=s

                while(s > 0 & y <= MaxLife+2){
                    s <- s-rbinom(1,s,m_hazard[y])
                    y <- y+1
                    l <- c(l,s)
                }
            
                c_mat[i,i:(i+length(l)-1)] <- l
            }

        return(colSums(c_mat))
    }


### sim_butterfly_count() simulates the number of adult butterfly present during the monitoring seasons for a univoltine or multivoltine species, from month 4 to 9 (April to September)

sim_butterfly_nbr <- function(d_season, CompltSeason=TRUE, GenNumb=1, PeakPos=c(25,75), sdPeak=c(1,2), sigE=0.15,
                                betE=3, TotalEmerg=100, MaxLife=10, ShapeA=0.5, ShapeB=0.2) {

        if(length(PeakPos) < GenNumb) {stop("For multivoltine species, you need to provide a vector of distinct peak positions [PeakPos] to cover each emergence \n")}
        
        sdPeak <- rep(sdPeak,GenNumb)
        sigE <- rep(sigE,GenNumb)
        betE <- rep(betE,GenNumb)
        TotalEmerg <- rep(TotalEmerg,GenNumb)
        MaxLife <- rep(MaxLife,GenNumb)

        d_season <- data.table::copy(d_season)

        if(isTRUE(CompltSeason)){
            sim_season=d_season[COMPLT_SEASON==1,unique(M_SEASON)]
        }else{
            sim_season=d_season[,unique(M_SEASON)]
        }

        GenSim <- 1
            
        while(GenSim<=GenNumb) {
            for (i in sim_season){
                if(i==0){next}
                    t_s <- seq_along(1:d_season[M_SEASON==i,.N])
                    emerg_curve <- sim_emerg_curve(t_s,PeakPos=PeakPos[GenSim],sdPeak=sdPeak[GenSim],sigE=sigE[GenSim],betE=betE[GenSim])
                    emerg_nbr <- sim_emerg_nbr(emerg_curve,TotalEmerg=TotalEmerg[GenSim])
                    adult_nbr <- sim_adult_nbr(emerg_nbr,MaxLife=MaxLife[GenSim],ShapeA, ShapeB)
                    d_season[M_SEASON==i,ADLT_NBR:=as.integer(adult_nbr[1:d_season[M_SEASON==i,.N]])]
            }
            if(GenSim==1){
                cumul_count <- d_season[,ADLT_NBR]
            }else{
            cumul_count <- cumul_count + d_season[,ADLT_NBR]
            }
            GenSim <- GenSim+1
        }

        d_season[,ADLT_NBR:=cumul_count]
        d_season[,SITE_ID:=1]
        d_season[,SPECIES:='sim']

    return(d_season)
}

# sim_monitoring_visit() simulates monitoring visits by volunteers, based on monitoring frequency set by the protocol c('weekly','fortnightly','monthly') or c('none') for all days within season

sim_monitoring_visit <- function(d_season, MonitoringFreq=c('none')){

        d_season[,WEEK:=data.table::isoweek(DATE)]

        if(MonitoringFreq=='weekly'){
            is_even <- sample(c(1,2),1) == 2
            monitoring_day <- d_season[M_SEASON!=0L & (DAY_SINCE %% 2 == 0L)==is_even,sample(DAY_SINCE,1),by=.(M_YEAR,WEEK)][,V1]    
        }

        if(MonitoringFreq=='fortnightly'){
            is_even <- sample(c(1,2),1) == 2
            monitoring_day <- d_season[M_SEASON!=0L & (WEEK %% 2 == 0L)==is_even,sample(DAY_SINCE,1),by=.(M_YEAR,WEEK)][,V1]
        }

        if(MonitoringFreq=='monthly'){
            is_even <- sample(c(1,2),1) == 2
            monitoring_day <- d_season[M_SEASON!=0L & (WEEK %% 2 == 0L)==is_even,sample(DAY_SINCE,1),by=.(M_YEAR,MONTH)][,V1]
        }

        if(MonitoringFreq=='none'){
            is_even <- sample(c(1,2),1) == 2
            monitoring_day <- d_season[M_SEASON!=0L,DAY_SINCE]
        }

    return(monitoring_day)
}

### sim_sites_count function to simulate butterfly count across sites, assuming a common flight curve,
###                     but with potential shift in the position of the Peak

sim_butterfly_count <- function(d_season,NbrSite=10,FullSeason=TRUE,GenNumb=1,PeakPos=c(25,75),sdPeak=c(1,2),sigE=0.15,betE=3,TotalEmerg=100,MaxLife=10,MonitoringFreq=c('none'),PerctSampled=100,DetectProb=1){

        site_count_list <- vector("list",NbrSite)
        site_visit_list <- vector("list",NbrSite)

        for(s in 1:(NbrSite)){
            d_season_count <- sim_butterfly_count(d_season,FullSeason=FullSeason,GenNumb=GenNumb,PeakPos=PeakPos,sdPeak=sdPeak,sigE=sigE,betE=betE,TotalEmerg=TotalEmerg,MaxLife=MaxLife)
            m_day <- sim_monitoring_visit(d_season,MonitoringFreq=MonitoringFreq)
            site_count <- d_season_count[DAY_SINCE %in% sample(m_day,round((length(m_day)*PerctSampled)/100),replace=FALSE),COUNT:=rbinom(1,ADLT_NBR,DetectProb)][,SITE_ID:=SITE_ID+(s-1)]
            site_count_list[[s]] <- site_count
        }
   
        m_count <- data.table::rbindlist(site_count_list)

    return(m_count)
}


### build_map_obj builds map for specific region, using function from the sf() simple feature package

build_map_obj <- function(region=c('Africa','Antarctica','Americas','Asia','Europe','Oceania','Russia'),GadmWorld=GadmWorld,level=0) {

    iso_code <- region[nchar(region)==3]

    region_name <- region[nchar(region)!=3]

    test_region <- !(region_name %in% c('Africa','Antarctica','America','Asia','Europe','Oceania','Russia'))
    if(sum(test_region)>0) {cat(paste(dQuote(region_name[test_region]),'is not a recognized region. It must be one of these:','\n',
                                        'Africa,','Antarctica,','Americas,','Asia,','Europe,','Oceania or','Russia','\n'))}


    gadm_set <- GadmWorld[(unregion2 %in% region_name | iso3 %in% iso_code) & !is.na(unregion2),]

    for (i in seq_along(unlist(gadm_set[,iso3]))){

        cat(unlist(gadm_set[,name_english][i]),"\n")

        country_sf <- sf::st_as_sf(raster::getData(name = "GADM", country = gadm_set[,iso3][i], level = level))

            if (i == 1) {
                combined_sf <- country_sf
            }else{ 
                combined_sf <- rbind(combined_sf,country_sf)
            }

        }

return(combined_sf)

}

