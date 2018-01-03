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
##      rbms_simulate: function for simulating butterfly count
##
##==========================================

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