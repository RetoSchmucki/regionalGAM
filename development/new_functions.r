##==========================================
## Name Convention
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##
##      Date: 20.07.2017
##
##==========================================

###  ts_date_seq  function to build a full time series sequence of date since initial year to an ending year
###         to subset a time-series template for specified monitoring season, independent of 
###         the year

ts_date_seq = function(InitYear=1970,LastYear=format(Sys.Date(),"%Y")) {

    init_date <- as.Date(paste((InitYear-1), "01-01", sep = "-"))
    last_date <- as.Date(paste((as.numeric(LastYear)+1), "12-31", sep = "-"))

    date_serie <- as.POSIXct(seq(from=init_date, to= last_date, by = "day"), format = "%Y-%m-%d")
    date_serie <- date_serie[!format(date_serie,'%Y') %in% c((InitYear-1),as.numeric(LastYear)+1)]

    return(date_serie)

    }

###  ts_dwmy_table  function to build a full time series sequence of days, iso weeks, week-days (1:monday, 7:sunday), 
###         months and years to subset a time-series template for specified monitoring season,
###         independent of the year

ts_dwmy_table = function(InitYear=1970,LastYear=format(Sys.Date(),"%Y"),WeekDay1='monday') {
    
        date_seq <- ts_date_seq(InitYear,LastYear)

        if(WeekDay1=='monday'){
            w <- c(7,1:6)
        }else{
            w <- c(1:7)
        }

        dt_iso_dwmy <- data.table::data.table(
                                date=data.table::as.IDate(date_seq),
                                day_since=seq_along(date_seq),
                                year=data.table::year(date_seq),
                                month=data.table::month(date_seq),
                                week=data.table::isoweek(date_seq),
                                week_day=w[data.table::wday(date_seq)],
                                month_day=data.table::mday(date_seq))

        return(dt_iso_dwmy)

    }

###  ts_monit_season  function to build a full time_series sequence of monitoring season, with specific starting end ending
###             months and days, the season start in year y and end in year y+1.

ts_monit_season = function(d_series,StartMonth=4,EndMonth=9,StartDay=1,EndDay=NULL,FullSeason=TRUE){
        
        d_series <- data.table::copy(d_series)

        ## NO leap year
        if(is.null(EndDay)) {
            EndDay <- max(data.table::mday(seq(from=as.Date(paste0('2017-',EndMonth,'-01')),by='day',length=31)))
        } 

        if (StartMonth < EndMonth){
            s_month <- c(StartMonth:EndMonth)
            s_md_out <- c(paste(StartMonth,c(0:(StartDay-1))),paste(EndMonth,c((EndDay+1):32)))
            y_series <- data.table::data.table(season_year=as.factor(data.table::year(d_series$date)),
                                               season=ifelse(((data.table::month(d_series$date)%in%(s_month))
                                                        & (!(paste(data.table::month(d_series$date),data.table::mday(d_series$date)) %in% s_md_out))),
                                                        1,0))                               
        }

        if (StartMonth > EndMonth){

            s_month <- c(StartMonth:12,1:EndMonth)
            s_md_out <- c(paste(StartMonth,c(0:(StartDay-1))),paste(EndMonth,c((EndDay+1):32)))
            y_series <- data.table::data.table(season_year=as.factor(ifelse(data.table::month(d_series$date)>=StartMonth,data.table::year(d_series$date),(data.table::year(d_series$date)-1))),
                                                season=ifelse(((data.table::month(d_series$date)%in%(s_month))
                                                        & (!(paste(data.table::month(d_series$date),data.table::mday(d_series$date))%in%s_md_out))),
                                                        1,0))                               
        }

        d_series <- d_series[,c("season_year","season") := y_series[,.(season_year,(as.numeric(season_year)*season))]]

        ## identify (trim) the full seasons
        if(isTRUE(FullSeason)){
            
            d_series[,start_end:=d_series[,ifelse(month==StartMonth & month_day==StartDay,1L,0L)]
                                            + d_series[,ifelse(month==EndMonth & month_day==EndDay,1L,0L)]]
            
            d_series[,full_season:=ifelse(season!= 0 & 
                                         season_year %in% (d_series[,sum(start_end),by=season_year][V1==2,season_year]),
                                         1,0)*season]

            d_series[,start_end:=NULL]
        }   

        return(d_series)
    }

### ts_site_visit function to initialize a time series with all visit and site with "zeros" while leaving all non visited day 
###                 with an <NA>, this can then be used to add the observed count for specific species
###                 only have time series for years when a site has been monitored.

ts_site_visit = function(m_visit, d_season) {
       
        data.table::setkey(d_season,year)
        data.table::setkey(m_visit,YEAR,SITENO)

        r_year <- d_season[,range(data.table::year(date))]

        site_l <- m_visit[data.table::year(VISIT_DATE)>=min(r_year) & 
                            data.table::year(VISIT_DATE)<=max(r_year),
                            .(site=.SD[,unique(SITENO)]),by=YEAR]

        data.table::setkey(site_l,YEAR,site)
        
        d_site <- merge(d_season,site_l,by.x="year",by.y="YEAR",allow.cartesian=TRUE)

        data.table::setkey(d_site,site,date)
        data.table::setkey(m_visit,SITENO,VISIT_DATE)

        d_site <- d_site[m_visit,count:=0L]

        return(d_site)

    }

### b_count_perday function to standardize the butterfly species count for one visit per day. Two options,
###                 use the mean rounded to the higher integer or delete site where count could
###                 not be unambiguously attributed to a single visit 
###                 (method: "average" or "delete")

b_count_perday = function(m_count,m_visit,UniMethod="average") {

        data.table::setkey(m_count,SITENO,SPECIES,VISIT_DATE)
        data.table::setkey(m_visit,SITENO,VISIT_DATE)
        
        t_m_count <- m_count[,.(total_count=sum(COUNT)),by=.(SITENO,SPECIES,VISIT_DATE)]
        t_d_count <- m_visit[,.(nbr_visit=.N),by=.(SITENO,VISIT_DATE)]

        data.table::setkey(t_m_count,SITENO,VISIT_DATE)
        data.table::setkey(t_d_count,SITENO,VISIT_DATE)

        if(UniMethod=="average") {t_m_count[t_d_count,total_day_count:=ceiling(total_count/nbr_visit)]}
        if(UniMethod=="delete") {t_m_count[nbr_visit==1,t_d_count,total_day_count:=total_count]}

        return(t_m_count)

    }

### ts_count_site_visit function to generate a full time series of observed count, including zeros and missing
###                     observation for one species for each day since a starting and ending years

ts_count_site_visit = function(d_series,c_t_day,sp=1) {

    d_series <- data.table::copy(d_series)

    data.table::setkey(c_t_day,SITENO,VISIT_DATE)
    data.table::setkey(d_series,site,date)
    c_site_series <- d_series[c_t_day[SPECIES %in% sp,],count:=as.integer(total_day_count)]
    c_site_series[,species:=sp]

    return(c_site_series)

    }


### SIMULATE DATA

### sim_emerg_curve() estimates an emergence curve shape following a logistic distribution along a time series [t_series], with peak position [peak_pos] relative along a  , 
### vector using the percentile and a standard deviation around the peak [sd_peak] in days, with two shape parameters
### sigma [sig] for left or right skewness (left when > 1, right when < 1, logistic when  = 1) and bet [bet] for a scale parameter. 

### return a vector of relative emergence along a vector of length t.

### Calabrese, J.M. (2012) How emergence and death assumptions affect count-based estimates of butterfly abundance and lifespan. Population Ecology, 54, 431–442.


sim_emerg_curve <- function (t_series, PeakPos=50, sdPeak=1, sig=0.15, bet=3) {
               
            # EMERGENCE 1
                
            # peak + sample 1 from normal distribution (0,sd)
            u1 <- round(PeakPos * length(t_series)/100) + rnorm(1, 0, sdPeak)

            # Logistic with 
            fE1 <- (sig * (exp((t_series - u1)/bet)))/(bet * (1 + exp((t_series - u1)/bet)))^(sig + 1)
                
            # Standardize to AUC 1
            sdfe <- fE1/sum(fE1)  

            return(sdfe)
            
            }

### sim_emerg_count simulates emergence of n adults [TotalEmerg] according an emergence curve [sdfe] using a Poisson process

sim_emerg_count <- function(sdfe, TotalEmerg=100) {

            n_emerg <- unlist(lapply(TotalEmerg*sdfe,function(x) {rpois(1,x)}))

            return(n_emerg)

            }

### sim_adult_count simulates the number of adults in a population, according to individual maximum life span [max_life] and daily
### mortality risk based on a beta distribution with ShapeA and ShapeB parameters

sim_adult_count <- function(n_emerg, MaxLife=15, ShapeA=0.5, ShapeB=0.2){

           c_mat <- matrix(0,nrow=length(n_emerg),ncol=length(n_emerg)+MaxLife+1)

           ## daily hazard and decay (beta distribution)

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


### sim_butterfly_count() simulates count data along monitoring seasons for a univoltine or multivoltine species, from month 4 to 9 (April to September)

sim_butterfly_count <- function(d_season, FullSeason=TRUE, GenNumb=1, PeakPos=c(25,75), sdPeak=c(1,2), sig=0.15,
                                bet=3, TotalEmerg=100, MaxLife=10) {

        if(length(PeakPos) < GenNumb) {stop("For multivoltine species, you need to provide a vector of distinct peak positions [PeakPos] to cover each emergence \n")}
        
        sdPeak <- rep(sdPeak,GenNumb)
        sig <- rep(sig,GenNumb)
        bet <- rep(bet,GenNumb)
        TotalEmerg <- rep(TotalEmerg,GenNumb)
        MaxLife <- rep(MaxLife,GenNumb)

        d_season <- data.table::copy(d_season)

        if(isTRUE(FullSeason)){
            sim_season=unique(d_season$full_season)
            }else{
            sim_season=unique(d_season$season)
        }

        GenSim <- 1

        for (i in sim_season){

            if(i==0){next}

                t_s <- seq_along(1:d_season[season==i,.N])
                emerg_curve <- sim_emerg_curve(t_s,PeakPos=PeakPos[GenSim],sdPeak=sdPeak[GenSim],sig=sig[GenSim],bet=bet[GenSim])
                emerg_count <- sim_emerg_count(emerg_curve,TotalEmerg=TotalEmerg[GenSim])
                adult_count <- sim_adult_count(emerg_count,MaxLife=MaxLife[GenSim])
                d_season[season==i,count:=adult_count[1:d_season[season==i,.N]]]
        }

        cumul_count <- d_season[,count]
        GenSim <- GenSim+1
            
        while(GenNumb>=GenSim) {

            for (i in sim_season){

                if(i==0){next}

                    t_s <- seq_along(1:d_season[season==i,.N])
                    emerg_curve <- sim_emerg_curve(t_s,PeakPos=PeakPos[GenSim],sdPeak=sdPeak[GenSim],sig=sig[GenSim],bet=bet[GenSim])
                    emerg_count <- sim_emerg_count(emerg_curve,TotalEmerg=TotalEmerg[GenSim])
                    adult_count <- sim_adult_count(emerg_count,MaxLife=MaxLife[GenSim])
                    d_season[season==i,count:=adult_count[1:d_season[season==i,.N]]]
            }
        
        cumul_count <- cumul_count + d_season[,count]
        GenSim <- GenSim+1

        }

    d_season[,count:=cumul_count]

    return(d_season)

}

# sim_monitoring_visit() simulates monitoring visits by volunteers, based on monitoring frequency set by the protocol c('weekly','fortnightly','monthly')

sim_monitoring_visit <- function(d_season, MonitoringFreq=c('weekly')){

    if(MonitoringFreq=='weekly'){

        is_even <- (d_season[season!=0,day_since][1]) %% 2 == 0
        monitoring_day <- d_season[season!=0 & (day_since %% 2 == 0)==is_even,sample(day_since,1),by=.(season_year,week)][,V1]
    
    }

    if(MonitoringFreq=='fortnightly'){
        
        is_even <- (d_season[season!=0,week][1]) %% 2 == 0
        monitoring_day <- d_season[season!=0 & (week %% 2 == 0)==is_even,sample(day_since,1),by=.(season_year,week)][,V1]
    
    }

    if(MonitoringFreq=='monthly'){

        monitoring_day <- d_season[season!=0,sample(day_since,1),by=.(season_year,month)][,V1]
    
    }

    return(monitoring_day)

}


sim_monitoring_visit(d_season,'monthly')

monitoring_day <- d_season[complete_season!=0,sample(day_since,1),by=.(season_year,week)]

d_season[day_since %in% monitoring_day$V1,monitored_count:=count]

## test SIMULATION functions

par(mfrow=c(2,2))

for (t in 1:4){
    e <- sim_emerg_curve(1:100,50,sig=1,bet=10)
plot(1:100,e)
    c <- sim_emerg_count(e,100)
plot(1:100,c,,main=paste(sum(c),'realized emergences'))
points(1:100,100*e,col='red',pch=19)
    a1 <- sim_adult_count(c,MaxLife=10)
plot(seq_along(a1),a1,main=paste(sum(a1),'cumulative adult counts (b1)'))

    e <- sim_emerg_curve(1:100,70,sig=0.5,bet=2)
    c <- sim_emerg_count(e,100)
    a2 <- sim_adult_count(c,MaxLife=10)
plot(seq_along(1:max(length(a1),length(a2))),a1+a2,main=paste(sum(c(a1,a2)),'cumulative adult counts b1+2'))

    }

# ## load monitoring visit and fixing dates

m_visit <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Visits Table.txt",header=TRUE)
m_visit[,VISIT_DATE:=data.table::as.IDate(as.Date(m_visit$VISIT_DATE,format="%d-%b-%y"))]

m_count <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK CountsTable.txt",header=TRUE)
m_count[,VISIT_DATE:=data.table::as.IDate(as.Date(m_count$VISITDATE,format="%d/%m/%Y"))][,VISITDATE:=NULL]

section_geo <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Sections Table.txt",header=TRUE)
section_geo_sp <- section_geo[!is.na(section_geo$EAST)]
sp::coordinates(section_geo_sp) <- ~ EAST + NORTH


## test full time series allocation with real data 

m_count[,max(data.table::year(VISIT_DATE))]

d <- ts_dwmy_table(2015)
d_season <- ts_monit_season(d,4,9)
d_series <- ts_site_visit(m_visit, d_season)
c_t_day <- b_count_perday(m_count,m_visit)
c_site_series <- ts_count_site_visit(sp=2,d_series,c_t_day)

c_site_series[year==1997 & count==0,site]

## object size check
sort(sapply(c('d_series'),function(x){object.size(get(x))}))




d_series <- ts_dwmy_table(2017)
d_season <- ts_monit_season(d_series,4,9)

par(mfrow=c(2,2))
for(sim in 1:4){
 d_season_count <- sim_butterfly_count(d_season,GenNumb=2,PeakPos=c(25,65),TotalEmerg=c(100,200))
 plot(d_season_count$day_since,d_season_count$count,type='l',col='magenta')

 m_day <- sim_monitoring_visit(d_season,'weekly')
 points(m_day,d_season_count$count[m_day],pch=19,col='red')

 m_day <- sim_monitoring_visit(d_season,'fortnightly')
 points(m_day,d_season_count$count[m_day],pch=19,col='blue')
}

cbind(d_season_count$day_since[m_day],m_day)

monitoring_day <- d_season[complete_season!=0,sample(day_since,1),by=.(season_year,week)]
d_season[day_since %in% monitoring_day$V1,monitored_count:=count]

par(mfrow=c(1,2))
plot(d_season[complete_season!=0,sum(count),by=season_year])
plot(d_season[complete_season!=0,sum(monitored_count,na.rm=TRUE),by=season_year])

plot(d_season$day_since,d_season$count,col='magenta',pch=19,type='l')
points(d_season$day_since,d_season$monitored_count,col='blue',pch=19)

monitoring_day <- d_season[season!=0,sample(day_since,1),by=.(season_year,week)]
points(monitoring_day$V1,rep(0,length(monitoring_day$V1)))



DT[,.SD[sample(.N,3)],by = a]

###########################################


### SF build map

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


## test package sf

###build initial data set that should be included in the package's data set
system.time(gadm28 <- sf::st_read("C:/Users/RETOSCHM/Downloads/gadm28_levels.shp/gadm28_adm0.shp"))

gadm_World <- data.table::data.table( name_english    = as.character(gadm28$NAME_ENGLI),
                        iso3        = as.character(gadm28$ISO),
                        iso2        = as.character(gadm28$ISO2),
                        unregion2   = as.character(gadm28$UNREGION2))

unique(gadm_World$unregion2)

### fix some region names and use Russia a one region

gadm_World[iso3=='RUS',c("unregion2"):='Russia']
gadm_World[unregion2=='Antartica',c("unregion2"):='Antarctica']

write.csv(gadm_World,"W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/reto_workfiles/ebms_database/data/gadm_world.csv",row.names=FALSE)

## build Europe

region=c('Africa','GBR','CHE')

build_map_obj <- function(region=c('Africa','Antarctica','Americas','Asia','Europe','Oceania','Russia'),gadm_World=gadm_World,level=0) {

    iso_code <- region[nchar(region)==3]

    region_name <- region[nchar(region)!=3]

    test_region <- !(region_name %in% c('Africa','Antarctica','America','Asia','Europe','Oceania','Russia'))
    if(sum(test_region)>0) {cat(paste(dQuote(region_name[test_region]),'is not a recognized region. It must be one of these:','\n',
                                        'Africa,','Antarctica,','Americas,','Asia,','Europe,','Oceania or','Russia','\n'))}


    gadm.set <- gadm_World[(unregion2 %in% region_name | iso3 %in% iso_code) & !is.na(unregion2),]

    for (i in seq_along(unlist(gadm.set[,iso3]))){

        cat(unlist(gadm.set[,name_english][i]),"\n")

        country_sf <- sf::st_as_sf(raster::getData(name = "GADM", country = gadm.set[,iso3][i], level = level))

            if (i == 1) {
                combined_sf <- country_sf
            }else{ 
                combined_sf <- rbind(combined_sf,country_sf)
            }

        }

return(combined_sf)

}