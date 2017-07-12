##==========================================
## Name Convention
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
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

ts_monit_season = function(d_series,StartMonth=4,EndMonth=9,StartDay=1,EndDay=NULL){
        
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

            ## trim time-series around complete season
            a <- as.Date(paste(min(data.table::year(d_series$date)),StartMonth,StartDay,sep='-'))
            b <- as.Date(paste(max(data.table::year(d_series$date)),EndMonth,EndDay,sep='-'))
            d <- d[d$date >= a & d$date <= b]

            s_month <- c(StartMonth:12,1:EndMonth)
            s_md_out <- c(paste(StartMonth,c(0:(StartDay-1))),paste(EndMonth,c((EndDay+1):32)))
            y_series <- data.table::data.table(season_year=as.factor(ifelse(data.table::month(d_series$date)>=StartMonth,data.table::year(d_series$date),(data.table::year(d_series$date)-1))),
                                                season=ifelse(((data.table::month(d_series$date)%in%(s_month))
                                                        & (!(paste(data.table::month(d_series$date),data.table::mday(d_series$date))%in%s_md_out))),
                                                        1,0))                               
        }

        d_series <- d_series[,c("season_year","season") := y_series[,.(season_year,(as.numeric(season_year)*season))]]

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

ts_count_site_visit = function(sp=1,d_series,c_t_day) {

    data.table::setkey(c_t_day,SITENO,VISIT_DATE)
    data.table::setkey(d_series,site,date)
    c_site_series <- d_series[c_t_day[SPECIES %in% sp,],count:=as.integer(total_day_count)]
    c_site_series[,species:=sp]

    return(c_site_series)

    }


t.d.count[nbr_visit==3,]
c.t_day[VISITDATE=="2014-05-16" & SITENO==1842]

## load monitoring visit and fixing dates

m.visit <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Visits Table.txt",header=TRUE)
m.visit[,VISIT_DATE:=data.table::as.IDate(as.Date(m.visit$VISIT_DATE,format="%d-%b-%y"))]

m.count <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK CountsTable.txt",header=TRUE)
m.count[,VISIT_DATE:=data.table::as.IDate(as.Date(m.count$VISITDATE,format="%d/%m/%Y"))][,VISITDATE:=NULL]

section.geo <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Sections Table.txt",header=TRUE)
section.geo.sp <-section.geo
sp::coordinates(section.geo.sp) <- ~EAST + NORTH



## test

m.count[,max(data.table::year(VISIT_DATE))]

d <- t.table(2015)
d.season <- m.season(d,4,9)
d.series <- d.site.series(m.visit, d.season)
c.t_day <- s.count.t_day(m.count,m.visit)
c.site_series <- s.count.site.series(sp=2,d.series,c.t_day)

c.site_series[year==1997 & count==0,site]

## object size check
sort( sapply(c('d.series','d.series2'),function(x){object.size(get(x))}))


### SF build map

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

### SIMULATE DATA

### e_curve() estimates an emergence curve shape following a logistic distribution along a time series [t_series], with peak position [peak_pos] relative along a  , 
### vector using the percentile and a standard deviation around the peak [sd_peak] in days, with two shape parameters
### sigma [sig] for left or right skewness (left when > 1, right when < 1, logistic when  = 1) and bet [bet] for a scale parameter. 

### return a vector of relative emergence along a vector of length t.

### Calabrese, J.M. (2012) How emergence and death assumptions affect count-based estimates of butterfly abundance and lifespan. Population Ecology, 54, 431–442.


e_curve <- function (t_series, peak_pos=50, sd_peak=1, sig=0.15, bet=3) {
               
            # EMERGENCE 1
                
            # peak + sample 1 from normal distribution (0,sd)
            u1 <- round(peak_pos * length(t_series)/100) + rnorm(1, 0, sd_peak)

            # Logistic with 
            fE1 <- (sig * (exp((t_series - u1)/bet)))/(bet * (1 + exp((t_series - u1)/bet)))^(sig + 1)
                
            # Standardize to AUC 1
            sdfe <- fE1/sum(fE1)  

            return(sdfe)
            
            }

### c_curve simulate emergence of n adults [tot_emerg] according a emergence curve [sdfe] using a Poisson process

c_curve <- function(sdfe, tot_emerg=100) {

            N <- unlist(lapply(tot_emerg*sdfe,function(x) {rpois(1,x)}))

            return(N)

            }

### a_curve simulate the number of adults in a population, according to individual maximum life span [max_life] and daily
### mortality risk based on a beta distribution with shape_a and shape_b parameters

a_curve <- function(N, max_life=15, shape_a=0.5, shape_b=0.2){

           c_mat <- matrix(0,nrow=length(N),ncol=length(N)+max_life+1)

           ## daily hazard and decay (beta distribution)

           M_haz <- c(0,diff(pbeta(seq(0,1,(1/max_life)),shape_a,shape_b)),1)
        
           for (i in seq_along(N)){

                s <- N[i]
                y=2
                l=s

                while(s > 0 & y <= max_life+2){
                    s <- s-rbinom(1,s,M_haz[y])
                    y <- y+1
                    l <- c(l,s)
                    }

                c_mat[i,i:(i+length(l)-1)] <- l

                }

            return(colSums(c_mat))

            }

## test

par(mfrow=c(2,2))

for (t in 1:4){
    e <- e_curve(1:100,50,sig=0.5,bet=3)
plot(1:100,e)
    c <- c_curve(e,100)
plot(1:100,c,,main=paste(sum(c),'realized emergences'))
points(1:100,100*e,col='red',pch=19)
    a1 <- a_curve(c,max_life=10)
plot(seq_along(a1),a1,main=paste(sum(a1),'cumulative adult counts (b1)'))

    e <- e_curve(1:100,60,sig=0.5,bet=2)
    c <- c_curve(e,100)
    a2 <- a_curve(c,max_life=10)
plot(seq_along(1:max(length(a1),length(a2))),a1+a2,main=paste(sum(c(a1,a2)),'cumulative adult counts b1+2'))

    }



###########################################

## test package sf

###build initial data set
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