    
###  t.seq  function to build a full sequence of date since initial year to an ending year
###         to subset a time-series template for specified monitoring season, independent of 
###         the year

t.seq = function(init.year=1970,last.year=format(Sys.Date(),"%Y")) {

    init.date <- as.Date(paste((init.year-1), "01-01", sep = "-"))
    last.date <- as.Date(paste((as.numeric(last.year)+1), "12-31", sep = "-"))

    date.serie <- as.POSIXct(seq(from=init.date, to= last.date, by = "day"), format = "%Y-%m-%d")
    date.serie <- date.serie[!format(date.serie,'%Y') %in% c((init.year-1),as.numeric(last.year)+1)]

    return(date.serie)

    }

###  t.table  function to build a full sequence of days, iso weeks, week-days (1:monday, 7:sunday), 
###         months and years to subset a time-series template for specified monitoring season,
###         independent of the year

t.table = function(init.year=1970,last.year=format(Sys.Date(),"%Y"),week_day1='monday') {
    
        d.s <- t.seq(init.year,last.year)

        if(week_day1=='monday'){
            w <- c(7,1:6)
        }else{
            w <- c(1:7)
        }

        date.series <- data.table::data.table(date=data.table::as.IDate(d.s),
                                 day.since=seq_along(d.s),
                                 year=data.table::year(d.s),
                                 month=data.table::month(d.s),
                                 week=data.table::isoweek(d.s),
                                 week_day=w[data.table::wday(d.s)],
                                 month_day=data.table::mday(d.s))

        return(date.series)

    }

###  m.season  function to build a full sequence of monitoring season, with specific starting end ending
###             months and days, the season start in year y and end in year y+1.

m.season = function(d.series,start.month=4,end.month=9,start.day=1,end.day=NULL){
        
        ## NO leap year
        if(is.null(end.day)) {
            end.day <- max(data.table::mday(seq(from=as.Date(paste0('2017-',end.month,'-01')),by='day',length=31)))
        } 

        if (start.month < end.month){
            s.month <- c(start.month:end.month)
            s.md.out <- c(paste(start.month,c(0:(start.day-1))),paste(end.month,c((end.day+1):32)))
            y.series <- data.table::data.table(season.year=as.factor(data.table::year(d.series$date)),
                                                season=ifelse(((data.table::month(d.series$date)%in%(s.month))
                                                        & (!(paste(data.table::month(d.series$date),data.table::mday(d.series$date))%in%s.md.out))),
                                                        1,0))                               
        }

        if (start.month > end.month){

            ## trim time-series around complete season
            a <- as.Date(paste(min(data.table::year(d.series$date)),start.month,start.day,sep='-'))
            b <- as.Date(paste(max(data.table::year(d.series$date)),end.month,end.day,sep='-'))
            d <- d[d$date >= a & d$date <= b]

            s.month <- c(start.month:12,1:end.month)
            s.md.out <- c(paste(start.month,c(0:(start.day-1))),paste(end.month,c((end.day+1):32)))
            y.series <- data.table::data.table(season.year=as.factor(ifelse(data.table::month(d.series$date)>=start.month,data.table::year(d.series$date),(data.table::year(d.series$date)-1))),
                                                season=ifelse(((data.table::month(d.series$date)%in%(s.month))
                                                        & (!(paste(data.table::month(d.series$date),data.table::mday(d.series$date))%in%s.md.out))),
                                                        1,0))                               
        }

        d.series <- d.series[,c("season.year","season") := y.series[,.(season.year,(as.numeric(season.year)*season))]]

        return(d.series)
    }

### d.site.series function to initialize all visit and site with "zeros" while leaving all non visited day 
###                 with an <NA>, this can then be used to add the observed count for specific species
###                 only have time series for years when a site has been monitored.

d.site.series = function(m.visit, d.season) {
       
        data.table::setkey(d.season,year)
        data.table::setkey(m.visit,YEAR,SITENO)

        r.year <- d.season[,range(data.table::year(date))]

        site.l <- m.visit[data.table::year(VISIT_DATE)>=min(r.year) & 
                            data.table::year(VISIT_DATE)<=max(r.year),
                                .(site=.SD[,unique(SITENO)]),by=YEAR]

        data.table::setkey(site.l,YEAR,site)
        
        d.site <- merge(d.season,site.l,by.x="year",by.y="YEAR",allow.cartesian=TRUE)

        data.table::setkey(d.site,site,date)
        data.table::setkey(m.visit,SITENO,VISIT_DATE)

        d.site <- d.site[m.visit,count:=0L]

        return(d.site)

    }

### s.count.t_day function to standardize the species count for one visit per day. Two options,
###                 use the mean rounded to the higher integer or delete site where count could
###                 not be unambiguously attributed to a single visit 
###                 (method: "average" or "delete")

s.count.t_day = function(m.count,m.visit,method="average") {

        data.table::setkey(m.count,SITENO,SPECIES,VISIT_DATE)
        data.table::setkey(m.visit,SITENO,VISIT_DATE)
        
        t.m.count <- m.count[,.(total_count=sum(COUNT)),by=.(SITENO,SPECIES,VISIT_DATE)]
        t.d.count <- m.visit[,.(nbr_visit=.N),by=.(SITENO,VISIT_DATE)]

        data.table::setkey(t.m.count,SITENO,VISIT_DATE)
        data.table::setkey(t.d.count,SITENO,VISIT_DATE)

        if(method=="average") {t.m.count[t.d.count,total_day_count:=ceiling(total_count/nbr_visit)]}
        if(method=="delete") {t.m.count[nbr_visit==1,t.d.count,total_day_count:=total_count]}

        return(t.m.count)

    }

### s.count.site.series function to generate a full series of observed count, including zeros and missing
###                     observation for one species for each day since a starting and ending years

s.count.site.series = function(sp=1,d.series,c.t_day) {

    data.table::setkey(c.t_day,SITENO,VISIT_DATE)
    data.table::setkey(d.series,site,date)
    c.site_series <- d.series[c.t_day[SPECIES %in% sp,],count:=as.integer(total_day_count)]
    c.site_series[,species:=sp]

    return(c.site_series)

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

           ## decay

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
    e <- e_curve(1:100)
    plot(1:100,e)
    c <- c_curve(e,100)
    plot(1:100,c)
    a <- a_curve(c)
    plot(seq_along(a),a,main=paste('time',t))
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