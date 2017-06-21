    
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
###                 - remove all monitoring seasons before establishment for a specific site.


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

        d.site <- d.site[m.visit,count:=0L]

        return(d.site)

    } 

### s.count.t_day function to standardize the species count for one visit per day. Two options,
###                 use the mean rounded to the higher integer or delete site where count could
###                 not be unambigously be attibutre to a single visit 
###                 (method: "average" or "delete")

s.count.t_day = function(m.count,m.visit,method="average") {

        data.table::setkey(m.count,SITENO,SPECIES,VISITDATE)
        data.table::setkey(m.visit,SITENO,VISIT_DATE)
        
        t.m.count <- m.count[,.(total_count=sum(COUNT)),by=.(SITENO,SPECIES,VISITDATE)]
        t.d.count <- m.visit[,.(nbr_visit=.N),by=.(SITENO,VISIT_DATE)]

        data.table::setkey(t.m.count,SITENO,VISITDATE)
        data.table::setkey(t.d.count,SITENO,VISIT_DATE)

        if(method=="average") {t.m.count[t.d.count,total_day_count:=ceiling(total_count/nbr_visit)]}
        if(method=="delete") {t.m.count[nbr_visit==1,t.d.count,total_day_count:=total_count]}

        return(t.m.count)

    } 

### s.count.site.series function to generate a full series of observed count, including zeros and missing
###                     observation for one species for each day since a starting and ending years

s.count.site.series = function(sp=1,d.series,c.t_day) {

    data.table::setkey(c.t_day,SITENO,VISITDATE)
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

## test

d <- t.table(1970)
d.season <- m.season(d,4,9)
d.series <- d.site.series(m.visit, d.season)
c.t_day <- s.count.t_day(m.count,m.visit)
c.site_series <- s.count.site.series(sp=54,d.series,c.t_day)



c.site_series[date=="2002-08-11" & count>=0,]


## object size check
sort( sapply(ls(),function(x){object.size(get(x))}))


year_day_func = function(sp_data) {


