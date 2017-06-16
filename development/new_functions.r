    
###  t.seq  function to build a full sequence of days, weeks, week-days, months and years
###         to subset a time-series template for specified monitoring season, independent of 
###         the year

t.seq = function(init.year=1970,last.year=format(Sys.Date(),"%Y")) {

    init.date <- as.Date(paste((init.year-1), "01-01", sep = "-"))
    last.date <- as.Date(paste((as.numeric(last.year)+1), "12-31", sep = "-"))

    date.serie <- as.POSIXct(seq(from=init.date, to= last.date, by = "day"), format = "%Y-%m-%d")
    date.serie <- date.serie[!format(date.serie,'%Y') %in% c((init.year-1),as.numeric(last.year)+1)]

    return(date.serie)

    }

t.table = function(init.year=1970,last.year=format(Sys.Date(),"%Y"),week_day1='monday') {
    
        d.s <- t.seq(init.year,last.year)

        if(week_day1=='monday'){
            w <- c(7,1:6)
        }else{
            w <- c(1:7)
        }

        date.serie <- data.table::data.table(date=data.table::as.IDate(d.s),
                                 day.since=seq_along(d.s),
                                 month=data.table::month(d.s),
                                 week=data.table::isoweek(d.s),
                                 week_day=w[data.table::wday(d.s)],
                                 month_day=data.table::mday(d.s))

        return(date.serie)

    }   


year_day_func = function(sp_data) {






a <- t.table(1999,2000)



a <- data.table::as.IDate(init.date)
b <- data.table::as.IDate(last.date) 

data.table::isoweek(a)[5]


c(7,1,2,3,4,5,6)[data.table::wday(a)[5]]
data.table::month(a)[32]



a[1]

install.packages("ISOweek")

ISOweek::date2ISOweek(a)[5]

lubridate::wday(a)[1]






b-a







    dayno <- as.numeric(julian(date.serie, origin = as.Date(origin.d)) + 1)
    month <- as.numeric(strftime(date.serie, format = "%m"))
    week <- as.numeric(strftime(date.serie, format = "%W"))
    week_day <- as.numeric(strftime(date.serie, format = "%u"))
    day <- as.numeric(strftime(date.serie, format = "%d"))

