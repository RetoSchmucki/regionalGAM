## sandbox for testing new development for regionalGAM
## Date: 21.07.2017

## NEED TO SOURCE new_functions.r

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

d_series <- ts_dwmy_table(2008,2015)
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

monitoring_day <- d_season_count[full_season!=0,sample(day_since,1),by=.(season_year,week)]
d_season_count[day_since %in% monitoring_day$V1,monitored_count:=count]

par(mfrow=c(1,2))
plot(d_season_count[full_season!=0,sum(count),by=season_year])
plot(d_season_count[full_season!=0,sum(monitored_count,na.rm=TRUE),by=season_year])

plot(d_season_count$day_since,d_season_count$count,col='magenta',pch=19,type='l')
points(d_season_count$day_since,d_season_count$monitored_count,col='blue',pch=19)


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


sim_monitoring_visit(d_season,'monthly')

monitoring_day <- d_season[complete_season!=0,sample(day_since,1),by=.(season_year,week)]

d_season[day_since %in% monitoring_day$V1,monitored_count:=count]

###########################################

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

GadmWorld <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/reto_workfiles/ebms_database/data/gadm_world.csv",header=TRUE)

## build Europe

region=c('Africa','GBR','CHE')

my_map <- build_map_obj(region,GadmWorld,level=0)

plot(my_map$geometry,col=seq_along(my_map$ISO))