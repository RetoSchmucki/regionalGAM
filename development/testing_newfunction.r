## sandbox for testing new development for regionalGAM
## Date: 21.07.2017

## NEED TO SOURCE new_functions.r
r --vanilla


# ## load monitoring visit and fixing dates

m_visit <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Visits Table.txt",header=TRUE)
data.table::setnames(m_visit,"SITENO","SITE_ID"); data.table::setnames(m_visit,"VISIT_DATE","DATE")
m_visit[,DATE:=data.table::as.IDate(as.Date(m_visit$DATE,format="%d-%b-%y"))]


m_count <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK CountsTable.txt",header=TRUE)
data.table::setnames(m_count,"SITENO","SITE_ID"); data.table::setnames(m_count,"VISITDATE","DATE")
m_count[,DATE:=data.table::as.IDate(as.Date(m_count$DATE,format="%d/%m/%Y"))]
m_count[,COUNT:=sum(COUNT),by=.(SITE_ID,SPECIES,DATE,DAY,MONTH,YEAR)][,SECTION:=NULL]

section_geo <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Sections Table.txt",header=TRUE)
section_geo_sp <- section_geo[!is.na(section_geo$EAST)]
sp::coordinates(section_geo_sp) <- ~ EAST + NORTH


### Test New Workflow with monitoring_year as main year indicator

## we need two tables; 1) visits and 2) counts
## here we will use count per transect where transect are used as SITE_ID
## =====
## 1. load and organize data sets
m_visit <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK Visits Table.txt",header=TRUE)
data.table::setnames(m_visit,"SITENO","SITE_ID"); data.table::setnames(m_visit,"VISIT_DATE","DATE")
m_visit[,DATE:=data.table::as.IDate(as.Date(m_visit$DATE,format="%d-%b-%y"))]

## standardize the number of visit per date for each site to one
nbr_visitperday <- m_visit[,.N,by=.(SITE_ID,DATE)]
data.table::setkey(nbr_visitperday,SITE_ID,DATE)
data.table::setkey(m_visit,SITE_ID,DATE)
m_visit <- unique(m_visit[,.(SITE_ID,DATE)])

m_count <- data.table::fread("W:/PYWELL_SHARED/Pywell Projects/BRC/BMS/eBMS/DATASETS/UK CountsTable.txt",header=TRUE)
data.table::setnames(m_count,"SITENO","SITE_ID"); data.table::setnames(m_count,"VISITDATE","DATE")
m_count[,DATE:=data.table::as.IDate(as.Date(m_count$DATE,format="%d/%m/%Y"))][,SECTION:=NULL]
m_count[,COUNT:=sum(COUNT),by=.(SITE_ID,SPECIES,DATE,DAY,MONTH,YEAR)]
data.table::setkey(m_count)
m_count <- unique(m_count)
m_count <- m_count[!is.na(SITE_ID) & !is.na(SPECIES) & !is.na(DATE) & !is.na(COUNT),]

## standardize the count for one visit per date and site, using the average if multiple visits 
data.table::setkey(m_count,SITE_ID,DATE)
m_count <- merge(m_count,nbr_visitperday,all.x=FALSE)        ## delete counts that are not included in the visits
m_count[,COUNT:=ceiling(COUNT/N)][,N:=1]


## build example data sets

v <- m_visit[SITE_ID %in% c(1:250) & DATE>='2000-01-01' & DATE<='2005-01-01']
c <-  m_count[SITE_ID %in% c(1:250) & DATE>='2000-01-01' & DATE<='2005-01-01' & SPECIES %in% c(2,4)][,N:=NULL]

data.table::fwrite(v,file="m_visit.csv",row.names=FALSE)
data.table::fwrite(c,file="m_count.csv",row.names=FALSE)

# FOLLOW THESE STEPS WITH THE DATA PROVIDED AND TRY WITH YOUR OWN DATA SET
# NOTE THAT THESE DATA CONTAIN A SUBSET OF BUTTERFLY COUNT FOR TWO SPECIES 2 AND 4
# THAT HAVE BEEN MONITORED BETWEEN 2000 AND 2004 IN 140 sites across the UK
# THE MONITORING SEASON WAS FROM APRIL TO SEPTEMBER (e.g. 4 to 9), BUT FOR DEMONSTRATION
# YOU CAN GIVE IT A TRY IT FOR NOVEMBER TO AUGUST (e.g. 11 to 8) AND IT SHOULD WORK.

source(new_regionalgam_function.r)

### load your data monitoring visit and butterfly count.
m_visit <- data.table::fread("m_visit.csv",header=TRUE)
m_count <- data.table::fread("m_count.csv",header=TRUE)

### IMPORTANT! all function are now build on the data.table framework, so if your data are not in this format (e.g. data.frame), you need to convert them first.
## this is only required if you do not used data.table to load your data. Here is how it should be done
m_visit <- data.table::data.table(m_visit)
m_count <- data.table::data.table(m_count)

## =====
## 2. build a long time-series to cover all days for the period we are interested in
## define the full time-series of your analysis, where InitYear is the first year and
## LasYear the last year of monitoring in the data set.

ts_date <- ts_dwmy_table(InitYear=2000,LastYear=2010,WeekDay1='monday')

## Define your monitoring season, with StartMonth and EndMonth, StartDay and EndDay, if EndDay is not defined, the last day of the month
## will be used. If CompltSeason is set to TRUE, only these year with full monitoring season will be used. Anchor are extra zeros set at the
## begining and the end of the season to help closing the curve (length and lag are defining the weight of the Anchor) 
ts_season <- ts_monit_season(ts_date,StartMonth=4,EndMonth=9,StartDay=1,EndDay=NULL,CompltSeason=TRUE,Anchor=TRUE,AnchorLength=7,AnchorLag=7)

## The following two step need to done in this order
m_visit <- df_visit_season(m_visit,ts_season)
ts_season_visit <- ts_monit_site(ts_season,m_visit)

  # check the species available in your data set
  m_count[order(SPECIES),unique(SPECIES)]
  
# you can choose between species 2 and 4 here.
ts_season_count <- ts_monit_count_site(ts_season_visit,m_count,sp=2)

## compute the flight curve for the selected species in the previous step (e.g. sp=4)
## The NbrSample is the maximum number of site used to compute the regionalGAM
## MinVisit is the minimum number of visit a site needs to be included in the model
## MinOccur is the minimum number of visit with count that a site needs to be included in the model
## MaxTrial is the maximum of trial to fit the model, if your number of site is less then NbrSample
## increasing this will not help if it did not fit in the 3 first trials.
## GamFamily is the distribution of the error term used in the GAM.
## CompltSeason is a logical limiting modelling to complete season only.

ts_flight_curve <- flight_curve(ts_season_count,NbrSample=100,MinVisit=3,MinOccur=2,MinNbrSite=1,MaxTrial=3,FcMethod='regionalGAM',GamFamily='poisson',CompltSeason=TRUE)

  ## plot the flight curves
  plot(ts_flight_curve[M_YEAR==2000,trimDAYNO],ts_flight_curve[M_YEAR==2000,NM],type='l',xlab='Monitoring Year Day',ylab='Relative Abundance')
  c <- 2
  for(y in 2001:2015){
      points(ts_flight_curve[M_YEAR==y,trimDAYNO],ts_flight_curve[M_YEAR==y,NM],type='l',col=c)
      c <- c + 1
  }

  
## Use the output of the GAM model (flight curve) to impute values for missing counts, using a GLM 
site_year_sp_count <- impute_count(ts_season_count,ts_flight_curve)

  ## plot the fitted values and the observed observation
  plot(site_year_sp_count[SITE_ID==1 & M_YEAR==2003,DATE],site_year_sp_count[SITE_ID==1 & M_YEAR==2003,FITTED],col='blue',type='l',main='Season 2003',xlab='Monitoring Month',ylab='Fitted Count')
  points(site_year_sp_count[SITE_ID==1 & M_YEAR==2003,DATE],site_year_sp_count[SITE_ID==1 & M_YEAR==2003,COUNT],col='red')

  
## Compute the total number of butterfly days - abundance index per site and year.
butterfly_index <- butterfly_day(site_year_sp_count)

  ## plot the computed abundance indices as butterfly days (area under the curve of accumulated butterfly count)
  for(site in butterfly_index[,unique(SITE_ID)]){
  plot(butterfly_index[SITE_ID==site,M_YEAR],butterfly_index[SITE_ID==site,BUTTERFLY_DAY],main=paste("SITE:",site),xlab='Monitoring Year',ylab='Butterfly Days')
  Sys.sleep(0.5)
  }

## NOTE MONITORING YEAR'S NAME IS BASED ON THE YEAR WHEN IT START.


## test full time series allocation with real data 

m_count[,max(data.table::year(VISIT_DATE))]

d <- ts_dwmy_table(2012,2015)
d_season <- ts_monit_season(d,8,5)
d_series <- ts_site_visit(d_season, m_visit,Anchor=TRUE,AnchorLength=7,AnchorLag=10) #augment with zeros and anchors
c_t_day <- fd_count_perday(m_count,m_visit)
c_site_series <- ts_count_site_visit(d_series,c_t_day,sp=2)
f_curve_series <- flight_curve(c_site_series,NbrSample=100,MinVisit=3,MinOccur=2,MaxTrial=3,GamFamily='poisson',FullSeason=TRUE)

plot(f_curve_series[SEASON_YEAR==2012,trimDAYNO],f_curve_series[SEASON_YEAR==2012,NM],col='red',type='l')
points(f_curve_series[SEASON_YEAR==2013,trimDAYNO],f_curve_series[SEASON_YEAR==2013,NM],col='red',type='l')
points(f_curve_series[SEASON_YEAR==2014,trimDAYNO],f_curve_series[SEASON_YEAR==2014,NM],col='red',type='l')


dev.new()
plot(f_curve_series[SEASON_YEAR==2010,trimDAYNO],f_curve_series[SEASON_YEAR==2010,NM],type='n',ylim=c(0,f_curve_series[,max(NM)]),ylab='Relative abundance',xlab='Julian day')
points(f_curve_series[SEASON_YEAR==2010,trimDAYNO],f_curve_series[SEASON_YEAR==2010,NM],type='l')
points(f_curve_series[SEASON_YEAR==2011,trimDAYNO],f_curve_series[SEASON_YEAR==2011,NM],col='red',type='l')
points(f_curve_series[SEASON_YEAR==2012,trimDAYNO],f_curve_series[SEASON_YEAR==2012,NM],col='blue',type='l')
points(f_curve_series[SEASON_YEAR==2013,trimDAYNO],f_curve_series[SEASON_YEAR==2013,NM],col='magenta',type='l')
points(f_curve_series[SEASON_YEAR==2015,trimDAYNO],f_curve_series[SEASON_YEAR==2015,NM],col='goldenrod',type='l')

## object size checkS
sort(sapply(ls(),function(x){object.size(get(x))}))

par(mfrow=c(2,2))
d <- ts_dwmy_table(2014,2015)
d_season <- ts_monit_season(d,4,9)
for(sim in 1:4){
 d_season_count <- sim_butterfly_count(d_season,GenNumb=2,PeakPos=c(10,55),TotalEmerg=c(100,200))
 plot(d_season_count$DAY_SINCE,d_season_count$TOTAL_DAY_COUNT,type='l',col='magenta')

 m_day <- sim_monitoring_visit(d_season,'weekly')
 points(m_day,d_season_count$TOTAL_DAY_COUNT[m_day],pch=19,col='red')


 m_day <- sim_monitoring_visit(d_season,'fortnightly')
 points(m_day,d_season_count$TOTAL_DAY_COUNT[m_day],pch=19,col='blue')

 m_day <- sim_monitoring_visit(d_season)
 points(m_day,d_season_count$TOTAL_DAY_COUNT[m_day],pch=19,col='green')
}

d <- ts_dwmy_table(2013,2015)
d_season <- ts_monit_season(d,10,6)

d_season_count <- sim_butterfly_count(d_season,GenNumb=1,PeakPos=c(40),TotalEmerg=c(100,200))
plot(d_season_count$DAY_SINCE,d_season_count$COUNT,type='l',col='magenta')

c_t_day <- sim_sites_count(d_season,NbrSite=100,GenNumb=1,PeakPos=c(40),TotalEmerg=c(100,200),MonitoringFreq=c("weekly"))
sim_count <- c_t_day[[1]]
sim_visit <- c_t_day[[2]]

d_series <- ts_site_visit(d_season,sim_visit,Anchor=TRUE,AnchorLength=7,AnchorLag=10)
sim_count_site_series <- ts_count_site_visit(d_series,sim_count,sp="sim")
f_curve_series <- flight_curve(sim_count_site_series,NbrSample=100,MinVisit=3,MinOccur=2,MaxTrial=3,GamFamily='poisson')

dev.new()
plot(f_curve_series[SEASON_YEAR==2013,trimDAYNO],f_curve_series[SEASON_YEAR==2013,NM],type='n',ylim=c(0,f_curve_series[,max(NM,na.rm=TRUE)]),ylab='Relative abundance',xlab='Julian day')
points(f_curve_series[SEASON_YEAR==2013,trimDAYNO],f_curve_series[SEASON_YEAR==2013,NM],col='red',type='l')
points(f_curve_series[SEASON_YEAR==2014,trimDAYNO],f_curve_series[SEASON_YEAR==2014,NM],col='magenta',type='l')
points(f_curve_series[SEASON_YEAR==2015,trimDAYNO],f_curve_series[SEASON_YEAR==2015,NM],col='blue',type='l')



par(mfrow=c(1,2))
plot(d_season_count[FULL_SEASON!=0,sum(COUNT),by=SEASON_YEAR])
plot(d_season_count[FULL_SEASON!=0,sum(MONITORED_COUNT,na.rm=TRUE),by=SEASON_YEAR])

plot(d_season_count$DAY_SINCE,d_season_count$COUNT,col='magenta',pch=19,type='l')
points(d_season_count$DAY_SINCE,d_season_count$MONITORED_COUNT,col='blue',pch=19)


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

monitoring_day <- d_season[complete_season!=0,sample(DAY_SINCE,1),by=.(SEASON_YEAR,WEEK)]

d_season[DAY_SINCE %in% monitoring_day$V1,MONITORED_COUNT:=COUNT]

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





fit_gam <- function(dataset_y, NbrSample=NbrSample, GamFamily=GamFamily, MaxTrial=MaxTrial){

        if (length(dataset_y[,unique(SITE)]) > NbrSample) {
            sp_data_all <- data.table::copy(dataset_y[SITE %in% sample(dataset_y[,unique(SITE)],NbrSample,replace=FALSE),])
        }else{
            sp_data_all <- data.table::copy(dataset_y)
        }

## fit GAM model ##
        print(paste("Fitting the flight curve-GAM for species",as.character(sp_data_all$SPECIES[1]),"in year",sp_data_all$YEAR[1],"with",length(sp_data_all[,unique(SITE)]),"sites :",Sys.time()))

        if(length(sp_data_all[,unique(SITE)])>1){
            gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE) -1, data=sp_data_all, family=GamFamily), silent = TRUE)
        }else {
            gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1, data=sp_data_all, family=GamFamily), silent = TRUE)
        }

## subsequent trials if not previous did not converge ##
        t <- 2
        while(class(gam_obj_site)[1] == "try-error" & t<=MaxTrial){

            if (length(dataset_y[,unique(SITE)]) > NbrSample) {
                sp_data_all <- data.table::copy(dataset_y[SITE %in% sample(dataset_y[,unique(SITE)],NbrSample,replace=FALSE),])
            }else{
                sp_data_all <- data.table::copy(dataset_y)
            }

            print(paste("Fitting the flight curve-GAM for species",as.character(sp_data_all$SPECIES[1]),"in year",sp_data_all$YEAR[1],"with",length(sp_data_all[,unique(SITE)]),"sites :",Sys.time(),"-> trial",t))

            if(length(sp_data_all[,unique(SITE)])>1){
                gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE) -1, data=sp_data_all, family=GamFamily), silent = TRUE)
            }else {
                gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1, data=sp_data_all, family=GamFamily), silent = TRUE)
            }

            t <- t+1
        }

## predict from fitted model ##
        if (class(gam_obj_site)[1] == "try-error") {
            print(paste("Error in fitting the flight period for",as.character(sp_data_all$SPECIES[1]),"in year", sp_data_all$YEAR[1],"; Model did not converge after",t,"trials"))
            sp_data_all[,c("FITTED","NM"):=.(NA,NA)]
        }else{
            sp_data_all[,FITTED:=mgcv::predict.gam(gam_obj_site, newdata = sp_data_all[,c("trimDAYNO", "SITE")], type = "response")]
            sp_data_all[SEASON==0L,FITTED:=0]

            if(sum(is.infinite(sp_data_all[,FITTED]))>0){
                sp_data_all[,c("FITTED","NM"):=.(NA,NA)]
            }else{
                sp_data_all[,SITE_SUM:=sum(FITTED),by=SITE]
                sp_data_all[,NM:=FITTED/SITE_SUM]
            }
        }

        f_curve <- sp_data_all[,.(SPECIES,SEASON,YEAR,MONTH,MONTH_DAY,WEEK,DAY_SINCE,trimDAYNO,NM)]
        data.table::setkey(f_curve)
        f_curve <- unique(f_curve)

    return(f_curve)
}