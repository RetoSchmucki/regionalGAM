# FOLLOW THESE STEPS WITH THE DATA PROVIDED AND TRY WITH YOUR OWN DATA SET
# NOTE THAT THESE DATA CONTAIN A SUBSET OF BUTTERFLY COUNT FOR TWO SPECIES 2 AND 4
# THAT HAVE BEEN MONITORED BETWEEN 2000 AND 2004 IN 140 sites across the UK
# THE MONITORING SEASON WAS FROM APRIL TO SEPTEMBER (e.g. 4 to 9), BUT FOR DEMONSTRATION
# YOU CAN GIVE IT A TRY IT FOR NOVEMBER TO AUGUST (e.g. 11 to 8) AND IT SHOULD WORK.

# set you working directory and source the function script
source(new_regionalgam_function.R)

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

ts_date <- ts_dwmy_table(InitYear=2000,LastYear=2000,WeekDay1='monday')

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
ts_season_count <- ts_monit_count_site(ts_season_visit,m_count,sp=4)

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
for(y in 2001:2005){
  points(ts_flight_curve[M_YEAR==y,trimDAYNO],ts_flight_curve[M_YEAR==y,NM],type='l',col=c)
  c <- c + 1
}


## Use the output of the GAM model (flight curve) to impute values for missing counts, using a GLM 
site_year_sp_count <- impute_count(ts_season_count,ts_flight_curve)

## plot the fitted values and the observed observation
plot(site_year_sp_count[SITE_ID==1 & M_YEAR==2000,DATE],site_year_sp_count[SITE_ID==1 & M_YEAR==2000,FITTED],col='blue',type='l',main='Season 2003',xlab='Monitoring Month',ylab='Fitted Count')
points(site_year_sp_count[SITE_ID==1 & M_YEAR==2000,DATE],site_year_sp_count[SITE_ID==1 & M_YEAR==2000,COUNT],col='red')


## Compute the total number of butterfly days - abundance index per site and year.
butterfly_index <- butterfly_day(site_year_sp_count)

## plot the computed abundance indices as butterfly days (area under the curve of accumulated butterfly count)
for(site in butterfly_index[,unique(SITE_ID)]){
  plot(butterfly_index[SITE_ID==site,M_YEAR],butterfly_index[SITE_ID==site,BUTTERFLY_DAY],main=paste("SITE:",site),xlab='Monitoring Year',ylab='Butterfly Days')
  Sys.sleep(0.5)
}

## NOTE MONITORING YEAR'S NAME IS BASED ON THE YEAR WHEN IT START.