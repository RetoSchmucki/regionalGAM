### NORA's data
source("new_regionalgam/new_regionalgam_functions.R")
m_visit <- data.table::fread("new_regionalgam/d_visit.csv",header=TRUE)
m_count <- data.table::fread("new_regionalgam/d_count.csv",header=TRUE)

ts_date <- ts_dwmy_table(InitYear=2010,LastYear=2017,WeekDay1='monday')

ts_season <- ts_monit_season(ts_date,StartMonth=10,EndMonth=9,StartDay=1,EndDay=NULL,CompltSeason=TRUE,Anchor=TRUE,AnchorLength=7,AnchorLag=7)

m_visit <- df_visit_season(m_visit,ts_season)
ts_season_visit <- ts_monit_site(m_visit,ts_season)

m_count[order(SPECIES),unique(SPECIES)]

ts_season_count <- ts_monit_count_site(ts_season_visit,m_count,sp='Vanessa cardui')

ts_flight_curve <- flight_curve(ts_season_count,NbrSample=100,MinVisit=3,MinOccur=2,MinNbrSite=1,MaxTrial=3,FcMethod='regionalGAM',GamFamily='poisson',CompltSeason=TRUE)

## plot the flight curves
plot(ts_flight_curve[M_YEAR==2015,trimDAYNO],ts_flight_curve[M_YEAR==2015,NM],type='l',xlab='Monitoring Year Day',ylab='Relative Abundance')
c <- 2
for(y in 2015:2016){
  points(ts_flight_curve[M_YEAR==y,trimDAYNO],ts_flight_curve[M_YEAR==y,NM],type='l',col=c)
  c <- c + 1
}