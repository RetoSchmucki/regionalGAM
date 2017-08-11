m_count <- data.table::fread("C:/Users/RETOSCHM/OneDrive - Natural Environment Research Council/regionalGAM/development/Vcardui_desert.csv",header=TRUE)
m_count[,DATE:=data.table::as.IDate(as.Date(m_count$DATE,format="%d/%m/%Y"))]
m_count[,COUNT:=sum(COUNT),by=.(DAY,MONTH,YEAR,SITE_ID)]
m_count <- unique(m_count)


d <- ts_dwmy_table(2015,2017)
d_season <- ts_monit_season(d,10,6)

m_visit <- m_count[,.(SITE_ID,DATE,YEAR)]
m_visit <- unique(m_visit)

d_series <- ts_site_visit(d_season, m_visit,Anchor=TRUE,AnchorLength=7,AnchorLag=10) #augment with zeros and anchors
c_site_series <- ts_count_site_visit(d_series,m_count,sp=44)

f_curve_series <- flight_curve(c_site_series,NbrSample=100,MinVisit=3,MinOccur=3,MaxTrial=3,GamFamily='poisson',FullSeason=TRUE)
plot(f_curve_series[SEASON_YEAR==2015,trimDAYNO],f_curve_series[SEASON_YEAR==2015,NM])