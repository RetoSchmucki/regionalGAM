##==========================================
## Name Convention
##
##      FUNCTION: ts_snake_function()
##      ARGUMENT: CamelNotation
##      OBJECT: object_snake_name
##      VARIABLE NAME: UPPER_CASE
##
##      Date:   20.07.2017
##
##      Update: 14.08.2017
##
##==========================================

### check_names function is a generic function to verify if the provided data.table contains the required column names

check_names <- function(x, y){
        dt_names <- y %in% names(x)
        if(sum(dt_names)!=length(y)) {
            stop(paste('For this function, you need to have a variable named -', paste(y[!dt_names],collapse=' & '), '- in table',deparse(substitute(x)),'\n'))
        }
    }


###  ts_date_seq  function to build a full time series sequence of date since initial year to an ending year
###         to subset a time-series template for specified monitoring season, independent of 
###         the year
##
##  This function is used within the ts_dwmy_table() function

ts_date_seq = function(InitYear=1970,LastYear=format(Sys.Date(),"%Y")) {

    init_date <- as.Date(paste((InitYear-1), "01-01", sep = "-"))
    last_date <- as.Date(paste((as.numeric(LastYear)+1), "12-31", sep = "-"))

    date_series <- as.POSIXct(seq(from=init_date, to= last_date, by = "day"), format = "%Y-%m-%d")
    date_series <- date_series[!format(date_series,'%Y') %in% c((InitYear-1),as.numeric(LastYear)+1)]

    return(date_series)

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
                                DATE=data.table::as.IDate(date_seq),
                                DAY_SINCE=seq_along(date_seq))#,
        #                       YEAR=data.table::year(date_seq),
        #                       MONTH=data.table::month(date_seq),
        #                       DAY=data.table::mday(date_seq),
        #                       WEEK=data.table::isoweek(date_seq),
        #                       WEEK_DAY=w[data.table::wday(date_seq)])
        #                       

        return(dt_iso_dwmy)

    }


### set_anchor internal function used by ts_monit_season to define the days where anchor should be set to "0" where
###             AnchorLength and AnchorLag is provided by the function ts_monit_season()

set_anchor <- function(FirstObs,LastObs,AnchorLength=7,AnchorLag=7){
        
        x <- FirstObs$V1-(AnchorLength+AnchorLag)
        y <- FirstObs$V1-(AnchorLag+1)
        before_anchor <- c(apply(as.matrix(cbind(x,y)),1,function(w) w[1]:w[2]))
        
        x <- LastObs$V1+(AnchorLength+AnchorLag)
        y <- LastObs$V1+(AnchorLag+1)
        after_anchor <- c(apply(as.matrix(cbind(x,y)),1,function(w) w[2]:w[1]))

        anchor_day <- c(before_anchor,after_anchor)

        return(anchor_day)

    }


###  ts_monit_season  function to build a full time_series sequence of monitoring season, with specific starting 
###             and ending months and days, the season can start in year y and end in year y+1.
##
##      This function need the object (data.table) created by the "ts_dwmy_table()" function 
##      NOTE: if the monitoring season goes over two years (i.e. winter, December-January) the SEASON YEAR is shifted to set the monitoring season within a continuous year'
##      named according to the year at the start. Here the monitoring season is set in the middle or year' to give some room to set ANCHORs  before and after the monitoring 
##      season.

ts_monit_season = function(d_series,StartMonth=4,EndMonth=9,StartDay=1,EndDay=NULL,CompltSeason=TRUE,Anchor=TRUE,AnchorLength=7,AnchorLag=7){
        
        d_series <- data.table::copy(d_series)

        names(d_series) <- toupper(names(d_series))
        check_names(d_series,c("DATE"))

        ## NO leap year
        if(is.null(EndDay)) {
            EndDay <- max(data.table::mday(seq(from=as.Date(paste0('2017-',EndMonth,'-01')),by='day',length=31)))
        } 

        if (StartMonth < EndMonth){
            s_month <- c(StartMonth:EndMonth)
            month_days_out <- c(paste(StartMonth,c(0:(StartDay-1))),paste(EndMonth,c((EndDay+1):32)))
            y_series <- data.table::data.table(M_YEAR=as.factor(data.table::year(d_series$DATE)),
                                               M_SEASON=ifelse(((data.table::month(d_series$DATE)%in%(s_month)) & (!(paste(data.table::month(d_series$DATE),data.table::mday(d_series$DATE)) %in% month_days_out))),1L,0L))                               
        }

        if (StartMonth > EndMonth){
            s_month <- c(StartMonth:12,1:EndMonth)
            month_days_out <- c(paste(StartMonth,c(0:(StartDay-1))),paste(EndMonth,c((EndDay+1):32)))
            y_series <- data.table::data.table(M_YEAR=as.factor(ifelse(data.table::month(d_series$DATE)>=StartMonth-floor((12-length(s_month))/2),data.table::year(d_series$DATE),(data.table::year(d_series$DATE)-1))),
                                               M_SEASON=ifelse(((data.table::month(d_series$DATE)%in%(s_month)) & (!(paste(data.table::month(d_series$DATE),data.table::mday(d_series$DATE))%in%month_days_out))),1L,0L))                               
        }

        d_series <- d_series[,c("M_YEAR","M_SEASON") := y_series[,.(M_YEAR,(as.numeric(M_YEAR)*M_SEASON))]]

        if(isTRUE(CompltSeason)){       
            d_series[,START_END:=d_series[,ifelse(data.table::month(DATE)==StartMonth & data.table::mday(DATE)==StartDay,1L,0L)] + d_series[,ifelse(data.table::month(DATE)==EndMonth & data.table::mday(DATE)==EndDay,1L,0L)]]
            d_series[,M_SEASON:=ifelse(M_YEAR %in% (d_series[,sum(START_END),by=M_YEAR][V1==2,M_YEAR]),1L,0L)*M_SEASON][,START_END:=NULL]
            d_series[M_SEASON>0L,M_SEASON:=M_SEASON-(min(d_series[M_SEASON>0L,M_SEASON])-1)]
            d_series[, COMPLT_SEASON:=ifelse(M_YEAR %in% d_series[M_SEASON>0L,unique(M_YEAR)],1L,0L)]
        }
   
        if(isTRUE(Anchor)){
            first_obs <- d_series[M_SEASON>0L,min(DAY_SINCE),by=.(M_YEAR)]
            last_obs <- d_series[M_SEASON>0L,max(DAY_SINCE),by=.(M_YEAR)]
            anchor_day <- set_anchor(FirstObs=first_obs,LastObs=last_obs,AnchorLength=AnchorLength,AnchorLag=AnchorLag)
            d_series <- d_series[,ANCHOR:=0L][DAY_SINCE %in% anchor_day,ANCHOR:=1L][DAY_SINCE %in% anchor_day,COUNT:=0L]
        } 

        return(d_series)
    }


## df_visit_season function identify the monitoring season for each visit, the M_YEAR should then be used as reference to 
##              identify the site that have been monitored in specific monitoring year

df_visit_season <- function(m_visit,ts_season,DateFormat="%Y-%m-%d"){
            
            names(m_visit) <- toupper(names(m_visit))
            check_names(m_visit,c("DATE","SITE_ID"))
            m_visit <- data.table::copy(m_visit[,.(DATE,SITE_ID)])
            m_visit[,DATE:=data.table::as.IDate(as.Date(m_visit$DATE,format=DateFormat))]

            names(ts_season) <- toupper(names(ts_season))
            check_names(ts_season,c("DATE","M_YEAR"))

            season_year <- ts_season[,.(DATE,M_YEAR)] 

            data.table::setkey(m_visit,DATE)
            data.table::setkey(season_year,DATE)

            m_visit <- merge(m_visit,season_year,all.x=FALSE)

        return(m_visit)

    }

### ts_monit_site function to augment the time series in m_season with all sites and visits with "zeros", leaving all non visited day 
###                 with and <NA>, this can then be used to add the observed count for specific species
###                 only have time series for years when a site has been monitored.

ts_monit_site = function(m_visit,ts_season,DateFormat="%Y-%m-%d") {

            names(ts_season) <- toupper(names(ts_season))
            check_names(ts_season,c("DATE","M_YEAR","M_SEASON"))

            m_visit <- df_visit_season(m_visit,ts_season,DateFormat=DateFormat)

            data.table::setkey(ts_season,DATE)
            data.table::setkey(m_visit,DATE)

            r_year <-ts_season[,range(data.table::year(DATE))]
            m_visit <- m_visit[DATE %in% ts_season[M_SEASON>0L,DATE],]
        
            monit_syl <- m_visit[data.table::year(DATE)>=min(r_year) & 
                        data.table::year(DATE)<=max(r_year),
                        .(SITE_ID=.SD[,unique(SITE_ID)]),by=M_YEAR]

            data.table::setkey(monit_syl,M_YEAR,SITE_ID)
            data.table::setkey(ts_season,M_YEAR)

            ts_season_site <- merge(ts_season,monit_syl,by.x="M_YEAR",by.y="M_YEAR",allow.cartesian=TRUE)

            data.table::setkey(ts_season_site,DATE,SITE_ID)
            data.table::setkey(m_visit,DATE,SITE_ID)

            ts_season_site <- ts_season_site[m_visit,COUNT:=0L][M_SEASON==0L & ANCHOR==0L,COUNT:=NA]

        return(ts_season_site)

    }
  

### ts_count_site_visit function to generate a full time series of observed count,
###                     for each day since a starting and ending years of the defined
###                     time series

ts_monit_count_site = function(m_season_site,m_count,sp=1,DateFormat="%Y-%m-%d") {

        names(m_season_site) <- toupper(names(m_season_site))
        check_names(m_season_site,c("DATE","SITE_ID","COUNT"))

        names(m_count) <- toupper(names(m_count))
        check_names(m_count,c("DATE","SITE_ID","SPECIES","COUNT"))
        m_count[,DATE:=data.table::as.IDate(as.Date(m_count$DATE,format=DateFormat))]

        if(!sp %in% m_count[,unique(SPECIES)]){
            stop(paste("Species",sp,"is not found in your dataset, check your \"sp\" argument."))
        }else{
            m_sp_count <- m_count[SPECIES %in% sp,]
            data.table::setkey(m_sp_count,DATE,SITE_ID,SPECIES)
            data.table::setkey(m_season_site,DATE,SITE_ID)
            spcount_site_series <- m_season_site[m_sp_count,COUNT:=m_sp_count[,as.integer(COUNT)]]
            spcount_site_series[,SPECIES:=sp]
        }

        return(spcount_site_series)

    }


### fit_gam function to fit a GAM on series of count data
###         where the user can select the maximum number of site to be used to fit the model
###         the family of the model (error distribution) and the maximum number of trial to
###         perform if convergence has not been met

fit_gam <- function(dataset_y,NbrSample=NbrSample,GamFamily=GamFamily,MaxTrial=MaxTrial){

        if (length(dataset_y[,unique(SITE_ID)]) > NbrSample) {
            sp_data_all <- data.table::copy(dataset_y[SITE_ID %in% sample(dataset_y[,unique(SITE_ID)],NbrSample,replace=FALSE),])
        }else{
            sp_data_all <- data.table::copy(dataset_y)
        }

        ## fit GAM model ##
        # print(paste("Fitting the RegionalGAM for species",as.character(sp_data_all$SPECIES[1]),"and year",sp_data_all$M_YEAR[1],"with",length(sp_data_all[,unique(SITE_ID)]),"sites :",Sys.time()))

        # if(length(sp_data_all[,unique(SITE_ID)])>1){
        #     gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE_ID) -1, data=sp_data_all, family=GamFamily), silent = TRUE)
        # }else {
        #     gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1, data=sp_data_all, family=GamFamily), silent = TRUE)
        # }

        ## fit GAM model with maximum trials if previous did not converge ##
        tr <- 1
        gam_obj_site <- c()

        while((tr==1 | class(gam_obj_site)[1] == "try-error") & tr<=MaxTrial){

            if (length(dataset_y[,unique(SITE_ID)]) > NbrSample) {
                sp_data_all <- data.table::copy(dataset_y[SITE_ID %in% sample(dataset_y[,unique(SITE_ID)],NbrSample,replace=FALSE),])
            }else{
                sp_data_all <- data.table::copy(dataset_y)
            }

            print(paste("Fitting the RegionalGAM for species",as.character(sp_data_all$SPECIES[1]),"and year",sp_data_all$M_YEAR[1],"with",length(sp_data_all[,unique(SITE_ID)]),"sites :",Sys.time(),"-> trial",tr))

            if(length(sp_data_all[,unique(SITE_ID)])>1){
                gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr") + as.factor(SITE_ID) -1, data=sp_data_all, family=GamFamily), silent = TRUE)
            }else {
                gam_obj_site <- try(mgcv::gam(COUNT ~ s(trimDAYNO, bs = "cr")  -1, data=sp_data_all, family=GamFamily), silent = TRUE)
            }

            tr <- tr+1
        }

        ## predict from fitted model ##
        if (class(gam_obj_site)[1] == "try-error") {
            print(paste("Error in fitting the RegionalGAM for species",as.character(sp_data_all$SPECIES[1]),"and year", sp_data_all$M_YEAR[1],"; Model did not converge after",tr,"trials"))
            sp_data_all[,c("FITTED","NM"):=.(NA,NA)]
        }else{
            sp_data_all[,FITTED:=mgcv::predict.gam(gam_obj_site, newdata = sp_data_all[,c("trimDAYNO", "SITE_ID")], type = "response")]
            sp_data_all[M_SEASON==0L,FITTED:=0]

            if(sum(is.infinite(sp_data_all[,FITTED]))>0){
                sp_data_all[,c("FITTED","NM"):=.(NA,NA)]
            }else{
                sp_data_all[,SITE_SUM:=sum(FITTED),by=SITE_ID]
                sp_data_all[,NM:=round(FITTED/SITE_SUM,5)]
            }
        }

        f_curve <- sp_data_all[,.(SPECIES,DATE,DAY_SINCE,M_YEAR,M_SEASON,trimDAYNO,NM)]
        data.table::setkey(f_curve)
        f_curve <- unique(f_curve)

    return(f_curve)
}

### flight_curve function to compute the flight curve from a GAM (fit_gam) from series of counts
###             where the user can define the criteria required for site to be included in the flight
###             curve computation. So far, only one method is available, namely the regionalGAM.

flight_curve <- function(ts_season_count,NbrSample=100,MinVisit=3,MinOccur=2,MinNbrSite=1,MaxTrial=3,FcMethod='regionalGAM',GamFamily='poisson',CompltSeason=TRUE) {

    names(ts_season_count) <- toupper(names(ts_season_count))
    check_names(ts_season_count,c("COMPLT_SEASON","M_YEAR","SITE_ID","SPECIES","DATE","DAY_SINCE","M_SEASON","COUNT","ANCHOR"))

    if(isTRUE(CompltSeason)){
        ts_season_count <- ts_season_count[COMPLT_SEASON==1]
    }

        if(exists("f_pheno")){rm(f_pheno)}

        for (y in ts_season_count[,unique(as.integer(M_YEAR))]) {

            dataset_y <- ts_season_count[as.integer(M_YEAR)==y, .(SPECIES,SITE_ID,DATE,DAY_SINCE,M_YEAR,M_SEASON,COUNT,ANCHOR)]
            dataset_y[,trimDAYNO:=DAY_SINCE-min(DAY_SINCE)+1]

            ## filter for site with at least 3 visits and 2 occurrences
            visit_occ_site <- merge(dataset_y[!is.na(COUNT) & ANCHOR==0L,.N,by=SITE_ID],dataset_y[!is.na(COUNT) & ANCHOR==0L & COUNT>0,.N,by=SITE_ID],by="SITE_ID",all=TRUE)
            dataset_y <- data.table::copy(dataset_y[SITE_ID %in% visit_occ_site[N.x>=MinVisit&N.y>=MinOccur,SITE_ID]])

            if(dataset_y[,.N]<=MinNbrSite){
                dataset_y <- ts_season_count[as.integer(M_YEAR)==y, .(SPECIES,DATE,DAY_SINCE,M_YEAR,M_SEASON)]
                dataset_y[,trimDAYNO:=DAY_SINCE-min(DAY_SINCE)+1]
                f_curve <- dataset_y[,NM:=NA]
                data.table::setkey(f_curve,SPECIES,DAY_SINCE)
                f_curve <- unique(f_curve)
                print(paste("You have not enough sites with observations for estimating the flight curve for species",as.character(dataset_y$SPECIES[1]),"in", dataset_y$M_YEAR[1]))
            }else{
                if(FcMethod=='regionalGAM'){
                    f_curve <- fit_gam(dataset_y,NbrSample,GamFamily,MaxTrial)
                }else{
                    print("ONLY the regionalGAM method is available so far!")
                }
            }

            if ("f_pheno" %in% ls()) {
                f_pheno <- rbind(f_pheno, f_curve)
            }else {
                f_pheno <- f_curve
            }    
        }

    return(f_pheno)
}



### impute_count function to compute the Abundance Index across sites and years from 
###                 your count dataset and the regional flight curve


impute_count <- function(ts_season_count,ts_flight_curve) {
   
    if(exists("sp_ts_season_count")) {rm("sp_ts_season_count")}
    
    sp_ts_season_count <- data.table::copy(ts_season_count)
    sp_ts_season_count[,SPECIES:=ts_flight_curve$SPECIES[1]]

    data.table::setkey(sp_ts_season_count,DATE)
    data.table::setkey(ts_flight_curve,DATE)

    sp_count_flight <- merge(sp_ts_season_count,ts_flight_curve[,.(DATE,trimDAYNO,NM)],all.x=TRUE)
    data.table::setkey(sp_count_flight,M_YEAR,DATE,SITE_ID)

    for(y in sp_count_flight[,unique(as.integer(M_YEAR))]){
        
        sp_count_flight_y <-  data.table::copy(sp_count_flight[as.integer(M_YEAR)==y,])
        print(paste("Computing abundance indices for species",sp_count_flight_y[1,SPECIES],"monitored in year", sp_count_flight_y[1,M_YEAR],"across",sp_count_flight_y[unique(SITE_ID),.N],"sites:",Sys.time()))
        
        if(sp_count_flight_y[is.na(NM),.N]>0){
            tr<-1
            search_op <- rep(1:5,rep(2,5))*c(-1,1)
            alt_flight <- sp_count_flight[as.integer(M_YEAR)==y,.(M_YEAR,trimDAYNO,NM)]
           
            while(alt_flight[is.na(NM),.N]>0 & tr<20){
                alt_flight <- unique(sp_count_flight[as.integer(M_YEAR)==y+search_op[tr],.(M_YEAR,trimDAYNO,NM)])
                tr<-tr+1
            }
           
            if(alt_flight[is.na(NM),.N]>0){
                next(paste("No reliable flight curve available within a 5 year horizon of",sp_count_flight_y[1,M_YEAR,]))
            }else{
                warning(paste("We used the flight curve of",alt_flight[1,M_YEAR],"to compute abundance indices for year",sp_count_flight_y[1,M_YEAR,]))
                data.table::setkey(sp_count_flight_y,trimDAYNO)
                data.table::setkey(alt_flight,trimDAYNO)
                sp_count_flight_y[alt_flight,NM:=alt_flight[,NM]]
            }
        }

        sp_count_flight_y[M_SEASON==0L,COUNT:=NA]
        sp_count_flight_y[M_SEASON!=0L & NM==0,NM:=0.000001]
        non_zero <- sp_count_flight_y[,sum(COUNT,na.rm=TRUE),by=(SITE_ID)][V1>0,SITE_ID]
        zero <- sp_count_flight_y[,sum(COUNT,na.rm=TRUE),by=(SITE_ID)][V1==0,SITE_ID]

        if(length(non_zero)>1){
            if(sp_count_flight_y[unique(SITE_ID),.N]>1){
                glm_obj_site <- try(glm(COUNT ~ factor(SITE_ID) + offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,],
                    family = quasipoisson(link = "log"), control = list(maxit = 100)),silent=TRUE)
            } else {
                glm_obj_site <- try(glm(COUNT ~ offset(log(NM)) -1,data=sp_count_flight_y[SITE_ID %in% non_zero,],
                    family = quasipoisson(link = "log"), control = list(maxit = 100)),silent=TRUE)
            }
        }

        if (class(glm_obj_site)[1] == "try-error") {
            sp_count_flight_y[SITE_ID %in% non_zero,c("FITTED","COUNT_IMPUTED"):=.(NA,NA)]
            print(paste("Computation of abundance indices for year",sp_count_flight_y[1,M_YEAR,],"failed with the RegionalGAM, verify the data you provided for that year"))
            next()
        }else{
            sp_count_flight_y[SITE_ID %in% non_zero,FITTED:= predict.glm(glm_obj_site,newdata=sp_count_flight_y[SITE_ID %in% non_zero,],type = "response")]
        }

        sp_count_flight_y[SITE_ID %in% zero,FITTED:=0]
        sp_count_flight_y[is.na(COUNT),COUNT_IMPUTED:=FITTED][!is.na(COUNT),COUNT_IMPUTED:=as.numeric(COUNT)][M_SEASON==0L,COUNT_IMPUTED:=0] 

        data.table::setkey(sp_ts_season_count,SITE_ID,DAY_SINCE)
        data.table::setkey(sp_count_flight_y,SITE_ID,DAY_SINCE)

        if("FITTED" %in% names(sp_ts_season_count)){
            sp_ts_season_count[sp_count_flight_y,':='(trimDAYNO=i.trimDAYNO,NM=i.NM,FITTED=i.FITTED,COUNT_IMPUTED=i.COUNT_IMPUTED)]
        }else{
            sp_ts_season_count <- merge(sp_ts_season_count, sp_count_flight_y[,.(DAY_SINCE,SITE_ID,trimDAYNO,NM,FITTED,COUNT_IMPUTED)], all.x=TRUE) 
        }
    }

    return(sp_ts_season_count)
} 


### butterfly_day function to count cumulative butterfly count observed over one monitoring season.

butterfly_day <- function(sp_ts_season_count){

            b_day <- sp_ts_season_count[,sum(COUNT_IMPUTED),by=.(SPECIES,M_YEAR,SITE_ID)]
            data.table::setnames(b_day,"V1","BUTTERFLY_DAY")
        
        return(b_day)
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
            sim_season=unique(d_season$FULL_SEASON)
            }else{
            sim_season=unique(d_season$SEASON)
        }

        GenSim <- 1

        for (i in sim_season){

            if(i==0){next}

                t_s <- seq_along(1:d_season[SEASON==i,.N])
                emerg_curve <- sim_emerg_curve(t_s,PeakPos=PeakPos[GenSim],sdPeak=sdPeak[GenSim],sig=sig[GenSim],bet=bet[GenSim])
                emerg_count <- sim_emerg_count(emerg_curve,TotalEmerg=TotalEmerg[GenSim])
                adult_count <- sim_adult_count(emerg_count,MaxLife=MaxLife[GenSim])
                d_season[SEASON==i,COUNT:=adult_count[1:d_season[SEASON==i,.N]]]
        }

        cumul_count <- d_season[,COUNT]
        GenSim <- GenSim+1
            
        while(GenNumb>=GenSim) {

            for (i in sim_season){

                if(i==0){next}

                    t_s <- seq_along(1:d_season[SEASON==i,.N])
                    emerg_curve <- sim_emerg_curve(t_s,PeakPos=PeakPos[GenSim],sdPeak=sdPeak[GenSim],sig=sig[GenSim],bet=bet[GenSim])
                    emerg_count <- sim_emerg_count(emerg_curve,TotalEmerg=TotalEmerg[GenSim])
                    adult_count <- sim_adult_count(emerg_count,MaxLife=MaxLife[GenSim])
                    d_season[SEASON==i,COUNT:=adult_count[1:d_season[SEASON==i,.N]]]
            }
        
        cumul_count <- cumul_count + d_season[,COUNT]
        GenSim <- GenSim+1

        }

    d_season[,COUNT:=cumul_count]
    d_season[,SITE_ID:=1]
    d_season[,SPECIES:='sim']

    return(d_season)

}

# sim_monitoring_visit() simulates monitoring visits by volunteers, based on monitoring frequency set by the protocol c('weekly','fortnightly','monthly') or c('none') for all days within season

sim_monitoring_visit <- function(d_season, MonitoringFreq=c('none')){

    if(MonitoringFreq=='weekly'){

        is_even <- sample(c(1,2),1) == 2
        monitoring_day <- d_season[SEASON!=0L & (DAY_SINCE %% 2 == 0L)==is_even,sample(DAY_SINCE,1),by=.(SEASON_YEAR,WEEK)][,V1]
    
    }

    if(MonitoringFreq=='fortnightly'){
        
        is_even <- sample(c(1,2),1) == 2
        monitoring_day <- d_season[SEASON!=0L & (WEEK %% 2 == 0L)==is_even,sample(DAY_SINCE,1),by=.(SEASON_YEAR,WEEK)][,V1]
    
    }

    if(MonitoringFreq=='monthly'){
        is_even <- sample(c(1,2),1) == 2
        monitoring_day <- d_season[SEASON!=0L & (WEEK %% 2 == 0L)==is_even,sample(DAY_SINCE,1),by=.(SEASON_YEAR,MONTH)][,V1]
    
    }

    if(MonitoringFreq=='none'){

        is_even <- sample(c(1,2),1) == 2
        monitoring_day <- d_season[SEASON!=0L,DAY_SINCE]
    
    }

    return(monitoring_day)

}

### sim_sites_count function to simulate butterfly count across sites, assuming a common flight curve,
###                     but with potential shift in the position of the Peak

sim_sites_count <- function(d_season,NbrSite=10,FullSeason=TRUE,GenNumb=1,PeakPos=c(25,75),sdPeak=c(1,2),sig=0.15,bet=3,TotalEmerg=100,MaxLife=10,MonitoringFreq=c('none'),PerctSampled=100){

    site_count_list <- vector("list",NbrSite)
    site_visit_list <- vector("list",NbrSite)

    for(s in 1:(NbrSite)){
        d_season_count <- sim_butterfly_count(d_season,FullSeason=FullSeason,GenNumb=GenNumb,PeakPos=PeakPos,sdPeak=sdPeak,sig=sig,bet=bet,TotalEmerg=TotalEmerg,MaxLife=MaxLife)
        m_day <- sim_monitoring_visit(d_season,MonitoringFreq=MonitoringFreq)
        site_visit <- d_season_count[DAY_SINCE %in% m_day,.(DATE,YEAR,SITE_ID)][,SITE_ID:=SITE_ID+(s-1)]
        site_count <- d_season_count[DAY_SINCE %in% sample(m_day,round((length(m_day)*PerctSampled)/100),replace=FALSE),][,SITE_ID:=SITE_ID+(s-1)]
        site_count_list[[s]] <- site_count
        site_visit_list[[s]] <- site_visit
    }
    
    multisite_count <- data.table::rbindlist(site_count_list)
    multisite_visit <- data.table::rbindlist(site_visit_list)
    return(list(multisite_count,multisite_visit))
}


### build_map_obj builds map for specific region, using function from the sf() simple feature package

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

