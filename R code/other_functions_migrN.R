
# Functions needed to run the mcmc sampler
##########################################

dXdt = function(times,x.input,theta.input){

  alpha.0.input = theta.input[1] # Log beta0
  alpha.2.input = theta.input[2]
  alpha.4.input = theta.input[3]
  alpha.5.input = theta.input[4]
  beta.env.input = theta.input[5]
  beta.wash = theta.input[6]
  gamma.input   = 0.2*7
  epsilon.input = theta.input[7]
  cont.input = 10^theta.input[8]
  decay.input = theta.input[9]
  water.half = 10^6
  lambda = theta.input[10]
  lambda_c = theta.input[11]

  beta.env.t = approxfun(x=c(0, 46, 118), y=c(beta.env.input, beta.env.input, beta.env.input-beta.wash))
  # 1st WASH intervention was achieved on october 3rd 2014 (46th time step)

  beta.t = exp(alpha.0.input)
  cont.env.t = x.input[4]/(x.input[4]+water.half)
  rain.amp.t = (1 + lambda*(rain_amp_t(times)))
  decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

  lam     = rep(0,7)
  
  lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.t(times)*cont.env.t*rain.amp.t*x.input[1] # contact 
  lam[2]  = eta(times)*s_capt*x.input[1]  # vaccination
  lam[3]  = gamma.input*x.input[2] # recover
  lam[4]  = epsilon.input*x.input[3] # immunity waning
  lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
  lam[6]  = decay.input.t*x.input[4] # environmental drivers
  lam[7]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.t(times)*cont.env.t*rain.amp.t*x.input[1] # contact 

  x.dot    = R%*%lam
  
  return( list( x.dot ) ) 
  
}

rk4 = function(t.input,y.input,h.input,f.input,theta.input){
  
  k1 = f.input(t.input,y.input,theta.input)
  k2 = f.input(t.input+h.input/2,y.input+h.input*k1/2,theta.input)
  k3 = f.input(t.input+h.input/2,y.input+h.input*k2/2,theta.input)
  k4 = f.input(t.input+h.input, y.input + h.input*k3, theta.input)
  
  return(y.input + (h.input/6)*(k1+2*k2+2*k3+k4))
  
}

log_post_theta.1 = function(theta.1.input,theta.3.input,I.int.input){
  
  tau.input = theta.3.input[1]
  disp.input = 10^(theta.3.input[2])
  
  loglike = sum(
              dnbinom(
                y,
                size= I.int.input*disp.input,
                mu=I.int.input*tau.input,
                log=T))
  
  theta.1.p = sum(dtnorm( theta.1.input,
                        mean  = theta.1.prior.mean,
                        lower = theta.1.lower,
                        upper = theta.1.upper,
                        sd    = sigma.theta.1,
                        log   = T ) )
  
  return(loglike+theta.1.p)
  
  
}

log_post_theta.2 = function(theta.2.input,theta.3.input,I.int.input){
  
  tau.input = theta.3.input[1]
  disp.input = 10^(theta.3.input[2])
  
  loglike = sum(
              dnbinom(
                y,
                size= I.int.input*disp.input,
                mu=I.int.input*tau.input,
                log=T))
  
  theta.2.p = sum(dtnorm( theta.2.input,
                          mean  = theta.2.lower,
                          lower = theta.2.lower,
                          upper = theta.2.upper,
                          sd    = sigma.theta.2,
                          log   = T ) )
  
  return(loglike+theta.2.p)
}


log_post_theta.3 = function(theta.3.input,I.int.input){
  
  tau.input = theta.3.input[1]
  disp.input = 10^(theta.3.input[2])
  
  loglike = sum(
              dnbinom(
                y,
                size= I.int.input*disp.input,
                mu=I.int.input*tau.input,
                log=T))

  theta.3.p = sum(dtnorm( theta.3.input,
                          mean  = theta.3.prior.mean,
                          lower = theta.3.lower,
                          upper = theta.3.upper,
                          sd    = sigma.theta.3,
                          log   = T ) )

  return(loglike+theta.3.p)
}

# Other functions
#################

rad_lag2 <- function(lag=0, model=FALSE){ 
  # The raw data were extracted from the VIIRS VNP46A1 product
  # Those data are freely available from EARTHDATA (https://search.earthdata.nasa.gov/search?portal=idn&q=VNP46A1_NRT_1)
  rad_smooth <- rad_kal %>%
                  dplyr::mutate(
                    day=ISOweek::ISOweekday(date),
                    week=lubridate::isoweek(date),
                    year=lubridate::isoyear(date)) %>%
                  dplyr::mutate(
                    week_decimal=week+(day-1)/7)

  rad_smooth$rad_mean <- NA
  rad_smooth$rad_mean[7:nrow(rad_smooth)] <- lapply(7:nrow(rad_smooth), function(x){
    return(
      weighted.mean(rad_smooth$rad[(x-6):x], rad_smooth$n[(x-6):x], na.rm=TRUE)
      )}) %>%
  do.call("c", .)
  rad_smooth$n_sum <- NA
  rad_smooth$n_sum[7:nrow(rad_smooth)] <- lapply(7:nrow(rad_smooth), function(x){
    return(
      sum(rad_smooth$n[(x-6):x], na.rm=TRUE)
      )}) %>%
  do.call("c", .)
  rad_smooth$index <- 1:nrow(rad_smooth)

  rad_d_spline <- rad_smooth %>%
                  dplyr::filter(n_sum>0 & !is.na(n_sum)) %>%
                  data.frame()

  rad_mod <- rad_smooth %>%
            dplyr::arrange(date) %>%
            dplyr::mutate(
              clock=as.numeric(difftime(date, lubridate::ymd("2012-01-01"), units="weeks")),
              index=as.numeric(difftime(date, lubridate::ymd("2012-01-01"), units="days"))+1)

  mod <- mgcv::gam(rad_mean ~ s(week_decimal, bs = "cc"),
          data=rad_mod)

  if(!model){
    rad_mod <- rad_mod %>%
                dplyr::filter(index>15 & index<2557) %>%
                dplyr::mutate(pred=predict(mod, newdata=.))
    temp <- rad_mod %>%
              dplyr::select(date, clock, pred)
    eps <- 1e-10
    temp$clock <- temp$clock+lag
    pred <- approxfun(temp$clock, temp$pred)
    pred_0 <- pred(seq(min(temp$clock), max(temp$clock), by=0.01))
    pred_1 <- pred((seq(min(temp$clock), max(temp$clock), by=0.01)+eps))
    d1 <- (pred_1-pred_0)/eps
    d_neg <- which(d1< -4)
    d_index <- which(diff(d_neg)>1)+1
    d1[d_neg[1:(d_index-1)]] <- approx(c((d_neg[1]-1), (d_neg[d_index-1]+1)), c(d1[(d_neg[1]-1)], d1[(d_neg[d_index-1]+1)]), xout=d_neg[1:(d_index-1)])$y
    d1[d_neg[d_index:length(d_neg)]] <- approx(c(d_neg[d_index]-1, d_neg[length(d_neg)]+1), c(d1[d_neg[d_index]-1], d1[d_neg[length(d_neg)]+1]), xout=d_neg[d_index:length(d_neg)])$y
    # The positive ones
    d_pos <- which(diff(d1)>0.5|diff(d1)< -0.6)
    d1[d_pos[1]:d_pos[2]] <- approx(c((d_pos[1]-1), (d_pos[length(d_pos)]+1)), c(d1[(d_pos[1]-1)], d1[(d_pos[length(d_pos)]+1)]), xout=d_pos[1]:d_pos[2])$y

    rad_mod <- rad_mod %>%
                  dplyr::filter(index>15 & index<2557) %>%
                  dplyr::mutate(d1=approx(x=seq(min(temp$clock), max(temp$clock), by=0.01), y=d1, xout=rad_mod$clock)$y) %>%
                  dplyr::select(date, year, week_decimal, pred, d1)
    return(rad_mod)
  }else{
    temp <- rad_mod %>%
              dplyr::mutate(d1=gratia::derivatives(mod, newdata=.)$derivative)  %>%
              dplyr::mutate(clock=as.numeric(difftime(date, lubridate::dmy("11/11/2013"), units="weeks")))
    return(approxfun(x=(temp$clock), y=temp$d1))
  }
}

rain_lag2 <- function(lag=0, normalize=FALSE, regular=TRUE, raw=FALSE){
  # The raw data come from "Development of a 50-yr high-resolution global dataset of meteorological forcings for land" from Sheffield et al
  # They are freely available on http://hydrology.princeton.edu/data/index.html
  prcp_data <- prcp_sat %>%
                  dplyr::group_by(date) %>%
                  dplyr::summarize(
                    prcp=mean(prcp),
                    n=dplyr::n()) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(
                    day=ISOweek::ISOweekday(date),
                    week=lubridate::isoweek(date),
                    year=lubridate::isoyear(date)) %>%
                  dplyr::mutate(week_decimal=week+(day-1)/7)

  prcp_data$prcp_mean <- NA
  prcp_data$prcp_mean[7:nrow(prcp_data)] <- lapply(7:nrow(prcp_data), function(x){
    return(
      weighted.mean(prcp_data$prcp[(x-6):x], prcp_data$n[(x-6):x], na.rm=TRUE)
      )}) %>%
  do.call("c", .)
  prcp_data$n_sum <- NA
  prcp_data$n_sum[7:nrow(prcp_data)] <- lapply(7:nrow(prcp_data), function(x){
    return(
      sum(prcp_data$n[(x-6):x], na.rm=TRUE)
      )}) %>%
  do.call("c", .)

  model <- prcp_data %>%
            mgcv::gam(prcp_mean~s(week_decimal, bs="cc"),
              weights=n_sum,
              data=.)

  prcp_mod <- prcp_data %>%
                dplyr::mutate(pred=predict(model, newdata=.)) %>% 
                dplyr::mutate(clock=as.numeric(difftime(date, lubridate::dmy("11/11/2013"), units="weeks"))) %>%
                dplyr::filter(clock>=-10 & clock<=119) %>%
                dplyr::arrange(date) %>%
                dplyr::mutate(
                  st=(pred-mean(pred))/sd(pred),
                  norm=(pred+abs(min(pred)))) %>%
                dplyr::mutate(
                  norm=norm/max(norm),
                  clock=clock+lag)

  prcp_irr <- rain_lag(lag=0, corr=FALSE, trunc=FALSE) %>%
                dplyr::mutate(
                  isoW=lubridate::isoweek(date),
                  isoY=lubridate::isoyear(date),
                  isoD=ISOweek::ISOweek2date(
                    paste(isoY,
                      paste0("W", stringr::str_pad(isoW, width=2, side="left", pad="0")),
                    7,
                    sep="-"))) %>%
                dplyr::group_by(isoD) %>%
                dplyr::summarise(cum_rain=sum(pred)) %>%
                dplyr::mutate(clock=as.numeric(difftime(isoD, lubridate::dmy("11/11/2013"), units="weeks"))) %>%
                dplyr::filter(clock>=-10 & clock<=119) %>%
                dplyr::arrange(isoD) %>%
                dplyr::mutate(
                  st=(cum_rain-mean(cum_rain))/sd(cum_rain),
                  norm=(cum_rain+abs(min(cum_rain)))) %>%
                dplyr::mutate(
                  norm=norm/max(norm),
                  clock=clock+lag)

  prcp_raw <- prcp_data
  prcp_raw$prcp_mean <- NA
  prcp_raw$prcp_mean[7:nrow(prcp_raw)] <- lapply(7:nrow(prcp_raw), function(x){
    return(
      weighted.mean(prcp_raw$prcp[(x-6):x], prcp_raw$n[(x-6):x], na.rm=TRUE)
      )}) %>%
  do.call("c", .)
  prcp_raw$n_sum <- NA
  prcp_raw$n_sum[7:nrow(prcp_raw)] <- lapply(7:nrow(prcp_raw), function(x){
    return(
      sum(prcp_raw$n[(x-6):x], na.rm=TRUE)
      )}) %>%
  do.call("c", .)

  prcp_raw <- prcp_raw %>%
                dplyr::mutate(clock=as.numeric(difftime(date, lubridate::dmy("11/11/2013"), units="weeks"))) %>%
                dplyr::filter(clock>=-10 & clock<=119) %>%
                dplyr::arrange(date) %>%
                dplyr::mutate(
                  norm=predict(loess(prcp_mean~clock, ., span=0.1))) %>%
                dplyr::mutate(
                  st=(norm-mean(norm))/sd(norm),
                  norm=(norm+abs(min(norm)))) %>%
                dplyr::mutate(
                  norm=norm/max(norm),
                  clock=clock+lag)

  if(normalize){
    if(raw){
      return(
        approxfun(x=prcp_raw$clock, y=prcp_raw$norm))
    }else{
      if(regular){
        return(
          approxfun(x=prcp_mod$clock, y=prcp_mod$norm))
      }else{
        return(
          approxfun(x=prcp_irr$clock, y=prcp_irr$norm))
      }
    }
  }else{
    if(raw){
      return(
        approxfun(x=prcp_raw$clock, y=prcp_raw$st))
    }else{
      if(regular){
        return(
          approxfun(x=prcp_mod$clock, y=prcp_mod$st))
      }else{
        return(
          approxfun(x=prcp_irr$clock, y=prcp_irr$st))
      }
    }
  }
}

sst_lag <- function(raw=FALSE){ 
  # The raw data come from AQUA MODIS
  # The raw data are freely available on the USGS website (https://earthexplorer.usgs.gov/)
	sst <- as.data.frame(temp) %>%
			dplyr::group_by(date) %>%
			dplyr::summarize(
				sst=mean(sst, na.rm=TRUE),
				n=dplyr::n()) %>%
			dplyr::ungroup() %>%
			dplyr::right_join(.,
				data.frame(
					date=seq.Date(min(as.data.frame(temp)$date), max(as.data.frame(temp)$date), by="days")),
				by="date") %>%
			dplyr::mutate(
				week=lubridate::isoweek(date),
				year=lubridate::isoyear(date)) %>%
			dplyr::arrange(date)

	sst$sst_mean <- NA
	sst$sst_mean[7:nrow(sst)] <- lapply(7:nrow(sst), function(x){
	return(
	  weighted.mean(sst$sst[(x-6):x], sst$n[(x-6):x], na.rm=TRUE)
	  )}) %>%
	do.call("c", .)
	sst$N <- NA
	sst$N[7:nrow(sst)] <- lapply(7:nrow(sst), function(x){
	return(
	  sum(sst$n[(x-6):x], na.rm=TRUE)
	  )}) %>%
	do.call("c", .)

	sst <-  sst %>%
	          dplyr::mutate(
	          	day=ISOweek::ISOweekday(date),
	          	week_decimal=week+(day-1)/7,
	            clock=as.numeric(difftime(date, lubridate::dmy("11/11/2013"), units="weeks")),
	            index=as.numeric(difftime(date, lubridate::ymd("2012-01-01"), units="days"))+1)

	model <- sst %>%
	     	   mgcv::gam(sst_mean~s(index, bs="cr")+s(week_decimal, bs="cc")+
	          ti(index, week_decimal, bs=c("cr", "cc")),
	          weights=N,
	          data=.)

	l_model <- sst %>%
				loess(sst_mean~clock, ., weights=N, span=0.05)

	sst <- sst %>%
				dplyr::mutate(
					pred=predict(model, newdata=.),
					st=predict(l_model, newdata=.)) %>%
              	dplyr::filter(clock>=-1 & clock<=119) %>%
				dplyr::mutate(
					pred=(pred-mean(pred))/sd(pred),
              		st=(st-mean(st))/sd(st))
    
    if(raw){
	    return(approxfun(x=(sst$clock), y=sst$st))
    }else{
	    return(approxfun(x=(sst$clock), y=sst$pred))
    }
}

chlor_lag <- function(raw=FALSE, expand=FALSE){
  # The raw data come from MODIS Chlorophyll-a Concentration
  # The raw data are freely available on the USGS website (https://earthexplorer.usgs.gov/)
	if(expand){
		chlor_temp <- chlor_30
	}else{
		chlor_temp <- chlor
	}
	chlor_a <- as.data.frame(chlor_temp) %>%
			dplyr::group_by(date) %>%
			dplyr::summarize(
				chlor_a=mean(chlor_a, na.rm=TRUE),
				n=dplyr::n()) %>%
			dplyr::ungroup() %>%
			dplyr::right_join(.,
				data.frame(
					date=seq.Date(min(as.data.frame(chlor_temp)$date), max(as.data.frame(chlor)$date), by="days")),
				by="date") %>%
			dplyr::mutate(
				week=lubridate::isoweek(date),
				year=lubridate::isoyear(date)) %>%
			dplyr::arrange(date)

	chlor_a$chlor_mean <- NA
	chlor_a$chlor_mean[7:nrow(chlor_a)] <- lapply(7:nrow(chlor_a), function(x){
	return(
	  weighted.mean(chlor_a$chlor_a[(x-6):x], chlor_a$n[(x-6):x], na.rm=TRUE)
	  )}) %>%
	do.call("c", .)
	chlor_a$N <- NA
	chlor_a$N[7:nrow(chlor_a)] <- lapply(7:nrow(chlor_a), function(x){
	return(
	  sum(chlor_a$n[(x-6):x], na.rm=TRUE)
	  )}) %>%
	do.call("c", .)

	chlor_a <-  chlor_a %>%
		          dplyr::mutate(
		          	day=ISOweek::ISOweekday(date),
		          	week_decimal=week+(day-1)/7,
					clock=as.numeric(difftime(date, lubridate::dmy("11/11/2013"), units="weeks")),
		            index=as.numeric(difftime(date, lubridate::ymd("2012-01-01"), units="days"))+1)

	ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
	model <- chlor_a %>%
	     	   mgcv::gamm(chlor_mean~s(index, bs="cc")+s(week_decimal, bs="cc"),
	          correlation = nlme::corARMA(form = ~ 1|year, p = 1),
	          weights=N,
	          data=.,
	          verbosePQL=FALSE)

	chlor_a$pred <- predict(model$gam, newdata=chlor_a)
	l_model <- loess(chlor_mean~clock, chlor_a[chlor_a$N>0 & !is.na(chlor_a$N),], weights=N, span=0.05)
	chlor_a$st <- predict(l_model, newdata=data.frame(clock=chlor_a$clock))

	chlor_a <- chlor_a %>%
					dplyr::filter(clock>=-1 & clock<=119) %>%
					dplyr::mutate(
						pred=(pred-mean(pred, na.rm=TRUE))/sd(pred, na.rm=TRUE),
	              		st=(st-mean(st, na.rm=TRUE))/sd(st, na.rm=TRUE))
    
	if(raw){
		return(approxfun(x=(chlor_a$clock), y=chlor_a$st))
	}else{
		return(approxfun(x=(chlor_a$clock), y=chlor_a$pred))
	}
}

# After running the sampler, build the chain using the chunks saved by the sampler
##################################################################################

# For the names you use in this function do not include the suffix with the numbers
# It is used to reconstruct the chain with the right order

build_sample2 <- function(name=c("theta1_NAME_OF_YOUR_FILE", "theta2_NAME_OF_YOUR_FILE", "theta3_NAME_OF_YOUR_FILE")){
  file_list <- list.files("PATH TO THE FILE WHERE THE FILES ARE SAVED"))
  file_list <- gsub("rda", "", file_list)
  file_list <- gsub("\\.", "", file_list)
  theta1_list <- file_list[grepl(name[1], file_list)]
  theta2_list <- file_list[grepl(name[2], file_list)]
  theta3_list <- file_list[grepl(name[3], file_list)]

  theta1_list <- gsub(name[1], "", theta1_list)
  theta2_list <- gsub(name[2], "", theta2_list)
  theta3_list <- gsub(name[3], "", theta3_list)

  theta1_list <- theta1_list[grepl("^[[:digit:]]+$", theta1_list)] %>%
          as.integer() %>%
          sort()
  theta2_list <- theta2_list[grepl("^[[:digit:]]+$", theta2_list)] %>%
          as.integer() %>%
          sort()
  theta3_list <- theta3_list[grepl("^[[:digit:]]+$", theta3_list)] %>%
          as.integer() %>%
          sort()

  lim1 <- max(theta1_list, na.rm=TRUE)
  lim2 <- max(theta2_list, na.rm=TRUE)
  lim3 <- max(theta3_list, na.rm=TRUE)
  
  theta1 <- data.frame()
  for(i in 1:lim1){
    load(file.path(path.root, paste0("PATH TO THE FILE WHERE THE FILES ARE SAVED", name[1], i, ".rda")))
    theta1 <- rbind(
          theta1,
          data.frame(theta.1.states))
  }
  theta2 <- data.frame()
  for(i in 1:lim2){
    load(file.path(path.root, paste0("PATH TO THE FILE WHERE THE FILES ARE SAVED", name[2], i, ".rda")))
    theta2 <- rbind(
          theta2,
          data.frame(theta.2.states))
  }
  theta3 <- data.frame()
  for(i in 1:lim3){
    load(file.path(path.root, paste0("PATH TO THE FILE WHERE THE FILES ARE SAVED", name[3], i, ".rda")))
    theta3 <- rbind(
          theta3,
          data.frame(theta.3.states))
  }

  return(
    list(
      theta1=theta1,
      theta2=theta2,
      theta3=theta3))
}

# Producing the estimates to assess the impact of the intervention (RUN THOSE ON A CLUSTER)
###########################################################################################

fit_impact_intervention <- function(){
  env_sample <- build_sample2(name=paste0(c("theta1_", "theta2_", "theta3_"), "SUFFIX YOU USED TO SAVE THE CHUNKS OF THE CHAIN"))

  theta1_sample <- env_sample$theta1 %>%
                  setNames(., c("alpha0", "alpha1", "alpha2",
                    "alpha4", "alpha5", "beta_env",
                    "beta_wash", "epsilon", "cont", "decay",
                    "lambda", "lambda_c"))

  theta2_sample <- env_sample$theta2 %>%
            setNames(., c("S0", "I0", "B0"))

  theta3_sample <- env_sample$theta3 %>%
            setNames(., c("tau", "disp"))
  
  dXdt = function(times,x.input,theta.input){
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.env.t = approxfun(x=c(0, 46, 118), y=c(beta.env.input, beta.env.input, beta.env.input-beta.wash))

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.t(times)*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = eta(times)*s_capt*x.input[1]  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.t(times)*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_v = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.env.t = approxfun(x=c(0, 46, 118), y=c(beta.env.input, beta.env.input, beta.env.input-beta.wash))

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.t(times)*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = 0  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.t(times)*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_w = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = eta(times)*s_capt*x.input[1]  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_w_v = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = 0  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  rad_t <- rad_lag2(model=TRUE)
  rain_t <- rain_lag2(lag=0, normalize=FALSE, regular=FALSE, raw=TRUE)
  rain_amp_t <- rain_lag2(lag=0, normalize=TRUE, regular=FALSE, raw=TRUE)
  sst_t <- sst_lag(raw=TRUE)
  chlor_t <- chlor_lag(expand=TRUE)

  # Adjust depending on your assumptions for the penalty accounting for the spatially tageted vaccination 
  s_capt=0.7
  # Adjust depending on your assumptions on the ratio between susceptible and infected 
  ratio_S_I=40

  total.time = 118 #total weeks 
  t.obs.seq  = seq(0,total.time,by=1)
  N          = 262963
  y          = cases

  # Change which vacci_step to use based on the assumption on vaccine effectiveness (see README.md)
  vacc.ag    = vacci_step3_s1$n

  ### VACCINATION RATE
  eta = stepfun( x = t.obs.seq, 
                 y = c(0,vacc.ag/N,0))
                   y = c(0,vacc.ag/N,0) )

  R1 = c(-1, 1, 0, 0, 0) # contact s to i 
  R2 = c(-1, 0, 1, 0, 0) # vax s to r
  R3 = c(0,-1, 1, 0, 0) # recover/remove i to r
  R4 = c(1, 0, -1, 0, 0) # immunity waning
  R5 = c(0, 0, 0, 1, 0) # fecal contamination
  R6 = c(0, 0, 0, -1, 0) # environmental drivers
  R7 = c(1, 0, 0, 0, 0) # migration of susceptibles
  R8 = c(0, 1, 0, 0, 0) # migration of infected
  R9 = c(0, 0, 0, 0, 1) # new infected

  R  = cbind(R1, R2, R3, R4, R5, R6, R7, R8, R9)

  dt            = 0.01 
  X.mesh        = round(seq(0, total.time, by = dt),2)
  t.obs.seq  = seq(0,total.time,by=1)

  nsim <- 10000
  set.seed(654321)
  index <- sample(1:nrow(theta1_sample), nsim, replace=TRUE)

  est <- function(x){
    theta.1.temp <- as.numeric(theta1_sample[index[x],])
    theta.2.temp <- as.numeric(theta2_sample[index[x],])
    theta.3.temp <- as.numeric(theta3_sample[index[x],])
    tau.temp = theta.3.temp[1]
    disp.temp = 10^(theta.3.temp[2])
    # No changes
    init.curr     = c(theta.2.temp[1:2], N-theta.2.temp[1]-theta.2.temp[2], 10^theta.2.temp[3], 0)
    X.ode.curr    = deSolve::ode(init.curr,X.mesh,dXdt,theta.1.temp,"rk4")[,-1]
    I.ode.curr    = X.ode.curr[-1,5]
    I.int.curr    = colSums( matrix(c(0, diff(I.ode.curr)),nrow=1/dt)) 
    # No vacci
    X.ode.vacci    = deSolve::ode(init.curr,X.mesh,dXdt_v,theta.1.temp,"rk4")[,-1]
    I.ode.vacci    = X.ode.vacci[-1,5]
    I.int.vacci    = colSums( matrix(c(0, diff(I.ode.vacci)),nrow=1/dt)) 
    # No WASH
    X.ode.wash    = deSolve::ode(init.curr,X.mesh,dXdt_w,theta.1.temp,"rk4")[,-1]
    I.ode.wash    = X.ode.wash[-1,5]
    I.int.wash    = colSums( matrix(c(0, diff(I.ode.wash)),nrow=1/dt)) 
    # No WASH/vacci
    X.ode.wash.vacci    = deSolve::ode(init.curr,X.mesh,dXdt_w_v,theta.1.temp,"rk4")[,-1]
    I.ode.wash.vacci    = X.ode.wash.vacci[-1,5]
    I.int.wash.vacci    = colSums( matrix(c(0, diff(I.ode.wash.vacci)),nrow=1/dt)) 

    rain.amp.t = (1 + theta.1.temp[11]*(rain_amp_t(0:118)))
    cont.env.t = X.ode.curr[X.mesh %in% c(0:118),4]/(X.ode.curr[X.mesh %in% c(0:118),4]+10^6)
    beta.env.t = approxfun(x=c(0, 46, 118), y=c(theta.1.temp[6], theta.1.temp[6], theta.1.temp[6]-theta.1.temp[7]))
    
    return(
      list(
        timeseries=data.frame(
          t=1:118,
          I=rnbinom(
            length(I.int.curr),
            size=I.int.curr*disp.temp,
            mu=I.int.curr*tau.temp),
          I_v=rnbinom(
            length(I.int.vacci),
            size=I.int.vacci*disp.temp,
            mu=I.int.vacci*tau.temp),
          I_w=rnbinom(
            length(I.int.wash),
            size=I.int.wash*disp.temp,
            mu=I.int.wash*tau.temp),
          I_w_v=rnbinom(
            length(I.int.wash.vacci),
            size=I.int.wash.vacci*disp.temp,
            mu=I.int.wash.vacci*tau.temp),
          i=I.int.curr,
          i_v=I.int.vacci,
          i_w=I.int.wash,
          i_w_v=I.int.wash.vacci,
          v=I.int.curr-I.int.vacci,
          w=I.int.curr-I.int.wash,
          w_v=I.int.curr-I.int.wash.vacci),
        Nt=data.frame(
          N=apply(matrix(X.ode.curr[,1:3], ncol=3),1 , sum)),
        prop=data.frame(
          t=X.mesh,
          s=X.ode.curr[,1]/apply(matrix(X.ode.curr[,1:3], ncol=3),1 , sum),
          i=X.ode.curr[,2]/apply(matrix(X.ode.curr[,1:3], ncol=3),1 , sum),
          s_v=X.ode.vacci[,1]/apply(matrix(X.ode.vacci[,1:3], ncol=3),1 , sum),
          i_v=X.ode.vacci[,2]/apply(matrix(X.ode.vacci[,1:3], ncol=3),1 , sum),
          s_w=X.ode.wash[,1]/apply(matrix(X.ode.wash[,1:3], ncol=3),1 , sum),
          i_w=X.ode.wash[,2]/apply(matrix(X.ode.wash[,1:3], ncol=3),1 , sum),
          s_w_v=X.ode.wash.vacci[,1]/apply(matrix(X.ode.wash.vacci[,1:3], ncol=3),1 , sum),
          i_w_v=X.ode.wash.vacci[,2]/apply(matrix(X.ode.wash.vacci[,1:3], ncol=3),1 , sum)),
        foi = data.frame(
          t=0:118,
          foi_e=beta.env.t(0:118)*cont.env.t*rain.amp.t,
          foi_h=exp(theta.1.temp[1]*(X.ode.curr[X.mesh %in% 0:118, 1]/apply(matrix(X.ode.curr[X.mesh %in% 0:118, 1:3], ncol=3),1 , sum))))
      )
    )
  }

  ##########################################
  # We recommend running this on a cluster #
  ##########################################
  numCores <- parallel::detectCores()-1
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  `%dopar%` <- foreach::`%dopar%`

  estimate <- foreach::foreach(i=1:nsim) %dopar% {
    est(i)
  }

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  save(
    estimate,
    file="PATH TO THE FOLDER WHERE YOU WANT TO SAVE THE ESTIMATES/estimate.rda")

}

migration_environment <- function(){
  env_sample <- build_sample2(name=paste0(c("theta1_", "theta2_", "theta3_"), "SUFFIX YOU USED TO SAVE THE CHUNKS OF THE CHAIN"))

  theta1_sample <- env_sample$theta1 %>%
                  setNames(., c("alpha0", "alpha1", "alpha2",
                    "alpha4", "alpha5", "beta_env",
                    "beta_wash", "epsilon", "cont", "decay",
                    "lambda", "lambda_c"))

  theta2_sample <- env_sample$theta2 %>%
            setNames(., c("S0", "I0", "B0"))

  theta3_sample <- env_sample$theta3 %>%
            setNames(., c("tau", "disp"))
  
  dXdt_w_v = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = 0  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_h = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = 0
    migr_I = 0
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = 0  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_e = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) # contact 
    lam[2]  = 0  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3])

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_c = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    # 1st WASH intervention was achieved on october 3rd 2014 (46th time step)

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = eta(times)*s_capt*x.input[1]  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = 0 # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  rad_t <- rad_lag2(model=TRUE)
  rain_t <- rain_lag2(lag=0, normalize=FALSE, regular=FALSE, raw=TRUE)
  rain_amp_t <- rain_lag2(lag=0, normalize=TRUE, regular=FALSE, raw=TRUE)
  sst_t <- sst_lag(raw=TRUE)
  chlor_t <- chlor_lag(expand=TRUE)

  # Adjust depending on your assumptions for the penalty accounting for the spatially tageted vaccination 
  s_capt=0.7
  # Adjust depending on your assumptions on the ratio between susceptible and infected 
  ratio_S_I=40

  total.time = 118 #total weeks 
  t.obs.seq  = seq(0,total.time,by=1)
  N          = 262963
  y          = cases

  # Change which vacci_step to use based on the assumption on vaccine effectiveness (see README.md)
  vacc.ag    = vacci_step3_s1$n

  ### VACCINATION RATE
  eta = stepfun( x = t.obs.seq, 
                 y = c(0,vacc.ag/N,0))
                   y = c(0,vacc.ag/N,0) )

  R1 = c(-1, 1, 0, 0, 0) # contact s to i 
  R2 = c(-1, 0, 1, 0, 0) # vax s to r
  R3 = c(0,-1, 1, 0, 0) # recover/remove i to r
  R4 = c(1, 0, -1, 0, 0) # immunity waning
  R5 = c(0, 0, 0, 1, 0) # fecal contamination
  R6 = c(0, 0, 0, -1, 0) # environmental drivers
  R7 = c(1, 0, 0, 0, 0) # migration of susceptibles
  R8 = c(0, 1, 0, 0, 0) # migration of infected
  R9 = c(0, 0, 0, 0, 1) # new infected

  R  = cbind(R1, R2, R3, R4, R5, R6, R7, R8, R9)

  dt            = 0.01 
  X.mesh        = round(seq(0, total.time, by = dt),2)
  t.obs.seq  = seq(0,total.time,by=1)

  cont <- function(x){
    theta.1.temp <- as.numeric(theta1_sample[index[x],])
    theta.2.temp <- as.numeric(theta2_sample[index[x],])
    theta.3.temp <- as.numeric(theta3_sample[index[x],])
    tau.temp = theta.3.temp[1]
    disp.temp = 10^(theta.3.temp[2])
    # No intervention
    init.curr     = c(theta.2.temp[1:2], N-theta.2.temp[1]-theta.2.temp[2], 10^theta.2.temp[3], 0)
    X.ode.curr    = deSolve::ode(init.curr,X.mesh,dXdt_w_v,theta.1.temp,"rk4")[,-1]
    I.ode.curr    = X.ode.curr[-1,5]
    I.int.curr    = colSums( matrix(c(0, diff(I.ode.curr)),nrow=1/dt))
    # No interhuman transmission
    X.ode.hum    = deSolve::ode(init.curr,X.mesh,dXdt_h,theta.1.temp,"rk4")[,-1]
    I.ode.hum    = X.ode.hum[-1,5]
    I.int.hum    = colSums( matrix(c(0, diff(I.ode.hum)),nrow=1/dt)) 
    # No environmental transmission
    X.ode.env    = deSolve::ode(init.curr,X.mesh,dXdt_e,theta.1.temp,"rk4")[,-1]
    I.ode.env    = X.ode.env[-1,5]
    I.int.env    = colSums( matrix(c(0, diff(I.ode.env)),nrow=1/dt)) 
    # No environmental contamination
    X.ode.cont    = deSolve::ode(init.curr,X.mesh,dXdt_c,theta.1.temp,"rk4")[,-1]
    I.ode.cont    = X.ode.cont[-1,5]
    I.int.cont    = colSums( matrix(c(0, diff(I.ode.cont)),nrow=1/dt)) 
    return(
      list(
        timeseries=data.frame(
          t=1:118,
          i=I.int.curr,
          i_h=I.int.hum,
          i_e=I.int.env,
          i_c=I.int.cont,
          h=I.int.curr-I.int.hum,
          e=I.int.curr-I.int.env,
          c=I.int.curr-I.int.cont),
        prop=data.frame(
          t=X.mesh,
          net_decay=theta.1.temp[10]-exp(theta.1.temp[3]+theta.1.temp[4]*sst_t(X.mesh) + theta.1.temp[5]*chlor_t(X.mesh)),
          s=X.ode.curr[,1]/apply(matrix(X.ode.curr[,1:3], ncol=3),1 , sum),
          i=X.ode.curr[,2]/apply(matrix(X.ode.curr[,1:3], ncol=3),1 , sum),
          b=X.ode.curr[,4],
          foi_h=exp(theta.1.temp[1])*X.ode.curr[,2]/apply(matrix(X.ode.curr[,1:3], ncol=3),1 , sum),
          foi_e=theta.1.temp[6]*(X.ode.curr[,4]/(X.ode.curr[,4]+10^6))*(1+theta.1.temp[11]*(rain_amp_t(X.mesh))),
          s_h=X.ode.hum[,1]/apply(matrix(X.ode.hum[,1:3], ncol=3),1 , sum),
          i_h=X.ode.hum[,2]/apply(matrix(X.ode.hum[,1:3], ncol=3),1 , sum),
          b_h=X.ode.hum[,4],
          s_e=X.ode.env[,1]/apply(matrix(X.ode.env[,1:3], ncol=3),1 , sum),
          i_e=X.ode.env[,2]/apply(matrix(X.ode.env[,1:3], ncol=3),1 , sum),
          b_e=X.ode.env[,4],
          s_c=X.ode.cont[,1]/apply(matrix(X.ode.cont[,1:3], ncol=3),1 , sum),
          i_c=X.ode.cont[,2]/apply(matrix(X.ode.cont[,1:3], ncol=3),1 , sum),
          b_c=X.ode.cont[,4]) %>%
          dplyr::filter(t %in% seq(0, 118, by=0.5))
      )
    )
  }

  numCores <- parallel::detectCores()-1
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)

  contrib <- foreach::foreach(
    i=1:nsim,
    .export='%>%') %dopar% {
    cont(i)
  }

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  save(
    contrib,
    file="PATH TO THE FOLDER WHERE YOU WANT TO SAVE THE ESTIMATES/contrib.rda")

}

alternative_vaccination_strategies <- function(){
  env_sample <- build_sample2(name=paste0(c("theta1_", "theta2_", "theta3_"), "SUFFIX YOU USED TO SAVE THE CHUNKS OF THE CHAIN"))

  theta1_sample <- env_sample$theta1 %>%
                  setNames(., c("alpha0", "alpha1", "alpha2",
                    "alpha4", "alpha5", "beta_env",
                    "beta_wash", "epsilon", "cont", "decay",
                    "lambda", "lambda_c"))

  theta2_sample <- env_sample$theta2 %>%
            setNames(., c("S0", "I0", "B0"))

  theta3_sample <- env_sample$theta3 %>%
            setNames(., c("tau", "disp"))
  
  dXdt_w_v = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = 0  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
    
  }

  dXdt_var = function(times,x.input,theta.input){
    
    alpha.0.input = theta.input[1] # Log beta0
    alpha.1.input = theta.input[2] # Time varying component
    alpha.2.input = theta.input[3]
    alpha.4.input = theta.input[4]
    alpha.5.input = theta.input[5]
    beta.env.input = theta.input[6]
    beta.wash = theta.input[7]
    gamma.input   = 0.2*7
    epsilon.input = theta.input[8]
    cont.input = 10^theta.input[9]
    decay.input = theta.input[10]
    water.half = 10^6
    lambda = theta.input[11]
    lambda_c = theta.input[12]

    beta.t = exp(alpha.0.input)
    migr_S = ratio_S_I*alpha.1.input*rad_t(times)
    migr_I = alpha.1.input*rad_t(times)
    cont.env.t = x.input[4]/(x.input[4]+water.half)
    rain.amp.t = (1 + lambda*(rain_amp_t(times)))
    decay.input.t = decay.input-exp(alpha.2.input+alpha.4.input*sst_t(times) + alpha.5.input*chlor_t(times))

    lam     = rep(0,9)
    
    lam[1]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1] # contact 
    lam[2]  = eta_var(times)*s_capt*x.input[1]  # vaccination
    lam[3]  = gamma.input*x.input[2] # recover
    lam[4]  = epsilon.input*x.input[3] # immunity waning
    lam[5]  = cont.input*(1 + lambda_c*(rain_amp_t(times)))*x.input[2] # fecal contamination
    lam[6]  = decay.input.t*x.input[4] # environmental drivers
    lam[7]  = migr_S # migrating susceptibles
    lam[8]  = migr_I # migrating infected
    lam[9]  = beta.t*prod(x.input[1:2])/sum(x.input[1:3]) + beta.env.input*cont.env.t*rain.amp.t*x.input[1]

    x.dot    = R%*%lam
    
    return( list( x.dot ) ) 
  }

  rad_t <- rad_lag2(model=TRUE)
  rain_t <- rain_lag2(lag=0, normalize=FALSE, regular=FALSE, raw=TRUE)
  rain_amp_t <- rain_lag2(lag=0, normalize=TRUE, regular=FALSE, raw=TRUE)
  sst_t <- sst_lag(raw=TRUE)
  chlor_t <- chlor_lag(expand=TRUE)

  # Adjust depending on your assumptions for the penalty accounting for the spatially tageted vaccination 
  s_capt=0.7
  # Adjust depending on your assumptions on the ratio between susceptible and infected 
  ratio_S_I=40

  total.time = 118 #total weeks 
  t.obs.seq  = seq(0,total.time,by=1)
  N          = 262963
  y          = cases

  # Change which vacci_step to use based on the assumption on vaccine effectiveness (see README.md)
  vacc.ag    = vacci_step3_s1$n

  ### VACCINATION RATE
  eta = stepfun( x = t.obs.seq, 
                 y = c(0,vacc.ag/N,0))
                   y = c(0,vacc.ag/N,0) )

  R1 = c(-1, 1, 0, 0, 0) # contact s to i 
  R2 = c(-1, 0, 1, 0, 0) # vax s to r
  R3 = c(0,-1, 1, 0, 0) # recover/remove i to r
  R4 = c(1, 0, -1, 0, 0) # immunity waning
  R5 = c(0, 0, 0, 1, 0) # fecal contamination
  R6 = c(0, 0, 0, -1, 0) # environmental drivers
  R7 = c(1, 0, 0, 0, 0) # migration of susceptibles
  R8 = c(0, 1, 0, 0, 0) # migration of infected
  R9 = c(0, 0, 0, 0, 1) # new infected

  R  = cbind(R1, R2, R3, R4, R5, R6, R7, R8, R9)

  dt            = 0.01 
  X.mesh        = round(seq(0, total.time, by = dt),2)
  t.obs.seq  = seq(0,total.time,by=1)

  vacci_res <- data.frame(
          t=comb$t,
          N=comb$N,
          avoid=rep(NA, length(nrow(comb))),
          l=rep(NA, length(nrow(comb))),
          u=rep(NA, length(nrow(comb))))

  # Function used to change the parameters of the vaccination strategy
  vacci_variation <- function(R=210370, t=1){
    `%>%` <- magrittr::`%>%`
    vacci <- data.frame(
          year=c(rep(2013, 6), rep(2014, 52), rep(2015, 53), rep(2016, 7)),
          week=c(47:52, 1:52, 1:53, 1:7),
          t=1:118,
          R=0)

    vacci$R[vacci$t==(t+1)] <- R

    return(
      vacci %>%
        dplyr::select(t, R))
  }

  set.seed(87654)
  index_v <- sample(1:nrow(theta1_sample), 500, replace=TRUE)
  comb <- expand.grid(
        N=seq(50000, 200000, length=6),
        t=seq(1, 118, length=14))

  pb = txtProgressBar(min = 0, max = nrow(comb), style = 3) 

  for(i in 1:nrow(comb)){
    # Redefining the step function based on the vaccination strategy
    eta_var = stepfun(
          x = 0:118, 
          y = c(0, vacci_variation(R=comb$N[i], t=comb$t[i])$R/N, 0))

    vacci_optim <- lapply(index_v, function(x){
      theta.1.temp <- as.numeric(theta1_sample[x,])
      theta.2.temp <- as.numeric(theta2_sample[x,])
      theta.3.temp <- as.numeric(theta3_sample[x,])
      init.curr     = c(theta.2.temp[1:2], N-theta.2.temp[1]-theta.2.temp[2], 10^theta.2.temp[3], 0)
      # No changes
      X.ode.curr    = deSolve::ode(init.curr,X.mesh,dXdt_w_v,theta.1.temp,"rk4")[,-1]
      I.ode.curr    = X.ode.curr[-1,5]
      I.int.curr    = colSums( matrix(c(0, diff(I.ode.curr)),nrow=1/dt)) 
      # Vacci changes
      X.ode.vacci    = deSolve::ode(init.curr,X.mesh,dXdt_var,theta.1.temp,"rk4")[,-1]
      I.ode.vacci    = X.ode.vacci[-1,5]
      I.int.vacci    = colSums( matrix(c(0, diff(I.ode.vacci)),nrow=1/dt)) 
      return(
        data.frame(
          samp=x,
          t=comb$t[i],
          imm=comb$N[i],
          v=sum(I.int.curr-I.int.vacci))
      )
    }) %>%
      do.call("rbind", .) %>%
      dplyr::group_by(t, imm) %>%
      dplyr::summarize(
        l=bayestestR::ci(v, method="HDI")$CI_low,
        u=bayestestR::ci(v, method="HDI")$CI_high,
        v=mean(v)) %>%
      dplyr::ungroup()

    vacci_res$avoid[i] <- vacci_optim$v
    vacci_res$l[i] <- vacci_optim$l
    vacci_res$u[i] <- vacci_optim$u

    setTxtProgressBar(pb,i)

  }
  close(pb)

  save(
    vacci_res,
    file="PATH TO THE FOLDER WHERE YOU WANT TO SAVE THE ESTIMATES/vacci_res.rda")

}
