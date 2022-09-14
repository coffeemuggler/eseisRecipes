## PART 0 - SCRIPT INFORMATION SECTION ----------------------------------------

#' The script can be used to pick seismic events as recorded by a network of 
#' sensors, to plot the waveforms, to locate the seismic source and to make 
#' comprehensive plots of signals, PSDs and location maps of the events.
#' 
#' The script is in a transient stage, it works but is not thoroughly 
#' documented, tested, debugged yet. So please use care when using it!
#' 
#' Author: Michael Dietze, Uni Goettingen (michael.dietze)
#' 
#' Script version: 0.2.0
#' 
#' Changes to previous versions
#'   - a
#' 
#' Dependencies: 
#'   - R 4.n
#'   - eseis 0.7.1
#'   
#'   Parameter descriptions:
#'   - opt, general options
#'     - service, service ID which must meet the service directory name where 
#'       the service's raw cube or centaur data are stored. Only directory 
#'       name is to be used, the path will be added, see below.
#'       
## PART 1 - SETUP PARAMETERS --------------------------------------------------

pars <- list(
  
  opt = list(
    pick = TRUE,
    pick_plot = TRUE,
    locate_plot = FALSE,
    pick_locate_plot = FALSE,
    locate = FALSE,
    station_skip = c(),
    pkgs = c("eseis", "raster")
  ),
  
  manual = list(
    enable = FALSE,
    start = as.POSIXct("2020-10-29 04:59:20 UTC"),
    duration = 5,
    snr_reject = FALSE
  ),
  
  path = list(
    station = paste0("/data/STORAGE/seismic/continent/",
                     "country/site/data/seismic/",
                     "station/station_info_site.txt"),
    sac = paste0("/data/STORAGE/Seismic/continent/country/",
                 "site/data/seismic/sac/"),
    picks = paste0("/data/DATA/username/site/picks/"),
    plot = paste0("/data/DATA/username/site/plot/"),
    dem = paste0("/data/STORAGE/Seismic/continent/country/",
                 "site/data/geodata/dem/dem_10.tif"),
    hs = paste0("/data/STORAGE/Seismic/continent/country/",
                "site/data/geodata/dem/hs_10.tif"),
    d_map = paste0("/data/STORAGE/Seismic/continent/country/",
                   "site/data/geodata/distance/D_20.rda"),
    aoi = paste0("/data/STORAGE/Seismic/continent/country/",
                 "site/data/geodata/aoi/aoi_main_20.tif"),
    coupling = paste0("/data/STORAGE/Seismic/continent/country/",
                      "site/data/seismic/coupling/coupling.rda")
  ),
  
  pick = list(
    pick_ID = "picks_01",
    t_process = as.POSIXct(x = c("2022-07-20 00:00:00",
                                 "2022-08-16 04:00:00"),
                           tz = "UTC"),
    t_snippet = 3600,
    t_buffer = c(0, 120),
    f = c(5, 40),
    sta = 0.5 * 200,
    lta = 3 * 60 * 200,
    on = 3,
    off = 1,
    freeze = TRUE,
    dur_min = 3,
    dur_max = 240,
    t_spacing = 5,
    t_coherent = 2,
    n_coherent = 3),
  
  locate = list(
    i_start = 1,
    buffer = c(5, 5),
    append = c(1, 1),
    f = c(5, 20),
    n_taper = 200,
    method = c("migrate", "amplitude")[1],
    v_wave = 800,
    q = 30,
    quantile = 0.99,
    snr_min = 3,
    snr_mean_min = 3,
    station_min = 3
  ),
  
  plot = list(
    f = c(1, 45),
    f2 = c(10, 20),
    p_taper = 0.1,
    col_map = colorRampPalette(colors = c("grey50", "darkblue", "blue",
                                          "green", "yellow")),
    width = 2000,
    height = 2000,
    res = 200)
)

## PART 2 - INITIATE PROCESSES ------------------------------------------------

## load required packages
x <- lapply(pars$opt$pkgs, require, character.only = TRUE)
rm(x)
Sys.setenv(TZ = "UTC")

## load station data
stations <- read.table(file = pars$path$station,
                       header = TRUE,
                       stringsAsFactors = FALSE)

## remove stations not to use
if(length(pars$opt$station_skip) > 0) {
  
  id_out <- which(stations$ID %in% pars$opt$station_skip)
  stations <- stations[-id_out,]
}

## create processing intervals
t <- seq(from = pars$pick$t_process[1],
         to = pars$pick$t_process[2],
         by = pars$pick$t_snippet)

## create dummy putput data set
picks_out <- data.frame(ID = NA,
                        start = Sys.time(),
                        duration = NA,
                        stalta = NA,
                        n_stations = NA)

## PART 3 - FIND PICKS --------------------------------------------------------
if(pars$opt$pick == TRUE) {
  
  print(paste0(Sys.time(), ", Start picking (n = ", length(t), ")"))
  
  for(i in 1:length(t)) {
    
    ## print progress
    print(paste0("  ", i))
    
    ## read data sets
    s <- lapply(X = stations$ID, FUN = function(x) {
      
      s_out <- try(aux_getevent(start = t[i] - pars$pick$t_buffer[1],
                                duration = pars$pick$t_snippet +
                                  pars$pick$t_buffer[1] +
                                  pars$pick$t_buffer[2],
                                station = x,
                                component = "BHZ",
                                dir = pars$path$sac), 
                   silent = TRUE)
      if(class(s_out)[1] != "try-error") {
        
        return(s_out)
      } else {
        
        return(NA)
      }
    })
    
    ## identify successful read attempts
    s_ok <- try(do.call(c, lapply(X = s, FUN = function(s) {
      
      !is.na(s)[1]
    })))
    
    ## remove unsuccessful signals
    s <- try(s[s_ok])
    
    ## demean data set
    s <- try(signal_demean(data = s),
             silent = TRUE)
    
    ## detrend data set
    s <- try(signal_detrend(data = s),
             silent = TRUE)
    
    ## filter data set
    s_f <- try(signal_filter(data = s,
                             f = pars$pick$f),
               silent = TRUE)
    
    ## detrend data set
    s_f <- try(signal_detrend(data = s_f),
               silent = TRUE)
    
    ## calculate signal envelope
    s_e <- try(signal_envelope(data = s_f),
               silent = TRUE)
    
    ## pick events at all stations
    p <- lapply(X = s_e, FUN = function(s_e, pars) {
      
      try(eseis::pick_stalta(data = s_e,
                             sta = pars$pick$sta,
                             lta = pars$pick$lta,
                             on = pars$pick$on,
                             off = pars$pick$off)$picks)
    }, pars)
    
    ## convert list to data frame
    p <- try(do.call(rbind, p))
    
    ## identify joint picks within time limits
    p_joint <- try(lapply(X = 1:nrow(p), FUN = function(ii, p, pars) {
      
      abs(p$start[ii] - p$start) < pars$pick$t_coherent
      
    }, p, pars))
    
    ## convert data structure
    p_joint <- try(do.call(cbind, p_joint))
    
    ## identify common picks
    n_stations <- try(rowSums(p_joint))
    
    ## append number of stations that picked signal
    p <- try(cbind(p, n_stations))
    
    ## isolate common picks
    p_joint <- try(p[n_stations >= pars$pick$n_coherent,])
    
    ## sort picks in ascending start order
    p_joint_asc <- try(p_joint[order(p_joint$start),])
    
    ## calculate differences between picks
    dt_asc <- try(c(100, diff(p_joint_asc$start)))
    
    ## isolate unique picks across stations
    p_joint_asc <- try(p_joint_asc[dt_asc > pars$pick$t_spacing,])
    
    ## remove picks outside time window of focus
    p_joint_asc <- try(p_joint_asc[p_joint_asc$start > t[i] &
                                     p_joint_asc$start < t[i] +
                                     pars$pick$t_snippet,])
    
    ## return output
    if(class(p_joint_asc) != "try-error") {
      
      picks_out <- rbind(picks_out, p_joint_asc)
    }
  }
  
  ## remove first dummy entry
  picks_out <- picks_out[-1,]
  
  ## remove NA cases
  picks_out <- picks_out[!is.na(picks_out$ID),]
  
  ## remove too long and too short events
  picks_out <- picks_out[picks_out$duration >= pars$pick$dur_min &
                           picks_out$duration <= pars$pick$dur_max,]
  
  ## print picked events
  print(picks_out)
  
  picks <- picks_out
  
  save(picks,
       pars,
       file = paste0(pars$path$picks,
                     pars$pick$pick_ID,
                     ".rda"))
} else {
  
  ## save pars
  pars_save <- pars
  
  ## load picked events
  load(file = paste0(pars$path$picks,
                     pars$pick$pick_ID,
                     ".rda"))
  
  ## restore pars
  pars <- pars_save
}

## PART 4 - LOCATE PICKS ------------------------------------------------------

if(pars$opt$locate == TRUE) {
  
  ## load spatial data sets
  try(dem <- raster(x = pars$path$dem), silent = TRUE)
  try(hs <- raster(x = pars$path$hs), silent = TRUE)
  try(aoi <- raster(x = pars$path$aoi), silent = TRUE)
  try(values(aoi) <- as.logical(values(aoi)), silent = TRUE)
  try(aoi_plot <- aoi, silent = TRUE)
  try(values(aoi_plot) <- ifelse(values(aoi_plot) == FALSE, 1, NA), 
      silent = TRUE)
  try(load(pars$path$d_map), silent = TRUE)
  try(load(pars$path$coupling), silent = TRUE)
  
  if(exists("D_02")) {D <- D_02}
  if(exists("D_10")) {D <- D_10}
  if(exists("D_20")) {D <- D_20}
  if(exists("D_30")) {D <- D_30}
  if(exists("D_50")) {D <- D_50}
  
  ## remove data sets of stations not used
  try(if(length(pars$opt$station_skip) > 0) {
    
    D$maps <- D$maps[-id_out]
    D$stations <- D$stations[-id_out, -id_out]
  }, silent = TRUE)
  try(coupling_save <- coupling)
}

## process all events
if(pars$opt$locate == TRUE | pars$opt$pick_plot == TRUE) {
  
  ## print progress info
  if(pars$manual$enable == TRUE) {
    
    print(paste0(Sys.time(),
                 ", Start processing picks (n = ",
                 "1",
                 ")"))
  } else {
    
    print(paste0(Sys.time(),
                 ", Start processing picks (n = ",
                 nrow(picks),
                 ")"))
  }
  
  ## process all picks
  if(pars$manual$enable == TRUE) {
    
    n_process <- 1
  } else {
    
    n_process <- nrow(picks)
  }
  
  for(i in pars$locate$i_start:n_process) {
    
    if(pars$manual$enable == TRUE) {
      
      ## print progress
      print(paste0("  process manually defined event"))
      
      ## read data sets
      s <- aux_getevent(start = pars$manual$start -
                          pars$locate$buffer[1],
                        duration = pars$manual$duration +
                          pars$locate$buffer[1] +
                          pars$locate$buffer[2],
                        station = stations$ID,
                        component = "BHZ",
                        dir = pars$path$sac,
                        try = TRUE)
    } else {
      
      ## print progress
      print(paste0("  ", i))
      
      ## read data sets
      s <- lapply(X = stations$ID, FUN = function(x) {
        
        s_out <- try(aux_getevent(start = picks$start[i] -
                                    pars$locate$buffer[1],
                                  duration = picks$duration[i] +
                                    pars$locate$buffer[1] +
                                    pars$locate$buffer[2],
                                  station = x,
                                  component = "BHZ",
                                  dir = pars$path$sac), 
                     silent = TRUE)
        if(class(s_out)[1] != "try-error") {
          
          if(1/s_out$meta$dt > 100) {
            s_out <- signal_aggregate(data = s_out, n = 2)
          }
          
          return(s_out)
        } else {
          
          return(NA)
        }
      })
    }
    
    ## provide list item names
    names(s) <- stations$ID
    
    ## identify successful read attempts
    s_ok <- do.call(c, lapply(X = s, FUN = function(s) {
      
      !is.na(s)[1]
    }))
    
    ## remove unsuccessful imports, correct other data sets, as well
    s <- try(s[s_ok])
    
    if(pars$opt$locate == TRUE) {
      
      D_ok <- try(D, silent = TRUE)
      D_ok$maps <- try(D_ok$maps[s_ok], silent = TRUE)
      D_ok$stations <- try(D_ok$stations[s_ok, s_ok], silent = TRUE)
      coupling <- try(coupling_save[s_ok], silent = TRUE)
    }
    
    ## deconvolve data set
    try(for(j in 1:length(s)) {
      
      s[[j]] <- try(signal_deconvolve(data = s[[j]],
                                      sensor = stations$sensor_type[j],
                                      logger = stations$logger_type[j],
                                      gain = stations$gain[j]),
                    silent = TRUE)
    })       
    
    ## demean data set
    s <- try(signal_demean(data = s),
             silent = TRUE)
    
    ## detrend data set
    s <- try(signal_detrend(data = s),
             silent = TRUE)
    
    ## filter data set
    s_f <- try(signal_filter(data = s,
                             f = pars$locate$f),
               silent = TRUE)
    
    ## detrend and taper data set
    s_f <- try(signal_detrend(data = s_f),
               silent = TRUE)
    
    s_f <- try(signal_taper(data = s_f,
                            n = pars$locate$n_taper))
    
    ## filter signals for plotting
    s_plot <- try(eseis::signal_filter(data = s,
                                       f = pars$plot$f))
    
    ## taper signals for plotting
    s_plot <- try(eseis::signal_taper(data = s_plot,
                                      n = pars$plot$p_taper))
    
    ## filter signals for plotting
    s_plot_2 <- try(eseis::signal_filter(data = s,
                                         f = pars$plot$f2))
    
    ## taper signals for plotting
    s_plot_2 <- try(eseis::signal_taper(data = s_plot_2,
                                        n = pars$plot$p_taper))
    
    ## remove buffer
    if(pars$manual$enable == TRUE) {
      
      s_f <- try(signal_clip(data = s_f,
                             limits = c(pars$manual$start -
                                          pars$locate$append[1],
                                        pars$manual$start +
                                          pars$manual$duration +
                                          pars$locate$append[2])))
    } else {
      
      s_f <- try(signal_clip(data = s_f,
                             limits = c(picks$start[i] -
                                          pars$locate$append[1],
                                        picks$start[i] +
                                          picks$duration[i] +
                                          pars$locate$append[2])))
    }
    
    if(pars$opt$locate == TRUE) {
      
      ## calculate signal envelope
      s_e <- try(signal_envelope(data = s_f),
                 silent = TRUE)
      
      ## calcuate snr
      snr <- try(signal_snr(data = s_e))
      snr <- try(do.call(c, lapply(X = snr, FUN = function(snr) {snr$snr})))
      
      ## remove further data based on SNR criteria
      if(pars$manual$enable == TRUE) {
        
        if(pars$manual$snr_reject == TRUE) {
          
          snr_ok <- try(snr > pars$locate$snr_min)
          
          try(if(mean(snr) < pars$locate$snr_mean_min) {
            
            snr_ok <- rep(FALSE, n = length(snr_ok))
          })
          
          try(if(sum(snr_ok) < pars$locate$station_min) {
            
            snr_ok <- rep(FALSE, n = length(snr_ok))
          })
        } else {
          
          snr_ok <- try(rep(TRUE, length(snr)))
        }
      } else {
        
        snr_ok <- try(snr > pars$locate$snr_min)
        
        try(if(mean(snr) < pars$locate$snr_mean_min) {
          
          snr_ok <- rep(FALSE, n = length(snr_ok))
        })
        
        try(if(sum(snr_ok) < pars$locate$station_min) {
          
          snr_ok <- rep(FALSE, n = length(snr_ok))
        })
      }
      
      ## remove unsuitable data
      s <- try(s[snr_ok])
      s_e <- try(s_e[snr_ok])
      D_ok$maps <- try(D_ok$maps[snr_ok])
      D_ok$stations <- try(D_ok$stations[snr_ok, snr_ok])
      coupling_ok <- try(coupling[snr_ok])
      
      ## locate seismic source
      if(pars$locate$method == "migrate") {
        
        l <- try(spatial_migrate(data = s_e,
                                 d_stations = D_ok$stations,
                                 d_map = D_ok$maps,
                                 v = pars$locate$v))
      }
      
      if(pars$locate$method == "amplitude") {
        
        l <- try(spatial_amplitude(data = s_e,
                                   coupling = coupling_ok,
                                   d_map = D_ok$maps,
                                   aoi = aoi,
                                   q = pars$locate$q,
                                   v = pars$locate$v,
                                   f = mean(pars$locate$f), 
                                   normalise = FALSE))
        print(max(values(l), na.rm = TRUE))
      }
      
      ## reduce to values above threshold
      l_p <- try(spatial_clip(data = l,
                              quantile = pars$locate$quantile))
      
      ## get P_max or pixels above a given threshold
      xy_max <- try(coordinates(l)[values(l) == max(values(l),
                                                    na.rm = TRUE)[1],])
    } else {
      
      l <- try(qualle/water, silent = TRUE)
    }
    
    ## optionally plot signals ------------------------------------------------
    if(pars$opt$pick_plot == TRUE) {
      
      if(pars$manual$enable == TRUE) {
        
        try(jpeg(filename = paste0(pars$path$plot,
                                   "event_signals_manual_",
                                   pars$manual$start,
                                   ".jpg"),
                 width = pars$plot$width,
                 height = pars$plot$height,
                 res = pars$plot$res))
      } else {
        
        try(jpeg(filename = paste0(pars$path$plot,
                                   "event_signals_",
                                   picks$start[i],
                                   ".jpg"),
                 width = pars$plot$width,
                 height = pars$plot$height,
                 res = pars$plot$res))
      }
      
      ## set up plot area
      try(par(mar = c(2, 3, 3, 1), mfcol = c(length(s_plot), 1)))
      
      ## find global signal range
      s_range <- try(lapply(X = s_plot, FUN = function(x) {
        
        range(x$signal, na.rm = TRUE)
      }))
      s_range <- try(range(do.call(c, s_range), na.rm = TRUE))
      
      ## assign station IDs to signal list
      names(s_ok) <- stations$ID
      
      for(j in 1:length(s_ok)) {
        
        if(s_ok[j] == TRUE) {
          
          i_plot <- (1:length(s_plot))[names(s_plot) ==
                                         names(s_ok)[j]]
          
          try(plot_signal(data = s_plot[[i_plot]],
                          main = paste(names(s_ok)[j], " | ",
                                       picks$start[i], " + ",
                                       round(picks$duration[i], 1), 
                                       " s (grey = filter_1 & global, ",
                                       " blue = filter_2 & individual)"),
                          col = "grey30",
                          lwd = 1,
                          ylim = s_range))
          
          try(par(new = TRUE))
          
          try(plot_signal(data = s_plot_2[[i_plot]],
                          ann = FALSE,
                          axes = FALSE,
                          col = adjustcolor(col = "darkblue",
                                            alpha.f = 0.7),
                          lwd = 0.8))
          
          try(axis(side = 4, col = "darkblue", col.ticks = "darkblue"))
        } else {
          
          plot(0, main = names(s_ok)[j])
        }
      }
      
      ## optionally close plto device
      if(pars$opt$pick_plot == TRUE) {
        
        try(dev.off())
      }
    }
    
    ## optionally plot signals and location estimate --------------------------
    if(pars$opt$pick_locate_plot == TRUE & class(l) != "try-error") {
      
      if(pars$opt$pick_locate_plot == TRUE) {
        
        if(pars$manual$enable == TRUE) {
          
          try(jpeg(filename = paste0(pars$path$plot,
                                     "event_location_manual_",
                                     pars$manual$start,
                                     ".jpg"),
                   width = pars$plot$width,
                   height = pars$plot$height,
                   res = pars$plot$res))
        } else {
          
          try(jpeg(filename = paste0(pars$path$plot,
                                     "event_location_",
                                     picks$start[i],
                                     ".jpg"),
                   width = pars$plot$width,
                   height = pars$plot$height,
                   res = pars$plot$res))
        }
      }
      
      try(layout(mat = cbind(1:length(s_ok),
                             rep(length(s_ok) + 1,
                                 length(s_ok)))))
      
      par(mar = c(2, 3, 3, 1))
      
      s_range <- try(lapply(X = s_plot, FUN = function(x) {
        
        range(x$signal, na.rm = TRUE)
      }))
      
      s_range <- try(range(do.call(c, s_range), na.rm = TRUE))
      
      names(s_ok) <- stations$ID
      
      for(j in 1:length(s_ok)) {
        
        if(s_ok[j] == TRUE) {
          
          i_plot <- (1:length(s_plot))[names(s_plot) ==
                                         names(s_ok)[j]]
          
          try(plot_signal(data = s_plot[[i_plot]],
                          main = paste(names(s_ok)[j], " | ",
                                       "SNR: ",
                                       paste(round(snr, 1),
                                             collapse = " - ")),
                          col = "grey30",
                          lwd = 1,
                          ylim = s_range))
          
          try(par(new = TRUE))
          
          try(plot_signal(data = s_plot_2[[i_plot]],
                          ann = FALSE,
                          axes = FALSE,
                          col = adjustcolor(col = "darkblue",
                                            alpha.f = 0.7),
                          lwd = 0.8))
          
          try(axis(side = 4, col = "darkblue", col.ticks = "darkblue"))
        } else {
          
          plot(0, main = names(s_ok)[j])
        }
      }
      
      if(pars$manual$enable == TRUE) {
        
        plot_title <- paste("Event manual", " | ",
                            pars$manual$start, " + ",
                            pars$manual$duration,
                            " s | ",
                            round(max(values(l), na.rm = TRUE), 3),
                            sep = "")
      } else {
        
        plot_title <- paste("Event ", i, " | ",
                            picks$start[i], " + ",
                            round(picks$duration[i], 1),
                            " s | ",
                            round(max(values(l), na.rm = TRUE), 3),
                            sep = "")
      }
      
      try(image(hs,
                col = grey.colors(200),
                main = plot_title))
      
      par(new = TRUE)
      
      try(image(aoi_plot,
                col = adjustcolor(1, 0.4),
                ann = FALSE,
                axes = FALSE))
      
      par(new = TRUE)
      
      try(image(l_p,
                col = adjustcolor(col = pars$plot$col_map(200),
                                  alpha.f = 0.5),
                ann = FALSE,
                axes = FALSE))
      
      try(points(x = xy_max[1],
                 y = xy_max[2],
                 col = "white"))
      
      try(points(x = stations$x,
                 y = stations$y,
                 col = "black",
                 pch = 10,
                 cex = 3))
      
      try(text(x = stations$x,
               y = stations$y,
               col = "red",
               labels = stations$ID,
               adj = c(0, 1),
               cex = 2))
      
      if(pars$opt$pick_locate_plot == TRUE) {
        
        try(dev.off())
      }
    }
  }
}
