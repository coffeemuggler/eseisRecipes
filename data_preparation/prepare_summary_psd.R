## PART 0 - SCRIPT INFORMATION SECTION ----------------------------------------

#' The script can be used to generate spectrograms of long time periods, 
#' possibly several years.
#' 
#' Author: Michael Dietze, Uni Goettingen (michael.dietze)
#' 
#' Script version: 0.4.0
#' 
#' Changes to previous versions
#'   - not documented
#' 
#' Dependencies: 
#'   - R 4.n
#'   - eseis 0.7.1
#'   
#'   Parameter descriptions:
#'   - opt, general options
#'     - ID
#'     

## PART 1 - SETTINGS SECTION --------------------------------------------------

pars <- list(
  
  opt = list(ID = c("STA1", "STA2"),
             cmp = c("BHZ", "BHZ"),
             sensor = c("PE6B", "PE6B"),
             logger = c("Cube3ext", "Cube3ext"),
             gain = c(32, 32),
             
             t_start = c("2022-07-20"),
             t_stop = c("2022-08-16"),
             h_ok = 0:23,
             duration = 1800,
             zlim = c(-160, -100),
             cpu = 5 / 11,
             
             n_try = 200,
             length_f = 1000,
             p_taper = 0.01,
             welch = TRUE,
             window_sub = 300,
             overlap_sub = 0.7,
             
             save_data = FALSE,
             save_jpeg = TRUE, 
             clear_data = TRUE,
             
             jpg_width = 6000, 
             jpg_height = 2000,
             jpg_res = 300,
             col_scale = c("jet", "inferno", "plasma", "viridis")[2],
             
             
             packages = "eseis", "parallel", "colorspace"),
  
  path = list(sac = paste0("/data/DATA/sac/"),
              save = paste0("/data/DATA/rda/"),
              jpeg = paste0("/data/DATA/plots/")))

## PART 2 - PROCESSING OF DATA ------------------------------------------------

## prepare environment
x <- lapply(pars$opt$packages, require, character.only = TRUE)
rm(x)
Sys.setenv(TZ = "UTC")


## process all data sets
for (i in 1:length(pars$opt$ID)) {
  
  ## print progress info
  print(paste0("PROCESSING STATION-CMP ", i))
  
  ## create time slice vector
  t <- seq(from = as.POSIXct(x = pars$opt$t_start, tz = "UTC"), 
           to = as.POSIXct(x = pars$opt$t_stop, tz = "UTC") - 0.1, 
           by = pars$opt$duration)
  
  ## extract hours of the dates
  t_hour <- as.numeric(format(x = t, 
                              format = "%H", 
                              tz = "UTC"))
  
  ## keep only hours specified in parameters
  t <- t[!is.na(match(x = t_hour, table = pars$opt$h_ok))]
  
  ## calculate test PSD to get frequency vector -------------------------------
  
  ## read dummy file
  j <- 1
  while(j <= pars$opt$n_try) {
    
    print(paste("  Frequency estimation, attempt", j))
    
    s <- try(eseis::aux_getevent(start = t[runif(n = 1,
                                                 min = 1,
                                                 max = length(t))],
                                 duration = pars$opt$duration,
                                 station = pars$opt$ID[i],
                                 component = pars$opt$cmp[i],
                                 dir = pars$path$sac),
             silent = TRUE)
    
    if(class(s) == "try-error") {
      
      j <- j + 1
    } else {
      
      j <- pars$opt$n_try + 1
    }
  }
  
  ## calculate dummy spectrogram
  p <- try(eseis::signal_spectrogram(data = s, 
                                     Welch = pars$opt$welch,
                                     window = pars$opt$duration - 1,
                                     overlap = 0, 
                                     window_sub = pars$opt$window_sub,
                                     overlap_sub = pars$opt$overlap_sub))
  
  ## extract frequency vector
  f <- try(p$PSD$f)
  
  ## append frequency vector length to parameter list
  pars$length_f <- length(f)
  
  ## create aggregation index vector
  i_agg <- floor(seq(from = 1, 
                     to = length(f), 
                     length.out = pars$length_f))
  
  ## initiate cluster
  cores <- parallel::detectCores()
  n_cpu <- floor(cores * pars$opt$cpu)
  cores <- ifelse(cores < n_cpu, cores, n_cpu)
  cl <- parallel::makeCluster(getOption("mc.cores", cores))
  
  ## print information
  print(paste("  Preparation finished.",
              length(t), 
              "PSDs to calculate in total, using", 
              cores, 
              "CPUs. Starting job NOW."))
  
  ## PART 3 - PROCESSING THE DATA -----------------------------------------------
  
  ## get time before processing
  t_sys_0 <- Sys.time()
  
  ## process all time slices
  psd_list <- parallel::parLapply(cl = cl, X = t, fun = function(t, pars, i) {
    
    ## read file
    s <- try(eseis::aux_getevent(start = t, 
                                 duration = pars$opt$duration, 
                                 station = pars$opt$ID[i], 
                                 component = pars$opt$cmp[i], 
                                 dir = pars$path$sac))
    
    ## deconvolve data
    s <- try(eseis::signal_deconvolve(data = s, 
                                      sensor = pars$opt$sensor[i], 
                                      logger = pars$opt$logger[i], 
                                      gain = pars$opt$gain[i]))
    
    ## demean and detrend data
    s <- try(eseis::signal_demean(data = s))
    s <- try(eseis::signal_detrend(data = s))
    
    ## taper signal
    s <- try(eseis::signal_taper(data = s, p = pars$opt$p_taper))
    
    ## calculate spectrogram
    p <- try(eseis::signal_spectrogram(data = s, 
                                       Welch = pars$opt$welch,
                                       window = pars$opt$duration - 1,
                                       overlap = 0, 
                                       window_sub = pars$opt$window_sub,
                                       overlap_sub = pars$opt$overlap_sub))
    
    ## extract spectrogram matrix
    p <- try(p$PSD$S[,1]) 
    
    ## return output
    if(class(p) == "try-error") {
      
      return(rep(x = NA, times = pars$length_f))
    } else {
      
      return(p)
    }
    
  }, pars, i)
  
  ## stop cluster
  try(parallel::stopCluster(cl = cl))
  
  ## get time after processing
  t_sys_1 <- Sys.time()
  
  ## get time difference
  dt <- as.numeric(difftime(time1 = t_sys_1, 
                            time2 = t_sys_0, 
                            units = "mins"))
  
  ## get number of successfull calculations
  p_ok <- try(do.call(c, lapply(X = psd_list, FUN = function(psd_list) {
    
    sum(is.na(psd_list)) == 0
  })))
  
  ## print processing information
  try(print(paste("  Processing finished. Time: ",
                  round(x = dt, digits = 1),
                  " min. ",
                  sum(p_ok),
                  " successfull calculations (",
                  round(x = sum(p_ok) / length(p_ok) * 100, digits = 1),
                  " %).", 
                  sep = "")))
  
  ## convert list to matrix
  psd <- try(do.call(rbind, psd_list))
  
  ## aggregate data
  f_agg <- f[i_agg]
  psd_agg <- psd[, i_agg]
  
  ## scale psd for plotting
  psd_plot <- psd_agg
  psd_plot[psd_plot < pars$opt$zlim[1]] <- pars$opt$zlim[1]
  psd_plot[psd_plot > pars$opt$zlim[2]] <- pars$opt$zlim[2]
  
  
  ## PART 4 - SAVING AND PLOTTING ---------------------------------------------
  
  ## save the raw data
  if(pars$opt$save_data == TRUE) {
    
    save(pars, psd, t, f, f_agg, psd_agg, 
         file = paste0(pars$path$save, "PSD_", 
                       pars$opt$ID[i], "_", 
                       pars$opt$cmp[i], ".rda"))
  }
  
  ## optionally, plot data set to jpg file
  if(pars$opt$save_jpeg == TRUE) {
    
    jpeg(filename = paste0(pars$path$save, "PSD_", 
                           pars$opt$ID[i], "_", 
                           pars$opt$cmp[i], ".jpg"),
         width = pars$opt$jpg_width, 
         height = pars$opt$jpg_height,
         res = pars$opt$jpg_res)
  }
  
  if(pars$opt$col_scale == "jet") {
    
    try(fields::image.plot(x = t, 
                           y = f_agg,
                           z = psd_plot, 
                           zlim = pars$opt$zlim))
  } else if(pars$opt$col_scale == "inferno") {
    
    try(fields::image.plot(x = t, 
                           y = f_agg,
                           z = psd_plot, 
                           zlim = pars$opt$zlim, 
                           col = colorspace::sequential_hcl(
                             200, palette = "Inferno")))
  } else if(pars$opt$col_scale == "plasma") {
    
    try(fields::image.plot(x = t, 
                           y = f_agg,
                           z = psd_plot, 
                           zlim = pars$opt$zlim, 
                           col = colorspace::sequential_hcl(
                             200, palette = "Plasma")))
  } else if(pars$opt$col_scale == "viridis") {
    
    try(fields::image.plot(x = t, 
                           y = f_agg,
                           z = psd_plot, 
                           zlim = pars$opt$zlim, 
                           col = colorspace::sequential_hcl(
                             200, palette = "Viridis")))
  }
  
  ## close plot device
  if(pars$opt$save_jpeg == TRUE) {
    dev.off()
  }
  
  ## optionally, remove data
  if(pars$opt$clear_data == TRUE) {
    
    rm(psd, t, f, f_agg, psd_agg)
  }
}

## optionally, remove data
if(pars$opt$clear_data == TRUE) {
  
  rm(list = ls())
}
