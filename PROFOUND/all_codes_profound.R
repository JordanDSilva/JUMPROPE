library(Rfits)
library(Rwcs)
library(ProFound)
library(ProPane)
library(celestial)
library(magicaxis)
library(foreach)
library(doParallel)
library(data.table)
library(Cairo)
library(stringr)
library(checkmate)
library(dplyr)
library(Highlander)
library(matrixStats)

source("./ProFound_settings.R")

jumprope_version = "1.3.5"

frame_info = function(ref_dir){
  
  message("Running frame_info")
  
  patch_stack_dir = paste0(ref_dir, "/Patch_Stacks/")
  
  mosaic_patch_dir = paste0(ref_dir, "/Mosaic_Stacks/Patch/")
  if(dir.exists(mosaic_patch_dir)){
    message("Found aligned mosaics in ", mosaic_patch_dir)
  }

  filelist = c(
    list.files(
      path = patch_stack_dir,
      pattern = glob2rx("*.fits"),
      full.names = T
    ),
    list.files(
      path = mosaic_patch_dir,
      pattern = glob2rx("*.fits"), 
      full.names = T)
  )
  
  filenames = c(
    list.files(
      path = patch_stack_dir,
      pattern = glob2rx("*.fits"),
      full.names = F
    ),
    list.files(
      path = mosaic_patch_dir,
      pattern = glob2rx("*.fits"), 
      full.names = F)
  )
  
  VID_list = unname(sapply(filenames, function(x){
    split_fname = str_split(x, "_")[[1]]
    split_fname = split_fname[!(grepl("patch|stack|short|long|F.+W|F.+M|F.+N|MIRIMAGE|.fits", split_fname))]
    vid = split_fname[1]
    if(is.na(vid)){
      return(split_fname[1])
    }else{
      return(vid)
    }
  }))
  
  MODULE_list = unname(sapply(filenames, function(x){
    split_fname = str_split(x, "_")[[1]]
    split_fname = split_fname[!(grepl("patch|stack|short|long|F.+W|F.+M|F.+N|MIRIMAGE|.fits", split_fname))]
    if(is.na(split_fname[2])){
      module = split_fname[1]
    }else{
      module = split_fname[2]
    }
    return(module)
  }))
  
  PIX_list = unname(sapply(filenames, function(x){
    split_fname = tail(str_split(x, "_")[[1]], 1)
    pixscale_ = str_split_1(split_fname, ".fits")[1]
    return(pixscale_)
  }))

  frame_info = data.frame(
    "filenames" = filenames,
    "VISIT_ID" = VID_list,
    "MODULE" = MODULE_list,
    "PIXSCALE" = PIX_list
  )
  
  stack_grid = unique(frame_info[,c("VISIT_ID", "MODULE", "PIXSCALE")])
  
  stack_grid$PROPOSAL_ID = substr(stack_grid$VISIT_ID, 1, 4)
  stack_grid$PROPOSAL_ID[stack_grid$MODULE == stack_grid$VISIT_ID] = stack_grid$VISIT_ID[stack_grid$MODULE == stack_grid$VISIT_ID]
  
  foo = foreach(i = 1:dim(stack_grid)[1], .combine = "rbind") %do% {
    wcs_info = propaneFrameFinder(
      rad = 2*pi*3600,
      filelist = filelist[
        grepl(stack_grid$VISIT_ID[i], filelist) & grepl(stack_grid$MODULE[i], filelist)
      ][1],
      plot = F
    )
    
    mid_info = data.frame(
      cbind("RA_MBL" = (wcs_info$centre_RA + wcs_info$corner_BL_RA)/2.0, "DEC_MBL" = (wcs_info$centre_Dec + wcs_info$corner_BL_Dec)/2.0),
      cbind("RA_MTL" = (wcs_info$centre_RA + wcs_info$corner_TL_RA)/2.0, "DEC_MTL" = (wcs_info$centre_Dec + wcs_info$corner_TL_Dec)/2.0),
      cbind("RA_MTR" = (wcs_info$centre_RA + wcs_info$corner_BR_RA)/2.0, "DEC_MTR" = (wcs_info$centre_Dec + wcs_info$corner_BR_Dec)/2.0),
      cbind("RA_MBR" = (wcs_info$centre_RA + wcs_info$corner_TR_RA)/2.0, "DEC_MBR" = (wcs_info$centre_Dec + wcs_info$corner_TR_Dec)/2.0)
    )
    
    mid_corner_info = data.frame(
      cbind( "RA_CL" = (wcs_info$corner_BL_RA + wcs_info$corner_TL_RA)/2.0,  "DEC_CL" = (wcs_info$corner_BL_Dec + wcs_info$corner_TL_Dec)/2.0),
      cbind( "RA_CT" = (wcs_info$corner_TL_RA + wcs_info$corner_TR_RA)/2.0,  "DEC_CT" = (wcs_info$corner_TL_Dec + wcs_info$corner_TR_Dec)/2.0),
      cbind( "RA_CR" = (wcs_info$corner_TR_RA + wcs_info$corner_BR_RA)/2.0,  "DEC_CR" = (wcs_info$corner_TR_Dec + wcs_info$corner_BR_Dec)/2.0),
      cbind( "RA_CB" = (wcs_info$corner_BL_RA + wcs_info$corner_BR_RA)/2.0,  "DEC_CB" = (wcs_info$corner_BL_Dec + wcs_info$corner_BR_Dec)/2.0)
    )
    
    wcs_info_trim = wcs_info[, c("centre_RA", "centre_Dec", 
                                 "corner_BL_RA", "corner_BL_Dec",
                                 "corner_TL_RA", "corner_TL_Dec",
                                 "corner_TR_RA", "corner_TR_Dec",
                                 "corner_BR_RA", "corner_BR_Dec")]
    names(wcs_info_trim) = c("RA_CEN", "DEC_CEN", 
                             "RA_BL", "DEC_BL", 
                             "RA_TL", "DEC_TL",
                             "RA_TR", "DEC_TR",
                             "RA_BR", "DEC_BR")
    
    ret = cbind(
      stack_grid[i, c("PROPOSAL_ID","VISIT_ID", "MODULE", "PIXSCALE")],
      
      mid_info,
      
      wcs_info_trim,
      
      mid_corner_info
    )
    return(ret)
  }
  
  foo[["jumprope_version"]] = rep(jumprope_version, dim(foo)[1])
  csv_stub = paste0(ref_dir, "/ProFound/warp_info.csv")
  fwrite(foo, csv_stub)
} ##<--Compute the long warp frame info needed for querying GAIA and HST via MAST

## Moving files for source detection codes
warp_short_to_long = function(input_args){
  
  ## Safest to avoid covariance between pixels
  
  message("Running warp_short_to_long")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  
  if(!(grepl("NRC", MODULE, fixed = T))){
    MODULE = paste0("NRC", MODULE)
  }
  ## Going to warp the short to the long 
  message("Start warp!")
  
  patch_stack_dir = paste0(ref_dir, "/Patch_Stacks/")
  warp_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/")
  
  if(dir.exists(warp_dir)){
    warp_propane_files = list.files(
      path = warp_dir,
      pattern = ".fits$",
      full.names = T
    )
    warp_propane_files = warp_propane_files[
      !grepl("hst", warp_propane_files, ignore.case = T)
    ]
    file.remove(
      warp_propane_files
    )
  }else{
    dir.create(warp_dir, recursive = T)
  }
  
  ## here we read in the Patched_Stacks
  filepaths <- list.files(patch_stack_dir, pattern=".fits", full.names=TRUE) #list all the .fits files in a directory
  fitsnames <- list.files(patch_stack_dir, pattern=".fits", full.names=FALSE) #list all the .fits files in a directory
  
  #native scales
  short_idx = grepl(VID, filepaths) & grepl(MODULE, filepaths)  & grepl("short", filepaths) & grepl("F070W|F090W|F115W|F150W|F200W|F140M|F162M|F182M|F210M|F164N|F187N|F212N", filepaths) 
  long_idx = grepl(VID, filepaths) & grepl(MODULE, filepaths)  & grepl("long", filepaths) & grepl("F277W|F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430|F460M|F480M|F323N|F405N|F466N|F470N", filepaths)
  
  message(paste0(fitsnames[short_idx], collapse = "\n"))
  message(paste0(fitsnames[long_idx], collapse = "\n"))
  
  frames_short = lapply(filepaths[short_idx], function(x) Rfits_read_all(filename = x, pointer = F)) #to be warped
  names(frames_short) = fitsnames[short_idx]
  frames_long = lapply(filepaths[long_idx], function(x) Rfits_read_all(filename = x, pointer = F)) #reference warp
  names(frames_long) = fitsnames[long_idx]
  
  frames_short_to_long = lapply(frames_short, function(x){
    warp = propaneWarpProPane(x, 
                              keyvalues_out = frames_long[[1]]$image$keyvalues, 
                              magzero_out = 23.9, magzero_in = 23.9)
    warp$info = x$info
    
    warp$image$imDat[is.infinite(warp$image$imDat)] = NA
    
    warp
  })
  
  
  names(frames_short_to_long) = str_replace(names(frames_short), "short", "long")
  for(i in 1:length(frames_short_to_long)){
    class(frames_short_to_long[[i]]) = "Rfits_list"
  }
  message(paste0("Finished warping: ", names(frames_short_to_long), sep = "\n"))
  N_frames_short = length(frames_short)
  
  for(i in 1:length(frames_short_to_long)){
    message(paste0("Saving ", names(frames_short_to_long)[i], " to ", warp_dir))
    Rfits_write(data = frames_short_to_long[[i]], filename = paste0(warp_dir, "warp_", names(frames_short_to_long)[i]))
  }
  for(i in 1:length(frames_long)){
    message(paste0("Copying ", names(frames_long)[i], " to ", warp_dir))
    file.copy(from = filepaths[long_idx][i], to = paste0(warp_dir, "warp_", names(frames_long)[i]), overwrite = T)
  }
  message("Done!")
  # return(c(frames_short_to_long, frames_long))
} ##<-- Downwarp the short wavelength, short pixel scale to long 
copy_frames = function(input_args){
  
  ## Copy the long pixel scale from the PROCESS dirs to the DATA dir
  ## Keep separete DATA dir incase the user puts their own frames in there
  ## Will be data redundancy but I'm assuming you will have enough disk space if you're using this code!
  
  message("Running copy")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  message("Start copy!")
  
  patch_stack_dir = paste0(ref_dir, "/Patch_Stacks/")
  mosaic_patch_dir = paste0(ref_dir, "/Mosaic_Stacks/Patch/")
  warp_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/")
  
  if(dir.exists(warp_dir)){
    warp_propane_files = list.files(
      path = warp_dir,
      pattern = ".fits$",
      full.names = T
    )
    warp_propane_files = warp_propane_files[
      !grepl("hst", warp_propane_files, ignore.case = T)
    ]
    file.remove(
      warp_propane_files
    )
  }else{
    dir.create(warp_dir, recursive = T)
  }
  
  ## here we read in the Patched_Stacks
  filepaths <- c(
    list.files(patch_stack_dir, pattern=".fits", full.names=TRUE), #list all the .fits files in a directory
    list.files(mosaic_patch_dir, pattern=".fits", full.names = TRUE)
  )
  filepaths = filepaths[!is.na(filepaths)]
  fitsnames <- c(
    list.files(patch_stack_dir, pattern=".fits", full.names=FALSE), #list all the .fits files in a directory
    list.files(mosaic_patch_dir, pattern=".fits", full.names = FALSE)
  )
  fitsnames = fitsnames[!is.na(fitsnames)]
  
  #native scales
  long_idx = grepl(VID, filepaths) & grepl(MODULE, filepaths)  & grepl(PIXSCALE, filepaths)
  
  #frames_long = lapply(filepaths[long_idx], function(x) Rfits_read_all(filename = x, pointer = F))
  names_frames_long = fitsnames[long_idx]
  
  for(i in 1:length(names_frames_long)){
    message(paste0("Copying ", names_frames_long[i], " to ", warp_dir))
    file.copy(from = filepaths[long_idx][i], to = paste0(warp_dir, "warp_", names_frames_long[i]), overwrite = T)
  }
  message("Done!")
  # return(c(frames_short_to_long, frames_long))
} ##<-- Copy the long pixelscale to Data dir. File redundancy but safer option for detects.

## Star mask codes
star_mask = function(input_args){
  
  ## Why did they have to make the PSF so darn complicated :(
  
  ## Testing purposes - FML
  # input_args = list(
  #   ref_dir = "/Volumes/BRAIDY/JWST/JP",
  #   VID = 1180018001,
  #   MODULE = "NRCB",
  #   PIXSCALE = "short"
  # )

  message("Running star mask")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  star_mask_dir = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/")
  if(dir.exists(star_mask_dir)){
    unlink(star_mask_dir, recursive = T)
  }
  
  dir.create(star_mask_dir, showWarnings = F, recursive = T)
  data_dir = paste0(ref_dir, "/Patch_Stacks/") #directory with the propanes
  gaia_dir = paste0(ref_dir, "/ProFound/GAIA_Cats/", VID, "/")
  file_list = list.files(data_dir, pattern=glob2rx(paste0("*", VID, "*", MODULE, "*", PIXSCALE, "*.fits")), full.names=TRUE) #list all the .fits files in a directory
  file_names = list.files(data_dir, pattern=glob2rx(paste0("*", VID, "*", MODULE, "*", PIXSCALE, "*.fits")), full.names=FALSE)
  
  #make star masks on the F277W filter, and if it's missing use the shortest filter
  file_idx = (
    grepl("F277W", file_list) & 
    grepl(VID, file_list) & 
    grepl(MODULE, file_list) &
    !grepl("hst", file_list)
    )
  
  if(sum(file_idx) == 0){
    file_list = file_list[grepl(VID, file_list) & 
                            grepl(MODULE, file_list) & 
                            !grepl("hst", file_list)][1]
    file_names = file_names[grepl(VID, file_names) & 
                              grepl(MODULE, file_names) & 
                              !grepl("hst", file_names)][1]
  }else{
    file_list = file_list[file_idx]
    file_names = file_names[file_idx]
  }
  
  message(paste("Using", file_names, "to build star mask"))
  
  ref = Rfits_read(filename = file_list, pointer = FALSE)
  imdim = dim(ref$image)
  im_ra_dec = Rwcs_p2s(x = imdim[1], y = imdim[2], keyvalues = ref$image$keyvalues) #get extent of frame in RA and DEC coords
  #read in GAIA catalogue
  gaia_files = list.files(path = gaia_dir, pattern = glob2rx(paste0("*", VID, "*", MODULE, "*.csv")), full.names = T, recursive = T) 
  gaia = data.frame(
    fread(
      file = gaia_files
    )
  )
  delta_time = as.double(substr(ref$image$keyvalues$DATE_END, 1, 4)) - gaia$ref_epoch
  gaia[, c("xpix", "ypix")] = Rwcs_s2p(RA = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time, ## Do proper motion correction
                                       Dec = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time, ## pmra[pmdec] = mas/yr
                                       keyvalues = ref$image$keyvalues)
  gaia[, c("xpix_nopm", "ypix_nopm")] = Rwcs_s2p(RA = gaia$ra, 
                                                 Dec = gaia$dec, 
                                                 keyvalues = ref$image$keyvalues)
  gaia$ra_fix = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time
  gaia$dec_fix = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time
  gaia_idx = gaia$xpix >= 0 & gaia$xpix <= imdim[1] & gaia$ypix >= 0 & gaia$ypix <= imdim[2] & gaia$classprob_dsc_combmod_star > 0.98 & gaia$phot_bp_rp_excess_factor < 2.5
  gaia_idx[is.na(gaia_idx)] = FALSE
  gaia_ra_dec = gaia[gaia_idx,]
  gaia_trim = gaia_ra_dec[order(gaia_ra_dec$phot_g_mean_mag), ]

  message(paste0("Building star mask for VID: ", VID, ", MODULE: ", MODULE))
  psf_mask_function = function(image, xcen=dim(image)[1]/2, ycen=dim(image)[2]/2, rad=dim(image)[2]/2){
    mask =
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad/6.0, axrat=1) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad/2, axrat=0.03, ang=90) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad, axrat=0.03, ang=60) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad, axrat=0.03, ang=-60) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad, axrat=0.03, ang=0)
    return(mask > 0)
  } ## Define the PSF mask
  psf_mask = psf_mask_function(image = ref$image$imDat)

  core_sersic <- function(R, p) {
    ## Pretty much a two component Sersic profile to capture breaks in the 1D PSF
    ## https://arxiv.org/abs/astro-ph/0306023
    ## https://arxiv.org/abs/astro-ph/0403659
    
    # p = c(0.0, 1.0, 1.0, 4.0, 0.0) ## For testing purposes
    
    Ib = 10^p[1]
    Rb = 10^p[2]
    Re = 10^p[3]
    n = p[4]
    gamma_ = p[5]
    alpha = 4.0
    
    bn = 1.9992*n - 0.3271
    
    # Iprime = Ib * 2^(-1*gamma_/alpha) * exp(bn * (2^(1/alpha) * Rb/Re)^(1/n))
    # term1 = (1 + (Rb/R)^(alpha))^(gamma_/alpha)
    # term2 = exp(-1*bn * ( (R^alpha + Rb^alpha)/Re^alpha )^(1/(alpha*n)))
    # 
    # tot_ = log10( Iprime * term1 * term2 )
    
    ## Compute log directly to avoid infinities in exponents later on
    Iprime_log = log10(Ib) + log10(2)*(-1*gamma_/alpha) + log10(exp(1))*(bn * (2^(1/alpha) * Rb/Re)^(1/n))
    term1_log = log10(1 + (Rb/R)^(alpha))*(gamma_/alpha)
    term2_log = log10(exp(1))*(-1*bn * ( (R^alpha + Rb^alpha)/Re^alpha )^(1/(alpha*n)))
    
    tot_log = ( Iprime_log + term1_log + term2_log ) ## Confirmed that this is equal to log calculation above
    return(
      tot_log
    )
  }

  surf_fun = function(x, p){
    core_sersic(x, p)
  }
  
  LL = function(p, Data){ 
    ## Sersic LogLikelihood
    loglike = sum(
      dnorm(
        x = surf_fun(Data$rr, p),
        mean = log10( Data$ff ),
        sd = Data$ff_err / (log(10) * Data$ff), 
        log = T
      ) * (1/(Data$rr)), 
      na.rm = T)
    return(loglike)
  }
  
  ## meat of the function
  star_masker = function(ref, gaia, pro, psf){
    box = ceiling(
      0.3 / pixscale(ref, unit = "asec")
    ) ## Use 0.3 arcsec aperture to calculate ray 

    segim_map = Rfits_create_image(data = pro$segim, keyvalues = ref$keyvalues)
    sky_map = Rfits_create_image(data = pro$sky, keyvalues = ref$keyvalues)
    skyRMS_map = Rfits_create_image(data = pro$skyRMS, keyvalues = ref$keyvalues)
    
    ## Find local background levels to truncate diffraction spike
    depth = 1.0 * abs( diff(quantile( ref$imDat[pro$objects_redo == 0], probs = c(0.50,0.16), na.rm = TRUE)) * 2 )
    
    star_rad_list = foreach(kk = 1:dim(gaia)[1], .combine = "c") %do% {
    # star_rad_list = foreach(kk = 1, .combine = "c") %do% {
      
      test_coords = c(gaia$xpix[kk], gaia$ypix[kk])
      
      dr_vec = 10^seq(0, log10(dim(ref)[1]), length.out = 100)
      flux_mat = matrix(0L, nrow = 6, ncol = length(dr_vec))
      ## Calculate flux profiles along 6 radial rays
      for(dr in 1:length(dr_vec)){
        count = 1
        for(theta in c(pi/6.0, pi/2.0, 5*pi/6)){
          for(sgn in c(1,-1)){
            new_coords = test_coords + (sgn * c(dr_vec[dr] * cos(theta), dr_vec[dr] * sin(theta)) )
            flux = mean(ref[new_coords[1], new_coords[2], box = box]$imDat, na.rm = FALSE)
            flux_mat[count, dr] = flux
            count = count + 1
          }
        }
      }
      
      combined_flux = colMedians(
        flux_mat, na.rm = TRUE
      )
      combined_flux_err = rowDiffs(colQuantiles(
        flux_mat,
        probs = c(0.16, 0.84),
        na.rm = TRUE
      ))/2.0
      sn_tot = combined_flux/combined_flux_err
      sn_tot[is.na(sn_tot)] = 0
      
      N_nonna_rays = colSums(
        1-is.na(flux_mat)
      )
      
      suppressWarnings(magplot(
        dr_vec, 
        combined_flux, 
        log = "xy",
        xlim = c(1, 6000),
        ylim = 10^c(log10(0.01*depth), 6), 
        xlab = "Radius [pixels]",
        ylab = "Flux [microJy]"
      ))
      legend(
        x = "topright",
        paste0(kk), 
        col = "darkorange"
      )
      suppressWarnings(magerr(
        dr_vec, 
        combined_flux, 
        ylo = sqrt( depth^2 + combined_flux_err^2 ),
        col = rgb(0,0,0,0.2)
      ))
      for(ii in 1:6){
        lines(
          dr_vec,
          flux_mat[ii, ],
          col = "darkgrey"
        )
      }

      rr = c(dr_vec[combined_flux > 0.5*depth & !is.na(combined_flux)])
      ff = c(combined_flux[combined_flux > 0.5*depth & !is.na(combined_flux)])
      ff_err = c(combined_flux_err[combined_flux > 0.5*depth & !is.na(combined_flux)])
      ff_err[ff_err == 0] = ff[ff_err == 0]
      sn = ff / ff_err
      
      points(
        rr, 
        ff, 
        pch = 16
      )
      
      ## Fit average profile with Highlander
      Data = list(
        "rr" = c(rr),
        "ff" = c(ff),
        "ff_err" = sqrt( rep(depth, length(ff))^2 + ff_err^2 )
      )
      
      if(length(Data$rr) > 6){
        highout = suppressMessages(Highlander(
          parm = c(0.0, 1.0, 1.0, 4.0, 0.0),
          Data = Data,
          likefunc = LL,
          likefunctype = "CMA", liketype = "max",
          lower = c(
            log10(depth)-1.0, ## log10(Ib)
            log10(1), ## log10(Rb)
            log10(1), ## log10(Re)
            1.0, ## n (sersic index)
            -0.1 ## gamma (power law slope)
          ), 
          upper = c(
            log10(depth)+5.0,
            log10(dim(ref)[1]),
            log10(dim(ref)[1]),
            20.0,
            3.5
          ),
          Niters = c(1000,1000),
          NfinalMCMC = 0,
          seed = 666, 
          parm.names = c("Ib", "Rb", "Re", "n", "gamma")
        ))
        
        fitparm = highout$parm
        sfbFun = function(x){
          10^surf_fun(x = x, fitparm)
        }
        rtest = 1:imdim[1]
        
        curve(
          sfbFun, 1, 5000, add = TRUE, lwd = 3, col = "red"
        )
        
        ridx = (abs(sfbFun(rtest) - depth) / depth) < 0.01
        ridx[is.na(ridx)] = FALSE
        
        if(sum(ridx)==0){
          opt_rad = 50
        }else{
          opt_rad = mean(rtest[ridx], na.rm = TRUE)
        }
        
        box_check = 5 / pixscale(ref, unit = "asec")
        check_flux = combined_flux
        check_flux[is.na(check_flux)] = 0
        check_NA_centre = sum(is.na(ref[test_coords[1], test_coords[2], box = box_check]$imDat)) / box_check^2
        check_rays_noise = median(check_flux[check_flux>0])/depth
        check_NA_rays = sum(N_nonna_rays == 0) / length(N_nonna_rays)
        
        if(check_NA_centre >= 0.5 & check_rays_noise <= 0.8){
          opt_rad = 50
        }
      }else{
        opt_rad = 50
      }

      star_rad = opt_rad * 2.5 ## Because the PSF mask spurs are rad/2.0 and increase by 50%

      abline(
        v = star_rad / 2.0, lty = 2, lwd = 2
      )
      
      abline(
        h = depth, lty = 2, lwd = 2
      )
      legend(
        x = "topleft", 
        pch = c(1, 16, NA, NA, NA, NA),
        lty = c(NA, NA, 1, 1, 2, 2),
        col = c("black", "black", "grey", "red", "black", "black"),
        legend = c("Median profile", "Fit data", "6 rays", "Fit", "Depth", "Star radius"), 
        ncol = 2
      )

      return(star_rad)
    }
    
    star_rad_list[is.na(star_rad_list)] = 1
    
    all_mask = profoundApplyMask(
      image = ref$imDat, mask = psf>0, 
      xcen = gaia$xpix, ycen = gaia$ypix,
      xsize = star_rad_list,
      ysize = star_rad_list
    )
    ## Extra dilation step to capture more flux if needed
    pro_new = profoundProFound(
      image = ref, 
      segim = all_mask$mask, 
      sky = pro$sky,
      magzero = 23.9, 
      rem_mask = TRUE, 
      cliptol = 50, 
      tolerance = Inf
    )
    mask = pro_new$objects_redo
    gaia$PSFMASKRAD = star_rad_list
    if(dim(gaia)[1] == dim(pro_new$segstats)[1]){
      gaia = cbind(gaia, pro_new$segstats)
    }
    return(list("star_mask" = mask, "new_gaia" = gaia))
  }
  
  ref$image$imDat[is.na(ref$image$imDat) | is.infinite(ref$image$imDat)] = NA ## For safefty just get rid of any Infs or NA...
  pro_stars = profoundProFound(
    image = ref$image, 
    magzero = 23.9, 
    rem_mask = TRUE, 
    cliptol = 50,
    boundstats = TRUE
    )
  
  match_gaia = coordmatch(
    coordref = gaia_trim[, c("ra_fix", "dec_fix")],
    coordcompare = pro_stars$segstats[, c("RAmax", "Decmax")]
  )
  
  if(length(match_gaia$bestmatch)>1){ ## Make sure that there are GAIA stars in frame, length==1 => NA and no matches (we should at least have refID and compareID)
    ## calculate radius of stars/diffraction spikes
    
    plot_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_", PIXSCALE, "_star_mask.pdf")
    CairoPDF(plot_stub, width = 10, height = 10)
    
    star_mask_save = star_masker(ref$image, gaia_trim, pro_stars, psf_mask)
    
    message("Star mask complete")
    star_mask_redo = list()
    star_mask_redo$mask = star_mask_save[["star_mask"]]
  
    par(mfrow = c(1,1), mar = rep(0,4), oma = rep(0,4))
    magimage(ref$image$imDat, flip = T, sparse = 1)
    magimage(star_mask_save[["star_mask"]], col = c(NA, rgb(1,0,0,0.4)), add = T)
    for(ii in 1:dim(star_mask_save$new_gaia)[1]){
      points(
        star_mask_save$new_gaia$xpix[ii]+2,
        star_mask_save$new_gaia$ypix[ii]+2,
        pch = paste0(ii), 
        cex = 2.0, 
        col = "darkorange"
      )
    }
    legend(x = "topleft", paste0(VID, "_", MODULE))
    dev.off()
    
    csv_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_", PIXSCALE, "_star_mask.csv")
    fwrite(data.frame(star_mask_save[["new_gaia"]]), file = csv_stub)
  }else{
    star_mask_redo = list(mask = matrix(0, nrow = dim(ref$image)[1], ncol = dim(ref$image)[2]))
  }

  message("Saving star mask")
  
  save_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_", PIXSCALE, "_star_mask.rds")
  saveRDS(object = star_mask_redo, save_stub)
  message("Done!")
} ##<-- Create an empirical star mask from GAIA sources 
star_mask_tile = function(input_args){
  ## Make star mask from input VID star masks
  ## Account for rotation
  
  message("Running Star Mask for Tiles. Making big mask from per VID star masks...")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  if(MODULE == VID){
    message("VID [", VID, "] is the same as MODULE [", MODULE, "]")
    message("Detected big mosaic!")
    
    star_mask_dir = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/")
    if(dir.exists(star_mask_dir)){
      unlink(star_mask_dir, recursive = T)
    }
    dir.create(
      star_mask_dir,
      recursive = T, showWarnings = F
    )
    
    data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") #directory with the propanes
    
    file_list = list.files(data_dir, pattern=paste0(PIXSCALE,".fits"), full.names=TRUE) #list all the .fits files in a directory
    file_names = list.files(data_dir, pattern=paste0(PIXSCALE,".fits"), full.names=FALSE)
    
    file_idx = (
      grepl("F277W", file_list) & 
        grepl(VID, file_list) & 
        grepl(MODULE, file_list) &
        !grepl("hst", file_list)
    )
    
    if(sum(file_idx) == 0){
      file_list = file_list[grepl(VID, file_list) & 
                              grepl(MODULE, file_list) & 
                              !grepl("hst", file_list)][1]
      file_names = file_names[grepl(VID, file_names) & 
                                grepl(MODULE, file_names) & 
                                !grepl("hst", file_names)][1]
    }else{
      file_list = file_list[file_idx]
      file_names = file_names[file_idx]
    }
    
    ref = Rfits_read_image(filename = file_list, ext = 1)
    keyvalues =ref$keyvalues
    input_info = Rfits_read_table(filename = file_list, ext = 'info')
    ## Check if the frames are on the disk
    new_star_dirs = c()
    unique_star_dirs = unique(input_info$path)
    print(head(input_info$full))
    for(i in 1:length(unique_star_dirs)){
      if(!dir.exists(unique_star_dirs[i])){
        message('Patch stack directory that I am reading is missing. Please enter a new one (no quotes):')
        new_star_dir = toString(readLines("stdin", n=1))
      }else{
        new_star_dir = unique_star_dirs[i]
      }
      input_info$path[input_info$path == unique_star_dirs[i]] = new_star_dir
    }
    input_info$full = paste0(input_info$path, "/", input_info$file)

    input_VID = na.exclude(
      str_extract(
        sapply(
          input_info$stub, function(x)str_split_1(x, pattern = "_")
        ),
        "\\b[0-9].{1,9}\\b"
      )
    )
    input_MODULE = na.exclude(
      str_extract(
        sapply(
          input_info$stub, function(x)str_split_1(x, pattern = "_")
        ),
        "\\bNRCA|NRCB\\b"
      )
    )
    
    ## produce star masks for input frames that where the star mask did not already exist
    star_mask_files = foreach(i = 1:length(input_VID), .combine = "c")%do%{
      paste0(ref_dir, "/ProFound/Star_Masks/", input_VID[i], "/", 
                            input_MODULE[i], "/", 
                            input_VID[i], "_", input_MODULE[i], "_", PIXSCALE, "_star_mask.rds")
      }
    check_star_masks_exist = file.exists(star_mask_files)
    check_files_idx = which(!check_star_masks_exist)
    check_grid = data.frame(input_VID[check_files_idx], input_MODULE[check_files_idx])
    temp_args = list(ref_dir = ref_dir)
    
    if(dim(check_grid)[1] > 0){
      for(i in 1:dim(check_grid)[1]){
        temp_args$VID = check_grid[i,1]
        temp_args$MODULE = check_grid[i,2]
        temp_args$PIXSCALE = PIXSCALE
        query_gaia(temp_args)
        star_mask(temp_args)
      }
    }
    
    read_star_masks = foreach(i = 1:length(input_VID))%do%{ 
        star_mask_file = paste0(ref_dir, "/ProFound/Star_Masks/", input_VID[i], "/", 
                                input_MODULE[i], "/", 
                                input_VID[i], "_", input_MODULE[i], "_", PIXSCALE, "_star_mask.rds")
        
        message("Using: ", star_mask_file)
        star_mask = readRDS(star_mask_file)
        
        input_file_keyvalues = Rfits_read_image(
          filename = input_info$full[i],
          ext = 1
        )$keyvalues
        
        image_in = Rfits_create_image(star_mask$mask, keyvalues = input_file_keyvalues)
    }
    class(read_star_masks) = "Rfits_list"
    
    temp_fits_fname = paste0(star_mask_dir, "/star_mask_temp.fits")
    Rfits_write(
      read_star_masks,
      filename = temp_fits_fname
    ) ## Write to disk 
    rm(read_star_masks)
    read_star_masks = Rfits_read(temp_fits_fname, pointer = TRUE)
    big_star_mask = propaneStackWarpInVar(
      image_list = read_star_masks, 
      keyvalues_out = keyvalues, 
      cores = 1, 
      cores_warp = 1, 
      direction = "backward", 
      magzero_in = 23.9,
      magzero_out = 23.9
    )$image$imDat
    big_star_mask[big_star_mask != 0] = 1
    unlink(temp_fits_fname, recursive = TRUE) 
    
    ## Storing large matrices in a list, wtf? Try something new^
    # read_star_masks = foreach(i = 1:length(input_VID))%do%{ 
    #   star_mask_file = paste0(ref_dir, "/ProFound/Star_Masks/", input_VID[i], "/", 
    #                           input_MODULE[i], "/", 
    #                           input_VID[i], "_", input_MODULE[i], "_", PIXSCALE, "_star_mask.rds")
    #   
    #   message("Using: ", star_mask_file)
    #   star_mask = readRDS(star_mask_file)
    #   
    #   input_file_keyvalues = Rfits_read_image(
    #     filename = input_info$full[i],
    #     ext = 1
    #   )$keyvalues
    #   
    #   temp_warp = propaneWarp(
    #     image_in = Rfits_create_image(star_mask$mask, keyvalues = input_file_keyvalues),
    #     keyvalues_out = keyvalues,
    #     doscale = F,
    #     cores = 1
    #   )
    #   star_mask_warp = temp_warp$imDat
    #   # star_mask_warp[star_mask_warp != 0 | !is.na(star_mask_warp)] = 1
    #   star_mask_warp[is.na(star_mask_warp)] = 0
    #   return(star_mask_warp)
    # }
    # 
    # big_star_mask = Reduce("+", read_star_masks)
    # big_star_mask[big_star_mask != 0] = 1
    
    message("Star mask complete")
    message("Saving star mask...")
    
    plot_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_", PIXSCALE, "_star_mask.png")
    png(plot_stub, width = dim(ref)[1], height = dim(ref)[2], res = 72)
    par(mfrow = c(1,1), mar = rep(0,4), oma = rep(0,4))
    magimage(ref$imDat, flip = T, sparse = 2)
    magimage(profoundDilate(big_star_mask, size = 21) - big_star_mask, col = c(NA, "magenta"), add = T, sparse = 2)
    legend(x = "topleft", paste0(VID, "_", MODULE))
    dev.off()
    
    save_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_", PIXSCALE, "_star_mask.rds")
    saveRDS(object = list(mask=big_star_mask), save_stub)
    message("Done!")
    
  }
} ## <-- Create star mask for a big mosaic using per VID star masks

## ProFound source detection codes
do_detect = function(input_args, detect_bands = detect_bands_load, profound_function = profound_detect_master){
  
  message("Running ProDetect")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  star_mask_dir = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/")
  data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") 
  detect_dir = paste0(ref_dir, "/ProFound/Detects/", VID, "/", MODULE, "/") 
  
  if(dir.exists(detect_dir)){
    unlink(detect_dir, recursive = T)
  }
  
  dir.create(detect_dir, recursive = T, showWarnings = F)
  
  message("Start detect!")
  file_list <- list.files(data_dir, pattern=paste0(PIXSCALE, ".fits"), full.names=TRUE) #list all the .fits files in a directory
  file_names <- list.files(data_dir, pattern=paste0(PIXSCALE, ".fits"), full.names=FALSE) #list all the .fits files in a directory
  
  message(paste0("Using detect_band=", detect_bands))
  if(detect_bands == "ALL"){
    idx = grepl(VID, file_list) & grepl(MODULE, file_list) & grepl("patch", file_list)
  }else{
    idx = grepl(VID, file_list) & grepl(MODULE, file_list)  & grepl("patch", file_list) & grepl(detect_bands, file_list)
  }
  
  propanes = lapply(file_list[idx], function(x) Rfits_read(filename = x))
  names(propanes) = file_names[idx]
  
  N_frames = length(propanes)
  message(paste0("Running ", N_frames, " frames from ", "VISIT: ", VID, ", module: ", MODULE))
  
  star_mask_file = list.files(star_mask_dir, pattern = paste0(PIXSCALE,"_star_mask.rds"), full.names = T)
  
  if(length(star_mask_file>0)){
    message(paste0("Load in star mask from: ", star_mask_file))
    star_mask = readRDS(star_mask_file)
  }else{
    star_mask = list()
    star_mask$mask = 0
  }
  
  message(paste0("Stacking detect band ", names(propanes), collapse="\n"))
  
  stack_image = propaneStackFlatInVar(image_list = lapply(propanes, function(x)x$image[,]), 
                                      skyRMS_list = lapply(propanes, function(x){1 / sqrt(x$inVar[,]$imDat)}),
                                      magzero_in = 23.9, magzero_out = 23.9)
  stack_image_fits = Rfits_create_image(
    data = stack_image$image,
    keyvalues = propanes[[1]]$image$keyvalues
  )
  
  med_stack = propaneStackFlatMed(
    image_list = lapply(propanes, function(x)x$image[,])
  )
  
  patch_stack_image = propanePatch(
    stack_image_fits,
    med_stack$image
  )

  message("Finding sources with ProFound...")
  profound = profound_function(frame = patch_stack_image$image, 
                               skyRMS = stack_image$skyRMS, 
                               star_mask = star_mask$mask)
  if(dim(profound$segstats)[1] == 0){
    message("No objects found!")
  }else{
    profound$segstats$MODULE = rep(MODULE, dim(profound$segstats)[1])
    profound$segstats$VID = rep(VID, dim(profound$segstats)[1])
    
    stack_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_", PIXSCALE, "_profound_stack.fits")
    profound_stack = list(
      stack = patch_stack_image$image[,],
      segim = profound$segim,
      segim_orig = profound$segim_orig
    )
    Rfits_write(data = profound_stack, filename = stack_stub)
    
    catalogue_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_", PIXSCALE, "_segstats.csv")
    fwrite(profound$segstats, file = catalogue_stub)
    
    profound_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_", PIXSCALE, "_profound.rds")
    saveRDS(profound, file = profound_stub)
    
    ## I have no idea why plot.profound spits out and error 
    ## Error in Rwcs_p2s(rep(xlo, leny), ylo:yhi, keyvalues = keyvalues, pixcen = "R",  : 
    ## Assertion on 'y' failed: Must have length 6432, but has length 6433
    ## So I just put plot profound into a tryCatch
    plot_profound = function(x){
      tryCatch(
        {plot(x, coord.type = "deg")
        },
        error = function(cond){
          magplot(
            NA, xlim = c(-1,1), ylim = c(-1,1)
          )
          text(0,0,labels="NA", cex = 3.0)
          legend(x = "topright", legend =  conditionMessage(cond))
        },
        warning = function(cond){
          NULL
        },
        finally = {
        }
      )
    }
    
    plot_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_", PIXSCALE, "_profound_plot.pdf")
    CairoPDF(plot_stub, width = 10, height = 10)
    plot_profound(profound)
    par(mfrow = c(1,2), mar = rep(0,4), oma = c(1.5, 1.5, 0.5, 0.5), mai = rep(0,4))
    profoundSegimPlot(image = profound$image, segim = profound$segim, header = profound$header, mask = profound$mask, sparse=1)
    legend(x = "bottomleft", legend = "segim")
    profoundSegimPlot(image = profound$image, segim = profound$segim_orig, header = profound$header, mask = profound$mask, sparse=1)
    legend(x = "bottomleft", legend = "segim_orig")
    par(mfrow = c(1,2), mar = rep(0,4), oma = c(1.5, 1.5, 0.5, 0.5), mai = rep(0,4))
    magimage(profound$objects)
    legend(x = "bottomleft", legend = "Objects")
    magimage(profound$objects_redo)
    legend(x = "bottomleft", legend = "Objects_redo")
    dev.off()
  }
  return(NULL)
}

## ProFound measurement codes
getCircle = function(r){
  d = 2*r
  m = matrix(1:d^2,d,d)
  g = expand.grid(1:d, 1:d)
  if(d%%2==0){c=r+0.5}
  if(d%%2==1){c=r+1}
  g$d2 = sqrt((g$Var1-c)^2 + (g$Var2-c)^2)
  g$inside = g$d2 <= r
  return(m[as.matrix(g[g$inside, c('Var1', 'Var2')])])
}
err_sampler = function(r, imdat, fname, mask, root=root_sample){ 
  use = getCircle(r)
  check = getCircle(r*1.5)
  set.seed(r)
  collect = c()
  n=1
  whilebreak_thresh = 10000
  whilebreak = 0
  dum_imdat = imdat
  dum_imdat[mask] = NA
  
  CairoPDF(file.path(root,paste0(fname,'_',r,'.pdf')))
  magimage(dum_imdat, qdiff=T)
  xy = dim(imdat)
  while(n <= 200){
    if(whilebreak>=whilebreak_thresh){
      message(paste0("Can't find a spot for the ", n, "th aperture in ", whilebreak_thresh, " tries :("))
      break
    }
    x = sample(1:xy[1], 1)
    y = sample(1:xy[2], 1)
    if(x-1.5*r<1 | x+1.5*r>xy[1] | y-1.5*r<1 | y+1.5*r>xy[2]){
      next
    }
    mask2d = mask[(x-1.5*r):(x+1.5*r-1),(y-1.5*r):(y+1.5*r-1)]
    
    cut2d = imdat[(x-r):(x+r-1), (y-r):(y+r-1)]
    out = cut2d[use]
    if( any(is.na(mask2d[check])) | any(is.nan(mask2d[check])) | any(mask2d[check]) | (any(is.na(out)) | any(is.nan(out))) ){
      whilebreak = whilebreak + 1
      next
    }else{
      n = n+1
      profoundDrawEllipse(x, y, r, col='red')
      collect = rbind(collect, data.frame(id = n, sum=sum(out)))                # options to output more stats mean=mean(out), sd=sd(out)
      mask[(x-r):(x+r-1), (y-r):(y+r-1)][use] = NA  # un-comment to avoid re-sampling the same pixels
    }
  }
  q = unname(quantile(collect$sum, probs=c(0.5, 0.16), na.rm=F))
  maghist(collect$sum, breaks = 50, verbose = F)
  abline(v=q, col=c('red','red'), lty=c(2,2))
  dev.off()
  return(collect)
}
error_scaling = function(flux_err, N100, m, c){
  flux_err[!is.na(flux_err) & (flux_err < N100^m * 10^c)] = N100[!is.na(flux_err) & (flux_err < N100^m * 10^c)]^m * 10^c
  return(flux_err)
}
do_measure = function(input_args){
  
  plot_profound = function(x){
    tryCatch(
      {plot(x, coord.type = "deg")
      },
      error = function(cond){
        magplot(
          NA, xlim = c(-1,1), ylim = c(-1,1)
        )
        text(0,0,labels="NA", cex = 3.0)
        legend(x = "topright", legend =  conditionMessage(cond))
      },
      warning = function(cond){
        NULL
      },
      finally = {
      }
    )
  }
  
  message("\n Running ProMeasure \n")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  sampling_cores = input_args$sampling_cores
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") #directory with the propanes
  
  detect_dir = paste0(ref_dir, "/ProFound/Detects/", VID, "/", MODULE)
  
  inspect_dir = paste0(ref_dir, "/ProFound/Inspect/", VID, "/", MODULE)
  
  sampling_dir = paste0(ref_dir, "/ProFound/Sampling/", VID, "/", MODULE)
  
  measurements_dir = paste0(ref_dir, "/ProFound/Measurements/", VID, "/", MODULE)
  
  ######## Load segim and mask ####################
  pro_path = list.files(detect_dir, pattern = glob2rx(paste0("*", VID, "*", MODULE, "*", PIXSCALE, "*rds")), full.names = T)
  assert(checkFileExists(pro_path))
  
  super_pro = readRDS(pro_path)
  segim = super_pro$segim
  segim_col = super_pro$segim_orig
  mask = super_pro$mask
  super_segstats = super_pro$segstats
  rm(super_pro)
  gc()
  
  if(dim(super_segstats)[1] == 0){
    message("No objects to measure!")
  }else{
    ####### Load data images ########################
    assert(checkDirectoryExists(data_dir))
    
    data.list = list.files(data_dir, 
                           pattern = paste0(PIXSCALE, ".fits$"), 
                           full.names = T, 
                           recursive = T)
    data.names = list.files(data_dir, 
                            pattern = paste0(PIXSCALE, ".fits$"), 
                            full.names = F, 
                            recursive = T)
    
    filter.names = toupper(str_extract(data.list, "F\\d{3,}[A-Z]{1,}|f\\d{3,}[a-z]{1,}"))
    filter.names[is.na(filter.names)] = sapply(data.list[is.na(filter.names)], function(x)str_split_1(x, "_")[6]) ## Filter name should hopefully always be in position 6
    
    bad_name = data.list[is.na(filter.names)]
    
    ext = sapply(data.list, function(x)Rfits_extname_to_ext(x, extname = "image"))
    images = Rfits_make_list(
      filelist = data.list,
      extlist = ext,
      pointer = T
    )
    names(images) = filter.names      
    
    inVar = lapply(data.list, function(x){
      
      ext = Rfits_extname_to_ext(filename = x, extname = "inVar")
      
      if(is.na(ext)){
        return(NULL)
      }else{
        inVar = Rfits_point(x, ext = ext)
        return(inVar)
      }
      
    })
    names(inVar) = filter.names
    # filter ordering is not important
    
    ####### Initiate for sampling error #############
    r = seq(2,10)
    master = data.frame(r=r, a=unlist(lapply(r, function(x) length(getCircle(x)))))
    
    if(!dir.exists(sampling_dir)){dir.create(sampling_dir, recursive = T)}else{
      message("Unlinking sample dir")
      unlink(sampling_dir, recursive = T)
      dir.create(sampling_dir, recursive = T)
    }                           # Create directory for the sampling stuff
    
    error_file = file.path(sampling_dir,paste(VID,MODULE,PIXSCALE,'error_fit.txt', sep='_'))
    cat('Filters','Slopes','Intercepts','\n', file = error_file, sep='\t')          # Create file for the fitted relation of error sampling
    
    ####### Initiate csv for saving output ##########
    csvout = data.frame(segID=super_segstats$segID)
    
    if(!dir.exists(measurements_dir)){dir.create(measurements_dir, recursive = T)}else{
      message("Unlinking measurement dir")
      unlink(measurements_dir, recursive = T)
      dir.create(measurements_dir, recursive = T)
    }                                                    # Create directory for the output stuff
    
    if(!dir.exists(inspect_dir)){dir.create(inspect_dir, recursive = T)}else{
      message("Unlinking inspect dir")
      unlink(inspect_dir, recursive = T)
      dir.create(inspect_dir, recursive = T)
    }                                                    # Create directory for inspection plots
    ###### Main #############
    for(ff in names(images)){
      message(paste0("Running ProMeasure on: ", ff, "\n"))
      filt = images[[ff]][,]
      filt_invar = inVar[[ff]][,]
      # filt[[ff]]$imDat[filt[[ff]]$imDat==0L] = NA
      
      dum_pro = measure_profound(filt, inVar = filt_invar, segim, mask)
      csvout[paste0(ff,'_fluxt')] = dum_pro$segstats$flux
      csvout[paste0(ff,'_fluxt_err')] = dum_pro$segstats$flux_err
      gc()
      
      csvout[paste0(ff,'_maskfrac')] = dum_pro$segstats$Nmask/dum_pro$segstats$Nedge
      
      img = filt$imDat - dum_pro$sky
      img[img == 0L] = NA
      img[mask] = NA
      sample_mask = (dum_pro$objects_redo | mask)
      
      message("\n ...Error sampling... \n")
      row=foreach(rr = r, .combine="c")%do%{
        dum = err_sampler(rr, img, ff, mask = sample_mask, root=sampling_dir)
        q = unname(quantile(dum$sum, probs=c(0.5, 0.16), na.rm=F))
        sig = q[1]-q[2]
        std = sd(dum$sum)
        # print(c(sig,std))
        sig
      }
      
      master[ff] = row
      Nsam = sample(length(dum_pro$segstats$segID), 200)
      fit_samp = lm(log10(master[[ff]]) ~ poly(log10(master$a), 1, raw = T))
      slope = unname(fit_samp$coefficients)[2]
      intercept = unname(fit_samp$coefficients)[1]
      cat(ff, slope, intercept, '\n', file = error_file, sep = '\t', append = T)      # Save the fitted relation
      
      csvout[paste0(ff,'_scaled_fluxt_err')] = error_scaling(dum_pro$segstats$flux_err, dum_pro$segstats$N100, slope, intercept)
      
      message(("\n ...Running colour photometry... \n"))
      dum_pro_col = measure_profound(filt, inVar = filt_invar, segim_col, mask, redosegim = F) #don't redilate segments e.g., colour photometry mode
      csvout[paste0(ff,'_fluxc')] = dum_pro_col$segstats$flux
      csvout[paste0(ff,'_fluxc_err')] = dum_pro_col$segstats$flux_err
      csvout[paste0(ff,'_scaled_fluxc_err')] = error_scaling(dum_pro_col$segstats$flux_err, dum_pro_col$segstats$N100, slope, intercept)
      rm(dum_pro_col)
      gc()
      saveRDS(dum_pro, file.path(measurements_dir,paste(VID,MODULE,PIXSCALE,ff,"results.rds",sep='_')))
      
      CairoPDF(file.path(inspect_dir,paste0("profound_",ff,"_inspect.pdf")), width = 10, height = 10 )
      plot_profound(dum_pro)
      par(mfrow = c(1,2), mar = rep(0,4), oma = rep(0,4))
      profoundSegimPlot(dum_pro, sparse=1, sky = dum_pro$sky, qdiff = T)
      legend(x="topleft", legend="Total photometry")
      dev.off()
      
      CairoPDF(file.path(sampling_dir,paste0(ff,"_sample_error.pdf")))                            # Plot the sampled relation
      magplot(master$a, master[[ff]], col='white', pch=1,
              log='xy', xlim = c(5,5000), ylim = c(0.0001,1.0),
              main = paste(ff,VID,MODULE,sep='_'), ylab = 'Error/micro Jy', xlab = 'Area/pix')
      points(dum_pro$segstats[Nsam, 'N100'], dum_pro$segstats[Nsam, 'flux_err'], col='darkgrey', cex=0.3, pch=3)
      scaled_error = error_scaling(dum_pro$segstats[Nsam, 'flux_err'], dum_pro$segstats[Nsam, 'N100'], slope, intercept)
      points(dum_pro$segstats[Nsam, 'N100'], scaled_error, col='black', cex=0.5, pch=1)
      points(master$a, master[[ff]], col='red')
      abline(intercept, slope, col='red', lty=2)
      legend('bottomright', legend = c('sample', 'profound', 'corrected'), col = c('red','darkgrey','black'), pch = c(1,3,1), cex=1.2)
      dev.off()
      
      gc()
    }
    write.csv(csvout, file = file.path(measurements_dir,paste(VID,MODULE,PIXSCALE,"photometry.csv",sep='_')))
  }
}

## HST codes
hst_warp_stack = function(input_args){
  
  message("Running hst_warp_stack")

  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  cores_stack = input_args$cores_stack
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") #directory with the propanes

  existing_hst_files = list.files(
    path = data_dir,
    pattern = ".fits",
    full.names = T
  )
  existing_hst_files = existing_hst_files[
    grepl("hst", existing_hst_files, ignore.case = T)
  ]
  file.remove(
    existing_hst_files
  )
  
  HST_cutout_dir = paste0(ref_dir, "/ProFound/HST_cutout/", VID, "/", MODULE, "/") ## Holding pen for the MAST Hubble MVMs / HAPs
  
  #load in the target reference frame
  jwst_ls = c(
    list.files(
      paste0(ref_dir, "/Patch_Stacks/"),
      pattern = glob2rx(paste0("*", VID, "*", MODULE, "*long.fits")),
      full.names = T
      ),
    list.files(
      paste0(ref_dir, "/Mosaic_Stacks/Patch/"),
      pattern = glob2rx(paste0("*", VID, "*", "*long.fits")),
      full.names = T
    )
  )
  
  target_path = jwst_ls[grepl(MODULE, jwst_ls) & grepl("F444W", jwst_ls)]
  if(length(target_path)==0){
    idx_path = length(jwst_list)
    target = Rfits_read(jwst_ls[idx_path])
  }else{
    target = Rfits_read(target_path)
  }
  
  ## check to see if there is a profound catalogue previously made
  dir_profound = paste0(ref_dir, "/ProFound/Detects/", VID, "/", MODULE, "/")
  pro_file = list.files(
    dir_profound, pattern = ".csv", full.names = T
  )
  if(length(pro_file)==1){
    message(paste0("Reading: ", pro_file))
    
    pro_ref = fread(pro_file)
    pro_ref = list(segstats=data.frame(
      pro_ref[pro_ref$mag <25, ]
      ))
  }else{
    pro_ref = profoundProFound(
      image = target$image[,],
      pixcut = 20,
      skycut = 10.0,
      cliptol = 100,
      tolerance = Inf,
      box = dim(target$image)/10.0,
      mask = is.infinite(target$image[,]$imDat) | is.na(target$image[,]$imDat),
      magzero = 23.9,
      rem_mask = T
    )
  }

  files_temp_temp = list.files(path = HST_cutout_dir, pattern = ".fits$", recursive = T, full.names = T)
  
  get_rid_of_combined = grepl("combined", files_temp_temp)
  HST_files = files_temp_temp[!get_rid_of_combined]
  
  extlist = foreach(i = 1:length(HST_files), .combine = "c")%do%{
    
    extname_try1 = Rfits_extname_to_ext(
      HST_files[i], "SCI"
    )
    extname_try2 = Rfits_extname_to_ext(
      HST_files[i], "image"
    )

    extname = c(extname_try1, extname_try2)
    if(sum(is.na(extname)==1)){
      extname = 1
    }
    
  }
  
  if(length(HST_files) == 0){
    message("No HST in directory!")
    return(NULL)
  }else{
    frame_find = propaneFrameFinder(filelist = HST_files, 
                                    extlist = extlist,
                                    RAcen = target$image$keyvalues$CRVAL1, 
                                    Deccen = target$image$keyvalues$CRVAL2, 
                                    rad = pixscale(target$image, unit="deg")*dim(target$image)[1],
                                    plot = F)
    files_temp = frame_find$full
    files_names = frame_find$file
    
    if(is.null(files_temp)){
      message("No HST to fold in!")
    }else{
      hst_key_scan = Rfits_key_scan(files_temp, keylist = c("FILTER", "DETECTOR", "INSTRUME"))
      
      if(sum(is.na(hst_key_scan[,c("FILTER", "DETECTOR", "INSTRUME")])) == 
         prod(dim(hst_key_scan[,c("FILTER","DETECTOR","INSTRUME")]))){
        do_sky_rem = F
        message("Skipping sky rem. No FILTER, DETECTOR, INSTRUME keyword. Big mosaic?")
      }else{
        do_sky_rem = F
      }
      
      filters_na = na.exclude(
        str_extract(
          sapply(files_names, function(x)str_split_1(x, pattern = "_")), 
          "\\bF\\d{3,}[A-Z]{1,}|f\\d{3,}[a-z]{1,}\\b"))
      detector_na = na.exclude(
        str_extract(
          sapply(files_names, function(x)str_split_1(x, pattern = "_")), 
          "\\IR|UVIS|WFC\\b"))
      instrume_na = na.exclude(
        str_extract(
          sapply(files_names, function(x)str_split_1(x, pattern = "_")), 
          "\\bACS|WFC3\\b"))
      
      hst_key_scan$FILTER[is.na(hst_key_scan$FILTER)] = filters_na[is.na(hst_key_scan$FILTER)]
      hst_key_scan$DETECTOR[is.na(hst_key_scan$DETECTOR)] = detector_na[is.na(hst_key_scan$DETECTOR)]
      hst_key_scan$INSTRUME[is.na(hst_key_scan$INSTRUME)] = instrume_na[is.na(hst_key_scan$INSTRUME)]
      
      hst_key_scan$FILTER = toupper(hst_key_scan$FILTER)
      hst_key_scan$DETECTOR = toupper(hst_key_scan$DETECTOR)
      hst_key_scan$INSTRUME = toupper(hst_key_scan$INSTRUME)
      
      hst_filters = unique(hst_key_scan$FILTER)
      NFilters = length(hst_filters)
      
      #loop through HST filters
      foreach(j = 1:NFilters)%dopar%{
        hst_info = hst_key_scan[hst_key_scan$FILTER==hst_filters[j], ]
        hst_filter = unique(hst_info$FILTER)
        hst_detector = unique(hst_info$DETECTOR)
        hst_instrume = unique(hst_info$INSTRUME)
        if(length(hst_instrume) > 1){
          # just get the acs filters 
          files_idx = grepl(hst_filters[j], files_temp, ignore.case = T) & grepl("acs", files_temp, ignore.case = T)
          files = files_temp[files_idx]
          # extlist = ifelse(
          #   sapply(files, Rfits_nhdu) == 1, 1, 2
          # )
          frames = Rfits_make_list(files, extlist = frame_find$ext[files_idx], pointer = F)
          hst_detector = hst_detector[grepl("wfc", hst_detector, ignore.case = T)]
          hst_instrume = hst_instrume[grepl("acs", hst_instrume, ignore.case = T)]
        }else{
          files_idx = grepl(hst_filters[j], files_temp, ignore.case = T)
          files = files_temp[files_idx]
          # extlist = ifelse(
          #   sapply(files, Rfits_nhdu) == 1, 1, 2
          # )
          frames = Rfits_make_list(files, extlist = frame_find$ext[files_idx], pointer = T)
          
        }
        
        #prepare input frames for stacking
        magzero_list = c()
        image_list = {}
        inVar_list = {}
        tweak_pars_df = {}
        idxs = {}
        for(i in 1:length(frames)){
          message("Running: ", files[i])
          
          hst_magzero = c(frames[[i]]$keyvalues$MAGZERO)
          if(is.null(hst_magzero)){
            hst_magzero = -2.5*log10(frames[[i]]$keyvalues$PHOTFLAM)-5*log10(frames[[i]]$keyvalues$PHOTPLAM)-2.408
           ## Assumes FLUXUNIT is ELECTRONS/SECOND
          }
          
          magzero_list = c(magzero_list, hst_magzero)
          
          temp = frames[[i]][,]
          temp$imDat[is.infinite(temp$imDat)] = NA
          
          tryCatch_profoundProFound = function(temp, hst_magzero){
            tryCatch(
              {profoundProFound(
                image = temp,
                pixcut = 20,
                skycut = 10.0,
                cliptol = 100,
                tolerance = Inf,
                box = dim(temp)/10.0,
                mask = is.infinite(temp$imDat) | is.na(temp$imDat),
                magzero = hst_magzero,
                rem_mask = T
              )
              },
              error = function(cond){
                NULL
              },
              warning = function(cond){
                NULL
              },
              finally = {
              }
            )
          }

          pro_test = tryCatch_profoundProFound(
            temp = temp,
            hst_magzero = hst_magzero
          )

          if(is.null(pro_test$segstats) | is.null(pro_test)){
            message("No sources in frame")
            next 
          }
          match_coord = coordmatch(
            coordref = pro_ref$segstats[, c("RAmax", "Decmax")],
            coordcompare = pro_test$segstats[, c("RAmax", "Decmax")]
          )
          
          if(length(match_coord$Nmatch)==1){
            message("Insufficient source density for Tweak. Skipping.")
            next
            # pro_tweak = list(
            #   par = c(0,0,0)
            # )
          }else{
            idxs = c(idxs, i)
            pro_tweak = propaneTweakCat(cat_ref = cbind(pro_ref$segstats$RAmax[match_coord$bestmatch$refID],
                                                              pro_ref$segstats$Decmax[match_coord$bestmatch$refID]),
                                              cat_pre_fix = cbind(pro_test$segstats$RAmax[match_coord$bestmatch$compareID],
                                                                  pro_test$segstats$Decmax[match_coord$bestmatch$compareID]),
                                              delta_max = c(20.0,1.0), mode = "coord", keyvalues_pre_fix = temp$keyvalues)
          }
            temp = propaneWCSmod(
            temp,
            pro_tweak$par[1], pro_tweak$par[2], pro_tweak$par[3]
          )
          
          tweak_pars_df = c(list(pro_tweak$par), tweak_pars_df)

          if(!do_sky_rem){
            image_list = c(
              list(
                temp
              ),
              image_list
            )
          }else{
            message(paste0("Profound: ", files[i])) 
            
            pro = profoundProFound(
            image = temp,
            box = 100,
            skycut = 2.0,
            pixcut = 5.0,
            tolerance = Inf,
            redoskysize = 21,
            roughpedestal = T,
            redosky = F,
            rem_mask = T
          )
            
          inVar_list = c(list(median(pro$skyRMS[pro$objects_redo == 0], na.rm = T)^-2), inVar_list)
          image_list = c(list(temp-pro$sky), image_list)
          }
        }
        
        if(length(image_list) > 0){
          #stack
          output_stack = propaneStackWarpInVar(image_list = image_list,
                                               # inVar_list = inVar_list,
                                               magzero_in = hst_magzero,
                                               magzero_out = 23.9,
                                               keyvalues_out = target$image$keyvalues,
                                               cores = cores_stack,
                                               cores_warp = 1)
          
          output_stack$image$imDat[output_stack$image$imDat == 0 | is.infinite(output_stack$image$imDat)] = NA
          
          tweak_pars_df = transpose(data.frame(tweak_pars_df))
          names(tweak_pars_df) = c('dx', 'dy', 'dphi')
          df_info = bind_cols(hst_info[idxs,], tweak_pars_df)
          
          final_output = list()
          final_output$image = output_stack$image
          final_output$image$keyvalues$EXTNAME = "image"
          final_output$info = df_info
          
          class(final_output) = "Rfits_list"
          
          file_name = paste0(data_dir, 
                             "/warp_hst_", 
                             hst_instrume, "_", hst_detector, "_", VID, "_", hst_filter, "_", 
                             MODULE, "_long.fits")
          Rfits_write(final_output, filename = file_name)
          
          rm(image_list)
          rm(output_stack)
        }

      }
    }
  }
}
copy_hst_for_tile = function(input_args){
  
  message("Running hst_warp_stack")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  cores_stack = input_args$cores_stack
  
  HST_cutout_dir = paste0(ref_dir, "/ProFound/HST_cutout/", VID, "/", MODULE, "/") ## Holding pen for the MAST Hubble MVMs / HAPs
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  if(VID == MODULE){
    message("VID [", VID, "] is the same as MODULE [", MODULE, "]")
    message("Detected big mosaic!")

    mosaic_files = list.files(paste0(ref_dir, "/Mosaic_Stacks/Patch/"), full.names=T, pattern = glob2rx(paste0("*", VID, "*.fits"))) 

    input_info = foreach(i = 1:length(mosaic_files), .combine = "bind_rows")%do%{
      info = Rfits_read_table(filename = mosaic_files[i], ext = 7)
      input_VID = na.exclude(
        str_extract(
          sapply(
            info$stub, function(x)str_split_1(x, pattern = "_")
          ),
          "\\b[0-9].{1,9}\\b"
        )
      )
      input_MODULE = na.exclude(
        str_extract(
          sapply(
            info$stub, function(x)str_split_1(x, pattern = "_")
          ),
          "\\bNRCA|NRCB\\b"
        )
      )
        
      df = data.frame(
        "VID" = input_VID,
        "MODULE" = input_MODULE
      )
    }
    input_info = unique(input_info)

    dir.create(HST_cutout_dir, recursive = T)
    
    previous_hst = foreach(i = 1:dim(input_info)[1], .combine = "c")%do%{
      list.files(
        paste0(ref_dir, "/ProFound/HST_cutout/", input_info$VID[i], "/", input_info$MODULE[i], "/"),
        pattern = glob2rx("*hst*.fits"),
        full.names = T
      )
    }
    
    previous_hst_names = foreach(i = 1:dim(input_info)[1], .combine = "c")%do%{
      list.files(
        paste0(ref_dir, "/ProFound/HST_cutout/", input_info$VID[i], "/", input_info$MODULE[i], "/"),
        pattern = glob2rx("*hst*.fits"),
        full.names = F
      )
    }
    
    message(
      "Copying: ",
      paste0(previous_hst_names, sep = "\n")
    )
    
    file.copy(
      from = previous_hst,
      to = paste0(HST_cutout_dir, "/", 1:length(previous_hst_names), "_", previous_hst_names),
      overwrite = T
    )
    do_sky_rem = F
  }
}


## Chop up the large mosaic
frame_chunker = function(input_args){
  
  ## For a given big mosaic, chop it up and save in Data and star mask
  
  ## For scripting use
  area_block = ifelse(is.null(input_args$chunk_size), 10, input_args$chunk_size) ## amin^2
  buffer_size = ifelse(is.null(input_args$buffer_size), 200, input_args$chunk_size) ## pixels
  
  message("Processing ~", area_block, "amin^2 chunks...")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  star_mask_dir = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/")
  data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") 
  
  file_list <- list.files(data_dir, pattern=paste0(PIXSCALE, ".fits"), full.names=TRUE) #list all the .fits files in a directory
  file_names <- list.files(data_dir, pattern=paste0(PIXSCALE, ".fits"), full.names=FALSE) #list all the .fits files in a directory
  
  get_middle <- function(vec) {
    n = length(vec)
    if (n < 2) {
      stop("Vector must have at least two elements")
    }
    mid_index = ceiling(n / 2)  # Find middle index (round up)
    return(c(vec[mid_index], vec[mid_index + (n %% 2 == 0)]))
  }
  
  N_blocks_save = 0
  chunk_cache = c()
  
  filter_list = unname(sapply(
    file_list, 
    function(x){
      load_info = Rfits_read_table(
        x, ext = 'info'
      )
      unique(load_info$FILT)
    }
  ))
  
  det_bands_idx = which( grepl(detect_bands_load, filter_list) )
  
  for(ii in 1:length(file_list)){
    
    message("Chopping up ", file_list[ii])
    frame = Rfits_read(
      file_list[ii],
      pointer = TRUE
    )
    
    tot_area = pixarea(frame$image, unit = "amin2") * prod(dim(frame$image))
    imdim = dim(frame$image)
    
    ## How many 10 arcmin^2 blocks will I need?
    N_blocks = ceiling( tot_area/area_block )
    if(ii == 1){
      N_blocks_save = N_blocks
    }
    div = seq_len(N_blocks)
    div_multiples = div[N_blocks %% div == 0]
    N_grid = get_middle(div_multiples) + 1
    
    xpos = seq(0, nrow(frame$image), length.out = ifelse(imdim[1] >= imdim[2], N_grid[2], N_grid[1])) ## Use the bigger number for width if width > height
    ypos = seq(0, ncol(frame$image), length.out = ifelse(imdim[2] >= imdim[1], N_grid[2], N_grid[1])) ## Use the bigger number for height if height > width
    xpos_cen = xpos[1:(length(xpos)-1)] + diff(xpos)/2.0
    ypos_cen = ypos[1:(length(ypos)-1)] + diff(ypos)/2.0
    
    coord_grid = expand.grid(
      xpos_cen, ypos_cen
    )
    
    grid_s2p = Rwcs_p2s(
      x = coord_grid[,1],
      y = coord_grid[,2],
      keyvalues = frame$image$keyvalues
    )
    
    boxsize = 2*ceiling(
      c(median(diff(xpos)/2.0),
        median(diff(ypos)/2.0))
    ) + buffer_size
    
    boxsize_no_buffer = 2*ceiling(
      c(median(diff(xpos)/2.0),
        median(diff(ypos)/2.0))
    )
    
    chunk_info = data.frame(
      'tile_xcen' = coord_grid[,1],
      'tile_ycen' = coord_grid[,2],
      'tile_racen' = grid_s2p[,1],
      'tile_deccen' = grid_s2p[,2],
      
      'buffer_size' =  rep(buffer_size, dim(grid_s2p)[1]),
      
      'boxsize_x' = rep(boxsize[1], dim(grid_s2p)[1]),
      'boxsize_y' = rep(boxsize[2], dim(grid_s2p)[1]),
      
      'boxsize_no_buffer_x' = rep(boxsize_no_buffer[1], dim(grid_s2p)[1]),
      'boxsize_no_buffer_y' = rep(boxsize_no_buffer[2], dim(grid_s2p)[1]),
      
      'tile_indices' = 1:dim(grid_s2p)[1]
    )
    
    for(jj in 1:dim(grid_s2p)[1]){
      
      fstub = file_names[ii]
      split_VID_name = str_split_1(fstub, VID)
      fstub_chunk = paste0(split_VID_name[1], VID, "chunk", jj, split_VID_name[2])
      
      if(VID == MODULE){
        data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "chunk", jj, "/", MODULE, "chunk", jj, "/")
      }else{
        data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "chunk", jj, "/", MODULE, jj, "/")
      }
      
      fname_chunk = paste0(data_dir, fstub_chunk)
      
      if(ii == 1){
        if(dir.exists(data_dir)){
          unlink(data_dir)
        }
      }
      
      message('Working on chunk ', jj)
      trim_frame_image = frame$image[grid_s2p[jj,1], grid_s2p[jj,2], box = boxsize, type = "coord"]
  
      dim_chunk = prod(dim(trim_frame_image))
      non_na = !is.na(trim_frame_image$imDat)
      dim_usable = sum(non_na)
    
      ## Some area must be usable
      if(dim_usable > 0){

        if(dim_usable/dim_chunk < 0.1){ ## Shuffle point to mid point of data points
          
          non_na_idx = colMedians(which(
            non_na, arr.ind = TRUE
          ))
          non_na_idx_coord=Rwcs_p2s(
            x = non_na_idx[1], y = non_na_idx[2],
            keyvalues = trim_frame_image$keyvalues
          )
          non_na_idx_orig_grid = Rwcs_s2p(
            RA = non_na_idx_coord[1], Dec = non_na_idx_coord[2],
            keyvalues = frame$image$keyvalues
          )
          chunk_info$tile_xcen[jj] = non_na_idx_orig_grid[1]
          chunk_info$tile_ycen[jj] = non_na_idx_orig_grid[2]
          chunk_info$tile_racen[jj] = non_na_idx_coord[1]
          chunk_info$tile_deccen[jj] = non_na_idx_coord[2]
          
          grid_s2p[jj, 1] = non_na_idx_coord[1]
          grid_s2p[jj, 2] = non_na_idx_coord[2]
          
          trim_frame_image = frame$image[grid_s2p[jj,1], grid_s2p[jj,2], box = boxsize, type = "coord"]
        }
        dir.create(data_dir, recursive = TRUE)

        trim_frame_inVar = frame$inVar[grid_s2p[jj,1], grid_s2p[jj,2], box = boxsize, type = "coord"]
        trim_frame_weight = frame$weight[grid_s2p[jj,1], grid_s2p[jj,2], box = boxsize, type = "coord"]

        chunk_info$tile_index = rep(jj, dim(grid_s2p)[1])

        chunk_Rfits = list(
          "image" = trim_frame_image,
          "inVar" = trim_frame_inVar,
          "weight" = trim_frame_weight,
          'chunk_info' = chunk_info,
          'info' = frame$info
        )
        class(chunk_Rfits) = "Rfits_list"

        Rfits_write(
          data = chunk_Rfits,
          filename = fname_chunk
        )
        
        if(ii %in% det_bands_idx){
          chunk_cache = c(chunk_cache, jj)
        }
      }else{
        message("Not enough data in chunk. Skipped...")
      }
    }
  }
  
  for(ii in unique(chunk_cache)){
    
    temp_args = list(ref_dir = ref_dir)
    
    ## By default, just run through all of the relevant steps
    temp_args = list(ref_dir = ref_dir)
    
    if(VID == MODULE){
      temp_args$VID = paste0(VID, "chunk", ii)
      temp_args$MODULE = paste0(MODULE, "chunk", ii)
      temp_args$PIXSCALE = PIXSCALE
      star_mask_tile(
        temp_args
      )
      do_detect(
        temp_args
      )
      do_measure(
        temp_args
      )
    }else{
      temp_args$VID = paste0(VID, "chunk", ii)
      temp_args$MODULE = MODULE
      temp_args$PIXSCALE = PIXSCALE
      star_mask(
        temp_args
      )
      do_detect(
        temp_args
      )
      do_measure(
        temp_args
      )
    }
  }
  
  combine_chunks(input_args) ## combine for a final catalogue at the end
  return(NULL)
}
## Combine chunk catalogue
combine_chunks = function(input_args){
  
  # Testing purposes - FML
  # input_args = list(
  #   ref_dir = "/Volumes/BRAIDY/JWST/JP",
  #   VID = "JADESGOODSS",
  #   MODULE = "JADESGOODSS",
  #   PIXSCALE = "short"
  # )
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  
  message(paste0("Combining chunks for ", VID, "-", MODULE, "-", PIXSCALE, "..."))
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") 
  detect_dir = paste0(ref_dir, "/ProFound/Detects/", VID, "/", MODULE, "/") 
  measure_dir = paste0(ref_dir, "/ProFound/Measurements/", VID, "/", MODULE, "/") 
  
  dir.create(
    detect_dir, recursive = TRUE, showWarnings = FALSE
  )
  dir.create(
    measure_dir, recursive = TRUE, showWarnings = FALSE
  )
  
  clean_det_stub = paste0(
    detect_dir, VID, "_", MODULE, "_", PIXSCALE, "_chunked_segstats.csv"
  )
  clean_measure_stub = paste0(
    measure_dir, VID, "_", MODULE, "_", PIXSCALE, "_chunked_photometry.csv"
  )
  
  ## Remove previous catalogues to avoid matching stuff up
  if(file.exists(clean_det_stub)){
    file.remove(clean_det_stub)
  }
  if(file.exists(clean_measure_stub)){
    file.remove(clean_measure_stub)
  }
  
  ## Combine detects
  list_detect_chunks = list.files(
    paste0(ref_dir, "/ProFound/Detects/"), 
    pattern = "chunk", recursive = TRUE, full.names = TRUE
  )
  list_detect_chunk_cats = grep(pattern = glob2rx(paste0("*", VID, "*", PIXSCALE, "*.csv")), list_detect_chunks, value = TRUE)
  list_detect_chunk_fits = grep(pattern = glob2rx(paste0("*", VID, "*", PIXSCALE, "*.fits")), list_detect_chunks, value = TRUE)

  ## Load in the big mosaic for WCS info
  file_list <- list.files(data_dir, pattern=paste0(PIXSCALE, ".fits"), full.names=TRUE) #list all the .fits files in a directory
  big_mosaicWCS = Rfits_point(
    file_list[1], ext = 'image'
  )$keyvalues
  
  ## Load in chunk info, should hopefully be the same table across all chunk frames
  chunk_list =  list.files(paste0(ref_dir, "/ProFound/Data/"), pattern = '.fits', full.names = TRUE, recursive = TRUE)
  chunk_list = grep(pattern = glob2rx(paste0("*", VID, "*chunk*", PIXSCALE, "*fits")), chunk_list, value = TRUE)
  chunk_info = Rfits_read_table(
    chunk_list[1], ext = 'chunk_info'
  )
  
  ## Load and combine the chunk catalogues
  combine_det_cat = rbindlist(
    lapply(list_detect_chunk_cats, fread), fill = TRUE
  )
  ## Make coordinates consistent with the large mosaic
  combine_det_cat$xmax_tile = combine_det_cat$xmax
  combine_det_cat$ymax_tile = combine_det_cat$ymax
  combine_det_cat$xcen_tile = combine_det_cat$xcen
  combine_det_cat$ycen_tile = combine_det_cat$ycen
  
  new_max = Rwcs_s2p(
    RA = combine_det_cat$RAmax,
    Dec = combine_det_cat$Decmax,
    keyvalues = big_mosaicWCS
  )
  combine_det_cat$xmax = new_max[,1]
  combine_det_cat$ymax = new_max[,2]
  
  new_cen = Rwcs_s2p(
    RA = combine_det_cat$RAcen,
    Dec = combine_det_cat$Deccen,
    keyvalues = big_mosaicWCS
  )
  combine_det_cat$xcen = new_cen[,1]
  combine_det_cat$ycen = new_cen[,2]

  ## Get the boundary objects since only those will be duplicated across tiles
  buffer_region_idx = unique(foreach(i = 1:dim(chunk_info)[1], .combine = "c") %do% {
    
    idx_outer=which(
      combine_det_cat$xmax >= (chunk_info$tile_xcen[i] - chunk_info$boxsize_x[i]/2.0) &
        combine_det_cat$xmax <= (chunk_info$tile_xcen[i] + chunk_info$boxsize_x[i]/2.0) &
        combine_det_cat$ymax >= (chunk_info$tile_ycen[i] - chunk_info$boxsize_y[i]/2.0) &
        combine_det_cat$ymax <= (chunk_info$tile_ycen[i] + chunk_info$boxsize_y[i]/2.0)
    )
    
    idx_inner=which(
      combine_det_cat$xmax >= (chunk_info$tile_xcen[i] - chunk_info$boxsize_no_buffer_x[i]/2.0) &
        combine_det_cat$xmax <= (chunk_info$tile_xcen[i] + chunk_info$boxsize_no_buffer_x[i]/2.0) &
        combine_det_cat$ymax >= (chunk_info$tile_ycen[i] - chunk_info$boxsize_no_buffer_y[i]/2.0) &
        combine_det_cat$ymax <= (chunk_info$tile_ycen[i] + chunk_info$boxsize_no_buffer_y[i]/2.0)
    )
    
    idx_edge = idx_outer[!(idx_outer %in% idx_inner)]
    return(
      idx_edge
    )
  })
  
  clean_cat_idx = internalclean(
    RA = combine_det_cat$RAmax[buffer_region_idx],
    Dec = combine_det_cat$Decmax[buffer_region_idx], 
    tiebreak = combine_det_cat$mag[buffer_region_idx], 
    rad = 0.01
  )
  clean_det_cat = data.frame(
    combine_det_cat[!(1:dim(combine_det_cat)[1] %in% buffer_region_idx[!(1:length(buffer_region_idx) %in% clean_cat_idx)]), ]
  )
  #duplicate_det_cat = combine_det_cat[(1:dim(combine_det_cat)[1] %in% buffer_region_idx[clean_cat_idx]), ]
  
  ## Now load the measure chunks
  list_measure_chunks = list.files(
    paste0(ref_dir, "/ProFound/Measurements/"), 
    pattern = "chunk", recursive = TRUE, full.names = TRUE
  )
  list_measure_chunk_cats = grep(pattern = glob2rx(paste0("*", VID, "*", PIXSCALE, "*.csv")), list_measure_chunks, value = TRUE)

  combine_measure_cat = rbindlist(
    lapply(list_measure_chunk_cats, fread), fill = TRUE
  )
  clean_measure_cat = data.frame(
    combine_measure_cat[!(1:dim(combine_measure_cat)[1] %in% buffer_region_idx[!(1:length(buffer_region_idx) %in% clean_cat_idx)]), ]
  )
  
  fwrite(
    clean_det_cat,
    clean_det_stub
  )
  fwrite(
    clean_measure_cat,
    clean_measure_stub
  )
  
}

## python codes. Super sketchy :P
query_gaia = function(input_args){
  
  message("Running query_gaia")
  
  VID = input_args$VID
  system( paste0("python3 query_gaia.py --VID ", VID ) )
  return(NULL)
}
query_hst = function(input_args){
  
  message("Running query_hst")
  
  VID = input_args$VID
  MODULE = input_args$MODULE
  PIXSCALE = input_args$PIXSCALE
  ref_dir = input_args$ref_dir
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  if(VID == MODULE){
    
    message("VID [", VID, "] is the same as MODULE [", MODULE, "]")
    message("Detected big mosaic!")
    
    mosaic_files = list.files(
      paste0(ref_dir, "/Mosaic_Stacks/Patch/"),
      full.names = T,
      pattern = glob2rx(paste0("*", VID, "*", PIXSCALE, "*.fits"))
    )
    input_info = foreach(i = 1:length(mosaic_files), .combine = "bind_rows")%do%{
      info = Rfits_read_table(filename = mosaic_files[i], ext = 7)
      input_VID = na.exclude(
        str_extract(
          sapply(
            info$stub, function(x)str_split_1(x, pattern = "_")
          ),
          "\\b[0-9].{1,9}\\b"
        )
      )
      input_MODULE = na.exclude(
        str_extract(
          sapply(
            info$stub, function(x)str_split_1(x, pattern = "_")
          ),
          "\\bNRCA|NRCB\\b"
        )
      )
      df = data.frame(
        "VID" = input_VID,
        "MODULE" = input_MODULE
      )
    }
    input_info = unique(input_info)
    
    for(j in 1:dim(input_info)[1]){
      system( paste0("python3 request_hap.py --VID ", input_info$VID[j], " --MODULE ", input_info$MODULE[j]) )
    }
    copy_hst_for_tile(input_args)
  }else{
    system( paste0("python3 request_hap.py --VID ", VID, " --MODULE ", MODULE) )
  }
  return(NULL)
}


copy_long = function(input_args){
  
  ## Copy the long pixel scale from the PROCESS dirs to the DATA dir
  ## Keep separate DATA dir incase the user puts their own frames in there
  ## Will be data redundancy but I'm assuming you will have enough disk space if you're using this code!
  
  message("Running copy_long")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  message("Start copy!")
  
  patch_stack_dir = paste0(ref_dir, "/Patch_Stacks/")
  mosaic_patch_dir = paste0(ref_dir, "/Mosaic_Stacks/Patch/")
  warp_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/")
  
  if(dir.exists(warp_dir)){
    warp_propane_files = list.files(
      path = warp_dir,
      pattern = ".fits$",
      full.names = T
    )
    warp_propane_files = warp_propane_files[
      !grepl("hst", warp_propane_files, ignore.case = T)
    ]
    file.remove(
      warp_propane_files
    )
  }else{
    dir.create(warp_dir, recursive = T)
  }
  
  ## here we read in the Patched_Stacks
  filepaths <- c(
    list.files(patch_stack_dir, pattern=".fits", full.names=TRUE), #list all the .fits files in a directory
    list.files(mosaic_patch_dir, pattern=".fits", full.names = TRUE)
  )
  filepaths = filepaths[!is.na(filepaths)]
  fitsnames <- c(
    list.files(patch_stack_dir, pattern=".fits", full.names=FALSE), #list all the .fits files in a directory
    list.files(mosaic_patch_dir, pattern=".fits", full.names = FALSE)
  )
  fitsnames = fitsnames[!is.na(fitsnames)]
  
  #native scales
  long_idx = grepl(VID, filepaths) & grepl(MODULE, filepaths)  & grepl("long", filepaths)
  
  #frames_long = lapply(filepaths[long_idx], function(x) Rfits_read_all(filename = x, pointer = F))
  names_frames_long = fitsnames[long_idx]
  
  for(i in 1:length(names_frames_long)){
    message(paste0("Copying ", names_frames_long[i], " to ", warp_dir))
    file.copy(from = filepaths[long_idx][i], to = paste0(warp_dir, "warp_", names_frames_long[i]), overwrite = T)
  }
  message("Done!")
  # return(c(frames_short_to_long, frames_long))
} ##<-- Copy the long pixelscale to Data dir. File redundancy but safer option for detects.
