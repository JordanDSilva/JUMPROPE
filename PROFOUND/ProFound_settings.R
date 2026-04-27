
detect_bands_load = "ALL" 
## Use "ALL" fo an all filter stack for source detection
## Use "F277W|F356W|F444W" notation to change to 3 filter long wavelength detect 

r_aperture_photometry = c(0.16, 0.32, 0.5, 0.6)
## radii in asec for multiband forced circular aperture photometry

profound_detect_master = function(frame, skyRMS, star_mask, pix_mask = NULL, segim = NULL){
  ## Inject this function into the detect runs
  ## frame is the stack with WCS
  ## stack is for the skyRMS matrix
  ## star_mask should be 1 or 0 matrix

  if(is.null(pix_mask)){
    mask = is.na(frame$imDat) | star_mask
  }else{
    mask = is.na(frame$imDat) | star_mask | pix_mask
  }
  pro = profoundProFound(
    image = frame,
    mask = mask | (frame$imDat == 0) | (is.na(frame$imDat) | (is.infinite(frame$imDat))),
    segim = segim,
    rem_mask = TRUE,
    magzero = 23.9,
    
    skyRMS = skyRMS,
    
    skycut = 1.5,
    pixcut = 7.0,
    ext = 1.0,
    
    smooth = TRUE,
    sigma = 0.8,
    
    tolerance = 3.0,
    reltol = -1.0,
    cliptol = 100,
    
    #size = 13,
    iters = 4,
    
    box = 100,
    grid = 100,
    boxiters = 0,
    
    roughpedestal = FALSE,
    pixelcov = FALSE,
    boundstats = TRUE,
    nearstats = TRUE, 
    redosky = TRUE,
    # redoskysize = 11,

    fluxtype = "microjansky",
    pixscale = pixscale(frame, unit = "asec"),
    
    verbose = FALSE,
  )
  message("Finished profound")
  return(pro)
}

profound_measure_master = function(super_img, inVar, segim, mask, redosegim = TRUE, sky = NULL, redosky = TRUE){
  ## Do multiband measurements with ProFound
  ## set sky = 0 and redosky = FALSE if frames are already sky subtracted and you don't want another sky model

  if(inherits(inVar, 'Rfits_image')){
    skyRMS = inVar$imDat^-0.5
  }else if(is.null(inVar)){
    skyRMS = NULL
  }else{
    skyRMS = inVar^-0.5
  }
  
  if(redosegim){
    iters = 2
  }else{
    iters = 0
  }
  
  super_pro = profoundProFound(
    super_img,
    magzero = 23.9,
    mask = mask | super_img$imDat == 0 | is.na(super_img$imDat) | is.infinite(super_img$imDat),
    # detection, segmentation and dilation
    segim = segim,
    redosegim = redosegim,
    
    size = 5,
    iters = iters,
    
    # sky estimate
    sky = NULL,        # Do we model the sky again?
    skyRMS = skyRMS,
    box = 100,
    grid = 100,
    boxiters = 2,
    # measurement mode
    dotot = TRUE,
    docol = TRUE,
    dogrp = TRUE,
    # stats
    boundstats = TRUE,
    segstats = TRUE,
    # misc
    pixelcov = TRUE,
    rem_mask = TRUE,
    redosky = redosky,
    redoskysize = 11,
    
    # image units
    fluxtype = "microjansky",
    pixscale = pixscale(super_img, unit = "asec"),

    verbose = FALSE,
  )
  return(super_pro)
}