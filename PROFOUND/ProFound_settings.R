
detect_bands_load = "F277W|F356W|F444W" 
## Use "ALL" fo an all filter stack for source detection
## Use "F277W|F356W|F444W" notation to change to 3 filter long wavelength detect 

profound_detect_master = function(frame, skyRMS, star_mask, pix_mask=NULL, segim = NULL){
  #Inject this function into the detect runs
  #frame is the stack with WCS
  #stack is for the skyRMS matrix
  #star_mask should be 1 or 0 matrix
  if(is.null(pix_mask)){
    mask = is.na(frame$imDat) | star_mask
  }else{
    mask = is.na(frame$imDat) | star_mask | pix_mask
  }
  pro = profoundProFound(
    image = frame,
    mask = mask | (frame$imDat == 0) | (is.na(frame$imDat) | (is.infinite(frame$imDat))),
    segim = segim,
    rem_mask = T,
    magzero = 23.9,
    
    skyRMS = skyRMS,
    
    skycut = 1.5,
    pixcut = 7.0,
    ext = 1.0,
    
    smooth = T,
    sigma = 0.8,
    
    tolerance = 3.0,
    reltol = -1.0,
    cliptol = 100,
    
    #size = 13,
    iters = 4,
    
    box = 100,
    grid = 100,
    boxiters = 0,
    
    roughpedestal = F,
    pixelcov = F,
    boundstats = T,
    nearstats = T, 
    redosky = T,
    # redoskysize = 11,
    
    verbose = F,
  )
  message("Finished profound")
  return(pro)
}

measure_profound = function(super_img = super_img, inVar = inVar, segim=segim, mask=mask, redosegim=T, sky=NULL, redosky=T){
  
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
    dotot = T,
    docol = T,
    dogrp = T,
    # stats
    boundstats = T,
    segstats = T,
    # misc
    pixelcov = T,
    rem_mask = T,
    verbose = F,
    redosky = redosky,
    redoskysize = 11
  )
  return(super_pro)
}
