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

source("./ProFound_settings.R")

jumprope_version = 2.0

######################
## for testing only ##
######################
input_args = list(
  ref_dir = "/Volumes/RAIDY/JWST/",
  VID = "1176241001",
  MODULE = "NRCA",
  cores_stack = 1
)
######################


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
      pattern = glob2rx("*long*.fits"),
      full.names = T
    ),
    list.files(
      path = mosaic_patch_dir,
      pattern = glob2rx("*long*.fits"), 
      full.names = T)
  )
  
  filenames = c(
    list.files(
      path = patch_stack_dir,
      pattern = glob2rx("*long*.fits"),
      full.names = F
    ),
    list.files(
      path = mosaic_patch_dir,
      pattern = glob2rx("*long*.fits"), 
      full.names = F)
  )
  VID_list = sapply(filenames, function(x){
    split_fname = str_split(x, "_")[[1]]
    split_fname = split_fname[!(grepl("patch|stack|long|F.+W|F.+M", split_fname))]
    vid = split_fname[1]
    if(is.na(vid)){
      return("")
    }else{
      return(vid)
    }
  })
  
  MODULE_list = sapply(filenames, function(x){
    split_fname = str_split(x, "_")[[1]]
    split_fname = split_fname[!(grepl("patch|stack|long|F.+W|F.+M", split_fname))]
    module = split_fname[2]
    if(is.na(module)){
      return(split_fname[1])
    }else{
      return(module)
    }
  })
  
  frame_info = data.frame(
    filenames = filenames,
    VISIT_ID = VID_list,
    MODULE = MODULE_list
  )
  
  stack_grid = unique(frame_info[,c("VISIT_ID", "MODULE")])
  
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
      stack_grid[i, c("PROPOSAL_ID","VISIT_ID", "MODULE")],
      
      mid_info,
      
      wcs_info_trim,
      
      mid_corner_info
    )
    return(ret)
  }
  
  foo[["jumprope_version"]] = rep(jumprope_version, dim(foo)[1])
  csv_stub = paste0(ref_dir, "/ProFound/long_warp_info.csv")
  fwrite(foo, csv_stub)
} ##<--Compute the long warp frame info needed for querying GAIA and HST via MAST


warp_short_to_long = function(input_args){
  
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
  short_idx = grepl(VID, filepaths) & grepl(MODULE, filepaths)  & grepl("short", filepaths) & grepl("F070W|F090W|F115W|F150W|F200W|F140M|F162M|F182M|F210M", filepaths) 
  long_idx = grepl(VID, filepaths) & grepl(MODULE, filepaths)  & grepl("long", filepaths) & grepl("F277W|F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430|F460M|F480M", filepaths)
  
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
} ##<-- Downwarp the short wavelength, short pixel scale to long [DEPRECATED
copy_long = function(input_args){
  
  ## Copy the long pixel scale from the PROCESS dirs to the DATA dir
  ## Keep separete DATA dir incase the user puts their own frames in there
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


star_mask = function(input_args){
  
  message("Running star mask")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  
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
  
  file_list = list.files(data_dir, pattern=glob2rx(paste0("*", VID, "*", MODULE, "*long.fits")), full.names=TRUE) #list all the .fits files in a directory
  file_names = list.files(data_dir, pattern=glob2rx(paste0("*", VID, "*", MODULE, "*long.fits")), full.names=FALSE)
  
  #make star masks on the F200W filter, and if it's missing use the shortest filter
  file_idx = (
    grepl("F200W", file_list) & 
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
  
  f200w_ref = Rfits_read(filename = file_list, pointer = F)
  imdim = dim(f200w_ref$image)
  im_ra_dec = Rwcs_p2s(x = imdim[1], y = imdim[2], keyvalues = f200w_ref$image$keyvalues) #get extent of frame in RA and DEC coords
  
  #read in GAIA catalogue
  gaia_files = list.files(path = gaia_dir, pattern = glob2rx(paste0("*", VID, "*", MODULE, "*.csv")), full.names = T, recursive = T) 
  
  gaia = data.frame(
    fread(
      file = gaia_files
    )
  )
  
  delta_time = as.double(substr(f200w_ref$image$keyvalues$DATE_END, 1, 4)) - gaia$ref_epoch
  gaia[, c("xpix", "ypix")] = Rwcs_s2p(RA = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time, ## Do proper motion correction
                                       Dec = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time, ## pmra[pmdec] = mas/yr
                                       keyvalues = f200w_ref$image$keyvalues)
  
  gaia[, c("xpix_nopm", "ypix_nopm")] = Rwcs_s2p(RA = gaia$ra, 
                                                 Dec = gaia$dec, 
                                                 keyvalues = f200w_ref$image$keyvalues)
  
  gaia$ra_fix = gaia$ra + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmra), 0, gaia$pmra) * delta_time
  gaia$dec_fix = gaia$dec + 1e-3 * 1/3600 * ifelse(is.na(gaia$pmdec), 0, gaia$pmdec) * delta_time
  
  gaia_idx = gaia$xpix >= 0 & gaia$xpix <= imdim[1] & gaia$ypix >= 0 & gaia$ypix <= imdim[2] & 
    gaia$classprob_dsc_combmod_star > 0.98 & gaia$phot_bp_rp_excess_factor < 2.5
  gaia_ra_dec = gaia[gaia_idx[!is.na(gaia_idx)],]
  gaia_trim = gaia_ra_dec[order(gaia_ra_dec$phot_g_mean_mag), ]

  message(paste0("Building star mask for VID: ", VID, ", MODULE: ", MODULE))
  psf_mask_function = function(image, xcen=dim(image)[1]/2, ycen=dim(image)[2]/2, rad=dim(image)[2]/2){
    mask =
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad/4.0, axrat=1) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad/2.0, axrat=0.05, ang=90) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad, axrat=0.07, ang=60) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad, axrat=0.07, ang=-60) +
      profoundEllipseSeg(image=image, xcen=xcen, ycen=ycen, rad=rad, axrat=0.07, ang=0)
    return(mask > 0)
  }
  
  pro_stars = profoundProFound(
    image = f200w_ref$image,
    skycut = 10,
    pixcut = 21,
    tolerance = 1,
    cliptol = 50,
    rem_mask = T
  )
  
  match_gaia = coordmatch(
    coordref = gaia_trim[, c("ra_fix", "dec_fix")],
    coordcompare = pro_stars$segstats[, c("RAmax", "Decmax")]
  )
  R100_list = rep(50, dim(gaia_trim)[1])
  R100_list[match_gaia$bestmatch$refID] = pro_stars$segstats$R100[match_gaia$bestmatch$compareID]/pixscale(f200w_ref$image)
  
  tweak_gaia = propaneTweakCat(
    cat_pre_fix = gaia_trim[match_gaia$bestmatch$refID, c("xpix", "ypix")],
    cat_ref = pro_stars$segstats[match_gaia$bestmatch$compareID, c("xmax", "ymax")]
  )
  
  temp = matrix(1, 1500, 1500)
  psf_mask = psf_mask_function(temp)
  
  all_mask = profoundApplyMask(
    image = f200w_ref$image$imDat, mask = psf_mask, 
    xcen = gaia_trim$xpix, ycen = gaia_trim$ypix,
    xsize= R100_list * 8, 
    ysize = R100_list * 8
  )

  message("Star mask complete")
  star_mask_redo = list()
  star_mask_redo$mask = profoundDilate(all_mask$mask > 0, size = 3)

  message("Saving star mask")
  
  plot_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_star_mask.pdf")
  CairoPDF(plot_stub, width = 10, height = 10)
  par(mfrow = c(1,1), mar = rep(0,4), oma = rep(0,4))
  magimage(f200w_ref$image$imDat, flip = T, sparse = 1)
  magimage(profoundDilate(star_mask_redo$mask, size = 21) - star_mask_redo$mask, col = c(NA, "magenta"), add = T)
  legend(x = "topleft", paste0(VID, "_", MODULE))
  dev.off()
  
  save_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_star_mask.rds")
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
  
  if(!(grepl("NRC", MODULE, fixed = T)) & MODULE != VID){
    MODULE = paste0("NRC", MODULE)
  }
  
  if(MODULE == VID){
    message("VID [", VID, "] is the same as MODULE [", MODULE, "]")
    message("Detected big mosaic!")
    
    data_dir = paste0(ref_dir, "/ProFound/Data/", VID, "/", MODULE, "/") #directory with the propanes
    
    file_list = list.files(data_dir, pattern=".fits", full.names=TRUE) #list all the .fits files in a directory
    file_names = list.files(data_dir, pattern=".fits", full.names=FALSE)
    
    file_idx = (
      grepl("F200W", file_list) & 
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
    
    f200w_ref = Rfits_read_image(filename = file_list, ext = 1)
    keyvalues =f200w_ref$keyvalues
    input_info = Rfits_read_table(filename = file_list, ext = 7)
    
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
    
    read_star_masks = foreach(i = 1:length(input_VID))%do%{
      star_mask_file = paste0(ref_dir, "/ProFound/Star_Masks/", input_VID[i], "/", 
                              input_MODULE[i], "/", 
                              input_VID[i], "_", input_MODULE[i], "_star_mask.rds")
      
      message("Using: ", star_mask_file)
      star_mask = readRDS(star_mask_file)
      
      input_file_keyvalues = Rfits_read_image(
        filename = input_info$full[i],
        ext = 1
      )$keyvalues
      
      temp_warp = propaneWarp(
        image_in = Rfits_create_image(star_mask$mask, keyvalues = input_file_keyvalues),
        keyvalues_out = keyvalues,
        doscale = F,
        cores = 1
      )
      star_mask_warp = temp_warp$imDat
      # star_mask_warp[star_mask_warp != 0 | !is.na(star_mask_warp)] = 1
      star_mask_warp[is.na(star_mask_warp)] = 0
      return(star_mask_warp)
    }
    
    big_star_mask = Reduce("+", read_star_masks)
    big_star_mask[big_star_mask != 0] = 1
    
    message("Star mask complete")
    
    dir.create(
      paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/"),
      recursive = T, showWarnings = F
    )
    
    message("Saving star mask...")
    plot_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_star_mask.png")
    png(plot_stub, width = dim(f200w_ref)[1], height = dim(f200w_ref)[2], res = 72)
    par(mfrow = c(1,1), mar = rep(0,4), oma = rep(0,4))
    magimage(f200w_ref$imDat, flip = T, sparse = 2)
    magimage(profoundDilate(big_star_mask, size = 21) - big_star_mask, col = c(NA, "magenta"), add = T, sparse = 2)
    legend(x = "topleft", paste0(VID, "_", MODULE))
    dev.off()
    
    save_stub = paste0(ref_dir, "/ProFound/Star_Masks/", VID, "/", MODULE, "/", VID, "_", MODULE, "_star_mask.rds")
    saveRDS(object = list(mask=big_star_mask), save_stub)
    message("Done!")
    
  }
} ## <-- Create star mask for a big mosaic using per VID star masks


## ProFound source detection codes
do_detect = function(input_args, detect_bands = "ALL", profound_function = profound_detect_master){
  
  message("Running ProDetect")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
  
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
  file_list <- list.files(data_dir, pattern=".fits", full.names=TRUE) #list all the .fits files in a directory
  file_names <- list.files(data_dir, pattern=".fits", full.names=FALSE) #list all the .fits files in a directory
  
  if(detect_bands == "ALL"){
    idx = grepl(VID, file_list) & grepl(MODULE, file_list) & grepl("patch", file_list)
  }else{
    idx = grepl(VID, file_list) & grepl(MODULE, file_list)  & grepl("patch", file_list) & grepl(detect_bands, file_list)
  }
  
  propanes = lapply(file_list[idx], function(x) Rfits_read(filename = x))
  names(propanes) = file_names[idx]
  
  N_frames = length(propanes)
  message(paste0("Running ", N_frames, " frames from ", "VISIT: ", VID, ", module: ", MODULE))
  
  star_mask_file = list.files(star_mask_dir, pattern = "mask.rds", full.names = T)
  
  if(length(star_mask_file>0)){
    message(paste0("Load in star mask from: ", star_mask_file))
    star_mask = readRDS(star_mask_file)
  }else{
    star_mask = list()
    star_mask$mask = 0
  }
  
  message(paste0("Stacking detect band ", names(propanes), collapse="\n"))
  
  #we need to do a stack image first because we cannot supply a skyRMS list to profoundMultiBand
  stack_image = propaneStackWarpInVar(image_list = lapply(propanes, function(x)x$image[,]),
                                      inVar_list = lapply(propanes, function(x)x$inVar[,]),
                                      keyvalues_out = propanes[[1]]$image$keyvalues,
                                      magzero_in = 23.9, magzero_out = 23.9, cores = 1, cores_warp = 1)
  
  pix_mask = propanes[[1]]$weight[,]$imDat
  pix_mask[pix_mask <= 1]= 0
  pix_mask[pix_mask > 0]=1
  pix_mask = 1 - pix_mask
  
  profound = profound_function(frame = stack_image$image[,], 
                               skyRMS = (stack_image$inVar[,]$imDat)^-0.5, 
                               star_mask = star_mask$mask, 
                               pix_mask=pix_mask)

  profound$segstats$MODULE = rep(MODULE, dim(profound$segstats)[1])
  profound$segstats$VID = rep(VID, dim(profound$segstats)[1])
  
  plot_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_profound_plot.pdf")
  CairoPDF(plot_stub, width = 10, height = 10)
  plot(profound)
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
  
  stack_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_profound_stack.fits")
  profound_stack = list(
    stack = stack_image$image[,],
    segim = profound$segim,
    segim_orig = profound$segim_orig
  )
  Rfits_write(data = profound_stack, filename = stack_stub)
  
  catalogue_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_segstats.csv")
  fwrite(profound$segstats, file = catalogue_stub)
  
  profound_stub = paste0(detect_dir, "/", VID, "_", MODULE, "_profound.rds")
  saveRDS(profound, file = profound_stub)
  
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
  
  message("Running ProMeasure")
  
  ref_dir = input_args$ref_dir
  VID = input_args$VID
  MODULE = input_args$MODULE
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
  pro_path = list.files(detect_dir, pattern = glob2rx(paste0("*", VID, "*", MODULE, "*rds")), full.names = T)
  assert(checkFileExists(pro_path))
  
  super_pro = readRDS(pro_path)
  segim = super_pro$segim
  mask = super_pro$mask
  
  ####### Load data images ########################
  assert(checkDirectoryExists(data_dir))
  
  # nirc.list = list.files(data_path, pattern = modl, full.names = T, recursive = T)
  data.list = list.files(data_dir, 
                         pattern = ".fits$", 
                         full.names = T, 
                         recursive = T)
  data.names = list.files(data_dir, 
                          pattern = ".fits$", 
                          full.names = F, 
                          recursive = T)
  
  # data.list = data.list[-grep('jwst_nirc', data.list)]                            # (should) still work even with only NIRCam data
  # data.list = c(data.list, nirc.list)
  
  filter.names = toupper(str_extract(data.list, "F\\d{3,}[A-Z]{1,}|f\\d{3,}[a-z]{1,}"))
  filter.names[is.na(filter.names)] = sapply(data.list[is.na(filter.names)], function(x)str_split_1(x, "_")[6]) ## Filter name should hopefully always be in position 6
  
  bad_name = data.list[is.na(filter.names)]
  str_extract(bad_name,"F+")
  
  images = lapply(data.list, function(x){
    ext = Rfits_extname_to_ext(x, extname = "image")
    message(paste("Using ext=", ext))
    Rfits_read_image(x, ext = ext)
  })
  
  names(images) = filter.names                                                    # filter ordering is not important
  
  ####### Initiate for sampling error #############
  r = seq(2,10)
  master = data.frame(r=r, a=unlist(lapply(r, function(x) length(getCircle(x)))))
  
  if(!dir.exists(sampling_dir)){dir.create(sampling_dir, recursive = T)}else{
    message("Unlinking sample dir")
    unlink(sampling_dir, recursive = T)
    dir.create(sampling_dir, recursive = T)
  }                           # Create directory for the sampling stuff
  
  error_file = file.path(sampling_dir,paste(VID,MODULE,'error_fit.txt', sep='_'))
  # error_file = paste0(sampling_dir, "/", VID, "_", MODULE, "_error_fit.txt")
  # error_file = "/Volumes/RAIDY/JWST/ProFound/Sampling/2738001001/NRCA/2738001001_NRCA_error_fit.txt"
  cat('Filters','Slopes','Intercepts','\n', file = error_file, sep='\t')          # Create file for the fitted relation of error sampling
  
  ####### Initiate csv for saving output ##########
  csvout = data.frame(segID=super_pro$segstats$segID)
  
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
    message(paste0("Running ProMeasure on: ", ff))
    filt = images[[ff]]
    # filt[[ff]]$imDat[filt[[ff]]$imDat==0L] = NA
    dum_pro = measure_profound(filt, segim, mask)
    dum_pro_col = measure_profound(filt, segim, mask, redosegim = F) #don't redilate segments e.g., colour photometry mode
    csvout[paste0(ff,'_fluxt')] = dum_pro$segstats$flux
    csvout[paste0(ff,'_fluxt_err')] = dum_pro$segstats$flux_err
    csvout[paste0(ff,'_fluxc')] = dum_pro_col$segstats$flux
    csvout[paste0(ff,'_fluxc_err')] = dum_pro_col$segstats$flux_err
    csvout[paste0(ff,'_maskfrac')] = dum_pro$segstats$Nmask/dum_pro$segstats$Nedge
    
    CairoPDF(file.path(inspect_dir,paste0("profound_",ff,"_inspect.pdf")), width = 10, height = 10 )
    plot(dum_pro)
    plot(dum_pro_col)
    par(mfrow = c(1,2), mar = rep(0,4), oma = rep(0,4))
    profoundSegimPlot(dum_pro, sparse=1)
    legend(x="topleft", legend="Total photometry")
    profoundSegimPlot(dum_pro_col, sparse=1)
    legend(x="topleft", legend="Colour photometry")
    dev.off()
    
    img = filt$imDat - dum_pro$sky
    img[img == 0L] = NA
    img[mask] = NA
    sample_mask = (dum_pro$objects_redo | mask)
    
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
    
    csvout[paste0(ff,'_scaled_fluxt_err')] = error_scaling(dum_pro$segstats$flux_err, dum_pro$segstats$N100, slope, intercept)
    csvout[paste0(ff,'_scaled_fluxc_err')] = error_scaling(dum_pro_col$segstats$flux_err, dum_pro_col$segstats$N100, slope, intercept)
    saveRDS(dum_pro, file.path(measurements_dir,paste(VID,MODULE,ff,"results.rds",sep='_')))
  }
  write.csv(csvout, file = file.path(measurements_dir,paste(VID,MODULE,"photometry.csv",sep='_')))
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
  
  target_path = jwst_ls[grepl(MODULE, jwst_ls) & grepl("F200W", jwst_ls)]
  if(length(target_path)==0){
    target = Rfits_read(jwst_ls[1])
  }else{
    target = Rfits_read(target_path)
  }
  
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
                                    rad = 2.0,
                                    plot = F)
    files_temp = frame_find$full
    files_names = frame_find$file
    
    if(is.null(files_temp)){
      message("No HST to fold in!")
    }else{
      hst_key_scan = Rfits_key_scan(files_temp, keylist = c("FILTER", "DETECTOR", "INSTRUME"))
      
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
      
      hst_filters = unique(hst_key_scan$FILTER)
      NFilters = length(hst_filters)
      
      #loop through HST filters
      foreach(j = 1:NFilters)%do%{
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
          frames = Rfits_make_list(files, extlist = frame_find$ext[files_idx], pointer = F)
          
        }
        
        #prepare input frames for stacking
        magzero_list = c()
        image_list = {}
        inVar_list = {}
        for(i in 1:length(frames)){
          
          hst_magzero = frames[[i]]$keyvalues$MAGZERO
          if(is.na(hst_magzero)){
            hst_magzero = -2.5*log10(frames[[i]]$keyvalues$PHOTFLAM)-5*log10(frames[[i]]$keyvalues$PHOTPLAM)-2.408
          }
          
          magzero_list = c(magzero_list, hst_magzero)
          
          temp = frames[[i]]
          
          if(!do_sky_rem){
            message("Running: ", files[i])
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
        
        #stack
        output_stack = propaneStackWarpInVar(image_list = image_list,
                                             # inVar_list = inVar_list,
                                             magzero_in = hst_magzero,
                                             magzero_out = 23.9,
                                             keyvalues_out = target$image$keyvalues,
                                             cores = cores_stack,
                                             cores_warp = 1)
        
  
          message("Tweaking ProFound source shift")
         
          pro_test = profoundProFound(
            image = output_stack$image,
            pixcut = 20,
            skycut = 10.0,
            cliptol = 100,
            tolerance = Inf,
            box = dim(output_stack$image)/10.0,
            mask = is.infinite(output_stack$image$imDat) | is.na(output_stack$image$imDat),
            magzero = 23.9,
            rem_mask = T
          )
          
          match_coord = coordmatch(
            coordref = pro_ref$segstats[, c("RAmax", "Decmax")],
            coordcompare = pro_test$segstats[, c("RAmax", "Decmax")]
          )
          
          propane_cat_tweak = propaneTweakCat(cat_ref = cbind(pro_ref$segstats$xcen[match_coord$bestmatch$refID], 
                                                              pro_ref$segstats$ycen[match_coord$bestmatch$refID]),
                                              cat_pre_fix = cbind(pro_test$segstats$xcen[match_coord$bestmatch$compareID], 
                                                                  pro_test$segstats$ycen[match_coord$bestmatch$compareID]), 
                                              delta_max = c(10,1.0), mode = "pix")
          
          tweaked_output_imDat = propaneTran(output_stack[["image"]][,]$imDat,
                                             delta_x = propane_cat_tweak$par[1],
                                             delta_y = propane_cat_tweak$par[2],
                                             delta_rot = propane_cat_tweak$par[3])
          
          tweaked_output = Rfits_create_image(data = tweaked_output_imDat, 
                                              keyvalues = target$image$keyvalues)
          
          final_output = list()
          final_output$image = tweaked_output
          final_output$image$keyvalues$EXTNAME = "image"
          final_output$tweak_sol = propane_cat_tweak$par
          
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
        paste0(ref_dir, "/ProFound/Data/", input_info$VID[i], "/", input_info$MODULE[i], "/"),
        pattern = glob2rx("*hst*.fits"),
        full.names = T
      )
    }
    
    previous_hst_names = foreach(i = 1:dim(input_info)[1], .combine = "c")%do%{
      list.files(
        paste0(ref_dir, "/ProFound/Data/", input_info$VID[i], "/", input_info$VID[i], "/"),
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
      to = paste0(HST_cutout_dir, "/", previous_hst_names),
      overwrite = T
    )
    do_sky_rem = F
  }
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
      pattern = glob2rx(paste0("*", VID, "*long.fits"))
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


