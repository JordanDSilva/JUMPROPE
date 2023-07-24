make_directory_structure = function(){
  message("Should I make the directory structure for you? (T/F): ")
  make_dirs = readLines("stdin", n = 1)
  if(make_dirs == "T"){
    message("Where should I make the directory (supply directory or nothing): ")
    ref_dir = readLines("stdin", n=1)
    if(ref_dir==""){
      ref_dir = getwd()
      message("No user input. Making in current working directory.")
    }
    print(ref_dir)
    dir.create(paste0(ref_dir, "/ProFound/Data/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/GAIA_Cats/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/Star_Masks/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/HST_cutout/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/Detects/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/Sampling/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/Inspect/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/ProFound/Measurements/"), recursive = T, showWarnings = F)
    return(ref_dir)
  }else{
    message("Continuing...")
    message("Where is the reference directory (supply directory or nothing): ")
    ref_dir = readLines("stdin", n=1)
    if(ref_dir==""){
      ref_dir = getwd()
      message("No user input. Assuming everything in current working directory.")
    }
    return(ref_dir)
  }
}

frame_info = function(ref_dir){
  
  frame_info = fread(
    paste0(ref_dir, "/Pro1oF/cal_sky_info.csv")
  )
  
  patch_stack_dir = paste0(ref_dir, "/Patch_Stacks/")
  filelist = list.files(
    path = patch_stack_dir,
    pattern = glob2rx("*long*.fits"),
    full.names = T
  )
  filenames = list.files(
    path = patch_stack_dir,
    pattern = glob2rx("*long*.fits"),
    full.names = F
  )
  VID_list = substr(filenames, 13, 22)
  MODULE_list = substr(filenames, 30, 33)
  
  frame_info = data.frame(
    filenames = filenames,
    VID = VID_list,
    MODULES = MODULE_list
  )
  
  stack_grid = unique(frame_info[,c("VID", "MODULES")])
  
  foo = foreach(i = 1:dim(stack_grid)[1], .combine = "rbind") %do% {
    wcs_info = propaneFrameFinder(
      rad = 2*pi*3600,
      filelist = filelist[
        grepl(stack_grid$VID[i], filelist) & grepl(stack_grid$MODULES[i], filelist)
          ][1],
      plot = F
    )
    
    mid_BL_RA = (wcs_info$centre_RA + wcs_info$corner_BL_RA)/2.0
    mid_TL_RA = (wcs_info$centre_RA + wcs_info$corner_TL_RA)/2.0
    mid_BR_RA = (wcs_info$centre_RA + wcs_info$corner_BR_RA)/2.0
    mid_TR_RA = (wcs_info$centre_RA + wcs_info$corner_TR_RA)/2.0
    
    mid_BL_Dec = (wcs_info$centre_Dec + wcs_info$corner_BL_Dec)/2.0
    mid_TL_Dec = (wcs_info$centre_Dec + wcs_info$corner_TL_Dec)/2.0
    mid_BR_Dec = (wcs_info$centre_Dec + wcs_info$corner_BR_Dec)/2.0
    mid_TR_Dec = (wcs_info$centre_Dec + wcs_info$corner_TR_Dec)/2.0
    
    ret = cbind(
      cbind("mid_BL_RA" = mid_BL_RA, "mid_BL_Dec" = mid_BL_Dec),
      cbind("mid_TL_RA" = mid_TL_RA, "mid_TL_Dec" = mid_TL_Dec),
      cbind("mid_TR_RA" = mid_TR_RA, "mid_TR_Dec" = mid_TR_Dec),
      cbind("mid_BR_RA" = mid_BR_RA, "mid_BR_Dec" = mid_BR_Dec),

      wcs_info[, c("centre_RA", "centre_Dec", 
                 "corner_BL_RA", "corner_BL_Dec",
                 "corner_TL_RA", "corner_TL_Dec",
                 "corner_TR_RA", "corner_TR_Dec",
                 "corner_BR_RA", "corner_BR_Dec")]
    )
    
    # magplot(unlist(ret[, grepl("RA", names(ret))]), unlist(ret[, grepl("Dec", names(ret))]))
    return(ret)
  }
  
}

warp_short_to_long = function(ref_dir, VID, module_label){
  ## Going to warp the short to the long 
  message("Start warp!")
  
  patch_stack_dir = paste0(ref_dir, "/Patch_Stacks/")
  warp_dir = paste0(ref_dir, "/ProFound/Data/")
  
  ## here we read in the Patched_Stacks
  filepaths <- list.files(patch_stack_dir, pattern=".fits", full.names=TRUE) #list all the .fits files in a directory
  fitsnames <- list.files(patch_stack_dir, pattern=".fits", full.names=FALSE) #list all the .fits files in a directory
  
  #native scales
  short_idx = grepl(VID, filepaths) & grepl(module_label, filepaths)  & grepl("short", filepaths) & grepl("F070W|F090W|F115W|F150W|F200W|F140M|F162M|F182M|F210M", filepaths) 
  long_idx = grepl(VID, filepaths) & grepl(module_label, filepaths)  & grepl("long", filepaths) & grepl("F277W|F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430|F460M|F480M", filepaths)
  
  message(paste0(fitsnames[short_idx], collapse = "\n"))
  message(paste0(fitsnames[long_idx], collapse = "\n"))
  
  frames_short = lapply(filepaths[short_idx], function(x) Rfits_read_all(filename = x, pointer = F)) #to be warped
  names(frames_short) = fitsnames[short_idx]
  framess_long = lapply(filepaths[long_idx], function(x) Rfits_read_all(filename = x, pointer = F)) #reference warp
  names(promos_long) = fitsnames[long_idx]
  
  frames_short_to_long = lapply(frames_short, function(x){
    warp = propaneWarpProPane(x, keyvalues_out = frames_long[[1]]$image$keyvalues)
    warp$info = x$info
    warp
  })
  
  names(frames_short_to_long) = str_replace(names(frames_short), "short", "long")
  for(i in 1:length(frames_short_to_long)){
    class(frames_short_to_long[[i]]) = "Rfits_list"
  }
  message(paste0("Finished warping: ", names(frames_short_to_long), sep = "\n"))
  N_frames_short = length(frames_short)
  
  for(i in 1:length(frames_short_to_long)){
    message(paste0("Saving ", names(frames_short_to_long)[i], " to ", dataDir))
    Rfits_write(data = frames_short_to_long[[i]], filename = paste0(warp_dir, "warp_", names(frames_short_to_long)[i]))
  }
  for(i in 1:length(frames_long)){
    message(paste0("Copying ", names(frames_long)[i], " to ", dataDir))
    file.copy(from = filepaths[long_idx][i], to = paste0(warp_dir, "warp_", names(frames_long)[i]), overwrite = T)
  }
  message("Done!")
  return(c(frames_short_to_long, frames_long))
}

