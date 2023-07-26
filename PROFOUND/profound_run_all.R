source("profound_all_codes.R")

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
    dir.create(paste0(ref_dir, "/Pro1oF/cal/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/Pro1oF/cal_sky/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/Pro1oF/cal_sky_renorm/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/sky_pro/sky_frames/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/sky_pro/sky_super/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/dump/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/InVar_Stacks/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/Median_Stacks/"), recursive = T, showWarnings = F)
    dir.create(paste0(ref_dir, "/Patch_Stacks/"), recursive = T, showWarnings = F)
    
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

main = function(){
  
  args = commandArgs(trailingOnly = T)
  if (length(args)==0) {
    message("Specify VID")
    VID = ""
  } else if (length(args)==1) {
    VID = toString(args[1])
  } else if (length(args)==2) {
    VID = toString(args[1])
    MODULE = toupper(toString(args[2]))
  } else {
    message(paste0("Error Nargs wrong: ", length(args)))
    q()
  }
  
  ref_dir = Sys.getenv(
    x = "JUMPROPE_REF_DIR"
  )
  
  if(ref_dir == ""){
    ref_dir = make_directory_structure()
  }
  
  frame_info(ref_dir)
  
  frame_info_file = data.frame(
    fread(
      paste0(ref_dir, "/ProFound/long_warp_info.csv")
    ) 
  )
  
  frame_grid = unique(frame_info_file[, c("VISIT_ID", "MODULE")])
  
  for(VID in grep(VID, frame_grid$VISIT_ID, value = T)){
    for(MODULE in grep(MODULE, frame_grid$MODULE, value = T)){
      # warp_short_to_long(ref_dir = ref_dir, VID = VID, MODULE = MODULE)
      
      system( paste0("python3 query_gaia.py --PID ", substr(VID, 1, 4)) )
      
      star_mask(ref_dir, VID, MODULE)
      
      do_detect(ref_dir, VID, MODULE)
      
    }
  }
  
}

if(sys.nframe()==0){
  main()
}
