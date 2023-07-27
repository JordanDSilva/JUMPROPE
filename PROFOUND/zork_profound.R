source("all_codes_profound.R")

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

select_code_func = function(){
  message("
          ############################################################
          ############################################################
          ## Thank you for choosing ZORK for your JWST processing!  ##
          ## Press:                                                 ##
          ## 1 = Warp short to long                                 ##
          ## 2 = Query GAIA                                         ## 
          ## 3 = Star masks                                         ##
          ## 4 = ProDetect                                          ##
          ## 5 = Query HST                                          ##
          ## 6 = Warp stack HST                                     ##
          ## 7 = ProMeasure                                         ##
          ##                                                        ##
          ## CONTROL + /\ to EXIT                                   ##
          ############################################################
          ############################################################"
  )
  select_code = toString(readLines("stdin", n=1))
  select_vector = as.integer(
    strsplit(
      select_code, ","
    )[[1]]
  )

  if( sum(select_vector %in% 1:7) == 0){
    message("Oops, I think you made a mistake. Trying again.")
    select_code_func()
  }
  else{
    return(c(1:7)[1:7 %in% select_vector])
  }
  
}


main = function(){
  
  args = commandArgs(trailingOnly = T)
  if (length(args)==0) {
    message("Specify VID")
    VID = ""
  } else if (length(args)==1) {
    VID = toString(args[1])
    MODULE = "NRCA|NRCB"
  } else if (length(args)==2) {
    VID = toString(args[1])
    MODULE = toupper(toString(args[2]))
  } else {
    message(paste0("Error Nargs wrong: ", length(args)))
    q()
  }
  
  env_var = Sys.getenv(
    x = c('JUMPROPE_RAW_DIR', 
          'JUMPROPE_do_NIRISS', 
          'JUMPROPE_cores_pro',
          'JUMPROPE_cores_stack',
          'JUMPROPE_tasks_stack')
  )
  
  if(any(env_var == "")){
    message("Please set the env variables.")
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
  frame_grid = frame_grid[
    grepl(VID, frame_grid$VISIT_ID) & grepl(MODULE, frame_grid$MODULE),
  ]
  
  select_code = select_code_func()

  code_organiser = list(
    'warp_short_to_long'=warp_short_to_long,
    'query_gaia'=query_gaia,
    'star_mask'=star_mask,
    'do_detect'=do_detect,
    'query_hst'=query_hst,
    'warp_stack_hst'=hst_warp_stack,
    'do_measure'=do_measure
    )
  
  for(VID in grep(VID, frame_grid$VISIT_ID, value = T)){
    for(MODULE in grep(MODULE, frame_grid$MODULE, value = T)){
      
      input_args = list(
        VID = VID,
        MODULE = MODULE,
        ref_dir = ref_dir,
        
        sampling_cores = env_var[["JUMPROPE_cores_stack"]],
        
        cores_stack = env_var[["JUMPROPE_cores_stack"]],
        profound_tweak = T
      )
      
      lapply(select_code,
             function(x){
               code_organiser[[x]](input_args)
             })
    }
  }

}

if(sys.nframe()==0){
  main()
}
