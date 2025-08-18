source("all_codes_profound.R")

main = function(){
  
  args = commandArgs(trailingOnly = T)
  if (length(args)==0) {
    message("Specify VID")
    VID = ""
    MODULE = ""
    PIXSCALE = "long"
  } else if (length(args)==1) {
    VID = toString(args[1])
    MODULE = ""
    PIXSCALE = "long"
  } else if (length(args)==2) {
    VID = toString(args[1])
    MODULE = toupper(toString(args[2]))
    PIXSCALE = "long"
  } else if (length(args)==3) {
    VID = toString(args[1])
    MODULE = toupper(toString(args[2]))
    PIXSCALE = toString(args[3])
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
      paste0(ref_dir, "/ProFound/warp_info.csv")
    ) 
  )
  
  frame_grid_unique = unique(frame_info_file[, c("PROPOSAL_ID", "VISIT_ID", "MODULE", "PIXSCALE")])
  
  frame_grid = frame_grid_unique[
    grepl(paste0("^",VID,"$"), frame_grid_unique$VISIT_ID) & grepl(paste0("^",MODULE,"$"), frame_grid_unique$MODULE, ignore.case = TRUE) & grepl(PIXSCALE, frame_grid_unique$PIXSCALE),
  ]
  if(dim(frame_grid)[1]==0){
    frame_grid = frame_grid_unique[
      grepl(paste0("^",VID,"$"), frame_grid_unique$PROPOSAL_ID) & grepl(paste0("^",MODULE,"$"), frame_grid_unique$MODULE, ignore.case = TRUE) & grepl(PIXSCALE, frame_grid_unique$PIXSCALE),
    ]
  }
  
  message("Using this frame grid:")
  print(
    frame_grid
  )
  
  code_organiser = list(
    'copy_frames'=copy_frames,
    'query_gaia'=query_gaia,
    'star_mask'=star_mask,
    'do_detect'=do_detect,
    'query_hst'=query_hst,
    'warp_stack_hst'=hst_warp_stack,
    'do_measure'=do_measure
  )
  
  for(i in 1:dim(frame_grid)[1]){
    input_args = list(
      VID = frame_grid$VISIT_ID[i],
      MODULE = frame_grid$MODULE[i],
      PIXSCALE = frame_grid$PIXSCALE[i],
      
      ref_dir = ref_dir,
      
      sampling_cores = as.numeric(env_var[["JUMPROPE_cores_stack"]]),
      
      cores_stack = as.numeric(env_var[["JUMPROPE_cores_stack"]]),
      
      profound_chunk = profound_chunk
    )
    
    lapply(1:length(code_organiser),
           function(x){
             code_organiser[[x]](input_args)
           })
  }
  
}

if(sys.nframe()==0){
  main()
}
