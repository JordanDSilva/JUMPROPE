source("all_codes_profound.R")

select_code_func = function(){
  message(" 
          ############################################################
          ############################################################
          ## Thank you for choosing ZORK for your JWST processing!  ##
          ## Press:                                                 ##
          ## 0 = Set up directory structure                         ##
          ## 1 = Warp short to long                                 ##
          ## 2 = Copy frames                                        ##
          ## 3 = Query GAIA                                         ## 
          ## 4 = Star masks                                         ##
          ## 5 = Star masks for big mosaic                          ##
          ## 6 = ProDetect                                          ##
          ## 7 = Query HST                                          ##
          ## 8 = Warp stack HST                                     ##
          ## 9 = Copy HST for big mosaic                            ##
          ## 10 = ProMeasure                                        ##
          ## 11 = Chop up frames                                    ##
          ## 12 = Combine chunks                                    ##
          ## 13 = Help                                              ##
          ##                                                        ##
          ## CONTROL + /\ to EXIT                                    ##
          ############################################################
          ## E.g., 2,3,4,6 for copying, stars, source detection     ##
          ############################################################
          "
  )
  select_code = toString(readLines("stdin", n=1))
  select_vector = as.integer(
    strsplit(
      select_code, ","
    )[[1]]
  )
  
  if( sum(select_vector %in% 0:13) == 0){
    message("Oops, I think you made a mistake. Trying again.")
    select_code_func()
  }
  return(
    select_vector+1 ## Plus one since can't index list with 0
  )
}


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
  
  select_code = select_code_func()
  
  make_dir_temp_func = function(input_args){make_directory_structure()}
  code_organiser = list(
    'make_directory_structure'=make_dir_temp_func, ## Hack to handle input_args
    'warp_short_to_long'=warp_short_to_long,
    'copy_frames'=copy_frames,
    'query_gaia'=query_gaia,
    'star_mask'=star_mask,
    'star_mask_tile'=star_mask_tile,
    'do_detect'=do_detect,
    'query_hst'=query_hst,
    'warp_stack_hst'=hst_warp_stack,
    'copy_hst_for_tile'=copy_hst_for_tile,
    'do_measure'=do_measure,
    'frame_chunker'=frame_chunker,
    'combine_chunks'=combine_chunks,
    'do_help' = do_help
  )
  
  for(i in 1:dim(frame_grid)[1]){
    input_args = list(
      VID = frame_grid$VISIT_ID[i],
      MODULE = frame_grid$MODULE[i],
      PIXSCALE = frame_grid$PIXSCALE[i],
      
      ref_dir = ref_dir,
      
      sampling_cores = as.numeric(env_var[["JUMPROPE_cores_stack"]]),
      
      cores_stack = as.numeric(env_var[["JUMPROPE_cores_stack"]])
      
    )
    if(14 %in% select_code){
      code_organiser[[14]](input_args)
      break
    }else{
      lapply(select_code,function(x){code_organiser[[x]](input_args)}) 
    }
  }
}

if(sys.nframe()==0){
  main()
}
