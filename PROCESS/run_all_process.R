#######################
## Script to run all ##
#######################

## Make sure to edit the proper things in the "main" function

library(magicaxis)
library(Rwcs)
library(Rfits)
library(ProPane)
library(ProFound)
library(stringr)
library(utils)
library(Cairo)
library(data.table)
library(foreach)
library(doParallel)
library(dplyr)

source("all_codes_process.R")
source("initialise_variables.R")

main = function(){
  args = commandArgs(trailingOnly = T)
  if (length(args)==0) {
    message("Specify VID")
    message("Running all")
    VID = ""
    FILT = ""
  } else if (length(args)==1) {
    VID = toString(args[1])
    FILT = ""
  } else if (length(args)==2) {
    VID = toString(args[1])
    FILT = toString(args[2])
  } else {
    message(paste0("Error Nargs wrong: ", length(args)))
    q()
  }
  
  additional_params_var = initialise_params()
  list2env(additional_params_var, .GlobalEnv)
  
  env_var = Sys.getenv(
    x = c('JUMPROPE_RAW_DIR', 
          'JUMPROPE_do_NIRISS', 
          'JUMPROPE_do_MIRI',
          'JUMPROPE_cores_pro',
          'JUMPROPE_cores_stack',
          'JUMPROPE_tasks_stack')
  )
  
  if(any(env_var == "")){
    message("Please set the env variables.")
    q()
  }
  
  do_NIRISS = ifelse(grepl("T|TRUE|True", env_var['JUMPROPE_do_NIRISS']), T, F) 
  do_MIRI = ifelse(grepl("T|TRUE|True", env_var['JUMPROPE_do_MIRI']), T, F) 
  
  cores_pro = as.integer(env_var['JUMPROPE_cores_pro'])
  cores_stack = as.integer(env_var['JUMPROPE_cores_stack'])
  tasks_stack = as.integer(env_var['JUMPROPE_tasks_stack'])
  
  dir_raw = env_var['JUMPROPE_RAW_DIR']
  
  ref_dir = Sys.getenv(x = 'JUMPROPE_REF_DIR')

  if(ref_dir == ""){
    ref_dir = make_directory_structure()
  }
  
  if(do_NIRISS){
      FILT = "CLEAR|NIS"
    }
  
  cal_sky_info_save_dir = paste0(ref_dir, "/Pro1oF/")
  Pro1oF_dir = paste0(ref_dir, "/Pro1oF/cal/")
  cal_sky_dir = paste0(ref_dir, "/Pro1oF/cal_sky/")
  cal_sky_renorm_dir = paste0(ref_dir, "/Pro1oF/cal_sky_renorm/")
  sky_pro_dir = paste0(ref_dir, "/sky_pro/")
  sky_frames_dir = paste0(ref_dir, "/sky_pro/sky_frames/")
  sky_super_dir = paste0(ref_dir, "/sky_pro/sky_super/")
  dump_dir = paste0(ref_dir, "/dump/")
  invar_dir = paste0(ref_dir, "/InVar_Stacks/")
  median_dir = paste0(ref_dir, "/Median_Stacks/")
  patch_dir = paste0(ref_dir, "/Patch_Stacks/")
  
  raw_files = load_raw_files(dir_raw = dir_raw)
  
  VID_list = ifelse(
    VID == "",
    c(""),
    unlist(strsplit(VID, "|", fixed = T))
  )
  
  for(VID in VID_list){
    
    input_args = list(
      filelist = raw_files,
      additional_params = additional_params,
      
      ref_dir = ref_dir,
      
      Pro1oF_dir = Pro1oF_dir,
      cal_sky_dir = cal_sky_dir,
      cal_sky_renorm_dir = cal_sky_renorm_dir,
      cal_sky_info_save_dir = cal_sky_info_save_dir,
      
      sky_frames_dir = sky_frames_dir,
      sky_super_dir = sky_super_dir,
      sky_pro_dir = sky_pro_dir,
      
      invar_dir = invar_dir, 
      median_dir = median_dir, 
      patch_dir = patch_dir,
      dump_dir = dump_dir, 
      
      magzero = 23.9,
      
      VID = VID,
      FILT = FILT,
      
      do_NIRISS = do_NIRISS,
      do_MIRI = do_MIRI,
      
      cores_pro = cores_pro,
      cores_stack = cores_stack,
      tasks_stack = tasks_stack,
  
      SIGMA_LO = NULL #keep blurring at wisp rem stage off by default
    )
   
    if(length(args) <= 1){
      input_args$FILT = 'F277W|F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430M|F460M|F480M|F323N|F405N|F466N|F470N' ## Process long wavelengths first

    }
    
    do_1of(input_args)
    do_cal_process(input_args)
    do_regen_sky_info(input_args)
    do_super_sky(input_args)
    do_apply_super_sky(input_args)
    do_modify_pedestal(input_args)
    do_cal_sky_info(input_args)
    do_gen_stack(input_args)
  
    if(input_args$do_NIRISS){
      input_args$FILT = "CLEAR"
      do_patch(input_args)
    }else if(input_args$do_MIRI){
      do_patch(input_args)
    }else{
      if(length(args) <= 1){
        input_args$FILT = "F070W|F090W|F115W|F150W|F200W|F140M|F162M|F182M|F210M|F164N|F187N|F212N"
      }
      
      ## Test for any short frames to process 
      # sky_info = fread(
      #   paste0(
      #     input_args$sky_pro, "/sky_info.csv"
      #   )
      # )
      
      frame_info = data.frame(Rfits_key_scan(
        filelist = raw_files, 
        keylist = c('FILTER'), 
        cores = cores_pro
      ))
      
      if( length(grep(input_args$FILT, frame_info$FILTER)) > 0 ){
        do_wisp_rem(input_args)
        
        wisp_fix_files = load_raw_files(dir_raw = dir_raw)
        input_args$filelist = wisp_fix_files
        
        do_1of(input_args)
        do_cal_process(input_args)
        do_regen_sky_info(input_args)
        do_super_sky(input_args)
        do_apply_super_sky(input_args)
        do_modify_pedestal(input_args)
        do_cal_sky_info(input_args)
        do_gen_stack(input_args)
      }
      
      input_args$FILT = ""
      do_patch(input_args)
      do_RGB(input_args)
    }
  }

}

if(sys.nframe()==0){
  main()
}
