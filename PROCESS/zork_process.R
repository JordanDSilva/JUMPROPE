############################################
## Script to run specific process code V2 ##
############################################

## Make sure to edit the proper things in the "zork" function

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

select_code_func = function(){
  message("
          ############################################################
          ############################################################
          ## Thank you for choosing ZORK for your JWST processing!  ##
          ## Press:                                                 ##
          ## 1 = Pro1oF                                             ##
          ## 2 = Create sky                                         ##
          ## 3 = Regen sky info                                     ##
          ## 4 = Make super sky                                     ##
          ## 5 = Remove super sky                                   ##
          ## 6 = Modify pedestal                                    ##
          ## 7 = Cal sky info                                       ##
          ## 8 = Stack                                              ##
          ## 9 = Patch                                              ##
          ## 10 = RGB                                               ##
          ## 11 = Wisp remove                                       ##
          ##                                                        ##
          ##  e.g., 1,2,3,4,5,6,7,8,11 for everything + wisp rem    ##
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
  
  if( sum(select_vector %in% 1:11) == 0){
    message("Oops, I think you made a mistake. Trying again.")
    select_code_func()
  }
  else{
    return(select_vector)
  }
  
}

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


load_raw_files = function(dir_raw){
  
  cal_files = c(
    list.files(dir_raw, pattern = ".fits$", full.names = T, recursive = T),
    NULL
  )
  
  return(cal_files)
}

zork = function(){
  
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
  
  select_code = select_code_func()

  code_organiser = list(
    'do_1oF'=do_1of,
    'do_cal_process'=do_cal_process,
    'do_regen_sky_info' = do_regen_sky_info,
    'do_super_sky'=do_super_sky,
    'do_apply_super_sky' = do_apply_super_sky,
    'do_modify_pedestal' = do_modify_pedestal,
    'do_cal_sky_info' = do_cal_sky_info,
    'do_gen_stack' = do_gen_stack,
    'do_patch' = do_patch,
    'do_rgb' = do_RGB,
    'do_wisp_rem' = do_wisp_rem
  )
  
  #VID_list = ifelse(VID != "", unlist(strsplit(VID, "|", fixed = T)), list(c("")))
  if(VID != ""){
   VID_list = unlist(strsplit(VID, "|", fixed = T)) 
  }else{
   VID_list = list(c("")) 
  }
  
  for(VID in VID_list){
    input_args = list(
      filelist = raw_files,
      additional_params = additional_params,
      
      Pro1oF_dir = Pro1oF_dir,
      sky_frames_dir = sky_frames_dir,
      sky_pro_dir = sky_pro_dir,
      cal_sky_dir = cal_sky_dir,
      cal_sky_renorm_dir = cal_sky_renorm_dir,
      cal_sky_info_save_dir = cal_sky_info_save_dir,
      ref_dir = ref_dir,
      invar_dir = invar_dir, 
      median_dir = median_dir, 
      patch_dir = patch_dir,
      
      magzero = 23.9,
      
      VID = VID,
      FILT = FILT,
      
      do_NIRISS = do_NIRISS,
      do_MIRI = do_MIRI,
      
      cores_pro = cores_pro,
      cores_stack = cores_stack,
      tasks_stack = tasks_stack,
      
      SIGMA_LO = NULL
    )
    
    lapply(select_code,
           function(x){
             code_organiser[[x]](input_args)
           })
  }
}

if(sys.nframe()==0){
  zork()
}
