##########################
## Script to run all V2 ##
##########################

## Make sure to edit the proper things in the "main" function

library(Rfits)
library(Rwcs)
library(ProPane)
library(ProFound)
library(magicaxis)
library(data.table)
library(stringr)
library(foreach)
library(doParallel)
library(utils)
library(Cairo)

source("all_codes.R")
source("initialise_variables.R")

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

main = function(){
  args = commandArgs(trailingOnly = T)
  if (length(args)==0) {
    message("Specify VID")
    message("Running all")
    VID = ""
  } else if (length(args)==1) {
    VID = toString(args[1])
  } else if (length(args)==2) {
    VID = toString(args[1])
    FILT = toString(args[2])
  } else {
    message(paste0("Error Nargs wrong: ", length(args)))
    q()
  }
  
  keep_trend_var = initialise_params()
  list2env(keep_trend_var, .GlobalEnv)
  
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
  
  do_NIRISS = ifelse(grepl("T|TRUE|True", env_var['JUMPROPE_do_NIRISS']), T, F) 
  
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
 
  if(length(args) <= 1){
    FILT = ""
  }
  
  do_1of(filelist = raw_files, keep_trend_data = keep_trend_data, Pro1oF_dir = Pro1oF_dir, VID = VID, FILT = FILT, cores = cores_pro)
  do_cal_process(Pro1oF_dir = Pro1oF_dir, sky_frames_dir = sky_frames_dir, VID = VID, FILT = FILT, cores = cores_pro, do_NIRISS = do_NIRISS)
  do_regen_sky_info(sky_pro_dir = sky_pro_dir, cores = cores_pro)
  do_super_sky(sky_pro_dir = sky_pro_dir, VID = VID, cores = cores_pro, do_NIRISS = do_NIRISS)
  do_apply_super_sky(Pro1oF_dir = Pro1oF_dir, cal_sky_dir = cal_sky_dir, sky_pro_dir = sky_pro_dir, VID = VID, FILT = FILT, cores = cores_pro)
  do_modify_pedestal(cal_sky_dir = cal_sky_dir, cal_sky_renorm_dir = cal_sky_renorm_dir, VID = VID, FILT = FILT, cores = cores_pro, do_NIRISS = do_NIRISS)
  do_cal_sky_info(cal_sky_renorm_dir = cal_sky_renorm_dir, cal_sky_info_save_dir = cal_sky_info_save_dir, cores = cores_pro)
  do_gen_stack(VID = VID, FILT = FILT, ref_dir = ref_dir, do_niriss = do_NIRISS, magzero_out = 23.9, cores = cores_stack, tasks = tasks_stack)

  if(do_NIRISS){
  do_patch(VID, 'CLEAR', invar_dir = invar_dir, median_dir = median_dir, patch_dir = patch_dir, cores = cores_stack)
  q()
  }

  if(length(args) <= 1){
    FILT = "F070W|F090W|F115W|F150W|F200W|F140M|F162M|F182M|F210M"
  }
  
  do_wisp_rem(filelist = raw_files, VID = VID, median_dir = median_dir, cores = cores_pro)
  wisp_fix_files = load_raw_files(dir_raw = dir_raw)
  do_1of(filelist = wisp_fix_files, keep_trend_data = keep_trend_data, Pro1oF_dir = Pro1oF_dir, VID = VID, FILT = FILT, cores = cores_pro)
  do_cal_process(Pro1oF_dir = Pro1oF_dir, sky_frames_dir = sky_frames_dir, VID = VID, FILT = FILT, cores = cores_pro, do_NIRISS = do_NIRISS)
  do_regen_sky_info(sky_pro_dir = sky_pro_dir, cores = cores_pro)
  do_super_sky(sky_pro_dir = sky_pro_dir, VID = VID, cores = cores_pro, do_NIRISS = do_NIRISS)
  do_apply_super_sky(Pro1oF_dir = Pro1oF_dir, cal_sky_dir = cal_sky_dir, sky_pro_dir = sky_pro_dir, VID = VID, FILT = FILT, cores = cores_pro)
  do_modify_pedestal(cal_sky_dir = cal_sky_dir, cal_sky_renorm_dir = cal_sky_renorm_dir, VID = VID, FILT = FILT, cores = cores_pro, do_NIRISS = do_NIRISS)
  do_cal_sky_info(cal_sky_renorm_dir = cal_sky_renorm_dir, cal_sky_info_save_dir = cal_sky_info_save_dir, cores = cores_pro)
  do_gen_stack(VID = VID, FILT = FILT, ref_dir = ref_dir, magzero_out = 23.9, cores = cores_stack, tasks = tasks_stack, do_niriss = do_NIRISS)

  FILT = ""
  do_patch(VID, FILT, invar_dir = invar_dir, median_dir = median_dir, patch_dir = patch_dir, cores = cores_stack)
  do_RGB(VID = VID, patch_dir = patch_dir, ref_dir = ref_dir)

}

if(sys.nframe()==0){
  main()
}
