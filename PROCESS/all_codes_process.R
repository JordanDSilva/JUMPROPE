## Required packages
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
library(imager)
library(celestial)
library(matrixStats)

pipe_version = "1.2.1" ## Change nominal from too high version 2.0 (1.0.0 being release on GitHub)

load_files = function(input_args, which_module, sky_info = NULL){
  ## Load the correct files for what ever task
  ## Correctly place files into whatever module
  ## Which module tells me which step to load files for
  
  VID = input_args$VID
  FILT = input_args$FILT
  
  Pro1oF_dir = input_args$Pro1oF_dir
  sky_pro_dir = input_args$sky_pro_dir
  cal_sky_dir = input_args$cal_sky_dir
  
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  
  if(which_module == "1oF"){
    files_1oF = input_args$filelist
    files_1oF = files_1oF[grepl(VID, files_1oF) & grepl(".fits$", files_1oF)]
    
    scan_1oF = Rfits_key_scan(filelist = files_1oF, keylist = c("FILTER", "PROGRAM", "VISIT_ID"))
    corr_pid = grep(VID, scan_1oF$VISIT_ID, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID] ## Make sure 4 digit PID is embedded in 10 digit VID
    pid_idx = scan_1oF$VISIT_ID %in% corr_pid
    
    files_1oF = files_1oF[pid_idx]
    files_1oF = files_1oF[grepl(FILT, scan_1oF$FILTER[pid_idx])]
    return(files_1oF)
  }
  
  if(which_module == "cal_process"){
    files_cal = c(
      list.files(Pro1oF_dir, 
                 full.names = TRUE, 
                 pattern = glob2rx(paste0("*",VID,"*.fits")))
    )
    if(do_NIRISS){
      files_cal = files_cal[!grepl('_mirimage_',files_cal) & grepl('_nis_',files_cal) & grepl(".fits$", files_cal)]
    }else if (do_MIRI){
      files_cal = files_cal[grepl('_mirimage_',files_cal) & !grepl('_nis_',files_cal) & grepl(".fits$", files_cal)]
    }else{
      files_cal = files_cal[!grepl('_mirimage_',files_cal) & !grepl('_nis_', files_cal) & grepl(".fits$", files_cal)]
    }
    
    scan_cal = Rfits_key_scan(filelist = files_cal,keylist = c("FILTER", "PROGRAM", "VISIT_ID"))
    corr_pid = grep(VID, scan_cal$VISIT_ID, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID] ## Make sure 4 digit PID is embedded in 10 digit VID
    pid_idx = scan_cal$VISIT_ID %in% corr_pid
    
    files_cal = files_cal[pid_idx]
    files_cal = files_cal[grepl(FILT, scan_cal$FILTER[pid_idx])]
    return(files_cal)
  }
  
  if(which_module == "super_sky"){
    sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
    sky_super_dir = paste0(sky_pro_dir, "/sky_super/")
    files_sky = list.files(sky_frames_dir, 
                           full.names = TRUE, 
                           pattern = glob2rx(paste0("*", VID, "*.fits")))

    if(do_NIRISS){
      files_sky = files_sky[!grepl('_mirimage_',files_sky) & grepl('_nis_',files_sky) & grepl(".fits$", files_sky)]
    }else if (do_MIRI){
      files_sky = files_sky[grepl('_mirimage_',files_sky) & !grepl('_nis_',files_sky) & grepl(".fits$", files_sky)]
    }else{
      files_sky = files_sky[!grepl('_mirimage_',files_sky) & !grepl('_nis_', files_sky) & grepl(".fits$", files_sky)]
    }
    
    scan_sky = Rfits_key_scan(filelist = files_sky, keylist = c("FILTER", "PROGRAM", "VISIT_ID"))
    corr_pid = grep(VID, scan_sky$VISIT_ID, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID] ## Make sure 4 digit PID is embedded in 10 digit VID
    pid_idx = scan_sky$VISIT_ID %in% corr_pid
    
    files_sky = files_sky[pid_idx]
    files_sky = files_sky[grepl(FILT, scan_sky$FILTER[pid_idx])]
    return(files_sky)
  }
  
  if(which_module == "apply_super"){
    
    corr_pid = grep(VID, sky_info$visit_id, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID] ## Make sure 4 digit PID is embedded in 10 digit VID
    pid_idx = paste0(sky_info$visit_id) %in% corr_pid
    
    sky_info = sky_info[grepl(VID, sky_info$fileim) & pid_idx & grepl(FILT, sky_info$filter), ]
    
    if(do_NIRISS){
      sky_info = sky_info[grepl("NIS", sky_info$detector), ]
    }else if (do_MIRI){
      sky_info = sky_info[grepl("MIRIMAGE", sky_info$detector), ]
    }
    
    return(list('sky_info' = sky_info, 'sky_filelist', sky_info$filesky))
  }
  
  if(which_module == "modify_pedestal"){
    files_cal_sky = c(
      list.files(cal_sky_dir, full.names = TRUE, pattern = glob2rx(paste0("*", VID, "*.fits"))) #If running Pro1/F outputs
    )
    
    if(do_NIRISS){
      files_cal_sky = files_cal_sky[!grepl('_mirimage_',files_cal_sky) & grepl('_nis_',files_cal_sky) & grepl(".fits$", files_cal_sky)]
    }else if (do_MIRI){
      files_cal_sky = files_cal_sky[grepl('_mirimage_',files_cal_sky) & !grepl('_nis_',files_cal_sky) & grepl(".fits$", files_cal_sky)]
    }else{
      files_cal_sky = files_cal_sky[!grepl('_mirimage_',files_cal_sky) & !grepl('_nis_', files_cal_sky) & grepl(".fits$", files_cal_sky)]
    }
    
    scan_cal_sky = Rfits_key_scan(filelist = files_cal_sky, keylist = c("FILTER", "PROGRAM", "VISIT_ID"))
    corr_pid = grep(VID, scan_cal_sky$VISIT_ID, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID] ## Make sure 4 digit PID is embedded in 10 digit VID
    pid_idx = scan_cal_sky$VISIT_ID %in% corr_pid
    
    files_cal_sky = files_cal_sky[pid_idx]
    files_cal_sky = files_cal_sky[grepl(FILT, scan_cal_sky$FILTER[pid_idx])]
    return(files_cal_sky)
  }
  
  if(which_module == "wisp_rem"){
    files_wisp = input_args$filelist
    files_wisp = files_wisp[grepl(VID, files_wisp) & grepl(".fits$", files_wisp)]
    
    scan_1oF = Rfits_key_scan(filelist = files_wisp,keylist = c("FILTER", "PROGRAM", "VISIT_ID"))
    corr_pid = grep(VID, scan_1oF$VISIT_ID, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID]
    pid_idx = scan_1oF$VISIT_ID %in% corr_pid
    
    files_wisp = files_wisp[pid_idx]
    return(files_wisp)
  }
  
}

do_1of = function(input_args){
  
  cat("\n")
  message("## Removing 1/f ##")
  cat("\n")
  Sys.sleep(time = 5)

  additional_params = input_args$additional_params
  filelist = input_args$filelist
  Pro1oF_dir = input_args$Pro1oF_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  
  do_MIRI = input_args$do_MIRI ## Because we need to alter the dimensions of the image for processing
  
  ### EXAMPLE BITS TO EDIT
  
  #very large sources
  # 1176341001 A
  # 1176341001 B
  trend_block_vlarge = 101 #good for larger sources
  
  #large sources
  # 1176211001 B
  # 1176241001 B
  # 2282010001 B
  # 2736001001 B
  trend_block_large = 501 # good for frames with fairly large sources
  
  ID_vlarge = additional_params$ID_vlarge
  ID_large = additional_params$ID_large
  
  ow_vlarge = additional_params$ow_vlarge
  ow_large = additional_params$ow_large
  
  ### BITS TO EDIT END ###
  
  registerDoParallel(cores=cores)
  
  filelist = load_files(input_args, which_module = "1oF")
  
  message("Showing first <10 files:")
  cat(head(filelist, 10), sep='\n')
  cat("...")
  cat('Processing',length(filelist),'files\n')
  
  lo_loop = 1
  hi_loop = length(filelist)
  
  dummy = foreach(i = lo_loop:hi_loop, .inorder = FALSE)%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    file_image = filelist[i]
    basename = strsplit(basename(file_image),'.fits$')[[1]]
    fullbase = paste(Pro1oF_dir,basename,sep='/')
    
    temp_image = Rfits_read(filelist[i], pointer=FALSE)
    
    if(do_MIRI){
      imdim = dim(temp_image$SCI)
      xpix = 360:imdim[1]
      # xpix = 1:imdim[1]
      ypix = 1:imdim[2]
      
      new_centre = suppressMessages(Rwcs_p2s(
        x = floor(median(xpix)),
        y = floor(median(ypix)),
        keyvalues = temp_image$SCI$keyvalues
      ))
      
      dq = temp_image$DQ$imDat
      
      temp_image$SCI = temp_image$SCI[xpix, ypix]
      temp_image$ERR$imDat = temp_image$ERR[xpix, ypix]$imDat
      temp_image$DQ$imDat = dq[xpix, ypix]
      temp_image$AREA$imDat = temp_image$AREA[xpix, ypix]$imDat
      temp_image$VAR_POISSON$imDat = temp_image$VAR_POISSON[xpix, ypix]$imDat
      temp_image$VAR_RNOISE$imDat = temp_image$VAR_RNOISE[xpix, ypix]$imDat
      temp_image$VAR_FLAT$imDat = temp_image$VAR_FLAT[xpix, ypix]$imDat
      
      for(ext in c("ERR", "DQ", "AREA", "AREA", "VAR_POISSON", "VAR_RNOISE", "VAR_FLAT")){
        temp_image[[ext]]$keyvalues$NAXIS1 = imdim[1] - 360
      }
      
      keep_trend = FALSE
      if(ow_large & ow_vlarge){
        message("Can't use both LARGE and VLARGE keep_trend \n Defaulting to keep_trend FALSE")
        keep_trend = FALSE
      }else{
        
        if(ow_vlarge){
          keep_trend = TRUE
          trend_block = trend_block_vlarge
          message("KT: TRUE", " TB ", trend_block_vlarge)
        }
        
        if(ow_large){
          keep_trend = TRUE
          trend_block = trend_block_large
          message("KT: TRUE", " TB ", trend_block_large)
        }
        
      }
      
      temp_mask = temp_image$DQ$imDat
      JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
      
      temp_zap = profoundSkyScan(image = temp_image$SCI$imDat,
                                 mask = (temp_image$SCI$imDat==0) | JWST_cal_mask,
                                 # mask = foo,
                                 clip = c(0.0,0.9),
                                 scan_block = c(nrow(temp_image$SCI$imDat), ncol(temp_image$SCI$imDat)),
                                 trend_block = trend_block,
                                 keep_trend = keep_trend)
      
      Rfits_write(
        data = temp_image,
        filename = paste0(fullbase,'_MIRI_trim.fits')
      )
      Rfits_write_pix(temp_zap$image_fix, paste0(fullbase,'_MIRI_trim.fits'), ext=2)
      
      check_Nhdu = Rfits_nhdu(paste0(fullbase,'_MIRI_trim.fits'))
      extloc = Rfits_extname_to_ext(paste0(fullbase,'_MIRI_trim.fits'), 'SKY_Pro1oF')
      
      if(is.na(extloc)){
        Rfits_write_image(temp_zap$row_map + temp_zap$col_map, filename=paste0(fullbase,'_MIRI_trim.fits'), create_file=FALSE)
        Rfits_write_key(filename=paste0(fullbase,'_MIRI_trim.fits'), ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SKY_Pro1oF', keycomment='extension name')
      }else{
        Rfits_write_pix(temp_zap$row_map + temp_zap$col_map, paste0(fullbase,'_MIRI_trim.fits'), ext=extloc)
      }
    }else{
      if(!(any(c(ow_vlarge, ow_large)))){
        if(any(temp_image[[1]]$keyvalues$VISIT_ID == ID_vlarge$VISIT_ID & temp_image[[1]]$keyvalues$MODULE == ID_vlarge$MODULE)){
          trend_block = trend_block_vlarge
          keep_trend = TRUE
          message(temp_image[[1]]$keyvalues$VISIT_ID,' ',temp_image[[1]]$keyvalues$MODULE, ' KT: TRUE', ' TB: ',trend_block)
        }else if(any(temp_image[[1]]$keyvalues$VISIT_ID == ID_large$VISIT_ID & temp_image[[1]]$keyvalues$MODULE == ID_large$MODULE)){
          trend_block = trend_block_large
          keep_trend = TRUE
          message(temp_image[[1]]$keyvalues$VISIT_ID,' ',temp_image[[1]]$keyvalues$MODULE, ' KT: TRUE', ' TB: ',trend_block)
        }else if(any(temp_image[[1]]$keyvalues$VISIT_ID == ID_large$VISIT_ID & temp_image[[1]]$keyvalues$INSTRUME == ID_large$MODULE)){
          trend_block = trend_block_large
          keep_trend = TRUE
          message(temp_image[[1]]$keyvalues$VISIT_ID,' ',temp_image[[1]]$keyvalues$INSTRUME, ' KT: TRUE', ' TB: ',trend_block)
        }else{
          keep_trend = FALSE
          message(temp_image[[1]]$keyvalues$VISIT_ID,' ',temp_image[[1]]$keyvalues$MODULE, ' KT: FALSE')
        }
      }else{
        ## set over write
        if(ow_large & ow_vlarge){
          message("Can't use both LARGE and VLARGE keep_trend \n Defaulting to keep_trend FALSE")
          keep_trend = FALSE
        }else{
          if(ow_vlarge){
            keep_trend = TRUE
            trend_block = trend_block_vlarge
            message("KT: TRUE", " TB ", trend_block_vlarge)
          }
          if(ow_large){
            keep_trend = TRUE
            trend_block = trend_block_large
            message("KT: TRUE", " TB ", trend_block_large)
          }
        }
      }
      
      temp_mask = temp_image$DQ$imDat
      JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
      
      if(nrow(temp_image$SCI$imDat) != 2048 | ncol(temp_image$SCI$imDat) != 2048){
        scan_block = c(
          nrow(temp_image$SCI$imDat),
          ncol(temp_image$SCI$imDat)
        )
      }else{
        scan_block = c(512, 2048)
      }
      
      temp_zap = profoundSkyScan(image = temp_image$SCI$imDat,
                                 mask = (temp_image$SCI$imDat==0) | JWST_cal_mask,
                                 clip = c(0.0,0.9),
                                 scan_block = scan_block,
                                 trend_block = trend_block,
                                 keep_trend = keep_trend)
      
      file.copy(filelist[i], paste0(fullbase,'_Pro1oF.fits'), overwrite=TRUE)
      Rfits_write_pix(temp_zap$image_fix, paste0(fullbase,'_Pro1oF.fits'), ext=2)

      check_Nhdu = Rfits_nhdu(paste0(fullbase,'_Pro1oF.fits'))
      extloc = Rfits_extname_to_ext(paste0(fullbase,'_Pro1oF.fits'), 'SKY_Pro1oF')
      
      if(is.na(extloc)){
        Rfits_write_image(temp_zap$row_map + temp_zap$col_map, filename=paste0(fullbase,'_Pro1oF.fits'), create_file=FALSE)
        Rfits_write_key(filename=paste0(fullbase,'_Pro1oF.fits'), ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SKY_Pro1oF', keycomment='extension name')
      }else{
        Rfits_write_pix(temp_zap$row_map + temp_zap$col_map, paste0(fullbase,'_Pro1oF.fits'), ext=extloc)
      }
    }
  
    return(NULL)
  }
}
do_cal_process = function(input_args){
  
  cat("\n")
  message("## Calculating sky statistics ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  Pro1oF_dir = input_args$Pro1oF_dir
  sky_frames_dir = input_args$sky_frames_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  
  filelist = load_files(input_args, which_module = "cal_process")
  
  message("Showing first <10 files:")
  cat(head(filelist, 10), sep='\n')
  cat("...")
  cat('Processing',length(filelist),'files...\n')

  #get main info:
  obs_info = Rfits_key_scan(filelist = filelist,
                            keylist = c('VISIT_ID',
                                        'OBS_ID',
                                        'EXPOSURE',
                                        'DETECTOR',
                                        'FILTER',
                                        'DATE-OBS',
                                        'EFFEXPTM',
                                        'CAL_VER',
                                        'CRDS_CTX'),
                            extlist = 1,
                            fileinfo = 'all',
                            cores = cores
  )
  
  registerDoParallel(cores=cores)
  
  lo_loop = 1
  hi_loop = dim(obs_info)[1]
  
  dummy = foreach(i = lo_loop:hi_loop, .inorder=FALSE)%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    #setup files
    file_image = as.character(obs_info[i,"full"])
    basename = strsplit(basename(file_image),'.fits',fixed=T)[[1]]
    fullbase = paste(sky_frames_dir,basename,sep='/')
    file_sky = paste0(fullbase,'_sky_',obs_info[i,"FILTER"],'.fits')
    
    #get cal data
    JWST_cal_image = Rfits_read_image(file_image, ext=2)
    temp_mask = Rfits_read_image(file_image, ext=4, header = FALSE)
    JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
    
    #basic info
    suppressMessages({ #do not really want to see lots of SIP warnings!
      RA_im = centre(JWST_cal_image)[1]
      Dec_im = centre(JWST_cal_image)[2]
    })
    NAXIS1 = JWST_cal_image$keyvalues$NAXIS1
    NAXIS2 = JWST_cal_image$keyvalues$NAXIS2
    Npix = NAXIS1*NAXIS2

    if(do_MIRI){
      box = 64
      redoskysize = 21
      sky_poly_deg = 2
    }else{
      box = 512
      redoskysize = 101
      sky_poly_deg = 2
    }
    
    ## TryCatch for ProFound in case there are not enough sky pixels to be sampled.
    ## Usually caused by zealous object dilation => lower redoskysize
    pro_func = function(x){
      tryCatch(
        {
          pro = profoundProFound(x$image, mask = x$mask, 
                                 skycut=2, 
                                 pixcut=5, 
                                 box=x$box, 
                                 redoskysize = x$redoskysize, 
                                 roughpedestal = TRUE, 
                                 tolerance=Inf)
          return(list(pro = pro, redoskysize=101))
        },
        error = function(cond){
          message("Adjusting redoskysize=", redoskysize, " from 101")
          redoskysize = 21
          pro = profoundProFound(x$image, mask = x$mask, 
                                 skycut=2, 
                                 pixcut=5, 
                                 box=x$box, 
                                 redoskysize = redoskysize, 
                                 roughpedestal = TRUE, 
                                 tolerance=Inf)
          return(list(pro=pro, redoskysize=redoskysize))
        },
        warning = function(cond){
          NULL
        },
        finally = {
          NULL
        }
      )
    }
    
    pro_out = pro_func(list(image = JWST_cal_image$imDat, mask = JWST_cal_mask, box = box, redoskysize = redoskysize))
    pro = pro_out$pro
    redoskysize = pro_out$redoskysize
    
    sky_redo_func = function(x){
        tryCatch(
          {
            sky_redo = profoundSkyPoly(x$image, objects=x$pro$objects_redo, degree=x$sky_poly_deg, quancut=0.99, mask=x$mask)
            return(sky_redo)
          },
          error = function(cond){
            message("Adjusting object mask for polynomial skies")
            sky_redo = profoundSkyPoly(x$image, objects=x$pro$objects, degree=x$sky_poly_deg, quancut=0.99, mask=x$mask)
            return(sky_redo)
          },
          warning = function(cond){
            NULL
          },
          finally = {
            NULL
          }
        )
    }
    
    if(do_MIRI){
      sky_redo = list(sky = pro$sky)
      pro_redo = pro
    }else{
      sky_redo = sky_redo_func(list(image = JWST_cal_image$imDat, pro = pro, mask = JWST_cal_mask, sky_poly_deg = sky_poly_deg))
      pro_redo = profoundProFound(JWST_cal_image$imDat, mask=JWST_cal_mask, skycut=2, pixcut=5, box=box, grid = box, redoskysize=redoskysize, sky=sky_redo$sky, redosky=FALSE, tolerance=Inf)
    }
    

    ## what to do if no objects in frame
    if(is.null(pro_redo$objects_redo)){
      pro_redo$objects_redo = matrix(0L, dim(JWST_cal_image$imDat)[1], dim(JWST_cal_image$imDat)[2])
    }
    sky_med = median(sky_redo$sky[JWST_cal_mask == 0 & pro_redo$objects_redo == 0], na.rm=TRUE)
    skyRMS_med = median(pro_redo$skyRMS[JWST_cal_mask == 0 & pro_redo$objects_redo == 0], na.rm=TRUE)
    
    if(is.null(pro_redo$skyChiSq)){
      pro_redo$skyChiSq = 1e6
    }else if(is.infinite(pro_redo$skyChiSq) | is.na(pro_redo$skyChiSq)){
      pro_redo$skyChiSq = 1e6 ## set to large value, Rfits no like Inf...
    }
    
    maskpix = sum(JWST_cal_mask!=0, na.rm=TRUE)/Npix
    objpix = sum(pro_redo$objects_redo!=0, na.rm=TRUE)/Npix
    goodpix = sum(JWST_cal_mask==0 & pro_redo$objects_redo==0, na.rm=TRUE)/Npix
    
    if(file.exists(file_sky)){
      message('Removing old sky file: ',file_sky)
      file.remove(file_sky)
    }
    
    Rfits_write_header(keyvalues = list(
      FILEIM = basename(file_image),
      PATHIM = dirname(file_image),
      FILESKY = basename(file_sky),
      PATHSKY = dirname(file_sky),
      VERSION = pipe_version,
      VISIT_ID = obs_info$VISIT_ID[i],
      OBS_ID = obs_info$OBS_ID[i],
      EXPOSURE = obs_info$EXPOSURE[i],
      DETECTOR = obs_info$DETECTOR[i],
      FILTER = obs_info$FILTER[i],
      DATE = obs_info$`DATE-OBS`[i],
      EXPTIME =  obs_info$EFFEXPTM[i],
      CAL_VER =  obs_info$CAL_VER[i],
      CRDS_CTX =  obs_info$CRDS_CTX[i],
      RA = RA_im,
      DEC = Dec_im,
      NAXIS1 = NAXIS1,
      NAXIS2 = NAXIS2,
      NPIX = Npix,
      SKY = ifelse(is.na(sky_med), -999, sky_med),
      SKYRMS = ifelse(is.na(skyRMS_med), -999, skyRMS_med),
      SKYCHI = pro_redo$skyChiSq,
      MASK = maskpix,
      OBJ = objpix,
      GOODPIX = goodpix
    ), filename=file_sky, create_file=TRUE, overwrite_file=TRUE)
    Rfits_write_key(filename=file_sky, ext=1, keyname='EXTNAME', keyvalue='INFO', keycomment='extension name')
    
    Rfits_write_image(sky_redo$sky, filename=file_sky, create_file=FALSE)
    Rfits_write_key(filename=file_sky, ext=2, keyname='EXTNAME', keyvalue='SKY', keycomment='extension name')
    
    Rfits_write_image(pro_redo$skyRMS, filename=file_sky, create_file=FALSE)
    Rfits_write_key(filename=file_sky, ext=3, keyname='EXTNAME', keyvalue='SKYRMS', keycomment='extension name')
    
    Rfits_write_image(JWST_cal_mask, filename=file_sky, create_file=FALSE)
    Rfits_write_key(filename=file_sky, ext=4, keyname='EXTNAME', keyvalue='PIXELMASK', keycomment='extension name')
    
    Rfits_write_image(pro_redo$objects_redo, filename=file_sky, create_file=FALSE)
    Rfits_write_key(filename=file_sky, ext=5, keyname='EXTNAME', keyvalue='OBJECTMASK', keycomment='extension name')
    
    return(NULL)
  }
}
do_regen_sky_info = function(input_args){
  
  cat("\n")
  message("## Generating sky statistics info ##")
  cat("\n")
  
  sky_pro_dir = input_args$sky_pro_dir
  cores = input_args$cores_pro
  
  registerDoParallel(cores=cores)
  
  sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
  filelist = list.files(sky_frames_dir, full.names = TRUE, pattern = glob2rx(paste0("*", input_args$VID, "*.fits")))
  filelist = grep('.fits$', filelist, value=TRUE)
  
  foo = list.files(sky_frames_dir, full.names = FALSE)
  vids = substr(grep('.fits$', foo, value=TRUE), 4, 13)
  
  cat('Processing',length(filelist),'files\n')
  
  sky_info = foreach(i = 1:length(filelist), .combine="rbind")%dopar%{
    temp_info = Rfits_read_header(filelist[i])$keyvalues
    temp_info[1:4] = NULL
    temp_info$EXTNAME = NULL
    return(temp_info)
  }
  if(length(filelist)==1){
    class(sky_info) = "list"
    sky_info = data.frame((sky_info))
  }
  colnames(sky_info) = tolower(colnames(sky_info))
  
  write.csv(sky_info, paste0(sky_pro_dir, '/sky_info.csv'), row.names=FALSE)
  
}
do_super_sky = function(input_args){
  
  cat("\n")
  message("## Making super skies ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  sky_pro_dir = input_args$sky_pro_dir
  VID = input_args$VID
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  
  sky_ChiSq_cut = 1.1
  good_pix_cut = 0.15
  
  registerDoParallel(cores=cores)
  
  sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
  sky_super_dir = paste0(sky_pro_dir, "/sky_super/")
  # filelist = list.files(sky_frames_dir, full.names = TRUE, pattern = ".fits$")
  filelist = load_files(input_args, which_module = "super_sky")
  sky_info = fread(paste0(sky_pro_dir, '/sky_info.csv'))
  sky_info = sky_info[gsub("//", "/", paste0(sky_info$pathsky, "/", sky_info$filesky))%in%gsub("//", "/",filelist), ]

  niriss_det = c("NIS")
  miri_det = c("MIRIMAGE")
  short_det = c("NRCA1", "NRCA2", "NRCA3", "NRCA4", "NRCB1", "NRCB2", "NRCB3", "NRCB4")
  long_det = c("NRCALONG", "NRCBLONG")
  
  short_filt = sort(unique(sky_info[detector %in% short_det,filter]))
  long_filt = sort(unique(sky_info[detector %in% long_det,filter]))
  nis_filt = sort(unique(sky_info[detector %in% niriss_det,filter]))
  miri_filt = sort(unique(sky_info[detector %in% miri_det,filter]))
  
  #short_filt = c("F090W", "F115W", "F150W", "F182M", "F200W", "F210M")
  #long_filt = c("F277W", "F300M", "F335M", "F356W", "F360M", "F410M", "F444W")
  #miri_filt = c("F1000W", "F1500W", "F1800W", "F770W") #we don't care about MIRI for now
  
  combine_grid_short = expand.grid(short_det, short_filt, stringsAsFactors=FALSE)
  combine_grid_long = expand.grid(long_det, long_filt, stringsAsFactors=FALSE)
  combine_grid_nis = expand.grid(niriss_det, nis_filt, stringsAsFactors=FALSE)
  combine_grid_miri = expand.grid(miri_det, miri_filt, stringsAsFactors=FALSE)
  #combine_grid_miri = expand.grid(miri_det, miri_filt, stringsAsFactors=FALSE) #we don't care about MIRI for now
  combine_grid = rbind(combine_grid_short, combine_grid_long, combine_grid_nis, combine_grid_miri)
  
  sky_filelist = list.files(sky_info[1,pathsky], pattern = ".fits$")

  # if(do_NIRISS){
  #   filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist) & grepl("nis", filelist)]
  #   sky_info = sky_info[grepl(VID, sky_info$visit_id) & grepl("NIS", sky_info$detector), ]
  #   niriss_det = c("NIS")
  #   nis_filt = sort(unique(sky_info[detector %in% niriss_det,filter]))
  #   combine_grid = expand.grid(niriss_det, nis_filt, stringsAsFactors=FALSE)
  #   sky_filelist = list.files(sky_info[1,pathsky])
  #   sky_filelist = sky_filelist[grepl('.fits$', sky_filelist) & grepl(VID, sky_filelist) & grepl("nis", sky_filelist)]
  # }
  
  dummy = foreach(i = 1:dim(combine_grid)[1], .inorder=FALSE)%dopar%{
    if(combine_grid[i,1] %in% c('NRCA1','NRCA2','NRCA3','NRCA4','NRCB1','NRCB2','NRCB3','NRCB4', 'MIRIMAGE')){
      temp_info = sky_info[detector == combine_grid[i,1] & filter == combine_grid[i,2] & skychi < sky_ChiSq_cut & goodpix >= good_pix_cut,]
    }else{
      temp_info = sky_info[detector == combine_grid[i,1] & filter == combine_grid[i,2] & skychi < sky_ChiSq_cut,]
    }
    
    Nsky = dim(temp_info)[1]
    
    message(combine_grid[i,1], ' ', combine_grid[i,2],': ', Nsky)
    
    if(Nsky > 0){
      sky_array = array(dim = c(temp_info[1,naxis1], temp_info[1,naxis2], Nsky)) #create empty array
      for(j in 1:dim(temp_info)[1]){
        #safe name (since sometimes the sky file name gets truncated due to 80 character FITS limits)
        file_sky = temp_info[j,filesky]
        file_sky = grep(temp_info[j,filesky], sky_filelist, value=TRUE)
        if(length(file_sky) == 0){
          stop('Missing',)
        }
        
        #fill array up with sky data for stacking
        sky_array[,,j] = Rfits_read_image(paste(temp_info[j,pathsky], file_sky, sep='/'), ext=2, header=FALSE) - temp_info[j,sky]
      }
      sky_mean = rowMeans(sky_array, na.rm = FALSE, dims = 2)
    }else{
      message("No usable data for ", combine_grid[i,1]," ",combine_grid[i,2])
      temp_file = sky_info[detector == combine_grid[i,1] & filter == combine_grid[i,2],fileim][1]
      NAXIS1 = sky_info[detector == combine_grid[i,1] & filter == combine_grid[i,2],naxis1][1]
      NAXIS2 = sky_info[detector == combine_grid[i,1] & filter == combine_grid[i,2],naxis2][1]
      if(is.na(NAXIS1) | is.na(NAXIS2)){
        sky_mean = NULL
      }else{
        sky_mean = matrix(0, NAXIS1, NAXIS2)
      }
    }
    
    if(!is.null(sky_mean)){
      CairoJPEG(filename = paste0(sky_super_dir,'/super_',combine_grid[i,1],'_',combine_grid[i,2],'.jpeg'), width=1000,height=1000)
      magimage(sky_mean, qdiff=T)
      dev.off()
      
      file_super_sky = paste0(sky_super_dir,'/super_',combine_grid[i,1],'_',combine_grid[i,2],'.fits')
      Rfits_write_image(sky_mean, filename=file_super_sky)
      
      header_super_sky = list(Nsky = Nsky)
      
      if(Nsky > 0){
        header_super_sky = c(header_super_sky, as.list(temp_info$filesky))
        names(header_super_sky) = c('NSKY', paste0('FILE',1:Nsky))
      }
      
      Rfits_write_header(file_super_sky, keyvalues=header_super_sky)
    }
  }
}
do_apply_super_sky = function(input_args){
  
  cat("\n")
  message("## Applying super skies ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  Pro1oF_dir = input_args$Pro1oF_dir
  cal_sky_dir = input_args$cal_sky_dir
  sky_pro_dir = input_args$sky_pro_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  
  registerDoParallel(cores=cores)
  
  sky_info = fread(paste0(sky_pro_dir, "/sky_info.csv"))
  sky_super_dir = paste0(sky_pro_dir, "/sky_super/")
  
  setcal_sky_remfile = function(file_image){
    path = paste0(strsplit(dirname(file_image), '/cal',fixed=T),'/cal_sky/')
    base = strsplit(basename(file_image),'.fits',fixed=T)[[1]]
    return(paste0(path,base,'_sky_rem.fits'))
  }
  
  sky_filelist = list.files(sky_info[1,pathsky])

    # sky_filelist = load_files(input_args, which_module = "apply_super", sky_info = sky_info)$sky_filelist
  sky_info = load_files(input_args, which_module = "apply_super", sky_info = sky_info)$sky_info
  sky_filelist = sky_info$filesky

  cat('Processing ',length(sky_filelist),'files...\n')
  
  lo_loop = 1
  hi_loop = length(sky_filelist)
  
  dummy = foreach(i = lo_loop:hi_loop, .inorder=FALSE)%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    file_image = paste(sky_info[i,pathim], sky_info[i,fileim], sep='/')
    temp_cal = Rfits_read(file_image, pointer=FALSE)
    temp_cal$DQ$keyvalues$BZERO = 0
    #safe name (since sometimes the sky file name gets truncated due to 80 character FITS limits)
    file_sky = sky_info[i,filesky]
    file_sky = grep(sky_info[i,filesky], sky_filelist, value=TRUE)
    if(length(file_sky) == 0){
      stop('Missing file_sky ',i)
    }
    file_sky = paste(sky_info[i,pathsky], file_sky, sep='/')
    
    temp_sky = Rfits_read(file_sky, pointer=FALSE)
    
    
    file_super_sky = paste0(sky_super_dir,'/super_',temp_cal[[1]]$keyvalues$DETECTOR,'_',temp_cal[[1]]$keyvalues$FILTER,'.fits')
    super_sky = Rfits_read_image(file_super_sky, header=FALSE)
    
    temp_lm = lm(temp_cal$SCI$imDat[temp_sky$PIXELMASK$imDat==0 & temp_sky$OBJECTMASK$imDat==0] ~ super_sky[temp_sky$PIXELMASK$imDat==0 & temp_sky$OBJECTMASK$imDat==0])$coefficients
    
    if(is.na(temp_lm[2])){
      temp_lm[2] = 0
    }
    
    if(temp_lm[2] < 0){ #force to be positive
      temp_lm[2] = 0
    }
    
    if(temp_lm[2] > 2){ #do not allow extreme values (expect to be near 1)
      temp_lm[2] = 2
    }
    
    sky_map = super_sky*temp_lm[2] + temp_lm[1]
    
    #maybe an extra sky stage here to clean things up per frame?
    
    rezero = which(temp_cal$SCI$imDat==0)
    temp_cal_sky = temp_cal$SCI$imDat - sky_map
    temp_cal_sky[!(temp_sky$PIXELMASK$imDat==0 & temp_sky$OBJECTMASK$imDat==0)] = NA
    
    #pedestal check - basically use the extremes of the frame to find a safe pedestal:
    
    pix_buffer = 10+1 ## Ignore first 10 pixels since these are often ratty
    imdim = dim(temp_cal_sky)

    edge_pix = (pix_buffer):(pix_buffer+4)
    
    LHS = median(temp_cal_sky[edge_pix,], na.rm=TRUE)
    RHS = median(temp_cal_sky[imdim[1]-edge_pix,], na.rm=TRUE)
    BHS = median(temp_cal_sky[,edge_pix], na.rm=TRUE)
    THS = median(temp_cal_sky[,imdim[2]-edge_pix], na.rm=TRUE)
    
    #corners (with masking, this gives about as many pixels as the edge above):
    corner_pix = (pix_buffer):(pix_buffer*10)
    BL = median(temp_cal_sky[corner_pix,corner_pix], na.rm=TRUE)
    TL = median(temp_cal_sky[corner_pix,imdim[2]-corner_pix], na.rm=TRUE)
    TR = median(temp_cal_sky[imdim[1]-corner_pix,imdim[2]-corner_pix], na.rm=TRUE)
    BR = median(temp_cal_sky[imdim[1]-corner_pix,corner_pix], na.rm=TRUE)
    
    #centre
    CC = median( magcutout(temp_cal_sky, box = c(100, 100))$image, na.rm=TRUE)
    
    # #edges (here we ignore the first 10 pixels since these are often ratty):
    # LHS = median(temp_cal_sky[11:15,], na.rm=TRUE)
    # RHS = median(temp_cal_sky[2049-11:15,], na.rm=TRUE)
    # BHS = median(temp_cal_sky[,11:15], na.rm=TRUE)
    # THS = median(temp_cal_sky[,2049-11:15], na.rm=TRUE)
    # #corners (with masking, this gives about as many pixels as the edge above):
    # BL = median(temp_cal_sky[11:110,11:110], na.rm=TRUE)
    # TL = median(temp_cal_sky[11:110,2049-11:110], na.rm=TRUE)
    # TR = median(temp_cal_sky[2049-11:110,2049-11:110], na.rm=TRUE)
    # BR = median(temp_cal_sky[2049-11:110,11:110], na.rm=TRUE)
    # #centre
    # CC = median(temp_cal_sky[974:1073,974:1073], na.rm=TRUE)
    
    pedestals = c(BL, LHS, TL, THS, TR, RHS, BR, BHS, CC)
    
    final_pedestal = quantile(pedestals, 0.5, na.rm=TRUE) #logic being a real issue due to a large object
    
    if(!is.finite(final_pedestal)){
      final_pedestal = 0
    }
    
    sky_map = sky_map + final_pedestal
    
    temp_cal$SCI$imDat = temp_cal$SCI$imDat - sky_map
    
    if(length(rezero) > 0){
      temp_cal$SCI$imDat[rezero] = 0
    }
    #temp_cal$SCI$keyvalues$EXTNAME = 'SCI_SKY'
    #temp_cal$SCI$keyvalues$SKY_M = temp_lm[2]
    #temp_cal$SCI$keyvalues$SKY_B = temp_lm[1]
    #temp_cal$SCI$keycomments$SKY_M = 'SKY M coef in SKY = M.SuperSky + B'
    #temp_cal$SCI$keycomments$SKY_B = 'SKY B coef in SKY = M.SuperSky + B'
    #temp_cal$SCI$keynames = names(temp_cal$SCI$keyvalues)
    
    #temp_cal$SKY = sky_map
    
    file_cal_sky = setcal_sky_remfile(file_image)
    
    if(file.exists(file_cal_sky)){
      message('Removing old cal_sky file: ',file_cal_sky)
      file.remove(file_cal_sky)
    }
    
    file.copy(file_image, file_cal_sky)
    Rfits_write_pix(temp_cal$SCI$imDat, file_cal_sky, ext=2)
    Rfits_write_key(file_cal_sky, keyname='SKY_M', keyvalue=temp_lm[2], keycomment='SKY M coef in SKY = M.SuperSky + B + P', ext=2)
    Rfits_write_key(file_cal_sky, keyname='SKY_B', keyvalue=temp_lm[1], keycomment='SKY B coef in SKY = M.SuperSky + B + P', ext=2)
    Rfits_write_key(file_cal_sky, keyname='SKY_P', keyvalue=final_pedestal, keycomment='SKY P coef in SKY = M.SuperSky + B + P', ext=2)
    #Rfits_write(temp_cal, file_cal_sky)
    
    # if(check_Ndhu == 10){
    #   Rfits_write_image(sky_map, file_cal_sky, create_file=FALSE, create_ext=TRUE)
    #   Rfits_write_key(filename=file_cal_sky, ext=check_Ndhu+1, keyname='EXTNAME', keyvalue='SKY_Super', keycomment='extension name')
    # }else if(check_Ndhu == 11){
    #   if(check_info$headers[[11]]$keyvalues$EXTNAME == 'SKY_Wisp'){
    #     #If we have made the SKY_Wisp extension.
    #     Rfits_write_image(sky_map, file_cal_sky, create_file=FALSE, create_ext=TRUE)
    #     Rfits_write_key(filename=file_cal_sky, ext=check_Ndhu+1, keyname='EXTNAME', keyvalue='SKY_Super', keycomment='extension name')
    #   }
    #   if(check_info$headers[[11]]$keyvalues$EXTNAME == 'SKY_Super'){
    #     #If we have previously made the SKY_Super extension
    #     Rfits_write_pix(sky_map, file_cal_sky, ext=check_Ndhu)
    #   }
    # }else if(check_Ndhu == 12){
    #   if(check_info$headers[[12]]$keyvalues$EXTNAME == 'SKY_Super'){
    #     #If we have previously made the SKY_Super extension. In this case we must have SKY_Pro1oF and SKY_Wisp extensions.
    #     Rfits_write_pix(sky_map, file_cal_sky, ext=check_Ndhu)
    #     Rfits_write_key(filename=file_cal_sky, ext=check_Nhdu, keyname='EXTNAME', keyvalue='SKY_Super', keycomment='extension name')
    #   }
    # }
    
    check_Ndhu = Rfits_nhdu(file_cal_sky)
    extloc = Rfits_extname_to_ext(file_cal_sky, 'SKY_Super')
    
    if(is.na(extloc)){
      Rfits_write_image(sky_map, file_cal_sky, create_file=FALSE, create_ext=TRUE)
      Rfits_write_key(filename=file_cal_sky, ext=check_Ndhu+1, keyname='EXTNAME', keyvalue='SKY_Super', keycomment='extension name')
    }else{
      Rfits_write_pix(sky_map, file_cal_sky, ext=extloc)
    }
    
    return(NULL)
  }
  
}
do_modify_pedestal = function(input_args){
  
  cat("\n")
  message("## Modifying pedestal ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  cal_sky_dir = input_args$cal_sky_dir
  cal_sky_renorm_dir = input_args$cal_sky_renorm_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  
  filelist = load_files(input_args, which_module = "modify_pedestal")
  
  cat('Processing',length(filelist),'file...\n')
  
  cal_sky_info_ext1 = Rfits_key_scan(filelist = filelist,
                                     keylist = c('VISIT_ID', 'OBS_ID', 'EXPOSURE', 'DETECTOR', 'MODULE', 'CHANNEL', 'FILTER'),
                                     extlist = 1,
                                     fileinfo = c('Full','File'),
                                     cores = cores
  )
  
  cal_sky_info_ext2 = Rfits_key_scan(filelist = filelist,
                                     keylist = c('SKY_M', 'SKY_B', 'SKY_P'),
                                     extlist = 2,
                                     fileinfo = NULL,
                                     cores = cores
  )
  
  registerDoParallel(cores=cores)
  
  cal_sky_info = as.data.table(cbind(cal_sky_info_ext1, cal_sky_info_ext2))
  
  ped_info = cal_sky_info[,list(ped_mean = mean(SKY_B + SKY_P), ped_med = median(SKY_B + SKY_P), ped_max = max(SKY_B + SKY_P), ped_min = min(SKY_B + SKY_P)), by=list(OBS_ID, EXPOSURE)]
  ped_info[,ped_diff := (ped_max - ped_min) / ped_mean]
  
  lo_loop = 1
  hi_loop = dim(cal_sky_info)[1]
  
  if(do_NIRISS){
    dummy = foreach(i = lo_loop:hi_loop)%dopar%{
      if(i %% 100 == 0){
        message('File ',i,' of ', hi_loop)
      }
      file_cal_sky_renorm = paste0(cal_sky_renorm_dir,'/',cal_sky_info[i,file])
      
      if(file.exists(file_cal_sky_renorm)){
        message('Removing old cal_sky file: ',file_cal_sky_renorm)
        file.remove(file_cal_sky_renorm)
      }
      file.copy(cal_sky_info[i,full], cal_sky_renorm_dir) #no pedestal adjustment for NIRISS
      return(NULL)
    }
  } else if(do_MIRI){
    dummy = foreach(i = lo_loop:hi_loop)%dopar%{
      if(i %% 100 == 0){
        message('File ',i,' of ', hi_loop)
      }
      file_cal_sky_renorm = paste0(cal_sky_renorm_dir,'/',cal_sky_info[i,file])

      if(file.exists(file_cal_sky_renorm)){
        message('Removing old cal_sky file: ',file_cal_sky_renorm)
        file.remove(file_cal_sky_renorm)
      }
      file.copy(cal_sky_info[i,full], cal_sky_renorm_dir) #no pedestal adjustment for NIRISS
      return(NULL)
    }
  }else{
    dummy = foreach(i = lo_loop:hi_loop)%dopar%{
      if(i %% 100 == 0){
        message('File ',i,' of ', hi_loop)
      }
      file_cal_sky_renorm = paste0(cal_sky_renorm_dir,'/',cal_sky_info[i,file])
      
      if(file.exists(file_cal_sky_renorm)){
        message('Removing old cal_sky file: ',file_cal_sky_renorm)
        file.remove(file_cal_sky_renorm)
      }
      file.copy(cal_sky_info[i,full], cal_sky_renorm_dir)
      
      if(cal_sky_info[i,CHANNEL] == 'SHORT'){
        temp_cal_sky = Rfits_read(cal_sky_info[i,full])
        ext_names = names(temp_cal_sky)
        
        message('Testing sky with ProFound: ',cal_sky_info[i,file])
        temp_mask = Rfits_read_image(cal_sky_info[i,full], ext = "DQ", header = F)
        JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
        
        pro_current_sky = profoundProFound(temp_cal_sky$SCI[,]$imDat, mask=JWST_cal_mask, sky=0, box=683, redosky = FALSE)

        if(pro_current_sky$skyChiSq < 0.9 | pro_current_sky$skyChiSq > 1.1){
          message('Bad Sky: ',cal_sky_info[i,file],'. Current: ',round(pro_current_sky$skyChiSq,2),' Trying to make better sky with ProFound')
          if(sum(pro_current_sky$objects) / prod(dim(temp_cal_sky$SCI)) < 0.2){
            pro_new_sky = profoundProFound(temp_cal_sky$SCI[,]$imDat, mask=JWST_cal_mask, box=683, roughpedestal = TRUE)
            if(pro_new_sky$skyChiSq > 0.9 & pro_new_sky$skyChiSq < 1.1){
              message('Using better ProFound Sky: ',cal_sky_info[i,file],' Old: ',round(pro_current_sky$skyChiSq,2),' New: ',round(pro_new_sky$skyChiSq,2))
              replace_SCI = temp_cal_sky$SCI[,]$imDat - pro_new_sky$sky
              replace_SKY_Super = temp_cal_sky$SKY_Super[,]$imDat + pro_new_sky$sky
              replace_SKY_P = temp_cal_sky$SCI$keyvalues$SKY_P + mean(pro_new_sky$sky, na.rm=TRUE)
              
              Rfits_write_pix(replace_SCI, file_cal_sky_renorm, ext=which(ext_names == 'SCI')[1])
              Rfits_write_pix(replace_SKY_Super, file_cal_sky_renorm, ext=which(ext_names == 'SKY_Super'))
              Rfits_write_key(file_cal_sky_renorm, keyname='SKY_P', keyvalue=replace_SKY_P, keycomment='SKY P coef in SKY = M.SuperSky + B + P', ext=2)
              Rfits_write_key(file_cal_sky_renorm, keyname='SKY_CHI', keyvalue=pro_new_sky$skyChiSq, keycomment='SKY Chi-Sq', ext=2)
            }else{
              Rfits_write_key(file_cal_sky_renorm, keyname='SKY_CHI', keyvalue=pro_current_sky$skyChiSq, keycomment='SKY Chi-Sq', ext=2)
            }
          }else{
            Rfits_write_key(file_cal_sky_renorm, keyname='SKY_CHI', keyvalue=pro_current_sky$skyChiSq, keycomment='SKY Chi-Sq', ext=2)
          }
        }
        #    }
      }
      return(NULL)
    }
  }
}
do_cal_sky_info = function(input_args){

  cat("\n")
  message("## Collecting cal_sky_renorm info for ProPane ##")
  cat("\n")
  
  cal_sky_renorm_dir = input_args$cal_sky_renorm_dir
  cal_sky_info_save_dir = input_args$cal_sky_info_save_dir
  cores = input_args$cores_pro
  
  cal_sky_info_WCS = Rfits_key_scan(dirlist=cal_sky_renorm_dir,
                                    keylist = c('PHOTMJSR','PIXAR_SR','CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'),
                                    extlist = 2,
                                    fileinfo = 'all',
                                    cores = cores
  )
  cal_sky_info = Rfits_key_scan(dirlist=cal_sky_renorm_dir,
                                keylist = c('VISIT_ID','DATE-OBS', 'DETECTOR', 'FILTER', 'EFFEXPTM', 'CRDS_CTX'),
                                extlist = 1,
                                fileinfo = 'none',
                                cores = cores
  )
  cal_sky_info = cbind(cal_sky_info_WCS, cal_sky_info)
  cal_sky_info$MAGZERO = -2.5*log10(cal_sky_info$PIXAR_SR*1e6) + 8.9 #note the -ve
  #cal_sky_info$magPHOTMJSR = -2.5*log10(cal_sky_info$PHOTMJSR) #note the -ve
  #cal_sky_info$MAGZERO_test = -2.5*log10(cal_sky_info$PIXAR_SR*1e6)+8.9
  cal_sky_info = as.data.table(cal_sky_info)
  
  pmap = cal_sky_info[,list(PHOTMJSR=mean(PHOTMJSR)), keyby=list(CRDS_CTX,DETECTOR,FILTER)]
  
  cal_sky_info$MAGZERO_FIX = 0.0
  
  for(i in 1:dim(cal_sky_info)[1]){
    #below is then Eqn [1] below
    temp_MAGZERO_FIX = cal_sky_info[i,MAGZERO] -2.5*log10(as.numeric(pmap[CRDS_CTX=='jwst_1179.pmap' & DETECTOR==cal_sky_info[i,'DETECTOR'] & FILTER==cal_sky_info[i,'FILTER'],PHOTMJSR]) / cal_sky_info[i,PHOTMJSR])
    if(length(temp_MAGZERO_FIX) == 1){
      cal_sky_info[i,MAGZERO_FIX := temp_MAGZERO_FIX] 
    }else{
      cal_sky_info[i,MAGZERO_FIX := cal_sky_info[i,MAGZERO]]
    }
  }
  
  # for(i in 1:dim(cal_sky_info)[1]){
  #   #below is then Eqn [1] below
  #   temp_MAGZERO_FIX = cal_sky_info[i,MAGZERO] -2.5*log10(as.numeric(pmap[CRDS_CTX=='jwst_0984.pmap' & DETECTOR==cal_sky_info[i,'DETECTOR'] & FILTER==cal_sky_info[i,'FILTER'],PHOTMJSR]) / cal_sky_info[i,PHOTMJSR])
  #   if(length(temp_MAGZERO_FIX) == 1){
  #     cal_sky_info[i,MAGZERO_FIX := temp_MAGZERO_FIX] 
  #   }else{
  #     cal_sky_info[i,MAGZERO_FIX := cal_sky_info[i,MAGZERO]]
  #   }
  # }
  
  write.csv(cal_sky_info, paste0(cal_sky_info_save_dir, '/cal_sky_info.csv'), row.names = FALSE)
}
do_gen_stack = function(input_args){
  
  cat("\n")
  message("## Making ProPane stacks ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  VID = input_args$VID
  FILT = input_args$FILT
  ref_dir = input_args$ref_dir
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  magzero = input_args$magzero
  cores = input_args$cores_stack
  tasks = input_args$tasks_stack
  
  ref_cat = input_args$additional_params$ref_cat
  if(is.null(ref_cat)){
    do_tweak = FALSE
  }else{
    do_tweak = TRUE
    ref_cat = data.frame(fread(ref_cat))
  }
  
  ## Optional editable params
  clip_tol = c(40,20) #default (80,40) hot and cold
  clip_sigma = 20
  clip_dilate = 3
  doclip = F
  
  ## Allow user to change these
  NAXIS_short = ifelse(
    is.null(input_args$additional_params$NAXIS_short), 
    6000, #good for ceers, need 6000 for cosmos web
    input_args$additional_params$NAXIS_short) 
  NAXIS_long = ifelse(
    is.null(input_args$additional_params$NAXIS_long), 
    3000, 
    input_args$additional_params$NAXIS_long) 
  
  
  invar_dir = paste0(ref_dir, "/InVar_Stacks/")
  median_dir = paste0(ref_dir, "/Median_Stacks/")
  sky_frames_dir = paste0(ref_dir, "/sky_pro/sky_frames/")
  dump_dir_stub = paste0(ref_dir, "/dump/")
  orig_cal_sky_info = fread(paste0(ref_dir, "/Pro1oF/cal_sky_info.csv"))
  
  unique_visits = unique(orig_cal_sky_info$VISIT_ID)
  
  message("Now running stack!")
  
  
  temp_vid = VID
  VID_list = grep(temp_vid, unique_visits, value = T)
  if(temp_vid == ""){
    not_pid_idx = 1:length(VID_list)
  }else{
    not_pid_idx = sapply(substr(VID_list,1,4), function(x) grepl(x, VID, fixed = T)) ## make sure PID is not embedded in string of VID
  }
  
  for(VID in VID_list[not_pid_idx]){
    
    if(do_NIRISS){
      cal_sky_info = orig_cal_sky_info[grepl(VID, orig_cal_sky_info$VISIT_ID) & grepl("NIS", orig_cal_sky_info$DETECTOR),]
    }else if(do_MIRI){
      cal_sky_info = orig_cal_sky_info[grepl(VID, orig_cal_sky_info$VISIT_ID) & grepl("MIRIMAGE", orig_cal_sky_info$DETECTOR),]
    }else{
      cal_sky_info = orig_cal_sky_info[grepl(VID, orig_cal_sky_info$VISIT_ID) & !grepl("NIS", orig_cal_sky_info$DETECTOR),]
    }
    cal_sky_info$MAGZERO_FIX[is.na(cal_sky_info$MAGZERO_FIX)] = cal_sky_info$MAGZERO[is.na(cal_sky_info$MAGZERO_FIX)]
    
    stack_grid = cal_sky_info[,list(FILTER=unique(FILTER)), keyby=VISIT_ID]
    stack_grid = stack_grid[grepl(FILT, stack_grid$FILTER), ]
    
    if(dim(stack_grid)[1] == 0){
      next ## If the FILTER/VID combo is null then skip 
    }else{
      print(stack_grid)
    }
    
    module_A_WCS_short = cal_sky_info[grepl('NRCA[1-4]', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    module_B_WCS_short = cal_sky_info[grepl('NRCB[1-4]', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    module_A_WCS_long = cal_sky_info[grepl('NRCALONG', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    module_B_WCS_long = cal_sky_info[grepl('NRCBLONG', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    if(dim(module_A_WCS_long)[1] == 0){
      module_A_WCS_long = module_A_WCS_short
      module_A_WCS_long$CD1_1 = module_A_WCS_long$CD1_1*2.0
      module_A_WCS_long$CD1_2 = module_A_WCS_long$CD1_2*2.0
    }
    
    if(dim(module_B_WCS_long)[1] == 0){
      module_B_WCS_long = module_B_WCS_short
      module_B_WCS_long$CD1_1 = module_B_WCS_long$CD1_1*2.0
      module_B_WCS_long$CD1_2 = module_B_WCS_long$CD1_2*2.0
    }
    
    if(dim(module_A_WCS_short)[1] == 0){
      module_A_WCS_short = module_A_WCS_long
      module_A_WCS_short$CD1_1 = module_A_WCS_short$CD1_1/2.0
      module_A_WCS_short$CD1_2 = module_A_WCS_short$CD1_2/2.0
    }
    
    if(dim(module_B_WCS_short)[1] == 0){
      module_B_WCS_short = module_B_WCS_long
      module_B_WCS_short$CD1_1 = module_B_WCS_short$CD1_1/2.0
      module_B_WCS_short$CD1_2 = module_B_WCS_short$CD1_2/2.0
    }
    
    module_NIS = cal_sky_info[grepl('NIS', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    module_MIRIMAGE = cal_sky_info[grepl('MIRIMAGE', DETECTOR),list(CRVAL1=median(CRVAL1), CRVAL2=median(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    module_idx = sapply(list(module_A_WCS_long, module_A_WCS_short, 
                             module_B_WCS_long, module_B_WCS_short,
                             module_NIS, module_MIRIMAGE), 
                        function(x){y = dim(x)[1]
                        y > 0})
    module_list = c('NRCA_short', 'NRCA_long', 'NRCB_short', 'NRCB_long', 'NIS', 'MIRIMAGE')[module_idx]
    
    registerDoParallel(cores=tasks)
    
    lo_loop = 1
    hi_loop = dim(stack_grid)[1]
    
    stack_grid = stack_grid[sample(hi_loop),] #randomise to spread out the CPU load (else all of the SMACS stacking end up in a single task and takes much longer)
    
    list_residual_temp_files = list.files(invar_dir, pattern = "temp_list_.+", full.names = T)
    if(length(list_residual_temp_files)>0){
      unlink(list_residual_temp_files, recursive = T)
    }
    
    dummy = foreach(i = lo_loop:hi_loop)%dopar%{
      message('Stacking ', i,' of ',hi_loop)
      for(j in module_list){
        message('  Processing ',stack_grid[i,VISIT_ID],' ',stack_grid[i,FILTER], ' ', j)
        
        if(j == 'MIRIMAGE'){
          input_info = cal_sky_info[VISIT_ID==stack_grid[i,VISIT_ID] & FILTER==stack_grid[i,FILTER] & grepl('MIRIMAGE', DETECTOR),]
          CRVAL1 = module_MIRIMAGE[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL1]
          CRVAL2 = module_MIRIMAGE[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL2]
          CD1_1 = module_MIRIMAGE[VISIT_ID == stack_grid[i,VISIT_ID], CD1_1]
          CD1_2 = module_MIRIMAGE[VISIT_ID == stack_grid[i,VISIT_ID], CD1_2]
          NAXIS = 1200
          CRPIX = 600
        }
        
        if(j == 'NIS'){
          input_info = cal_sky_info[VISIT_ID==stack_grid[i,VISIT_ID] & FILTER==stack_grid[i,FILTER] & grepl('NIS', DETECTOR),]
          CRVAL1 = module_NIS[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL1]
          CRVAL2 = module_NIS[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL2]
          CD1_1 = module_NIS[VISIT_ID == stack_grid[i,VISIT_ID], CD1_1]
          CD1_2 = module_NIS[VISIT_ID == stack_grid[i,VISIT_ID], CD1_2]
          NAXIS = NAXIS_long
          CRPIX = NAXIS_long/2.0
        }
        
        if(j == 'NRCA_short'){
          input_info = cal_sky_info[VISIT_ID==stack_grid[i,VISIT_ID] & FILTER==stack_grid[i,FILTER] & grepl('NRCA', DETECTOR),]
          CRVAL1 = module_A_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL1]
          CRVAL2 = module_A_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL2]
          CD1_1 = module_A_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CD1_1]
          CD1_2 = module_A_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CD1_2]
          NAXIS = NAXIS_short
          CRPIX = NAXIS_short/2.0
        }
        
        if(j == 'NRCB_short'){
          input_info = cal_sky_info[VISIT_ID==stack_grid[i,VISIT_ID] & FILTER==stack_grid[i,FILTER] & grepl('NRCB', DETECTOR),]
          CRVAL1 = module_B_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL1]
          CRVAL2 = module_B_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL2]
          CD1_1 = module_B_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CD1_1]
          CD1_2 = module_B_WCS_short[VISIT_ID == stack_grid[i,VISIT_ID], CD1_2]
          NAXIS = NAXIS_short
          CRPIX = NAXIS_short/2.0
        }
        
        if(j == 'NRCA_long'){
          input_info = cal_sky_info[VISIT_ID==stack_grid[i,VISIT_ID] & FILTER==stack_grid[i,FILTER] & grepl('NRCA', DETECTOR),]
          CRVAL1 = module_A_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL1]
          CRVAL2 = module_A_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL2]
          CD1_1 = module_A_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CD1_1]
          CD1_2 = module_A_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CD1_2]
          NAXIS = NAXIS_long
          CRPIX = NAXIS_long/2.0
        }
        
        if(j == 'NRCB_long'){
          input_info = cal_sky_info[VISIT_ID==stack_grid[i,VISIT_ID] & FILTER==stack_grid[i,FILTER] & grepl('NRCB', DETECTOR),]
          CRVAL1 = module_B_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL1]
          CRVAL2 = module_B_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CRVAL2]
          CD1_1 = module_B_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CD1_1]
          CD1_2 = module_B_WCS_long[VISIT_ID == stack_grid[i,VISIT_ID], CD1_2]
          NAXIS = NAXIS_long
          CRPIX = NAXIS_long/2.0
        }
        
        temp_proj = Rwcs_keypass(CRVAL1 = CRVAL1,
                                 CRVAL2 = CRVAL2,
                                 CRPIX1 = CRPIX,
                                 CRPIX2 = CRPIX,
                                 CD1_1 = CD1_1,
                                 CD1_2 = CD1_2,
                                 CD2_1 = CD1_2,
                                 CD2_2 = -CD1_1
        )
        temp_proj$NAXIS1 = NAXIS
        temp_proj$NAXIS2 = NAXIS
        temp_proj[is.na(temp_proj)] = NULL
        
        temp_proj$DATE_BEG = as.character(min(as.POSIXlt.Date(input_info$`DATE-OBS`), na.rm=TRUE))
        temp_proj$DATE_END = as.character(max(as.POSIXlt.Date(input_info$`DATE-OBS`), na.rm=TRUE))
        
        image_list = {}
        inVar_list = {}
        
        for(k in 1:dim(input_info)[1]){
          temp_image = Rfits_read_image(input_info[k,full], ext=2)
          
          temp_mask = Rfits_read_image(input_info[k,full], ext=4, header=FALSE)
          
          temp_image = propaneBadPix(
            image = temp_image,
            patch = T
          )
          
          JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
          temp_image$imDat[JWST_cal_mask==1] = NA
          
          if(do_tweak){
            pro_test = suppressMessages(profoundProFound(
              image = temp_image, 
              mask = JWST_cal_mask == 1,
              sky = 0,
              redosky = F,
              pixcut = 10,
              skycut = 3.0,
              cliptol = 100,
              tolerance = 1,
              iters = 0,
              redosegim = F,
              box = dim(temp_image)[1]/10.0,
              rem_mask = T
            ))
            
            match_cat_idx = coordmatch(
              coordref = ref_cat[, c("RA", "Dec")],
              coordcompare = pro_test$segstats[,c("RAmax", "Decmax")],
              rad = 2.0
            )
            tweak_solution = suppressMessages(propaneTweakCat(
              cat_ref = ref_cat[match_cat_idx$bestmatch$refID, c("RA", "Dec")], 
              cat_pre_fix = pro_test$segstats[match_cat_idx$bestmatch$compareID,c("RAmax", "Decmax")], 
              mode = "coord", keyvalues_pre_fix = temp_image$keyvalues, delta_max = c(100, 100)
            ))
            temp_image = propaneWCSmod(
              temp_image, 
              delta_x = tweak_solution$par[1], 
              delta_y = tweak_solution$par[2],
              delta_rot = tweak_solution$par[3]
            )
          }
          
          rm(temp_mask)
          rm(JWST_cal_mask)
          image_list = c(image_list, list(temp_image))
          
          sky_file = paste0(sky_frames_dir, sub('_rem', '', input_info[k,stub]),'_',input_info[k,FILTER],'.fits')
          inVar_list = c(inVar_list,
                         Rfits::Rfits_read_key(sky_file, 'SKYRMS')^-2
          )
        }
        
        temp_file = paste0(invar_dir, '/temp_list_',i,'.fits')
        Rfits_write(image_list, temp_file)
        rm(image_list)
        image_list_point = Rfits_read(temp_file, pointer=TRUE)
        
        dump_stub = paste0(dump_dir_stub, stack_grid[i,VISIT_ID], "/", stack_grid$FILTER[i], "/", j, "/")
        if(dir.exists(dump_stub)){
          unlink(dump_stub, recursive = T)
        }
        dir.create(dump_stub, recursive = T)
        
        output_stack = propaneStackWarpInVar(
          # image_list = image_list,
          image_list = image_list_point,
          inVar_list = inVar_list,
          exp_list = input_info$EFFEXPTM,
          magzero_in = input_info$MAGZERO_FIX,
          magzero_out = magzero,
          keyvalues_out = temp_proj,
          cores = cores, 
          direction = 'backward',
          doclip = doclip,
          clip_tol = clip_tol,
          clip_dilate = clip_dilate,
          clip_sigma = clip_sigma,
          keep_extreme_pix = TRUE,
          dump_frames = T,
          dump_dir = dump_stub,
          multitype="cluster"
        )
        
        output_stack$image$keyvalues$CLIP_MIN = clip_tol[1]
        output_stack$image$keyvalues$CLIP_MAX = clip_tol[2]
        output_stack$image$keyvalues$CLIP_DIL = clip_dilate
        
        exp_min = min(output_stack$exp$imDat[output_stack$weight$imDat > 0], na.rm=TRUE)
        exp_max = max(output_stack$exp$imDat[output_stack$weight$imDat > 0], na.rm=TRUE)
        
        mask_frc = sum(output_stack$weight$imDat==0 & output_stack$exp$imDat >= exp_min, na.rm=TRUE) / length(which(output_stack$weight$imDat > 0))
        clip_frc = sum(output_stack$clip$imDat>0, na.rm=TRUE) / length(which(output_stack$weight$imDat > 0))
        
        output_stack$image$keyvalues$CLIP_FRC = clip_frc
        output_stack$image$keyvalues$MASK_FRC = mask_frc
        output_stack$image$keyvalues$EXP_MIN = exp_min
        output_stack$image$keyvalues$EXP_MAX = exp_max
        output_stack$image$keyvalues$MAGZERO = magzero
        output_stack$image$keyvalues$BUNIT = "microJy"
        
        rms_min = 1/sqrt(max(output_stack$inVar$imDat[output_stack$weight$imDat > 0], na.rm=TRUE))
        rms_max = 1/sqrt(min(output_stack$inVar$imDat[output_stack$weight$imDat > 0], na.rm=TRUE))
        
        output_stack$image$keyvalues$RMS_MIN = rms_min
        output_stack$image$keyvalues$RMS_MAX = rms_max
        
        output_stack = c(output_stack, list(info = input_info[,c(1,5:20)]))
        
        output_stack$image$keynames = unique(names(output_stack$image$keyvalues))
        output_stack$image$keycomments = rep(list(''), length(unique(names(output_stack$image$keyvalues))))
        names(output_stack$image$keycomments) = output_stack$image$keynames
        
        median_stack = propaneStackWarpMed(dirlist = dump_stub,
                                           pattern = glob2rx('*image_warp*'), 
                                           keyvalues_out = output_stack$image$keyvalues, 
                                           cores = cores,
                                           multitype = "cluster")
        
        file.remove(temp_file)
        
        filestub = paste0(invar_dir, '/stack_',stack_grid[i,VISIT_ID],'_',stack_grid[i,FILTER],'_',j)
        
        Rfits_write(data = output_stack,
                    filename = paste0(filestub,'.fits')
        )
        
        filestub_med = paste0(median_dir, '/med_',stack_grid[i,VISIT_ID],'_',stack_grid[i,FILTER],'_',j)
        
        Rfits_write(data = median_stack,
                    filename = paste0(filestub_med,'.fits')
        )
        gc()
      }
    }
  }
  
}
do_wisp_rem = function(input_args){
  
  cat("\n")
  message("## Removing wisps ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  filelist = input_args$filelist
  VID = input_args$VID
  median_dir = input_args$median_dir
  cores = input_args$cores_pro
  
  additional_params  = input_args$additional_params

  SIGMA_LO = input_args$SIGMA_LO
  message(paste0("Using ", ifelse(SIGMA_LO)))
  
  wisp_poly = Rfits_read("wisp_poly.fits")
  
  filelist = load_files(input_args, which_module = "wisp_rem")
  
  info = Rfits_key_scan(filelist = filelist,
                        keylist=c('DETECTOR', 'MODULE', 'FILTER', 'VISIT_ID'), cores = cores)
  
  
  wcs_info = Rfits_key_scan(filelist = filelist, extlist = 2,
                            keylist = c("CRVAL1", "CRVAL2")) ## Find longest overlapping wisp frame
  info = bind_cols(info, wcs_info[,c("CRVAL1", "CRVAL2")])
  
  if(additional_params$do_claws){
    info_wisp = info[DETECTOR %in% c('NRCA1', 'NRCA2', 'NRCA3','NRCA4', 'NRCB1','NRCB2', 'NRCB3','NRCB4'),] ## Try do all short wavelength chips
  }else{
    info_wisp = info[DETECTOR %in% c('NRCA3','NRCA4', 'NRCB3','NRCB4'),] 
  }
  mod_visit_grid = info_wisp[,c("MODULE", "VISIT_ID", "CRVAL1", "CRVAL2", "DETECTOR")]
  
  cat('Processing',dim(info_wisp)[1],'files\n')
  
  ref_im_list = {}
  for(ii in 1:dim(mod_visit_grid)[1]){
    
    ## Look for long channel in the same VID
    ref_files = list.files(
      median_dir, 
      pattern = glob2rx(paste0(
        "*", mod_visit_grid$VISIT_ID[ii], "*", mod_visit_grid$MODULE[ii], "*", "short", "*.fits"
      )),
      full.names = T)
    ref_files = ref_files[!grepl("F150W2|F070W|F090W|F115W|F150W|F200W|F140M|F162M|F182M|F210M|F164N|F187N|F212N", ref_files)] ## Remove the short wavelength filters
    
    if(length(ref_files) == 0){
      ## If no long channel exists, then perform a spatial search and use any overlapping long channel for Wisp Rem
      ref_files = propaneFrameFinder(
        dirlist = median_dir, 
        RAcen = mod_visit_grid$CRVAL1[ii],
        Deccen = mod_visit_grid$CRVAL2[ii],
        rad = 0.1/3600, ## match within 0.1 asec
        plot = FALSE, 
        cores = cores
      )$full ## Now do a frame finder to try and find long channel NOT in VIISITID (but could be in PROGRAM)
      ref_files = ref_files[grepl("short.+fits$", ref_files)]
    }
    
    filter_long = c(
      str_match(ref_files, "F\\s*(.*?)\\s*M")[,2],
      str_match(ref_files, "F\\s*(.*?)\\s*W")[,2]
      )
    filter_long[is.na(filter_long)] = -1
    long_num = max(filter_long, na.rm = T)
    
    ref_file_long = ref_files[ which(grepl(long_num, ref_files))[1] ] #Get the longest filter
    print(ref_file_long)
    
    ref_im_list = c(ref_im_list, list(Rfits_point(ref_file_long))) ## Allow loop in wisp rem to read long wavelength reference
    message(paste("Loading reference for:", paste0(mod_visit_grid$VISIT_ID[ii]), mod_visit_grid$DETECTOR[ii] ))
  }
  
  registerDoParallel(cores = cores)  
  temp = foreach(ii = 1:dim(info_wisp)[1], .errorhandling = "pass")%dopar%{
    ## copy original data to directory where we keep the wisps
    vid = paste0(info_wisp$VISIT_ID[ii])
    modl = paste0("NRC", info_wisp$MODULE[ii])
    message(paste0("Running: ", info_wisp$file[ii]))
    
    # file.copy(info_wisp$full[ii], paste0(keep_wisp_stub, "/", vid, "/", info_wisp$file[ii]), overwrite = T)
    wisp_frame = Rfits_read_image(info_wisp$full[ii], ext = 2)
    
    check_Nhdu = Rfits_nhdu(info_wisp$full[ii])
    extloc = Rfits_extname_to_ext(info_wisp$full[ii], 'SCI_ORIG')
    
    if(is.na(extloc)){ ## only make the SCI_ORIG if it doesn't already exist, otherwise we risk overwriting with something that isn't the original science frame!
      Rfits_write_image(wisp_frame, filename=info_wisp$full[ii], create_file=FALSE)
      Rfits_write_key(filename=info_wisp$full[ii], ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SCI_ORIG', keycomment='No wisp corr')
    }
    
    # Rfits_write_image(data = wisp_frame, paste0(keep_wisp_stub, "/", vid, "/", info_wisp$file[ii]))
    
    ref_im = ref_im_list[[ii]] ## read in long wavelength reference 
    
    if(any( (vid == additional_params$ID_vlarge$VISIT_ID & modl == additional_params$ID_vlarge$MODULE) | 
            (vid == additional_params$ID_large$VISIT_ID & modl == additional_params$ID_large$MODULE)) ){
      sigma_lo = NULL
    }else{
      sigma_lo = SIGMA_LO
    }
    #poly = wisp_ploy[[ info_wisp$DETECTOR[ii] ]]
    poly = NULL    
    
    wisp_fix = wispFixer(wisp_im = wisp_frame, ref_im = ref_im, poly = poly, sigma_lo = sigma_lo)
    Rfits_write_pix(data = wisp_fix$wisp_fix$imDat, filename = info_wisp$full[ii], ext = 2)
    
    check_Nhdu = Rfits_nhdu(info_wisp$full[ii])
    extloc = Rfits_extname_to_ext(info_wisp$full[ii], 'REF_WISP_WARP')
    if(is.na(extloc)){ 
      Rfits_write_image(wisp_fix$ref_im_warp, filename=info_wisp$full[ii], create_file=FALSE)
      Rfits_write_key(filename=info_wisp$full[ii], ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='REF_WISP_WARP', keycomment='Reference long wavelength for wisp corr')
    }else{
      Rfits_write_pix(wisp_fix$ref_im_warp$imDat, filename = info_wisp$full[ii], ext = extloc)
    }
    
    check_Nhdu = Rfits_nhdu(info_wisp$full[ii])
    extloc = Rfits_extname_to_ext(info_wisp$full[ii], 'WISP_TEMPLATE')
    if(is.na(extloc)){
      Rfits_write_image(wisp_fix$wisp_template, filename=info_wisp$full[ii], create_file=FALSE)
      Rfits_write_key(filename=info_wisp$full[ii], ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='WISP_TEMPLATE', keycomment='Wisp template')
    }else{
      Rfits_write_pix(wisp_fix$wisp_template, filename = info_wisp$full[ii], ext = extloc)
    }
    return(NULL)
  }
}
do_miri_bkgnd = function(input_args){
  ## Remove interference substructure in MIRI images 
  ## algo inspired by https://arxiv.org/pdf/2405.15972
  
  cat("\n")
  message("## Removing MIRI backgrounds ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  filelist = load_files(input_args = input_args, which_module = "modify_pedestal") ## cal_sky
  VID = input_args$VID
  FILT = input_args$FILT
  median_dir = input_args$median_dir
  cores = input_args$cores_pro
  
  registerDoParallel(cores = cores)
  
  cal_info_raw = Rfits_key_scan(
    filelist = filelist,
    extlist = 1,
    keylist = c("FILTER", "VISIT_ID")
  )
  
  cal_info = cal_info_raw[
    grepl(FILT, cal_info_raw$FILTER) & grepl(VID, cal_info_raw$VISIT_ID), 
  ]
  
  for(VID in unique(cal_info$VISIT_ID)){
    for(FILT in unique(cal_info$FILTER)){
      
      bkgnd_grid = cal_info[
        cal_info$VISIT_ID == VID & cal_info$FILTER == FILT, 
      ]
      message("Running VISIT: ", VID, ", FILTER: ", FILT)
      message("...Processing ", dim(bkgnd_grid)[1], " files...")
      if(dim(bkgnd_grid)[1] == 0){
        next
      }
      ff_stack = Rfits_point(
        filename = paste0(median_dir, "/med_", VID, "_", FILT, "_MIRIMAGE.fits")
      )
      
      frames = lapply(
        bkgnd_grid$full, 
        function(x){
          ff = Rfits_read(x)
          
          ff_stackw = suppressMessages(propaneWarp(
            image_in = ff_stack,
            keyvalues_out = ff$SCI$keyvalues,
            magzero_in = 23.9, magzero_out = 23.9
          ))
          pro_stackw = suppressMessages(profoundProFound(
            image = ff_stackw, skycut = 10.0, pixcut = 9, magzero = 23.9, rem_mask = TRUE
          ))
          
          temp_cal = ff$SCI$imDat
          JWST_mask = profoundDilate(ff$DQ$imDat %% 2 == 1, size = 3)
          
          temp_cal[JWST_mask == 1L | pro_stackw$objects == 1] = NA
          temp_cal = temp_cal - median(temp_cal, na.rm = TRUE)
          
          return(temp_cal)
        }
      )
      
      ## Make median stack with sources clipped out
      median_stack = propaneStackFlatMed(
        image_list = frames, 
        na.rm = TRUE
      )
      
      ## Blur the median stack - interpolate
      med_blur = profoundImBlur(
        image = median_stack$image, 
        sigma = 1.0
      )

      sel = which(is.na(median_stack$image))
      median_stack$image[sel] = med_blur[sel]
      
      temp = foreach(i = 1:dim(bkgnd_grid)[1], .inorder = FALSE) %dopar% {
        
        message(
          "Adjusting MIRI background: ", bkgnd_grid$stub[i]
        )
        
        fname = bkgnd_grid$full[i]
        check_Nhdu = Rfits_nhdu(fname)
        extloc = Rfits_extname_to_ext(fname, 'SCI_ORIG')
        
        ff = Rfits_read(fname)
        temp_cal = ff$SCI
        temp_dq = Rfits_read_image(fname, ext = 4, header = FALSE)
        
        ff_stackw = suppressMessages(propaneWarp(
          image_in = ff_stack,
          keyvalues_out = ff$SCI$keyvalues,
          magzero_in = 23.9, magzero_out = 23.9
        ))
        pro_stackw = suppressMessages(profoundProFound(
          image = ff_stackw, magzero = 23.9, rem_mask = TRUE
        ))
        
        ## Save original extension 2
        if(is.na(extloc)){ ## only make the SCI_ORIG if it doesn't already exist, otherwise we risk overwriting with something that isn't the original science frame!
          Rfits_write_image(temp_cal, filename=fname, create_file=FALSE)
          Rfits_write_key(filename=fname, ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SCI_ORIG', keycomment='No background corr')
        }
        
        ## Remove median stack sky with ProFound
        pro_new = profoundProFound(
          temp_cal$imDat, 
          sky = median_stack$image,
          segim = pro_stackw$segim,
          objects = pro_stackw$objects_redo,
          mask = pro_stackw$objects_redo == 1 | profoundDilate(temp_dq %% 2 == 1 , size = 3)
        )
        
        bkgnd_fix = pro_new$image - pro_new$sky
        
        ## Remove residual 2D structure
        pro_tram = profoundSkyScan(
          image = bkgnd_fix,
          mask = profoundDilate(temp_dq %% 2 == 1, size = 3.0),
          scan_block = c(nrow(bkgnd_fix), ncol(bkgnd_fix)), 
          clip = c(0.0, 0.9), 
          keep_trend = FALSE,
          trend_block = 101
        )
        ## Overwrite extension 2
        Rfits_write_pix(data = pro_tram$image_fix, filename = fname, ext = 2)
        
        return( NULL )
      }
    }
  }
  
}
do_patch = function(input_args){
  
  cat("\n")
  message("## Patching ProPane stacks ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  VID = input_args$VID
  FILT = input_args$FILT
  invar_dir = input_args$invar_dir
  median_dir = input_args$median_dir
  patch_dir = input_args$patch_dir
  cores = input_args$cores_stack
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  
  registerDoParallel(cores = cores)
  patch_stub = patch_dir
  
  if(do_NIRISS){
    pixscale_list = c("NIS")
  }else if(do_MIRI){
    pixscale_list = c("MIRIMAGE")
  }else{
    pixscale_list = c("short", "long")
  }
  for(j in pixscale_list){
    invar_files = list.files(path=invar_dir, pattern = paste0(j, ".fits$"), full.names = T)
    median_files = list.files(path=median_dir, pattern = paste0(j, ".fits$"), full.names = T)
    invar_files=invar_files[grepl(VID, invar_files) & grepl(FILT, invar_files)]
    median_files=median_files[grepl(VID, median_files) & grepl(FILT, median_files)]
    
    if(length(median_files)==0){
      message("No frames") 
    }else{
      file_names = list.files(path=invar_dir, pattern = paste0(j, ".fits$"), full.names = F)
      file_names=file_names[grepl(VID, file_names) & grepl(FILT, file_names)]
      
      foreach(i = 1:length(invar_files))%dopar%{
        message(paste0("Patching: ", invar_files[i]))
        load_invar = Rfits_read(invar_files[i], pointer=F)
        load_median = Rfits_read(median_files[i], pointer=F)
        
        temp = load_invar$image
        fstub = strsplit(file_names[i], "stack")[[1]][2]
        med_patch=propanePatch(image_inVar = temp, 
                               image_med = load_median$image)
        
        ## compute the bounding box
        edge_mask = matrix(0L, nrow = dim(temp)[1], ncol = dim(temp)[2])
        rough_mask = is.na(temp$imDat)
        mask_label = as.matrix(imager::label(imager::as.cimg(rough_mask)))
        tally_pixels = tabulate(mask_label)
        sel = which(tally_pixels > 0)
        mask_loc = which(matrix(mask_label %in% sel, dim(temp)[1], dim(temp)[2]), arr.ind=TRUE)
        edge_mask[mask_loc] = 1
        edge_mask = 1- edge_mask
        
        big_NA_mask = matrix(0L, nrow = dim(temp)[1], ncol = dim(temp)[2])
        rough_mask = is.na(temp$imDat)
        mask_label = as.matrix(imager::label(imager::as.cimg(rough_mask)))
        tally_pixels = tabulate(mask_label)
        sel = which(tally_pixels > 150 & tally_pixels < max(tally_pixels))
        mask_loc = which(matrix(mask_label %in% sel, dim(temp)[1], dim(temp)[2]), arr.ind=TRUE)
        big_NA_mask[mask_loc] = 1
        # magimage(big_NA_mask)
        
        propane_mask_fix = propanePatchPix(
          image = med_patch$image, 
          mask = (edge_mask + big_NA_mask) == 1,
        )
        
        patch = list()
        patch$image = propane_mask_fix
        patch$patch = temp - propane_mask_fix
        patch$inVar = load_invar$inVar
        patch$weight = load_invar$weight
        patch$image$keyvalues$EXTNAME = "image"
        patch$patch$keyvalues$EXTNAME = "patch"
        patch$inVar$keyvalues$EXTNAME = "inVar"
        patch$weight$keyvalues$EXTNAME = "weight"
        
        Rfits_write_all(patch, paste0(patch_stub, "/stack_patch", fstub))
        
      }
      stopImplicitCluster()
    }
  }
}
do_RGB = function(input_args){
  
  cat("\n")
  message("## Making RGB image ##")
  cat("\n")
  Sys.sleep(time = 5)
  
  VID = input_args$VID
  ref_dir = input_args$ref_dir
  patch_dir = input_args$patch_dir
  
  blue_filters = "F070W|F090W|F115W|F150W|F140M|F162M"
  green_filters = "F200W|F277W|F182M|F210M"
  red_filters = "F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430|F460M|F480M|F560W|F770W|F1000W|F1130W|F1280W|F1500W|F1800W|F2100W|F2550W"
  locut = 1e-6
  hicut = 0.05
  
  cal_sky_info = fread(paste0(ref_dir, "/Pro1oF/cal_sky_info.csv"))
  vid_temp = VID
  unique_visits = grep(vid_temp, unique(cal_sky_info$VISIT_ID), value = T)
  
  for(VID in unique_visits){
    
    RGB_dir = paste0(ref_dir, "/RGB/", VID, "/")
    dir.create(RGB_dir, showWarnings = F, recursive = T)
    
    file_list = list.files(patch_dir, 
                           pattern = glob2rx(paste0("*", VID, "*short*fits")),
                           full.names = T)
    
    frame_info = Rfits_key_scan(filelist = file_list,
                                keylist = c("CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "NAXIS1", "NAXIS2"))
    temp_proj = Rwcs_keypass(CRVAL1 = mean(frame_info$CRVAL1),
                             CRVAL2 = mean(frame_info$CRVAL2),
                             CD1_1 = min(frame_info$CD1_1),
                             CD1_2 = min(frame_info$CD1_2),
                             CD2_1 = min(frame_info$CD1_2),
                             CD2_2 = -min(frame_info$CD1_1)
    )
    
    module_list = sapply(strsplit(file_list, "NRC"), function(x)(x[2]))
    
    if(length(unique(module_list))==2){
      temp_proj$NAXIS1 = 2.2 * max(frame_info$NAXIS1) #Modules in NIRCam always colinear
    }else{
      temp_proj$NAXIS1 = 1.2 * max(frame_info$NAXIS1)
    }
    temp_proj$NAXIS2 = 1.2 * max(frame_info$NAXIS1)
    temp_proj$CRPIX1 = temp_proj$NAXIS1/2.0
    temp_proj$CRPIX2 = temp_proj$NAXIS2/2.0
    temp_proj[is.na(temp_proj)] = NULL
    
    message("Making RGBs from")
    cat(file_list, sep="\n")
    blue_files = file_list[grepl(blue_filters, file_list)]
    green_files = file_list[grepl(green_filters, file_list)]
    red_files = file_list[grepl(red_filters, file_list)]
    
    if(length(red_files>0)){
      red_images = lapply(red_files, function(x){Rfits_read(x)$image[,]})
    }else{
      red_images = list(Rfits_create_image(matrix(0L, 
                                                  nrow=temp_proj$NAXIS1, 
                                                  ncol = temp_proj$NAXIS2), 
                                           keyvalues = temp_proj))
    }
    if(length(green_files>0)){
      green_images = lapply(green_files, function(x){Rfits_read(x)$image[,]})
    }else{
      green_images = list(Rfits_create_image(matrix(0L, 
                                                    nrow=temp_proj$NAXIS1, 
                                                    ncol = temp_proj$NAXIS2), 
                                             keyvalues = temp_proj))
    }
    if(length(blue_files>0)){
      blue_images = lapply(blue_files, function(x){Rfits_read(x)$image[,]})
    }else{
      blue_images = list(Rfits_create_image(matrix(0L, 
                                                   nrow=temp_proj$NAXIS1, 
                                                   ncol = temp_proj$NAXIS2), 
                                            keyvalues = temp_proj))
    }
    
    R_stack = propaneStackWarpInVar(red_images, keyvalues_out = temp_proj, magzero_in = 23.9, magzero_out = 23.9, cores=1)
    G_stack = propaneStackWarpInVar(green_images, keyvalues_out = temp_proj, magzero_in = 23.9, magzero_out = 23.9, cores=1)
    B_stack = propaneStackWarpInVar(blue_images, keyvalues_out = temp_proj, magzero_in = 23.9, magzero_out = 23.9, cores=1)
    
    temp = R_stack$image
    edge_mask = matrix(1L, nrow = dim(temp)[1], ncol = dim(temp)[2])
    rough_mask = is.na(temp$imDat)
    mask_label = as.matrix(imager::label(imager::as.cimg(rough_mask)))
    tally_pixels = tabulate(mask_label)
    sel = which(tally_pixels > 0)
    mask_loc = which(matrix(mask_label %in% sel, dim(temp)[1], dim(temp)[2]), arr.ind=TRUE)
    edge_mask[mask_loc] = 0L
    
    R_patch = propanePatchPix(R_stack$image[,], mask = edge_mask)
    G_patch = propanePatchPix(G_stack$image[,], mask = edge_mask)
    B_patch = propanePatchPix(B_stack$image[,], mask = edge_mask)
    
    Rfits_write_image(R_patch, paste0(RGB_dir, "/", VID, "_R.fits"))
    Rfits_write_image(G_patch, paste0(RGB_dir, "/", VID, "_G.fits"))
    Rfits_write_image(B_patch, paste0(RGB_dir, "/", VID, "_B.fits"))
    
    im_width = 12 #inches
    im_height = im_width * (temp_proj$NAXIS2 / temp_proj$NAXIS1)
    
    filestub = paste0(patch_dir, '/patch_RGB_stack_',VID, "_all", "_clear.png")
    CairoPNG(filename = filestub, width=im_width, height=im_height, 
             units='in', res=400, quality=100)
    par(mar=c(0,0,0,0), oma=rep(0,4))
    Rwcs_imageRGB(
      R = R_patch,
      G = G_patch,
      B = B_patch,
      locut=locut,
      hicut=hicut,
      type='num',
      sparse=1,
      decorate=F)
    dev.off()
  }

}

wispFixer = function(wisp_im, ref_im,
                     source_threshold=0.95, #Threshold on reference image for definitely real sources
                     scale_threshold=0.95, #When considering the +ve ratio map pixels, quantile to define very blue things
                     sky_threshold=0.5, #Safe threshold to find the correct sky level
                     clip_threshold=0.997, #We avoid very bright pixels in the wisp template
                     sigma_hi=2,
                     sigma_lo=20,
                     poly=NULL,
                     cores=1){
  
  ref_im$imdat[is.infinite(ref_im$imdat)] = NA #remove infinities
  wisp_im$imDat[is.infinite(wisp_im$imDat)] = NA
  
  #Step 1
  ref_im_warp_untweak = propaneWarp(ref_im, keyvalues_out=wisp_im$keyvalues, cores=cores, direction = "forward")
  
  ## 1i Align ref to wisp_im using propaneTweak
  ## Use only 1 core in case resources run out
  ref_im_warp = propaneTweakImage(
    image_ref = wisp_im[,], image_pre_fix = ref_im_warp_untweak[,], 
    delta_max = c(2, 0.1), quan_cut = c(0.995, 0.9999),
    shift_int = F, quick = T, verbose = F, cores = 1 
  )$image_post_fix
  
  
  #Step X (not in paper, but put here for sky sub reasons)
  real_source = quantile(ref_im_warp$imDat, source_threshold, na.rm=TRUE)
  sky_level = quantile(wisp_im$imDat[ref_im_warp$imDat < real_source], sky_threshold, na.rm=TRUE)
  wisp_im$imDat = (wisp_im$imDat - sky_level)
  
  #Step 2i
  relflux = (wisp_im$imDat / ref_im_warp$imDat)
  
  #Step 2ii
  #real_source = quantile(ref_im_warp$imDat, source_threshold, na.rm=TRUE) paper location
  relflux = relflux[ref_im_warp$imDat > real_source] #just looks at real sources
  
  #step 2iii
  relflux = relflux[relflux > 0] #ignore -ve noise
  
  #Step 2iv
  scale = quantile(relflux, scale_threshold, na.rm=T) #find blue things
  
  #Step 3
  wisp_template = wisp_im$imDat - ref_im_warp$imDat*scale
  
  #Step 4
  wisp_template[wisp_template > quantile(wisp_template[wisp_template > 0], clip_threshold, na.rm=T)] = NA
  wisp_template[wisp_template < quantile(wisp_template[wisp_template < 0], 1 - clip_threshold, na.rm=T)] = NA
  
  #Step 5
  wisp_template = profoundImBlur(wisp_template, sigma=sigma_hi)
  
  #Step 6
  wisp_template[is.na(wisp_im$imDat) | is.na(ref_im_warp$imDat)] = NA
  wisp_template[wisp_template < 0] = NA
  
  if(!is.null(sigma_lo)){
    #Step 7
    wisp_template_lo = profoundImBlur(wisp_template, sigma=sigma_lo)
    
    #Step 8
    sel = which(is.na(wisp_template) | wisp_template_lo > wisp_template)
    wisp_template[sel] = wisp_template_lo[sel]
  }else{
    wisp_template[is.na(wisp_template)] = 0
  }
  
  if(!is.null(poly)){
    wisp_mask = profoundDrawMask(wisp_im$imDat, poly=poly, mode='apply')$mask
    # wisp_template = magmap(wisp_template * profoundImBlur(image = wisp_mask, sigma = 50), range = c(0,1))$map
    wisp_template[wisp_mask==0L] = 0
  }
  
  #Just make sure we don't create non-zero values where we have missing data
  #Not sure where that's better or worse?
  #wisp_template[is.na(wisp_im$imDat) | is.na(ref_im_warp$imDat)] = 0
  
  #Step 9
  wisp_im$imDat  = wisp_im$imDat - wisp_template + sky_level
  
  return(list(wisp_fix=wisp_im, wisp_template=wisp_template, ref_im_warp=ref_im_warp))
}

