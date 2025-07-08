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
library(checkmate)

pipe_version = "1.4.1" 

load_files = function(input_args, which_module, sky_info = NULL){
  ## Load the correct files for what ever task
  ## Correctly place files into whatever module
  ## Which module tells me which step to load files for
  
  VID = input_args$VID
  FILT = input_args$FILT
  
  Pro1oF_dir = input_args$Pro1oF_dir
  sky_pro_dir = input_args$sky_pro_dir
  sky_frames_dir = input_args$sky_frames_dir
  sky_super_dir = input_args$sky_super_dir
  cal_sky_dir = input_args$cal_sky_dir
  
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  
  if(which_module == "1oF"){
    files_1oF = input_args$filelist
    files_1oF = files_1oF[grepl(VID, files_1oF) & grepl(".fits$", files_1oF)]
    if(do_NIRISS){
      files_1oF = files_1oF[!grepl('_mirimage_',files_1oF) & grepl('_nis_',files_1oF) & grepl(".fits$", files_1oF)]
    }else if (do_MIRI){
      files_1oF = files_1oF[grepl('_mirimage_',files_1oF) & !grepl('_nis_',files_1oF) & grepl(".fits$", files_1oF)]
    }else{
      files_1oF = files_1oF[!grepl('_mirimage_',files_1oF) & !grepl('_nis_', files_1oF) & grepl(".fits$", files_1oF)]
    }
    
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
    
    scan_cal = Rfits_key_scan(filelist = files_cal, keylist = c("FILTER", "PROGRAM", "VISIT_ID"))
    corr_pid = grep(VID, scan_cal$VISIT_ID, fixed = T, value = T)
    corr_pid = corr_pid[substring(corr_pid, 1, nchar(VID)) == VID] ## Make sure 4 digit PID is embedded in 10 digit VID
    pid_idx = scan_cal$VISIT_ID %in% corr_pid
    
    files_cal = files_cal[pid_idx]
    files_cal = files_cal[grepl(FILT, scan_cal$FILTER[pid_idx])]
    return(files_cal)
  }
  
  if(which_module == "super_sky"){
    # sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
    # sky_super_dir = paste0(sky_pro_dir, "/sky_super/")
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
    }else{
      sky_info = sky_info[!grepl("MIRIMAGE", sky_info$detector) & !grepl("NIS", sky_info$detector), ]
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

## Processing codes
do_1of = function(input_args){
  
  cat("\n")
  message("## Removing 1/f ##")
  cat("\n")

  additional_params = input_args$additional_params
  Pro1oF_dir = input_args$Pro1oF_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  parallel_type = input_args$additional_params$parallel_type
  
  do_MIRI = input_args$do_MIRI ## Because we need to alter the dimensions of the image for processing
  trend_block_vlarge = 101 #good for larger sources
  trend_block_large = 501 # good for frames with fairly large sources
  
  skip_completed_files = additional_params$skip_completed_files
  ID_vlarge = additional_params$ID_vlarge
  ID_large = additional_params$ID_large
  
  ow_vlarge = additional_params$ow_vlarge
  ow_large = additional_params$ow_large
  
  filelist = load_files(input_args, which_module = "1oF")
  
  message("Showing first <10 files:")
  cat(head(filelist, 10), sep='\n')
  cat("...")
  cat('Processing',length(filelist),'files\n')
  message("Skip completed files is ", ifelse(skip_completed_files, "ON", "OFF"))
  
  lo_loop = 1
  hi_loop = length(filelist)
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = parallel_type)
    registerDoParallel(cl)
  }

  dummy = foreach(i = lo_loop:hi_loop, .inorder = FALSE, .packages = c("Rfits", "Rwcs", "ProFound"))%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    file_image = filelist[i]
    basename = strsplit(basename(file_image),'.fits$')[[1]]
    fullbase = paste(Pro1oF_dir,basename,sep='/')
    
    if( !(skip_completed_files & file.exists(paste0(fullbase,'_Pro1oF.fits'))) ){
      temp_image = Rfits_read(filelist[i], pointer=FALSE)
  
      if(do_MIRI){
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
          filename = paste0(fullbase,'_MIRI.fits')
        )
        Rfits_write_pix(temp_zap$image_fix, paste0(fullbase,'_MIRI.fits'), ext=2)
  
        check_Nhdu = Rfits_nhdu(paste0(fullbase,'_MIRI.fits'))
        extloc = Rfits_extname_to_ext(paste0(fullbase,'_MIRI.fits'), 'SKY_Pro1oF')
  
        if(is.na(extloc)){
          Rfits_write_image(temp_zap$row_map + temp_zap$col_map, filename=paste0(fullbase,'_MIRI.fits'), create_file=FALSE)
          Rfits_write_key(filename=paste0(fullbase,'_MIRI.fits'), ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SKY_Pro1oF', keycomment='extension name')
        }else{
          Rfits_write_pix(temp_zap$row_map + temp_zap$col_map, paste0(fullbase,'_MIRI.fits'), ext=extloc)
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
                                   mask = (temp_image$SCI$imDat==0) | JWST_cal_mask | is.na(temp_image$SCI$imDat),
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
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
}
do_cal_process = function(input_args, filelist = NULL){
  
  cat("\n")
  message("## Calculating sky statistics ##")
  cat("\n")

  Pro1oF_dir = input_args$Pro1oF_dir
  sky_frames_dir = input_args$sky_frames_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  parallel_type = input_args$additional_params$parallel_type
  skip_completed_files = additional_params$skip_completed_files
  
  if(is.null(filelist)){
    filelist = (load_files(input_args, which_module = "cal_process"))
  }
  
  message("Showing first <10 files:")
  cat(head(filelist, 10), sep='\n')
  cat("...")
  cat('Processing',length(filelist),'files...\n')
  message("Skip completed files is ", ifelse(skip_completed_files, "ON", "OFF"))
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = parallel_type)
    registerDoParallel(cl)
  }
  
  #get main info:
  obs_info = Rfits_key_scan(
    filelist = filelist,
    keylist = c('VISIT_ID','OBS_ID','EXPOSURE','DETECTOR','FILTER','DATE-OBS','EFFEXPTM','CAL_VER','CRDS_CTX'),
    extlist = 1,
    fileinfo = 'all'
  )
  obs_info = data.frame(obs_info)
  
  lo_loop = 1
  hi_loop = dim(obs_info)[1]
  
  dummy = foreach(i = lo_loop:hi_loop, .inorder=FALSE, .packages = c("Rfits", "Rwcs", "ProFound"), .export = c("pipe_version"))%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    #setup files
    file_image = as.character(obs_info[i,"full"])
    basename = strsplit(basename(file_image),'.fits',fixed=T)[[1]]
    fullbase = paste(sky_frames_dir,basename,sep='/')
    file_sky = paste0(fullbase,'_sky_',obs_info[i,"FILTER"],'.fits')
    
    if( !(skip_completed_files & file.exists(file_sky)) ){
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
                                   tolerance=Inf 
            )
            return(list(pro = pro, redoskysize=101))
          },
          error = function(cond){
            redoskysize = 21
            message("Adjusting redoskysize=", redoskysize, " from 101: ", basename)
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
            message("Adjusting object mask for polynomial skies")
            sky_redo = profoundSkyPoly(x$image, objects=x$pro$objects, degree=x$sky_poly_deg, quancut=0.99, mask=x$mask)
            return(sky_redo)
          },
          finally = {
            NULL
          }
        )
      }
      
      if(do_MIRI){
        ## Do not do polynomial sky for MIRI - too hard
        sky_redo = list(sky = pro$sky)
        pro_redo = pro
      }else{
        sky_redo = sky_redo_func(list(image = JWST_cal_image$imDat, pro = pro, mask = JWST_cal_mask, sky_poly_deg = sky_poly_deg))
        pro_redo = profoundProFound(JWST_cal_image$imDat, mask=JWST_cal_mask, skycut=2, pixcut=5, box=box, grid = box, redoskysize=redoskysize, sky=sky_redo$sky, redosky=FALSE, tolerance=Inf)
      }
      
      if(is.null(sky_redo$sky)){
        ## if the polynomial sky is NULL, replace with the original profound sky
        ## If NULL it probably means that the polynomial sky had some problem as per the try-catch
        sky_redo$sky = pro$sky
      }
      if(is.null(pro_redo$skyRMS)){
        ## if the polynomial skyRMS is NULL, replace with the original profound skyRMS
        pro_redo$skyRMS = pro$skyRMS 
      }
      if(is.null(pro_redo$objects_redo)){
        pro_redo$objects_redo = pro$objects_redo 
      }
      
      if(is.null(pro_redo$skyChiSq)){
        pro_redo$skyChiSq = 1e6
      }
      if(is.na(pro_redo$skyChiSq)){
        pro_redo$skyChiSq = 1e6
      }
      if(is.infinite(pro_redo$skyChiSq)){
        pro_redo$skyChiSq = 1e6
      }
      
      if(pro_redo$skyChiSq >= 1e6){
        if(sum(pro_redo$objects_redo == 1) / prod(dim(pro_redo$objects_redo)) > 0.9){
          pro_redo$objects_redo = pro$objects 
        }else{
          pro_redo$objects_redo = pro$objects_redo
        }
      }
      
      # if(sum(pro_redo$objects_redo == 1)/prod(dim(pro_redo$objects_redo)) > 0.9){ ## Need to at least have 10% sky pixels?
      #   pro_redo$objects_redo = pro_redo$objects ## polynomial skies way below RMS of the true sky, whole frame is objects. Open clusters for example...
      # }
      
      sky_med = median(sky_redo$sky[JWST_cal_mask == 0 & pro_redo$objects_redo == 0], na.rm=TRUE)
      skyRMS_med = median(pro_redo$skyRMS[JWST_cal_mask == 0 & pro_redo$objects_redo == 0], na.rm=TRUE)
      
      ## Sky statistics being NULL/NA/Inf likely problem with the array slicing above, e.g., if the whole frame has objects 
      if(is.null(sky_med)){
        sky_med = 0.0 
      }
      if(is.na(sky_med)){
        sky_med = 0.0 
      }
      if(is.infinite(sky_med)){
        sky_med = 0.0 
      }
      
      if(is.null(skyRMS_med)){
        skyRMS_med = 1.0 
      }
      if(is.na(skyRMS_med)){
        skyRMS_med = 1.0 
      }
      if(is.infinite(skyRMS_med)){
        skyRMS_med = 1.0 
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
        SKY = sky_med,
        SKYRMS = skyRMS_med,
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
    }else{
      return(NULL)
    }
    
  }
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
}
do_regen_sky_info = function(input_args){
  
  cat("\n")
  message("## Generating sky statistics info ##")
  cat("\n")
  
  sky_pro_dir = input_args$sky_pro_dir
  sky_frames_dir = input_args$sky_frames_dir
  
  cores = input_args$cores_pro
  
  registerDoParallel(cores=cores)
  
  # sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
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

  sky_pro_dir = input_args$sky_pro_dir
  sky_frames_dir = input_args$sky_frames_dir
  sky_super_dir = input_args$sky_super_dir
  VID = input_args$VID
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  parallel_type = input_args$additional_params$parallel_type

  sky_ChiSq_cut = 1.1
  good_pix_cut = 0.15
  
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
  
  combine_grid_short = expand.grid(short_det, short_filt, stringsAsFactors=FALSE)
  combine_grid_long = expand.grid(long_det, long_filt, stringsAsFactors=FALSE)
  combine_grid_nis = expand.grid(niriss_det, nis_filt, stringsAsFactors=FALSE)
  combine_grid_miri = expand.grid(miri_det, miri_filt, stringsAsFactors=FALSE)
  combine_grid = rbind(combine_grid_short, combine_grid_long, combine_grid_nis, combine_grid_miri)
  
  sky_filelist = list.files(sky_info[1,pathsky], pattern = ".fits$")
  
  sky_info = data.frame(sky_info)
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = parallel_type)
    registerDoParallel(cl)
  }  
  
  dummy = foreach(i = 1:dim(combine_grid)[1], .inorder=FALSE, .packages = c("Rfits", "Rwcs", "ProFound"))%dopar%{
    if(combine_grid[i,1] %in% c('NRCA1','NRCA2','NRCA3','NRCA4','NRCB1','NRCB2','NRCB3','NRCB4')){
      temp_info = sky_info[sky_info$detector == combine_grid[i,1] & sky_info$filter == combine_grid[i,2] & sky_info$skychi < sky_ChiSq_cut & sky_info$goodpix >= good_pix_cut, ]
    }else if (combine_grid[i,1] %in% c('MIRIMAGE')){
      temp_info = sky_info[sky_info$detector == combine_grid[i,1] & sky_info$filter == combine_grid[i,2] & sky_info$skychi < sky_ChiSq_cut*1.5, ] ## Allow a little bit more legroom for MIRI?
    }else{
      temp_info = sky_info[sky_info$detector == combine_grid[i,1] & sky_info$filter == combine_grid[i,2] & sky_info$skychi < sky_ChiSq_cut, ]
    }
    
    Nsky = dim(temp_info)[1]
    
    message(combine_grid[i,1], ' ', combine_grid[i,2],': ', Nsky)
    
    if(Nsky > 0){
      sky_array = array(dim = c(temp_info[1,'naxis1'], temp_info[1,'naxis2'], Nsky)) #create empty array
      for(j in 1:dim(temp_info)[1]){
        #safe name (since sometimes the sky file name gets truncated due to 80 character FITS limits)
        file_sky = temp_info[j,'filesky']
        file_sky = grep(temp_info[j,'filesky'], sky_filelist, value=TRUE)
        if(length(file_sky) == 0){
          stop('Missing',)
        }
        
        #fill array up with sky data for stacking
        sky_array[,,j] = Rfits_read_image(paste(temp_info[j,'pathsky'], file_sky, sep='/'), ext=2, header=FALSE) - temp_info[j,'sky']
      }
      sky_mean = rowMeans(sky_array, na.rm = FALSE, dims = 2)
    }else{
      message("No usable data for ", combine_grid[i,1]," ",combine_grid[i,2])
      temp_file = sky_info[sky_info$detector == combine_grid[i,1] & sky_info$filter == combine_grid[i,2],'fileim'][1]
      NAXIS1 = sky_info[sky_info$detector == combine_grid[i,1] & sky_info$filter == combine_grid[i,2],'naxis1'][1]
      NAXIS2 = sky_info[sky_info$detector == combine_grid[i,1] & sky_info$filter == combine_grid[i,2],'naxis2'][1]
      if(is.na(NAXIS1) | is.na(NAXIS2)){
        sky_mean = NULL
      }else{
        sky_mean = matrix(0, NAXIS1, NAXIS2)
      }
    }
    
    if(!is.null(sky_mean)){
      
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
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
}
do_apply_super_sky = function(input_args){
  
  cat("\n")
  message("## Applying super skies ##")
  cat("\n")

  Pro1oF_dir = input_args$Pro1oF_dir
  cal_sky_dir = input_args$cal_sky_dir
  sky_pro_dir = input_args$sky_pro_dir
  sky_super_dir = input_args$sky_super_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  parallel_type = input_args$additional_params$parallel_type
  skip_completed_files = additional_params$skip_completed_files
  
  sky_info = fread(paste0(sky_pro_dir, "/sky_info.csv"))
  # sky_super_dir = paste0(sky_pro_dir, "/sky_super/")
  
  setcal_sky_remfile = function(file_image, cal_sky_dir){
    path = paste0(strsplit(dirname(file_image), '/cal',fixed=T),'/cal_sky/')
    if(path != cal_sky_dir){
      path = cal_sky_dir ## For scripting mode, normally JP will have pre-setup directories
    }
    base = strsplit(basename(file_image),'.fits',fixed=T)[[1]]
    return(paste0(path,base,'_sky_rem.fits'))
  }
  
  sky_filelist = list.files(sky_info[1,pathsky])
  sky_info = load_files(input_args, which_module = "apply_super", sky_info = sky_info)$sky_info
  sky_filelist = sky_info$filesky

  sky_info = data.frame(sky_info)
  
  cat('Processing',length(sky_filelist),'files...\n')
  message("Skip completed files is ", ifelse(skip_completed_files, "ON", "OFF"))
  
  lo_loop = 1
  hi_loop = length(sky_filelist)
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = parallel_type)
    registerDoParallel(cl)
  }  
  
  dummy = foreach(i = lo_loop:hi_loop, .inorder=FALSE, .packages = c("Rfits", "Rwcs", "ProFound"))%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    file_image = paste(sky_info[i,'pathim'], sky_info[i,'fileim'], sep='/')
    file_cal_sky = setcal_sky_remfile(file_image, cal_sky_dir)
    
    if( !(skip_completed_files & file.exists(file_cal_sky)) ){
      temp_cal = Rfits_read(file_image, pointer=FALSE)
      temp_cal$DQ$keyvalues$BZERO = 0
      #safe name (since sometimes the sky file name gets truncated due to 80 character FITS limits)
      file_sky = sky_info[i,'filesky']
      file_sky = grep(sky_info[i,'filesky'], sky_filelist, value=TRUE)
      if(length(file_sky) == 0){
        stop('Missing file_sky ',i)
      }
      
      file_sky = paste(sky_info[i,'pathsky'], file_sky, sep='/')
      temp_sky = Rfits_read(file_sky, pointer=FALSE)
      
      file_super_sky = paste0(sky_super_dir,'/super_',temp_cal[[1]]$keyvalues$DETECTOR,'_',temp_cal[[1]]$keyvalues$FILTER,'.fits')
      super_sky = Rfits_read_image(file_super_sky, header=FALSE)
      
      obj_idx = temp_sky$PIXELMASK$imDat==0 & temp_sky$OBJECTMASK$imDat==0
      if(sum(obj_idx) == 0){ ## otherwise the linear model is going to fail
        obj_idx = temp_sky$PIXELMASK$imDat==0 & temp_cal$SCI$imDat <= quantile(temp_cal$SCI$imDat, probs =  0.5, na.rm = TRUE) ## minimally select stuff less than the median
      }
      
      ## linear model to each pixel 
      temp_lm = lm(temp_cal$SCI$imDat[obj_idx] ~ super_sky[obj_idx])$coefficients
      
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
      temp_cal_sky[!(obj_idx)] = NA
      
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
      
      if(file.exists(file_cal_sky)){
        message('Removing old cal_sky file: ',file_cal_sky)
        file.remove(file_cal_sky)
      }
      
      file.copy(file_image, file_cal_sky)
      Rfits_write_pix(temp_cal$SCI$imDat, file_cal_sky, ext=2)
      Rfits_write_key(file_cal_sky, keyname='SKY_M', keyvalue=temp_lm[2], keycomment='SKY M coef in SKY = M.SuperSky + B + P', ext=2)
      Rfits_write_key(file_cal_sky, keyname='SKY_B', keyvalue=temp_lm[1], keycomment='SKY B coef in SKY = M.SuperSky + B + P', ext=2)
      Rfits_write_key(file_cal_sky, keyname='SKY_P', keyvalue=final_pedestal, keycomment='SKY P coef in SKY = M.SuperSky + B + P', ext=2)
      
      check_Ndhu = Rfits_nhdu(file_cal_sky)
      extloc = Rfits_extname_to_ext(file_cal_sky, 'SKY_Super')
      
      if(is.na(extloc)){
        Rfits_write_image(sky_map, file_cal_sky, create_file=FALSE, create_ext=TRUE)
        Rfits_write_key(filename=file_cal_sky, ext=check_Ndhu+1, keyname='EXTNAME', keyvalue='SKY_Super', keycomment='extension name')
      }else{
        Rfits_write_pix(sky_map, file_cal_sky, ext=extloc)
      }
      
      return(NULL)
    }else{
      return(NULL)
    }
    
  }
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
}
do_modify_pedestal = function(input_args){
  
  cat("\n")
  message("## Modifying pedestal ##")
  cat("\n")

  cal_sky_dir = input_args$cal_sky_dir
  cal_sky_renorm_dir = input_args$cal_sky_renorm_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  parallel_type = input_args$additional_params$parallel_type
  skip_completed_files = additional_params$skip_completed_files
  
  filelist = load_files(input_args, which_module = "modify_pedestal")
  
  cat('Processing',length(filelist),'files...\n')
  message("Skip completed files is ", ifelse(skip_completed_files, "ON", "OFF"))
  
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
  
  cal_sky_info = as.data.table(cbind(cal_sky_info_ext1, cal_sky_info_ext2))
  
  ped_info = cal_sky_info[,list(ped_mean = mean(SKY_B + SKY_P), ped_med = median(SKY_B + SKY_P), ped_max = max(SKY_B + SKY_P), ped_min = min(SKY_B + SKY_P)), by=list(OBS_ID, EXPOSURE)]
  ped_info[,ped_diff := (ped_max - ped_min) / ped_mean]
  
  cal_sky_info = data.frame(cal_sky_info)
  
  lo_loop = 1
  hi_loop = dim(cal_sky_info)[1]
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = parallel_type)
    registerDoParallel(cl)
  } 
  
  if(do_NIRISS){
    dummy = foreach(i = lo_loop:hi_loop, .packages = c("Rfits", "Rwcs", "ProFound"))%dopar%{
      if(i %% 100 == 0){
        message('File ',i,' of ', hi_loop)
      }
      file_cal_sky_renorm = paste0(cal_sky_renorm_dir,'/',cal_sky_info[i,'file'])
      
      if( !(skip_completed_files & file.exists(file_cal_sky_renorm)) ){
        if(file.exists(file_cal_sky_renorm)){
          message('Removing old cal_sky file: ',file_cal_sky_renorm)
          file.remove(file_cal_sky_renorm)
        }
        file.copy(cal_sky_info[i,'full'], cal_sky_renorm_dir) #no pedestal adjustment for NIRISS
        return(NULL)
      }else{
        return(NULL)
      }
    }
  }else if(do_MIRI){
    dummy = foreach(i = lo_loop:hi_loop, .packages = c("Rfits", "Rwcs", "ProFound"))%dopar%{
      if(i %% 100 == 0){
        message('File ',i,' of ', hi_loop)
      }
      file_cal_sky_renorm = paste0(cal_sky_renorm_dir,'/',cal_sky_info[i,'file'])

      if( !(skip_completed_files & file.exists(file_cal_sky_renorm)) ){
        if(file.exists(file_cal_sky_renorm)){
          message('Removing old cal_sky file: ',file_cal_sky_renorm)
          file.remove(file_cal_sky_renorm)
        }
        file.copy(cal_sky_info[i,'full'], cal_sky_renorm_dir) 
        return(NULL)
      }else{
        return(NULL)
      }
    }
  }else{
    dummy = foreach(i = lo_loop:hi_loop, .packages = c("Rfits", "Rwcs", "ProFound"))%dopar%{
      if(i %% 100 == 0){
        message('File ',i,' of ', hi_loop)
      }
      file_cal_sky_renorm = paste0(cal_sky_renorm_dir,'/',cal_sky_info[i,'file'])
      
      if( !(skip_completed_files & file.exists(file_cal_sky_renorm)) ){
        if(file.exists(file_cal_sky_renorm)){
          message('Removing old cal_sky file: ',file_cal_sky_renorm)
          file.remove(file_cal_sky_renorm)
        }
        file.copy(cal_sky_info[i,'full'], cal_sky_renorm_dir)
        
        if(cal_sky_info[i,'CHANNEL'] == 'SHORT'){
          temp_cal_sky = Rfits_read(cal_sky_info[i,'full'])
          ext_names = names(temp_cal_sky)
          
          temp_mask = Rfits_read_image(cal_sky_info[i,'full'], ext = "DQ", header = F)
          JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
          
          pro_current_sky = profoundProFound(temp_cal_sky$SCI[,]$imDat, mask=JWST_cal_mask, sky=0, box=683, redosky = FALSE)
          
          if(is.null(pro_current_sky$skyChiSq)){
            pro_current_sky$skyChiSq = 1e6
          }
          if(is.na(pro_current_sky$skyChiSq)){
            pro_current_sky$skyChiSq = 1e6
          }
          if(is.infinite(pro_current_sky$skyChiSq)){
            pro_current_sky$skyChiSq = 1e6
          }
          
          if(pro_current_sky$skyChiSq < 0.9 | pro_current_sky$skyChiSq > 1.1){
            message('Bad Sky: ',cal_sky_info[i,'file'],'. Current: ',round(pro_current_sky$skyChiSq,2),' Trying to make better sky with ProFound')
            # if(sum(pro_current_sky$objects) / prod(dim(temp_cal_sky$SCI)) < 0.2){ ## What was this for?
            
            pro_new_sky = profoundProFound(temp_cal_sky$SCI[,]$imDat, mask=JWST_cal_mask, box=683, roughpedestal = TRUE)
            
            if(is.null(pro_new_sky$skyChiSq)){
              pro_new_sky$skyChiSq = 1e6
            }
            if(is.na(pro_new_sky$skyChiSq)){
              pro_new_sky$skyChiSq = 1e6
            }
            if(is.infinite(pro_new_sky$skyChiSq)){
              pro_new_sky$skyChiSq = 1e6
            }
            
            if(pro_new_sky$skyChiSq > 0.9 & pro_new_sky$skyChiSq < 1.1){
              message('Using better ProFound Sky: ',cal_sky_info[i,'file'],' Old: ',round(pro_current_sky$skyChiSq,2),' New: ',round(pro_new_sky$skyChiSq,2))
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
          }
        }
        return(NULL)
      }else{
        return(NULL)
      }
    }
  }
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
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
  cal_sky_info = as.data.table(cal_sky_info)
  
  pmap = cal_sky_info[,list(PHOTMJSR=mean(PHOTMJSR)), keyby=list(CRDS_CTX,DETECTOR,FILTER)]
  
  cal_sky_info$MAGZERO_FIX = 0.0
  
  for(i in 1:dim(cal_sky_info)[1]){
    #below is then Eqn [1] below
    ## Tie to pmap 1179 for now
    temp_MAGZERO_FIX = cal_sky_info[i,MAGZERO] -2.5*log10(as.numeric(pmap[CRDS_CTX=='jwst_1179.pmap' & DETECTOR==cal_sky_info[i,'DETECTOR'] & FILTER==cal_sky_info[i,'FILTER'],PHOTMJSR]) / cal_sky_info[i,PHOTMJSR])
    if(length(temp_MAGZERO_FIX) == 1){
      cal_sky_info[i,MAGZERO_FIX := temp_MAGZERO_FIX] 
    }else{
      cal_sky_info[i,MAGZERO_FIX := cal_sky_info[i,MAGZERO]]
    }
  }
  write.csv(cal_sky_info, paste0(cal_sky_info_save_dir, '/cal_sky_info.csv'), row.names = FALSE)
}
do_gen_stack = function(input_args){
  
  cat("\n")
  message("## Making ProPane stacks ##")
  cat("\n")

  VID = input_args$VID
  FILT = input_args$FILT
  ref_dir = input_args$ref_dir
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  magzero = input_args$magzero
  cores = input_args$cores_stack
  tasks = input_args$tasks_stack
  
  skip_completed_files = additional_params$skip_completed_files
  parallel_type = input_args$additional_params$parallel_type
  
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
    input_args$additional_params$NAXIS_short
  ) 
  NAXIS_long = ifelse(
    is.null(input_args$additional_params$NAXIS_long), 
    3000, 
    input_args$additional_params$NAXIS_long
  ) 
  
  ## Better for scripting mode
  invar_dir = input_args$invar_dir
  median_dir = input_args$median_dir
  sky_frames_dir = input_args$sky_frames_dir
  dump_dir_stub = input_args$dump_dir
  orig_cal_sky_info = fread(paste0(input_args$cal_sky_info_save_dir, '/cal_sky_info.csv'))
  
  unique_visits = unique(orig_cal_sky_info$VISIT_ID)
  
  message("Now running stack!")
  message("Skip completed files is ", ifelse(skip_completed_files, "ON", "OFF"))
  
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
      cal_sky_info = orig_cal_sky_info[grepl(VID, orig_cal_sky_info$VISIT_ID) & !grepl("NIS|MIRIMAGE", orig_cal_sky_info$DETECTOR),]
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
    
    if(!is.null(input_args$additional_params$module_list)){
      module_list_temp = input_args$additional_params$module_list
      module_list = module_list[module_list %in% module_list_temp]
    }

    lo_loop = 1
    hi_loop = dim(stack_grid)[1]
    
    stack_grid = stack_grid[sample(hi_loop),] #randomise to spread out the CPU load (else all of the SMACS stacking end up in a single task and takes much longer)
    
    stack_grid = data.frame(stack_grid)
    cal_sky_info = data.frame(cal_sky_info)
    
    list_residual_temp_files = list.files(invar_dir, pattern = "temp_list_.+", full.names = T)
    if(length(list_residual_temp_files)>0){
      unlink(list_residual_temp_files, recursive = T)
    }
    
    if(is.null(parallel_type)){
      registerDoParallel(cores = tasks)
    }else{
      cl <- makeCluster(spec = tasks, type = parallel_type)
      registerDoParallel(cl)
    }
    
   dummy = foreach(i = lo_loop:hi_loop, .packages = c('Rwcs', 'Rfits', 'ProPane', 'ProFound'))%dopar%{
      message('Stacking ', i,' of ',hi_loop)
      for(j in module_list){
        if( !(skip_completed_files & file.exists(paste0(invar_dir, '/stack_',stack_grid[i,'VISIT_ID'],'_',stack_grid[i,'FILTER'],'_',j,'.fits'))) ){
          message('  Processing ',stack_grid[i,'VISIT_ID'],' ',stack_grid[i,'FILTER'], ' ', j)
          
          if(j == 'MIRIMAGE'){
            input_info = cal_sky_info[cal_sky_info[['VISIT_ID']]==stack_grid[i,'VISIT_ID'] & cal_sky_info[['FILTER']]==stack_grid[i,'FILTER'] & grepl('MIRIMAGE', cal_sky_info[['DETECTOR']]),]
            CRVAL1 = module_MIRIMAGE[module_MIRIMAGE[['VISIT_ID']] == stack_grid[i,'VISIT_ID'], 'CRVAL1']
            CRVAL2 = module_MIRIMAGE[module_MIRIMAGE[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CRVAL2']
            CD1_1 = module_MIRIMAGE[module_MIRIMAGE[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_1']
            CD1_2 = module_MIRIMAGE[module_MIRIMAGE[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_2']
            NAXIS = NAXIS_long
            CRPIX = NAXIS_long/2.0
          }
          
          if(j == 'NIS'){
            input_info = cal_sky_info[cal_sky_info[['VISIT_ID']]==stack_grid[i,'VISIT_ID'] & cal_sky_info[['FILTER']]==stack_grid[i,'FILTER'] & grepl('NIS', cal_sky_info[['DETECTOR']]),]
            CRVAL1 = module_NIS[module_NIS[['VISIT_ID']] == stack_grid[i,'VISIT_ID'], 'CRVAL1']
            CRVAL2 = module_NIS[module_NIS[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CRVAL2']
            CD1_1 = module_NIS[module_NIS[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_1']
            CD1_2 = module_NIS[module_NIS[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_2']
            NAXIS = NAXIS_long
            CRPIX = NAXIS_long/2.0
          }
          
          if(j == 'NRCA_short'){
            input_info = cal_sky_info[cal_sky_info[['VISIT_ID']]==stack_grid[i,'VISIT_ID'] & cal_sky_info[['FILTER']]==stack_grid[i,'FILTER'] & grepl('NRCA', cal_sky_info[['DETECTOR']]),]
            CRVAL1 = module_A_WCS_short[module_A_WCS_short[['VISIT_ID']] == stack_grid[i,'VISIT_ID'], 'CRVAL1']
            CRVAL2 = module_A_WCS_short[module_A_WCS_short[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CRVAL2']
            CD1_1 = module_A_WCS_short[module_A_WCS_short[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_1']
            CD1_2 = module_A_WCS_short[module_A_WCS_short[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_2']
            NAXIS = NAXIS_short
            CRPIX = NAXIS_short/2.0
          }
          
          if(j == 'NRCB_short'){
            input_info = cal_sky_info[cal_sky_info[['VISIT_ID']]==stack_grid[i,'VISIT_ID'] & cal_sky_info[['FILTER']]==stack_grid[i,'FILTER'] & grepl('NRCB', cal_sky_info[['DETECTOR']]),]
            CRVAL1 = module_B_WCS_short[module_B_WCS_short[['VISIT_ID']] == stack_grid[i,'VISIT_ID'], 'CRVAL1']
            CRVAL2 = module_B_WCS_short[module_B_WCS_short[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CRVAL2']
            CD1_1 = module_B_WCS_short[module_B_WCS_short[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_1']
            CD1_2 = module_B_WCS_short[module_B_WCS_short[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_2']
            NAXIS = NAXIS_short
            CRPIX = NAXIS_short/2.0
          }
          
          if(j == 'NRCA_long'){
            input_info = cal_sky_info[cal_sky_info[['VISIT_ID']]==stack_grid[i,'VISIT_ID'] & cal_sky_info[['FILTER']]==stack_grid[i,'FILTER'] & grepl('NRCA', cal_sky_info[['DETECTOR']]),]
            CRVAL1 = module_A_WCS_long[module_A_WCS_long[['VISIT_ID']] == stack_grid[i,'VISIT_ID'], 'CRVAL1']
            CRVAL2 = module_A_WCS_long[module_A_WCS_long[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CRVAL2']
            CD1_1 = module_A_WCS_long[module_A_WCS_long[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_1']
            CD1_2 = module_A_WCS_long[module_A_WCS_long[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_2']
            NAXIS = NAXIS_long
            CRPIX = NAXIS_long/2.0
          }
          
          if(j == 'NRCB_long'){
            input_info = cal_sky_info[cal_sky_info[['VISIT_ID']]==stack_grid[i,'VISIT_ID'] & cal_sky_info[['FILTER']]==stack_grid[i,'FILTER'] & grepl('NRCB', cal_sky_info[['DETECTOR']]),]
            CRVAL1 = module_B_WCS_long[module_B_WCS_long[['VISIT_ID']] == stack_grid[i,'VISIT_ID'], 'CRVAL1']
            CRVAL2 = module_B_WCS_long[module_B_WCS_long[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CRVAL2']
            CD1_1 = module_B_WCS_long[module_B_WCS_long[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_1']
            CD1_2 = module_B_WCS_long[module_B_WCS_long[['VISIT_ID']]  == stack_grid[i,'VISIT_ID'], 'CD1_2']
            NAXIS = NAXIS_long
            CRPIX = NAXIS_long/2.0
          }
          
          image_list = {}
          inVar_list = {}
          
          for(k in 1:dim(input_info)[1]){
            temp_image = Rfits_read_image(input_info[k,'full'], ext=2)
            
            temp_mask = Rfits_read_image(input_info[k,'full'], ext=4, header=FALSE)
            
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
            
            sky_file = paste0(sky_frames_dir, sub('_rem', '', input_info[k,'stub']),'_',input_info[k,'FILTER'],'.fits')
            inVar_list = c(inVar_list, Rfits::Rfits_read_key(sky_file, 'SKYRMS')^-2)
          }
          
          temp_file = paste0(invar_dir, '/temp_list_',i,'.fits')
          Rfits_write(image_list, temp_file)
          rm(image_list)
          image_list_point = Rfits_read(temp_file, pointer=TRUE)
          
          dump_stub = paste0(dump_dir_stub, stack_grid[i,'VISIT_ID'], "/", stack_grid[i,'FILTER'], "/", j, "/")
          if(dir.exists(dump_stub)){
            unlink(dump_stub, recursive = T)
          }
          dir.create(dump_stub, recursive = T)
          
          temp_proj = propaneGenWCS(
            image_list = image_list_point,
            CRVAL1 = as.numeric(CRVAL1),
            CRVAL2 = as.numeric(CRVAL2),
            CRPIX1 = as.numeric(CRPIX),
            CRPIX2 = as.numeric(CRPIX),
            CD1_1 = as.numeric(CD1_1),
            CD1_2 = as.numeric(CD1_2),
            CD2_1 = as.numeric(CD1_2),
            CD2_2 = -as.numeric(CD1_1),
            NAXIS1 = as.numeric(NAXIS),
            NAXIS2 = as.numeric(NAXIS)
          )
          
          temp_proj$DATE_BEG = as.character(min(as.POSIXlt.Date(input_info[["DATE.OBS"]]), na.rm=TRUE))
          temp_proj$DATE_END = as.character(max(as.POSIXlt.Date(input_info[["DATE.OBS"]]), na.rm=TRUE))
          
          output_stack = propaneStackWarpInVar(
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
            keep_extreme_pix = FALSE,
            dump_frames = TRUE,
            dump_dir = dump_stub,
            multitype="cluster" ## Just keep this default since this parallelisation occurs in the ProPane package (and the underlying C libraries)
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
          
          filestub = paste0(invar_dir, '/stack_',stack_grid[i,'VISIT_ID'],'_',stack_grid[i,'FILTER'],'_',j)
          
          Rfits_write(data = output_stack,
                      filename = paste0(filestub,'.fits')
          )
          
          filestub_med = paste0(median_dir, '/med_',stack_grid[i,'VISIT_ID'],'_',stack_grid[i,'FILTER'],'_',j)
          
          Rfits_write(data = median_stack,
                      filename = paste0(filestub_med,'.fits')
          )
          gc()
        }
        return(NULL)
      }

    }
   if(!is.null(parallel_type)){
     stopCluster(cl)
     registerDoSEQ()
   }else{
     stopImplicitCluster()
     registerDoSEQ()
   }
   gc()
  }
}
do_wisp_rem = function(input_args){
  
  cat("\n")
  message("## Removing wisps ##")
  cat("\n")

  filelist = input_args$filelist
  VID = input_args$VID
  median_dir = input_args$median_dir
  cores = input_args$cores_pro
  additional_params  = input_args$additional_params
  SIGMA_LO = input_args$SIGMA_LO
  
  skip_completed_files = additional_params$skip_completed_files
  parallel_type = additional_params$parallel_type
  
  if(is.null(SIGMA_LO)){
    message("No smoothing on wisp removal")
  }else{
    message(paste0("Using smoothing kernel SIGMA_LO = ", SIGMA_LO))
  }
  
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
  mod_visit_grid = info_wisp[,c("stub", "MODULE", "VISIT_ID", "CRVAL1", "CRVAL2", "DETECTOR")]
  
  cat('Processing',dim(info_wisp)[1],'files\n')
  message("Skip completed files is ", ifelse(skip_completed_files, "ON", "OFF"))
  
  idx_wisp_skip = !(unname(sapply( info_wisp$full, function(x){any(na.omit(Rfits_extnames(x)) == "SCI_ORIG")} )) & skip_completed_files)
  
  ref_im_list = {}
  for(ii in 1:dim(mod_visit_grid)[1]){
    if(idx_wisp_skip[ii]){
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
        str_match(ref_files, "F\\s*(.*?)\\s*[WMN]")[,2]
      )
      
      filter_long[is.na(filter_long)] = -1
      long_num = max(filter_long, na.rm = T)
      
      ref_file_long = ref_files[ which(grepl(long_num, ref_files))[1] ] #Get the longest filter
      
      ref_im_list = c(ref_im_list, list(Rfits_point(ref_file_long))) ## Allow loop in wisp rem to read long wavelength reference
      message(paste("Loading reference for:", paste0(mod_visit_grid$VISIT_ID[ii]), mod_visit_grid$DETECTOR[ii], mod_visit_grid$stub[ii] ))
    }
  }
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = 'PSOCK')
    registerDoParallel(cl)
  }
  
  temp = foreach(ii = 1:dim(info_wisp)[1], .errorhandling = "pass", .packages = c("Rwcs", "Rfits", "ProPane", "ProFound"))%dopar%{
    if(idx_wisp_skip[ii]){
      ## copy original data to directory where we keep the wisps
      vid = paste0(info_wisp$VISIT_ID[ii])
      modl = paste0("NRC", info_wisp$MODULE[ii])
      
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
    }else{
      return(NULL)
    }
  }
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
}
do_patch = function(input_args){
  
  cat("\n")
  message("## Patching ProPane stacks ##")
  cat("\n")

  VID = input_args$VID
  FILT = input_args$FILT
  invar_dir = input_args$invar_dir
  median_dir = input_args$median_dir
  patch_dir = input_args$patch_dir
  cores = input_args$cores_stack
  do_NIRISS = input_args$do_NIRISS
  do_MIRI = input_args$do_MIRI
  skip_completed_files = additional_params$skip_completed_files
  parallel_type = input_args$additional_params$parallel_type
  
  patch_stub = patch_dir
  
  if(do_NIRISS){
    pixscale_list = c("NIS")
  }else if(do_MIRI){
    pixscale_list = c("MIRIMAGE")
  }else{
    pixscale_list = c("short", "long")
  }
  
  for(j in pixscale_list){
    invar_files = list.files(path=invar_dir, pattern = glob2rx(paste0("*stack*", j, "*.fits$")), full.names = T)
    median_files = list.files(path=median_dir, pattern = glob2rx(paste0("*med*", j, "*.fits$")), full.names = T)
    invar_files = invar_files[grepl(VID, invar_files) & grepl(FILT, invar_files) & !grepl("patch", invar_files)]
    median_files = median_files[grepl(VID, median_files) & grepl(FILT, median_files) & !grepl("patch", median_files)]
    
    if(length(median_files)==0){
      message("No frames") 
    }else{
      if(is.null(parallel_type)){
        registerDoParallel(cores = cores)
      }else{
        cl <- makeCluster(spec = cores, type = parallel_type)
        registerDoParallel(cl)
      }
      
      dummy = foreach(i = 1:length(invar_files))%dopar%{
        
        message(paste0("Patching: ", invar_files[i]))
        info_load_invar = Rfits_read(invar_files[i], pointer=TRUE)
        info_temp = info_load_invar$image
        
        fstub = str_remove( str_remove(string = info_temp$filename, pattern = dirname(info_temp$filename)), ".+stack" )
        
        if( !(skip_completed_files & file.exists(paste0(patch_stub, "patch", fstub))) ){
          
          load_invar = Rfits_read(invar_files[i], pointer=FALSE)
          load_median = Rfits_read(median_files[i], pointer=FALSE)
          temp = load_invar$image
          
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
          
          patch_mask = temp - propane_mask_fix
          
          patch = list()
          
          patch$image = propane_mask_fix
          patch$patch = patch_mask
          patch$inVar = load_invar$inVar
          patch$weight = load_invar$weight
          
          patch$image$keyvalues$EXTNAME = "image"
          patch$patch$keyvalues$EXTNAME = "patch"
          patch$inVar$keyvalues$EXTNAME = "inVar"
          patch$weight$keyvalues$EXTNAME = "weight"
          
          class(patch) = "Rfits_list"
          
          Rfits_write_all(
            data = patch, 
            filename = paste0(patch_stub, "patch", fstub)
          )
          return(NULL)
        }else{
          return(NULL)
        }
      }
      
      if(!is.null(parallel_type)){
        stopCluster(cl)
        registerDoSEQ()
      }else{
        stopImplicitCluster()
        registerDoSEQ()
      }
      gc()
    }
  }
}
do_RGB = function(input_args){
  
  cat("\n")
  message("## Making RGB image ##")
  cat("\n")

  VID = input_args$VID
  ref_dir = input_args$ref_dir
  RGB_dir = input_args$RGB_dir
  patch_dir = input_args$patch_dir
  cal_sky_info_save_dir = input_args$cal_sky_info_save_dir
  
  blue_filters = "F070W|F090W|F115W|F150W|F140M|F162M|F164N"
  green_filters = "F200W|F277W|F182M|F210M|F187N|F212N"
  red_filters = "F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430M|F460M|F480M|F323N|F405N|F466N|F470N|F560W|F770W|F1000W|F1130W|F1280W|F1500W|F1800W|F2100W|F2550W"
  locut = 0.6
  hicut = 0.9999
  
  # cal_sky_info = fread(paste0(ref_dir, "/Pro1oF/cal_sky_info.csv"))
  cal_sky_info = fread(paste0(cal_sky_info_save_dir, "/cal_sky_info.csv"))
  vid_temp = VID
  unique_visits = grep(vid_temp, unique(cal_sky_info$VISIT_ID), value = T)
  
  for(VID in unique_visits){
    
    if(is.null(RGB_dir)){
      RGB_dir = paste0(ref_dir, "/RGB/", VID, "/")
      if(is.null(ref_dir)){
        stop('Please provide a parent directory to put the RGB folder')
      }else{
        dir.create(RGB_dir, showWarnings = F, recursive = T)
      }
    }
    
    file_list = list.files(patch_dir, 
                           pattern = glob2rx(paste0("*patch*", VID, "*short*fits")),
                           full.names = T)
    
    frame_info = Rfits_key_scan(filelist = file_list,
                                keylist = c("CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "NAXIS1", "NAXIS2"))
    # temp_proj = Rwcs_keypass(CRVAL1 = mean(frame_info$CRVAL1),
    #                          CRVAL2 = mean(frame_info$CRVAL2),
    #                          CD1_1 = min(frame_info$CD1_1),
    #                          CD1_2 = min(frame_info$CD1_2),
    #                          CD2_1 = min(frame_info$CD1_2),
    #                          CD2_2 = -min(frame_info$CD1_1)
    # )
    temp_proj = propaneGenWCS(
      filelist = file_list, 
      rotation = 'get'
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
    
    filestub = paste0(patch_dir, '/RGB_patch_image_',VID, "_all", "_clear.png")
    CairoPNG(filename = filestub, width=im_width, height=im_height, 
             units='in', res=400, quality=100)
    par(mar=c(0,0,0,0), oma=rep(0,4))
    Rwcs_imageRGB(
      R = R_patch,
      G = G_patch,
      B = B_patch,
      locut=locut,
      hicut=hicut,
      type='quan',
      stretch = 'log',
      sparse=1,
      decorate=F
    )
    dev.off()
  }

}
do_wisp_reverse = function(input_args){
  
  cat("\n")
  message("## Reversing wisp correction ##")
  cat("\n")
  
  filelist_all = input_args$filelist ## raw CAL files - straight from the calibration pipeline
  VID = input_args$VID
  cores = input_args$cores
  do_claws = input_args$additional_params$do_claws
  parallel_type = input_args$additional_params$parallel_type
  
  if(VID != ""){
    scan_wisp_rem = Rfits_key_scan(filelist = filelist_all, keylist = c("FILTER", "PROGRAM"))
    not_pid_idx = sapply(scan_wisp_rem$PROGRAM, function(x)grepl(x, VID, fixed = T)) ## make sure PID is not embedded in string of VID
    filelist = filelist_all[not_pid_idx]
  }else{
    filelist = filelist_all
  }
  ## stop edit 
  
  info = Rfits_key_scan(filelist = filelist,
                        keylist=c('DETECTOR', 'MODULE', 'FILTER', 'VID'), cores=1)
  
  if(do_claws){
    info_wisp = info[DETECTOR %in% c('NRCA1', 'NRCA2', 'NRCA3','NRCA4', 'NRCB1','NRCB2', 'NRCB3','NRCB4', 'MIRIMAGE'),] 
  }else{
    info_wisp = info[DETECTOR %in% c('NRCA3','NRCA4', 'NRCB3','NRCB4', 'MIRIMAGE'),] 
  }
  
  ## WISP remove reverseal for NIRCam (step 11)
  module_list = paste0("NRC", unique(info$MODULE))
  unique_visits = unique(info_wisp$VID)
  
  message("Reversing wisp correction on short wavelength detectors:")
  cat('Processing',dim(info_wisp)[1],'files\n')
  
  message("Showing <10 files")
  print(
    head(info_wisp[, c("stub", "DETECTOR", "MODULE")], 10)
  )
  message("...")
  
  if(is.null(parallel_type)){
    registerDoParallel(cores = cores)
  }else{
    cl <- makeCluster(spec = cores, type = parallel_type)
    registerDoParallel(cl)
  }
  
  temp = foreach(ii = 1:dim(info_wisp)[1], .errorhandling = "stop")%dopar%{
    ## copy original data to directory where we keep the wisps
    vid = paste0(info_wisp$VID[ii])
    modl = paste0("NRC", info_wisp$MODULE[ii])
    wisp_frame = Rfits_read(info_wisp$full[ii])
    if(!is.null(wisp_frame$SCI_ORIG)){
      Rfits_write_pix(data = wisp_frame$SCI_ORIG[,]$imDat, filename = info_wisp$full[ii], ext = 2) 
    }
    return(NULL)
  }
  
  if(!is.null(parallel_type)){
    stopCluster(cl)
    registerDoSEQ()
  }else{
    stopImplicitCluster()
    registerDoSEQ()
  }
  gc()
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

do_help = function(input_args){
  
  cat("\n")
  message("## Help ##")
  cat("\n")
  
  message("Author: Jordan C. J. D'Silva")
  message("Date: ", Sys.Date())
  
  message("\n")
  
  message(
    "Running zork: Rscript zork_process.R <VID, e.g., '2738|1176261001'> <FILT, e.g., 'F090W|F444W'> "
  )
  
  message("\n")
  
  message("It is advised to separate the files into organised directories to avoid unwanted clashes occurring during the sub-routines.")
  message("Recommended directory structure is, for example: ")
  message(
    'JP/
  ├─ Pro1oF/
  │  ├─ cal/
  │  ├─ cal_sky/
  │  ├─ cal_sky_renorm/
  │  ├─ cal_sky_info.csv
  ├─ sky_pro/
  │  ├─ sky_frames/
  │  ├─ sky_super/
  │  ├─ sky_info.csv
  ├─ InVar_Stacks/
  ├─ Median_Stacks/
  ├─ Patch_Stacks/
  ├─ dump/
  ├─ RGB/'
  )
  
  message("\n")
  
  message("Input arguments for each of the sub-routines in JP are as follows:")
  message(
    'input_args = list(
      filelist = <path to the cal files produced from stage2 of the calibration pipeline>, ## <-- vector of strings
      additional_params = additional_params, ## <-- extra 1/f settings, wisp removal settings, stacking settings as set in initialise_variables.R, list
      
      ref_dir = ref_dir, ## <-- Reference directory, where you want to put your JUMPROPE files, string
            
      Pro1oF_dir = Pro1oF_dir, ## <-- 1/f corrected file directory, string
      cal_sky_dir = cal_sky_dir, ## <-- Sky removed files directory, string
      cal_sky_renorm_dir = cal_sky_renorm_dir, ## <-- Sky removed and pedestal modified files directory, the input for making mosaics, string
      cal_sky_info_save_dir = cal_sky_info_save_dir, ## <-- Usually the parent directory of above cal stuff, directory to store a CSV file with info about the frames, string
      
      sky_frames_dir = sky_frames_dir, ## <-- sky frames from ProFound directory, string
      sky_super_dir = sky_super_dir, ## <-- super super, combined stack of sky files, directory, string
      sky_pro_dir = sky_pro_dir, ## <-- Usually the parent directory of above sky stuff to keep the sky things together, string

      invar_dir = invar_dir, ## <-- ProPane inverse variance stack directory, string
      median_dir = median_dir, ## <-- Median stack directory, string
      patch_dir = patch_dir, ## <-- Patched stack directory, string
      dump_dir = dump_dir, ## <-- Dump frames directory, save the warp fields from mosaicking, string
      
      RGB_dir = RGB_dir, ## <-- RGB image directory, string
      
      magzero = 23.9, ## <-- Magnitude zero-point for ProPane mosaics, 8.9/16.4/23.9 is Jy,milliJy/microJy, numeric
      
      VID = VID, ## <-- 4 digit program ID or 10 digit visit ID, string, e.g., "1234567890"
      FILT = FILT, ## <-- Set of filters to process, string separated by "|", e.g., "F090W|F150W" to do those 2 filters
      
      do_NIRISS = do_NIRISS, ## <-- Process NIRISS, boolean
      do_MIRI = do_MIRI, ## <-- Process MIRI, boolean
      
      cores_pro = cores_pro, ## <-- Cores for processing, numeric
      cores_stack = cores_stack, ## <-- Cores for stacking with ProPane, numeric
      tasks_stack = tasks_stack, ## <-- Cores for parallel stacking jobs, numeric
  
      SIGMA_LO = NULL #keep blurring at wisp rem stage off by default, otherwise numeric, standard deviation of Gaussian kernel to blur derived wisp template
    )'
  )
  
  message("\n")
  
  message("Contents of initialise_variables.R for some extra parameters to use in JUMPROPE:")
  message(
    "additional_params = list(
    ## Keep trend data determines the aggressiveness of the 1/f. If nothing in this data frame, then default 1/f proceeds. 
    ## Default should be good for most blank fields.
    
    ## vlarge for big objects that fill frame e.g., VV191. Least aggressive 1/f.
    ID_vlarge = data.frame(
      VISIT_ID = c(
        2561001003 #<-- Put 10 digit VISITIDs here, e.g., 1176341001, 1176361001
      ),
      MODULE = c(
       'A'   #<-- Put 'A'/'B' here 
      )
    ),
    
    ## large for crowded fields or wispy fields e.g., SMACS Cluster Module B. Less aggressive 1/f.
    ID_large = data.frame(
      VISIT_ID = c(
        2561001003
      ),
      MODULE = c(
       'B' 
      )
    ),
    
    ## overwrite and use one 1/f setting for everything. 
    ## so you don't have to laboriously type out every single VISITID and MODULE combination :D
    ow_vlarge = FALSE,
    ow_large = FALSE,
    
    ## Claws removal mode
    ## I.e., perform 'wisp removal' algorithm on all NIRCam short wavelength detectors
    ## and not only wisp affected [A3,A4,B3,B4]
    do_claws = FALSE, 
    
    ## Path to reference astrometric catalogue - to make the ProPane stacks
    tweak_catalogue = NULL, ## NULL will have no internal tweak
    
    NAXIS_long = NULL, ## Size of the long pixel scales mosaic, keep NULL for default (3000 pixels) 
    NAXIS_short = NULL, ## Size of the short pixel scales mosaic, keep NULL for default (6000 pixels) 
    module_list = NULL, ## What modules should we stack, options are ('NRCA_short', 'NRCA_long', 'NRCB_short', 'NRCB_long', 'NIS', 'MIRIMAGE')
    
    parallel_type = NULL ## type options for makeCluster in parallel for stacking and wisp removal, type='PSOCK' might be more stable on Linux systems "
  )
}

## Helper functions to set up JP directories
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
  }else if(make_dirs == "F"){
    message("Continuing...")
    message("Where is the reference directory (supply directory or nothing): ")
    ref_dir = readLines("stdin", n=1)
    if(ref_dir==""){
      ref_dir = getwd()
      message("No user input. Assuming everything in current working directory.")
    }
    return(ref_dir)
  }else{
    message("Must be 'T/F'. Please try again:")
    make_directory_structure()
  }
}
load_raw_files = function(dir_raw){
  
  cal_files = c(
    list.files(dir_raw, pattern = ".fits$", full.names = T, recursive = T),
    NULL
  )
  
  return(cal_files)
}
