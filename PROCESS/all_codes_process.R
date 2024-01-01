library(Rfits)
library(Rwcs)
library(ProFound)
library(data.table)
library(foreach)
library(doParallel)
library(Cairo)
library(ProPane)

pipe_version = 2.0

do_1of = function(input_args){
  
  keep_trend_data = input_args$keep_trend_data
  filelist = input_args$filelist
  Pro1oF_dir = input_args$Pro1oF_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  
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
  
  ID_vlarge = keep_trend_data$ID_vlarge
  ID_large = keep_trend_data$ID_large
  
  ### BITS TO EDIT END ###
  
  registerDoParallel(cores=cores)
  
  scan_filt = Rfits_key_scan(filelist = filelist, keylist = c("FILTER"))
  
  filelist = filelist[!grepl('miri',filelist)] #we don't care about MIRI for now
  filelist = filelist[grepl('.fits$',filelist) & grepl(VID, filelist) & grepl(FILT, scan_filt$FILTER)]
  
  cat(filelist, sep='\n')
  cat('Processing',length(filelist),'files\n')
  
  lo_loop = 1
  hi_loop = length(filelist)
  
  dummy = foreach(i = lo_loop:hi_loop, .inorder=FALSE)%dopar%{
    if(i %% 100 == 0){
      message('File ',i,' of ', hi_loop)
    }
    file_image = filelist[i]
    basename = strsplit(basename(file_image),'.fits$')[[1]]
    fullbase = paste(Pro1oF_dir,basename,sep='/')
    
    temp_image = Rfits_read(filelist[i], pointer=FALSE)
    
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
    
    temp_mask = temp_image$DQ$imDat
    JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
    # if(sum(temp_mask!=0)/(prod(dim(temp_mask))) > 0.8){
    #   JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
    # }else{
    #   JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1 | temp_mask == 4, size=3)
    # }
    
    
    temp_zap = profoundSkyScan(image = temp_image$SCI$imDat,
                               mask = (temp_image$SCI$imDat==0) | JWST_cal_mask,
                               # mask = foo,
                               clip = c(0.0,0.9),
                               scan_block = c(512, 2048),
                               trend_block = trend_block,
                               keep_trend = keep_trend)
    
    CairoJPEG(filename = paste0(fullbase,'_orig.jpeg'), width=1000,height=1000)
    par(mar=c(0.1,0.1,0.1,0.1))
    magimage(temp_image$SCI$imDat, axes=FALSE, rem_med = T)
    legend('topleft', legend=paste(basename, 'orig'))
    dev.off()
    
    CairoJPEG(filename = paste0(fullbase,'_Pro1oF.jpeg'), width=1000,height=1000)
    par(mar=c(0.1,0.1,0.1,0.1))
    magimage(temp_zap$image_fix, axes=FALSE, rem_med = T)
    legend('topleft', legend=paste(basename, 'Pro1oF'))
    dev.off()
    
    file.copy(filelist[i], paste0(fullbase,'_Pro1oF.fits'), overwrite=TRUE)
    Rfits_write_pix(temp_zap$image_fix, paste0(fullbase,'_Pro1oF.fits'), ext=2)
    
    # check_Nhdu = Rfits_nhdu(paste0(fullbase,'_Pro1oF.fits'))
    # if(check_Nhdu == 9){
    #   Rfits_write_image(temp_zap$row_map + temp_zap$col_map, filename=paste0(fullbase,'_Pro1oF.fits'), create_file=FALSE)
    #   Rfits_write_key(filename=paste0(fullbase,'_Pro1oF.fits'), ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SKY_Pro1oF', keycomment='extension name')
    # }else{
    #   Rfits_write_pix(temp_zap$row_map + temp_zap$col_map, paste0(fullbase,'_Pro1oF.fits'), ext=check_Nhdu)
    #   Rfits_write_key(filename=paste0(fullbase,'_Pro1oF.fits'), ext=check_Nhdu, keyname='EXTNAME', keyvalue='SKY_Pro1oF', keycomment='extension name')
    # }
    
    check_Nhdu = Rfits_nhdu(paste0(fullbase,'_Pro1oF.fits'))
    extloc = Rfits_extname_to_ext(paste0(fullbase,'_Pro1oF.fits'), 'SKY_Pro1oF')
    
    if(is.na(extloc)){
      Rfits_write_image(temp_zap$row_map + temp_zap$col_map, filename=paste0(fullbase,'_Pro1oF.fits'), create_file=FALSE)
      Rfits_write_key(filename=paste0(fullbase,'_Pro1oF.fits'), ext=check_Nhdu+1, keyname='EXTNAME', keyvalue='SKY_Pro1oF', keycomment='extension name')
    }else{
      Rfits_write_pix(temp_zap$row_map + temp_zap$col_map, paste0(fullbase,'_Pro1oF.fits'), ext=extloc)
    }
    
    return(NULL)
  }
}
do_cal_process = function(input_args){
  
  Pro1oF_dir = input_args$Pro1oF_dir
  sky_frames_dir = input_args$sky_frames_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  
  filelist = c(
    list.files(Pro1oF_dir, full.names = TRUE) #If running Pro1/F outputs
  )
  
  if(do_NIRISS){
    filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist) & grepl('_nis_',filelist)]
  }else{
    filelist = filelist[!grepl('_miri_',filelist)] #we don't care about MIRI for now
    filelist = filelist[!grepl('_nis_',filelist)] #we don't care about MIRI for now
    filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist)]
  }
  cat(filelist, sep ="\n")
  
  #get main info:
  scan_filt = Rfits_key_scan(filelist = filelist,
                             keylist = c('FILTER'),
                             extlist = 1,
                             cores = 6
  )
  
  filelist = filelist[grepl(FILT, scan_filt$FILTER)] #just get the FITS files to be safe
  # cat('Processing',length(filelist),'files\n')
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
    
    # if(sum(temp_mask!=0)/(prod(dim(temp_mask))) > 0.8){
    #   JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
    # }else{
    #   JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1 | temp_mask == 4, size=3)
    # }
    
    #basic info
    suppressMessages({ #do not really want to see lots of SIP warnings!
      RA_im = centre(JWST_cal_image)[1]
      Dec_im = centre(JWST_cal_image)[2]
    })
    NAXIS1 = JWST_cal_image$keyvalues$NAXIS1
    NAXIS2 = JWST_cal_image$keyvalues$NAXIS2
    Npix = NAXIS1*NAXIS2
    
    #run pro sky routines
    #pro = profoundProFound(JWST_cal_image$imDat, mask=JWST_cal_mask, box=256, redoskysize=101, roughpedestal=TRUE, tolerance=Inf, skycut=3, boxiters = 3)
    pro = profoundProFound(JWST_cal_image$imDat, mask=JWST_cal_mask,
                           skycut=2, pixcut=5, box=512, redoskysize=101,
                           roughpedestal=TRUE, tolerance=Inf)

    #sky_redo = profoundMakeSkyGrid(JWST_cal_image$imDat, objects=pro$objects_redo, mask=JWST_cal_mask, sky=pro$sky, box=256, grid=128, boxiters=3)
    sky_redo = profoundSkyPoly(JWST_cal_image$imDat, objects=pro$objects_redo, degree=2, quancut=0.99, mask=JWST_cal_mask)
    pro_redo = profoundProFound(JWST_cal_image$imDat, mask=JWST_cal_mask, skycut=2, pixcut=5,
                                box=512, redoskysize=101, sky=sky_redo$sky, redosky=FALSE, tolerance=Inf)
    
    
    sky_med = median(sky_redo$sky[JWST_cal_mask == 0 & pro_redo$objects_redo == 0], na.rm=TRUE)
    skyRMS_med = median(pro_redo$skyRMS[JWST_cal_mask == 0 & pro_redo$objects_redo == 0], na.rm=TRUE)
    
    CairoJPEG(filename = paste0(fullbase,'_sky_',obs_info[i,"FILTER"],'.jpeg'), width=1000,height=1000)
    layout(matrix(1:4,2))
    par(mar=c(0.1,0.1,0.1,0.1))
    magimage(JWST_cal_image$imDat, qdiff=T, rem_med=T, axes=FALSE)
    profoundSegimPlot(pro_redo)
    magimage(sky_redo$sky, qdiff=T, rem_med=T, axes=FALSE)
    magimage(pro_redo$skyRMS, axes=FALSE)
    dev.off()
    
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
    
    # return(c(
    #   sky = sky_med,
    #   skyRMS = skyRMS_med,
    #   skyChiSq = pro_redo$skyChiSq,
    #   maskpix = maskpix,
    #   objpix = objpix,
    #   goodpix = goodpix
    #   ))
    return(NULL)
  }
  
}
do_regen_sky_info = function(input_args){
  
  sky_pro_dir = input_args$sky_pro_dir
  cores = input_args$cores_pro
  
  registerDoParallel(cores=cores)
  
  sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
  filelist = list.files(sky_frames_dir, full.names = TRUE)
  filelist = grep('.fits$', filelist, value=TRUE)
  
  foo = list.files(sky_frames_dir, full.names = FALSE)
  vids = substr(grep('.fits$', foo, value=TRUE), 4, 13)
  
  cat('Processing',length(filelist),'files\n')
  
  sky_info = foreach(i = 1:length(filelist), .combine='rbind')%dopar%{
    temp_info = Rfits_read_header(filelist[i])$keyvalues
    temp_info[1:4] = NULL
    temp_info$EXTNAME = NULL
    return(temp_info)
  }
  colnames(sky_info) = tolower(colnames(sky_info))
  
  write.csv(sky_info, paste0(sky_pro_dir, '/sky_info.csv'), row.names=FALSE)
  
}
do_super_sky = function(input_args){
  
  sky_pro_dir = input_args$sky_pro_dir
  VID = input_args$VID
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  
  sky_ChiSq_cut = 1.1
  good_pix_cut = 0.15
  
  registerDoParallel(cores=cores)
  
  sky_frames_dir = paste0(sky_pro_dir, "/sky_frames/")
  sky_super_dir = paste0(sky_pro_dir, "/sky_super/")
  filelist = list.files(sky_frames_dir, full.names = TRUE)
  sky_info = fread(paste0(sky_pro_dir, '/sky_info.csv'))
  
  filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist)]
  sky_info = sky_info[grepl(VID, sky_info$visit_id), ]
  short_det = c("NRCA1", "NRCA2", "NRCA3", "NRCA4", "NRCB1", "NRCB2", "NRCB3", "NRCB4")
  long_det = c("NRCALONG", "NRCBLONG")
  #miri_det = 'MIRIMAGE' #we don't care about MIRI for now
  
  short_filt = sort(unique(sky_info[detector %in% short_det,filter]))
  long_filt = sort(unique(sky_info[detector %in% long_det,filter]))
  
  #short_filt = c("F090W", "F115W", "F150W", "F182M", "F200W", "F210M")
  #long_filt = c("F277W", "F300M", "F335M", "F356W", "F360M", "F410M", "F444W")
  #miri_filt = c("F1000W", "F1500W", "F1800W", "F770W") #we don't care about MIRI for now
  
  combine_grid_short = expand.grid(short_det, short_filt, stringsAsFactors=FALSE)
  combine_grid_long = expand.grid(long_det, long_filt, stringsAsFactors=FALSE)
  #combine_grid_miri = expand.grid(miri_det, miri_filt, stringsAsFactors=FALSE) #we don't care about MIRI for now
  combine_grid = rbind(combine_grid_short, combine_grid_long)
  
  sky_filelist = list.files(sky_info[1,pathsky])
  sky_filelist = sky_filelist[grepl('.fits$', sky_filelist) & grepl(VID, sky_filelist)]
  # sky_filelist = grep('.fits$', sky_filelist, value=TRUE) #just get the FITS files to be safe
  if(do_NIRISS){
    filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist) & grepl("nis", filelist)]
    sky_info = sky_info[grepl(VID, sky_info$visit_id) & grepl("NIS", sky_info$detector), ]
    niriss_det = c("NIS")
    nis_filt = sort(unique(sky_info[detector %in% niriss_det,filter]))
    combine_grid = expand.grid(niriss_det, nis_filt, stringsAsFactors=FALSE)
    sky_filelist = list.files(sky_info[1,pathsky])
    sky_filelist = sky_filelist[grepl('.fits$', sky_filelist) & grepl(VID, sky_filelist) & grepl("nis", sky_filelist)]
  }
  
  dummy = foreach(i = 1:dim(combine_grid)[1], .inorder=FALSE)%dopar%{
    if(combine_grid[i,1] %in% c('NRCA1','NRCA2','NRCA3','NRCA4','NRCB1','NRCB2','NRCB3','NRCB4')){
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
  
  sky_info = sky_info[grepl(VID, sky_info$visit_id), ]
  sky_filelist = list.files(sky_info[1,pathsky])
  sky_filelist = sky_filelist[grepl('.fits$', sky_filelist) & grepl(VID, sky_filelist)] #just get the FITS files to be safe
  
  scan_filt = Rfits_key_scan(paste0(sky_info$pathim, "/", sky_info$fileim), keylist = c("FILTER"))
  sky_info = sky_info[grepl(FILT, scan_filt$FILTER), ]
  sky_filelist = sky_filelist[grepl(FILT, scan_filt$FILTER)]
  cat('Processing',length(sky_filelist),'files\n')
  
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
    
    #edges (here we ignore the first 10 pixels since these are often ratty):
    LHS = median(temp_cal_sky[11:15,], na.rm=TRUE)
    RHS = median(temp_cal_sky[2049-11:15,], na.rm=TRUE)
    BHS = median(temp_cal_sky[,11:15], na.rm=TRUE)
    THS = median(temp_cal_sky[,2049-11:15], na.rm=TRUE)
    #corners (with masking, this gives about as many pixels as the edge above):
    BL = median(temp_cal_sky[11:110,11:110], na.rm=TRUE)
    TL = median(temp_cal_sky[11:110,2049-11:110], na.rm=TRUE)
    TR = median(temp_cal_sky[2049-11:110,2049-11:110], na.rm=TRUE)
    BR = median(temp_cal_sky[2049-11:110,11:110], na.rm=TRUE)
    #centre
    CC = median(temp_cal_sky[974:1073,974:1073], na.rm=TRUE)
    
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
  
  cal_sky_dir = input_args$cal_sky_dir
  cal_sky_renorm_dir = input_args$cal_sky_renorm_dir
  VID = input_args$VID
  FILT = input_args$FILT
  cores = input_args$cores_pro
  do_NIRISS = input_args$do_NIRISS
  
  filelist = c(
    list.files(cal_sky_dir, full.names = TRUE) #If running Pro1/F outputs
  )
  
  if(do_NIRISS){
    filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist) & grepl("nis", filelist)]
  }else{
    filelist = filelist[!grepl('_miri_',filelist)] #we don't care about MIRI for now
    filelist = filelist[!grepl('_nis_',filelist)] #we don't care about MIRI for now
    filelist = filelist[grepl('.fits$', filelist) & grepl(VID, filelist)] #just get the FITS files to be safe
  }

  scan_filt = Rfits_key_scan(filelist = filelist, keylist = c("FILTER"))
  filelist = filelist[grepl(FILT, scan_filt$FILTER)]
  
  cat('Processing',length(filelist),'files\n')
  
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
        
        # current_ped = temp_cal_sky$SCI$keyvalues$SKY_B + temp_cal_sky$SCI$keyvalues$SKY_P
        # new_ped = ped_info[OBS_ID == temp_cal_sky[[1]]$keyvalues$OBS_ID & EXPOSURE == temp_cal_sky[[1]]$keyvalues$EXPOSURE, ped_med]
        # 
        # if((current_ped - new_ped)/new_ped > 0.1){ #if the pedestal differs from the median by more than 20% then fix it
        #   message('Renorm pedestal: ',cal_sky_info[i,file])
        #   
        #   replace_SCI = temp_cal_sky$SCI[,]$imDat + current_ped - new_ped
        #   replace_SKY_Super = temp_cal_sky$SKY_Super[,]$imDat - current_ped + new_ped
        #   replace_SKY_P = temp_cal_sky$SCI$keyvalues$SKY_P - current_ped + new_ped
        #   
        #   message('Testing sky with ProFound: ',cal_sky_info[i,file])
        #   pro_current_sky = profoundProFound(temp_cal_sky$SCI[,]$imDat, mask=temp_cal_sky$DQ[,]$imDat != 0, sky=0, redosky = FALSE)
        #   if(pro_current_sky$skyChiSq < 0.95 | pro_current_sky$skyChiSq > 1.05){
        #     message('Bad Sky: ',cal_sky_info[i,file],'. Current: ',round(pro_current_sky$skyChiSq,2),' Trying to make better sky with ProFound')
        #     if(sum(pro_current_sky$objects) / prod(dim(replace_SCI)) < 0.2){
        #       pro_new_sky = profoundProFound(replace_SCI, mask=temp_cal_sky$DQ[,]$imDat != 0, box=512, roughpedestal = TRUE)
        #       if(pro_new_sky$skyChiSq > 0.95 & pro_new_sky$skyChiSq < 1.05){
        #         message('Using better ProFound Sky: ',cal_sky_info[i,file],' Old: ',round(pro_current_sky$skyChiSq,2),' New: ',round(pro_new_sky$skyChiSq,2))
        #         replace_SCI = replace_SCI - pro_new_sky$sky
        #         replace_SKY_Super = replace_SKY_Super + pro_new_sky$sky
        #         replace_SKY_P = replace_SKY_P + mean(pro_new_sky$sky, na.rm=TRUE)
        #       }
        #     }
        #   }
        #   Rfits_write_pix(replace_SCI, file_cal_sky_renorm, ext=which(ext_names == 'SCI'))
        #   Rfits_write_pix(replace_SKY_Super, file_cal_sky_renorm, ext=which(ext_names == 'SKY_Super'))
        #   Rfits_write_key(file_cal_sky_renorm, keyname='SKY_P', keyvalue=replace_SKY_P, keycomment='SKY P coef in SKY = M.SuperSky + B + P', ext=2)
        # }else{
        message('Testing sky with ProFound: ',cal_sky_info[i,file])
        
        temp_mask = temp_cal_sky$DQ[,]$imDat
        JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
        # if(sum(temp_mask!=0)/(prod(dim(temp_mask))) > 0.8){
        #   JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
        # }else{
        #   JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1 | temp_mask == 4, size=3) 
        # }
        
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
              
              Rfits_write_pix(replace_SCI, file_cal_sky_renorm, ext=which(ext_names == 'SCI'))
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
    temp_MAGZERO_FIX = cal_sky_info[i,MAGZERO] -2.5*log10(as.numeric(pmap[CRDS_CTX=='jwst_1084.pmap' & DETECTOR==cal_sky_info[i,'DETECTOR'] & FILTER==cal_sky_info[i,'FILTER'],PHOTMJSR]) / cal_sky_info[i,PHOTMJSR])
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
  
  VID = input_args$VID
  FILT = input_args$FILT
  ref_dir = input_args$ref_dir
  do_NIRISS = input_args$do_NIRISS
  magzero = input_args$magzero
  cores = input_args$cores_stack
  tasks = input_args$tasks_stack
  
  ## Optional editable params
  clip_tol = c(40,20) #default (80,40) hot and cold
  clip_sigma = 20
  clip_dilate = 3
  doclip = F
  NAXIS_short = 6000 #good for ceers, need 6000 for cosmos web
  NAXIS_long = 3000
  ## Keep these constant so don't use in function args
  
  invar_dir = paste0(ref_dir, "/InVar_Stacks/")
  median_dir = paste0(ref_dir, "/Median_Stacks/")
  sky_frames_dir = paste0(ref_dir, "/sky_pro/sky_frames/")
  dump_dir_stub = paste0(ref_dir, "/dump/")
  orig_cal_sky_info = fread(paste0(ref_dir, "/Pro1oF/cal_sky_info.csv"))
  
  unique_visits = unique(orig_cal_sky_info$VISIT_ID)
  
  temp_vid = VID
  message("Now running stack!")
  for(VID in grep(temp_vid, unique_visits, value = T)){
    
    if(do_NIRISS){
      cal_sky_info = orig_cal_sky_info[grepl(VID, orig_cal_sky_info$VISIT_ID) & grepl("NIS", orig_cal_sky_info$DETECTOR),]
    }else{
      cal_sky_info = orig_cal_sky_info[grepl(VID, orig_cal_sky_info$VISIT_ID) & !grepl("NIS", orig_cal_sky_info$DETECTOR),]
    }
    cal_sky_info$MAGZERO_FIX[is.na(cal_sky_info$MAGZERO_FIX)] = cal_sky_info$MAGZERO[is.na(cal_sky_info$MAGZERO_FIX)]
    
    stack_grid = cal_sky_info[,list(FILTER=unique(FILTER)), keyby=VISIT_ID]
    stack_grid = stack_grid[grepl(FILT, stack_grid$FILTER), ]
    
    print(stack_grid)
    
    module_A_WCS_short = cal_sky_info[grepl('NRCA[1-4]', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    module_B_WCS_short = cal_sky_info[grepl('NRCB[1-4]', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    module_A_WCS_long = cal_sky_info[grepl('NRCALONG', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    module_B_WCS_long = cal_sky_info[grepl('NRCBLONG', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    module_NIS = cal_sky_info[grepl('NIS', DETECTOR),list(CRVAL1=mean(CRVAL1), CRVAL2=mean(CRVAL2), CD1_1=mean(CD1_1), CD1_2=mean(CD1_2)), keyby=VISIT_ID]
    
    module_idx = sapply(list(module_A_WCS_long, module_A_WCS_short, 
                             module_B_WCS_long, module_B_WCS_short,
                             module_NIS), 
                        function(x){y = dim(x)[1]
                        y > 0})
    module_list = c('NRCA_short', 'NRCA_long', 'NRCB_short', 'NRCB_long', 'NIS')[module_idx]
    
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
          
          # sn_mask = matrix(0L, nrow = dim(temp_mask)[1], ncol = dim(temp_mask)[2])
          # rough_mask = temp_mask %% 4 == 0 & temp_mask > 0
          # mask_label = as.matrix(imager::label(imager::as.cimg(rough_mask)))
          # tally_pixels = tabulate(mask_label)
          # sel = which(tally_pixels >= 150) #find the big snowballs
          # mask_loc = which(matrix(mask_label %in% sel, dim(temp_mask)[1], dim(temp_mask)[2]), arr.ind=TRUE)
          # sn_mask[mask_loc] = 1L
          # snowball_mask_dilate = profoundDilate(sn_mask, size = 21) #extra dilation on big snowballs
          # naxis = dim(sn_mask)[1]
          # snowball_mask_dilate[1:20, ] = 0
          # snowball_mask_dilate[, 1:20] = 0
          # snowball_mask_dilate[(naxis):(naxis-20), ] = 0
          # snowball_mask_dilate[, (naxis):(naxis-20)] = 0
          # temp_image$imDat[snowball_mask_dilate==1] = NA
          
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
                                           keyvalues_out = output_stack$image$keyvalues)
        
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
  
  filelist = input_args$filelist
  VID = input_args$VID
  median_dir = input_args$median_dir
  cores = input_args$cores_pro
  
  keep_trend_data  = input_args$keep_trend_data

  SIGMA_LO = input_args$SIGMA_LO
  
  wisp_poly = Rfits_read("wisp_poly.fits")
  
  filelist = grep(VID, filelist, value = T)
  cat(filelist, sep = "\n")
  
  info = Rfits_key_scan(filelist = filelist,
                        keylist=c('DETECTOR', 'MODULE', 'FILTER', 'VISIT_ID'), cores = cores)
  info_wisp = info[DETECTOR %in% c('NRCA3','NRCA4','NRCB3','NRCB4'),]
  mod_visit_grid = unique(info_wisp[,c("MODULE", "VISIT_ID")])
  
  ref_im_list = {}
  for(ii in 1:dim(mod_visit_grid)[1]){
    ref_files = list.files(median_dir,
                           pattern = glob2rx(paste0(
                             "*", mod_visit_grid$VISIT_ID[ii], "*", mod_visit_grid$MODULE[ii], "*", "short", "*.fits"
                           )),
                           full.names = T)
    filter_long = c(
      str_match(ref_files, "F\\s*(.*?)\\s*M")[,2],
      str_match(ref_files, "F\\s*(.*?)\\s*W")[,2]
      )
      ref_file_long = ref_files[ which.max(filter_long[!is.na(filter_long)]) ] #Get the longest filter
    
    ref_im_list[[paste0(mod_visit_grid$VISIT_ID[ii])]][[paste0("NRC", mod_visit_grid$MODULE[ii])]] = Rfits_point(ref_file_long)
    message(paste("Loading reference for:", paste0(mod_visit_grid$VISIT_ID[ii]), paste0("NRC", mod_visit_grid$MODULE[ii])))
  }
  
  registerDoParallel(cores = cores)  
  temp = foreach(ii = 1:dim(info_wisp)[1])%dopar%{
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
    
    ref_im = ref_im_list[[vid]][[modl]]
    
    if(any( (vid == keep_trend_data$ID_vlarge$VISIT_ID & modl == keep_trend_data$ID_vlarge$MODULE) | 
            (vid == keep_trend_data$ID_large$VISIT_ID & modl == keep_trend_data$ID_large$MODULE)) ){
      sigma_lo = NULL
    }else{
      sigma_lo = SIGMA_LO
    }
    
    message("I am using this sigma_lo = ", ifelse(is.null(sigma_lo), "NULL", sigma_log))
    Sys.sleep(time = 2)

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
do_patch = function(input_args){
  
  VID = input_args$VID
  FILT = input_args$FILT
  invar_dir = input_args$invar_dir
  median_dir = input_args$median_dir
  patch_dir = input_args$patch_dir
  cores = input_args$cores_stack
  do_NIRISS = input_args$do_NIRISS
  
  registerDoParallel(cores = cores)
  patch_stub = patch_dir
  
  if(do_NIRISS){
    pixscale_list = c("NIS")
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
  
  VID = input_args$VID
  ref_dir = input_args$ref_dir
  patch_dir = input_args$patch_dir
  
  blue_filters = "F070W|F090W|F115W|F150W|F140M|F162"
  green_filters = "F200W|F277W|F182M|F210M"
  red_filters = "F356W|F444W|F250M|F300M|F335M|F360M|F410M|F430|F460M|F480M"
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
  ref_im_warp = propaneWarp(ref_im, keyvalues_out=wisp_im$keyvalues, cores=cores)
  
  
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

