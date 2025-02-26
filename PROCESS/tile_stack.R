## Deep stacker + Alignment

library(Rfits)
library(Rwcs)
library(ProPane)
library(ProFound)
library(celestial)
library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(Cairo)

## User defined inputs here
ref_cat = fread("put astrometric catalogue here") #example for HST catalogue in EGS field

names(ref_cat) = c("RA", "Dec")

input_args = list(
  ref_dir = "", ## Directory containing the ProPane/ProFound/JUMPROPE stuff
  super_name = "", ## Name of mosaic. Advise against including "_"  

  ref_cat = ref_cat, ## Reference catalogue for alignment
  RA = colMeans(ref_cat)[1], ## Coords to search for frames
  Dec = colMeans(ref_cat)[2],    
  
  calstack = FALSE, ## Make a mosaic from the input cal_sky_renorm files. Invalidates 'grid_size' argument, sets it to "". 
  grid_size = "short", ## short (0.03 asec/pixel) or long (0.06 asec/pixel) VISITID propane stack files to fetch 
  
  mosaic_size = NULL, ## Set the size of the search radius/mosaic size
  rotation = "North", ## Rotation argument for ProPaneGenWCS, e.g., "North" for North-Up alignment
  pixelscale = NULL, ## If we want to define our own pixelscale in asec/pixel
  magzero_out = 23.9, ## Magzero of final mosaic. Default means matrix will have flux units of microJy

  program_id = NULL,  # Maybe we we just want to stack a single program e.g., Primer (1837) and not COSMOS Web (1727)
  cores = 1 ## If NULL then use half the number of cores in the system
) 


deep_stacker = function(input_args){
  
  ## Find aligning frames in search radius, align with ProPaneTweakCat and stack with ProPane
  
  inVar_dir = paste0(input_args$ref_dir, "/InVar_Stacks")
  median_dir = paste0(input_args$ref_dir, "/Median_Stacks")
  patch_dir = paste0(input_args$ref_dir, "/Patch_Stacks")
  dump_dir = paste0(input_args$ref_dir, "/dump")
  cal_sky_renorm_dir = paste0(input_args$ref_dir, "/Pro1oF/cal_sky_renorm/")
  sky_frames_dir = paste0(input_args$ref_dir, "/sky_pro/sky_frames/")
  
  mosaic_dir = paste0(input_args$ref_dir, "/Mosaic_Stacks/")
  mosaic_invar = paste0(mosaic_dir, "/inVar/")
  mosaic_med = paste0(mosaic_dir, "/Median/")
  mosaic_patch = paste0(mosaic_dir, "/Patch/")
  
  dir.create(mosaic_dir, recursive = T)
  dir.create(mosaic_invar, recursive = T)
  dir.create(mosaic_med, recursive = T)
  dir.create(mosaic_patch, recursive = T)
  
  cal_sky_info = fread(
   paste0(
     input_args$ref_dir, "/Pro1oF/cal_sky_info.csv"
   )
  )

  VIDS = paste0(cal_sky_info$VISIT_ID)
  
  if(is.null(input_args$mosaic_size)){
    search_rad = 2.0
  }else{
    search_rad = input_args$mosaic_size
  }
  if(is.null(input_args$program_id)){
    program_id = ""
  }else{
    program_id = input_args$program_id
  }
  
  if(input_args$calstack){
    ## Pull the files form the cal_sky_renorm directory. Must have the cal_sky_info.csv
    
    list_files = cal_sky_info$full[
      grepl(program_id, cal_sky_info$VISIT_ID)
    ]
    list_files_names = cal_sky_info$stub[
      grepl(program_id, cal_sky_info$VISIT_ID)
    ]
    files_pids = substr(list_files_names, 4, 4+9) ## just get the pid
    not_pid_idx = sapply(files_pids, function(x)grepl(x, program_id, fixed = T))
    if(is.null(input_args$program_id)){
      not_pid_idx = 1:length(list_files)
    }
    message("Searching for frames in ", input_args$RA, " ", input_args$Dec, " rad=", search_rad, "deg")
    find_frames = propaneFrameFinder(
      filelist = list_files[not_pid_idx],
      RAcen = input_args$RA,
      Deccen = input_args$Dec,
      rad = search_rad,
      cores = 1, ## safe :)
      plot=F, 
      extlist = 2 ## WCS for cal_sky_renorm
    )
    
    find_frames$FILT = Rfits_key_scan(
      find_frames$full, 
      keylist = c("FILTER"),
      extlist = 1
    )$FILTER
    find_frames$MAGZERO_FIX = cal_sky_info$MAGZERO_FIX[cal_sky_info$full %in% find_frames$full]
    
  }else{
    ## Pull the files form the stack_patch directory
    list_files = list.files(
      path = patch_dir, 
      pattern = glob2rx(paste0("*stack_patch*", input_args$program_id, "*.fits")),
      full.names = T
    )
    
    list_files_names = list.files(
      path = patch_dir, 
      pattern = glob2rx(paste0("*stack_patch*", input_args$program_id, "*.fits")),
      full.names = F
    )
    
    files_pids = substr(list_files_names, 1+12, 4+12) ## just get the pid
    
    not_pid_idx = sapply(files_pids, function(x)grepl(x, program_id, fixed = T))
    if(is.null(input_args$program_id)){
      not_pid_idx = 1:length(list_files)
    }
    
    message("Searching for frames in ", input_args$RA, " ", input_args$Dec, " rad=", search_rad, "deg")
    find_frames = propaneFrameFinder(
      filelist = list_files[not_pid_idx],
      RAcen = input_args$RA,
      Deccen = input_args$Dec,
      rad = search_rad,
      cores = 1, ## safe :)
      plot=F
    )
    filters_w2 = str_extract(find_frames$file, pattern = regex("F(\\d+).W2")) ## string match F and W with only digits in between 
    filters_all = str_extract(find_frames$file, pattern = regex("F(\\d+).W|F(\\d+).M"))
    filters_all[which(filters_w2 == "F150W2")] = "F150W2"
    filters_all = na.omit(filters_all)
    
    find_frames$FILT = filters_all
  }
  
  UNIQUE_VIDS = VIDS[sapply(VIDS, function(x){any(grepl(x,find_frames$file))} )==TRUE]

  stack_grid = cal_sky_info[paste0(cal_sky_info$VISIT_ID) %in% UNIQUE_VIDS, c("VISIT_ID", "FILTER", "DETECTOR")]
 
  if(input_args$calstack){
    input_args$grid_size = c("") ## So grep all of the cal_sky_renorm files 
  }
  for(grid_size in c(input_args$grid_size)){
    
    file_info = find_frames[grepl(paste0(grid_size, "|MIRI"), find_frames$file), ] ## Force find MIRI
    
    if(input_args$calstack){
      wcs = propaneGenWCS(
        filelist = file_info$full, 
        rotation = input_args$rotation, 
        pixscale = input_args$pixelscale, 
        extlist = 2
        )
    }else{
      wcs = propaneGenWCS(
        filelist = file_info$full, 
        rotation = input_args$rotation, 
        pixscale = input_args$pixelscale, 
        extlist = 1
        )
    }
    
    ref_cat = input_args$ref_cat
    
    UNIQUE_FILTERS = unique(stack_grid$FILTER)
    UNIQUE_FILTERS = UNIQUE_FILTERS[!grepl("CLEAR", UNIQUE_FILTERS) & UNIQUE_FILTERS %in% file_info$FILT]
    for(i in 1:length(UNIQUE_FILTERS)){
      
      file_info_filt = file_info[file_info$FILT == paste0(UNIQUE_FILTERS[i]),]
      filenames = file_info_filt$full
 
      if(input_args$calstack){

        image_list = lapply(1:dim(file_info_filt)[1], function(kk){
          ## Same masking as in gen_stack
          temp_image = Rfits_read_image(file_info_filt[k,full], ext=2)
          temp_mask = Rfits_read_image(file_info_filt[k,full], ext=4, header=FALSE)
          temp_image = propaneBadPix(
            image = temp_image,
            patch = T
          )
          JWST_cal_mask = profoundDilate(temp_mask %% 2 == 1, size=3)
          temp_image$imDat[JWST_cal_mask==1] = NA
          return(temp_image)
        })
        
        inVar_list = lapply(1:dim(file_info_filt)[1], function(kk){
          sky_file = paste0(sky_frames_dir, sub('_rem', '', file_info_filt$stub[kk]),'_',file_info_filt$FILT[kk],'.fits')
          Rfits::Rfits_read_key(sky_file, 'SKYRMS')^-2
        }) ## Same logic that is in the gen_stack module in all_codes_process.R
        
        weight_list = lapply(1:dim(file_info_filt)[1], function(x){
        matrix(1L, file_info_filt$dim_1[x], file_info_filt$dim_2[x]) ## Might need this for MIRI
        })
        
        magzero_list = file_info_filt$MAGZERO_FIX
      }else{
        image_extlist = sapply(filenames, function(x)Rfits_extname_to_ext(x, extname = "image"))
        inVar_extlist = sapply(filenames, function(x)Rfits_extname_to_ext(x, extname = "inVar"))
        weight_extlist = sapply(filenames, function(x)Rfits_extname_to_ext(x, extname = "weight"))
        
        image_list = Rfits_make_list(
          filelist = filenames,
          extlist = image_extlist,
          pointer = T
        )
        weight_list = Rfits_make_list(
          filelist = filenames,
          extlist = weight_extlist,
          pointer = T,
          header = F
        )
        inVar_list = Rfits_make_list(
          filelist = filenames,
          extlist = inVar_extlist,
          pointer = T,
          header = F
        )
        magzero_list = rep(23.9, length(filenames))
      }
      
      if(i == 1){
        CairoPNG(paste0(mosaic_dir,input_args$super_name, "_footprints.png"))
        for(k in 1:length(image_list)){
          if(k == 1){
            Rwcs_overlap(
              keyvalues_test = image_list[[k]]$keyvalues, keyvalues_ref = wcs, plot = T
            )
          }else{
            Rwcs_overlap(
              keyvalues_test = image_list[[k]]$keyvalues, keyvalues_ref = wcs, plot = T, add = T
            )
          }
        }
        dev.off()
      }
      
      for(k in 1:length(image_list)){
        
        message(
          paste0("Tweaking ", filenames[k])
        )
        ## Default ProFound source detection, ends up being ~1.5 sigma source detect level
        pro_test = profoundProFound( 
          image = image_list[[k]][,],
          magzero = magzero_list[k],
          mask = is.infinite(image_list[[k]][,]$imDat),
          rem_mask = T
        )
        
        mag = ifelse(is.na(pro_test$segstats$mag), 35, pro_test$segstats$mag)
        match_cat = coordmatch(
          coordref = ref_cat,
          coordcompare = pro_test$segstats[mag < quantile(mag[mag<35], 0.5), c("RAmax", "Decmax")] ## tweak on sources brighter than 50th percentile
        )
        
        if(any(is.na(match_cat$bestmatch))){
          tweak_cat = list(par = c(0,0,0)) ## No matches => nothing to tweak on unfortunately...
        }else{
          pix_cat = Rwcs_s2p(RA = ref_cat$RA, Dec = ref_cat$Dec, keyvalues = image_list[[k]]$keyvalues)
          tweak_cat = propaneTweakCat(cat_ref = pix_cat[match_cat$bestmatch$refID, ], 
                                      cat_pre_fix = pro_test$segstats[match_cat$bestmatch$compareID,c("xmax","ymax")],
                                      delta_max = c(20,1),
                                      mode = "pix")
        }
        
        tempim = image_list[[k]][,]
        tempim$imDat[is.infinite(tempim$imDat)] = NA
        image_list[[k]] = propaneWCSmod(
          tempim,
          delta_x = tweak_cat$par[1],
          delta_y = tweak_cat$par[2],
          delta_rot = tweak_cat$par[3]
        )
        rm(tempim)
        
        if(!input_args$calstack){
          ## Have to adjust the WCS for Rfits objects when doing patch_stack mosaicking
          weight_list[[k]] = propaneWCSmod(
            weight_list[[k]][,],
            delta_x = tweak_cat$par[1],
            delta_y = tweak_cat$par[2],
            delta_rot = tweak_cat$par[3]
          )
          
          inVar_list[[k]] = propaneWCSmod(
            inVar_list[[k]][,],
            delta_x = tweak_cat$par[1],
            delta_y = tweak_cat$par[2],
            delta_rot = tweak_cat$par[3]
          )
        }
      }
      
      cores = ifelse(is.null(input_args$cores), floor(detectCores()/2), input_args$cores)
      deep_frame = propaneStackWarpInVar(
        image_list = image_list,
        inVar_list = inVar_list,
        weight_list = weight_list,
        keyvalues_out = wcs,
        magzero_in = magzero_list, ## Use appropriate magzero ins
        magzero_out = input_args$magzero_out,
        direction="backward",
        dump_frames = T,
        dump_dir = paste0(dump_dir, "/", input_args$super_name, "/", grid_size, "/", UNIQUE_FILTERS[i], "/"),
        multitype = "cluster",
        cores = cores,
        cores_warp = 1
      )
      
      deep_frame$image$imDat[is.infinite(deep_frame$image$imDat)] = NA
      deep_frame$info = file_info_filt
      
      Rfits_write(
        deep_frame,
        paste0(mosaic_invar, "/stack_", input_args$super_name, "_", UNIQUE_FILTERS[i], "_", grid_size, ".fits")
      )
      
      med_stack = propaneStackWarpMed(
        dirlist = deep_frame$dump_dir,
        pattern = glob2rx("*image*.fits"),
        keyvalues_out = deep_frame$image$keyvalues, 
        cores = cores, 
        multitype = "cluster"
      )
      
      Rfits_write(
        med_stack,
        paste0(mosaic_med, "/med_", input_args$super_name, "_", UNIQUE_FILTERS[i], "_", grid_size, ".fits"),
      )
      
      patch_stack = propanePatch(
        image_inVar = deep_frame$image,
        image_med = med_stack$image
      )
      temp = deep_frame
      temp$image = patch_stack$image
      
      Rfits_write(
        temp,
        paste0(mosaic_patch, "/patch_", input_args$super_name, "_", UNIQUE_FILTERS[i], "_", grid_size, ".fits")
      )
      
      rm(deep_frame)
      rm(med_stack)
      rm(patch_stack)
      rm(temp)
      gc()
    }
    
  }
  
}

deep_stacker(input_args = input_args)

