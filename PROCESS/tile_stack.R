## Deep stacker + Alignment

input_args = list(
  ref_dir = "/Volumes/RAIDY/JWST/",
  RA = 265.0355633,
  Dec = 68.974166,
  mosaic_size_deg = 10/60,
  super_name = "IDF"
)

propane_tweak_widget = function(load_frames, wcs){
  #chain up tweak solutions
  
  N_frames = length(load_frames)
  tweak_pars = {}
  counter = 0
  while(counter < N_frames){
    message("Counter = ", counter)
    if(counter == 0){
      message("Initialise reference frame.")
      ref_frame = propaneWarp(
        image_in = load_frames[[1]],
        keyvalues_out = wcs
      )
      tweak_pars = c(tweak_pars, list(c(0,0,0)))
      counter = counter + 1
      # store_frames = c(list(ref_frame), store_frames)
    }else{
      message("Tweak test frame.")
      
      message("Warp test to target WCS")
      warp_frame = propaneWarp(
        load_frames[[counter+1]],
        keyvalues_out = wcs
      )
      
      message("Make initial stack for weight selection")
      init_stack = propaneStackWarpInVar(
        image_list = list(ref_frame, warp_frame),
        cores = 1,
        keyvalues_out = wcs,
        magzero_in = 23.9,
        magzero_out = 23.9
      )
      
      array_idx = which(init_stack$weight$imDat >= 2, arr.ind = T)
      
      xmid_cent = quantile(array_idx[,1], 0.5, na.rm = T)
      ymid_cent = quantile(array_idx[,2], 0.5, na.rm = T)
      
      box = 2000
      
      if(
        sum(
          is.na(ref_frame[xmid_cent, ymid_cent, box = box]$imDat)
        )/prod(
          dim(
            ref_frame[xmid_cent, ymid_cent, box = box]$imDat
          )
        ) >= 0.8
      ){
        xmid_cent = quantile(array_idx[,1], 0.25, na.rm = T)
        ymid_cent = quantile(array_idx[,2], 0.25, na.rm = T) 
      }
      
      message("Tweaking for integer shifts...")
      tweak_soln = propaneTweak(
        image_pre_fix = warp_frame[xmid_cent, ymid_cent, box = box],
        image_ref = ref_frame[xmid_cent, ymid_cent, box = box],
        WCS_match = T,
        cutcheck = F,
        shift_int = T, 
        delta_max = c(10, 0.2),
        cores = 8
      )
      
      tweak_vals_translate = tweak_soln$optim_out$par
      # tweak_vals = c(tweak_soln$optim_out$par)
      
      new_frame_translate = propaneWCSmod(
        warp_frame,
        delta_x = tweak_vals_translate[1], 
        delta_y = tweak_vals_translate[2],
        delta_rot = 0
      )
      
      temp_ref = ref_frame
      temp_ref$imDat[((xmid_cent-box):(xmid_cent+box)), ((ymid_cent-box):(ymid_cent+box))] = NA
      
      temp_new_frame_translate = new_frame_translate
      temp_new_frame_translate$imDat[((xmid_cent-box):(xmid_cent+box)), ((ymid_cent-box):(ymid_cent+box))] = NA
      
      message("Reference source finding...")
      pro_ref = profoundProFound(
        image = temp_ref,
        skycut = 5.0,
        rem_mask = T,
        sky = 0,
        redosky = F, 
        box = 1000,
        grid = 1000,
        verbose = F,
        iters = 0
      )
      
      message("PreFix source finding")
      pro_test = profoundProFound(
        image = temp_new_frame_translate,
        skycut = 5.0,
        rem_mask = T,
        sky = 0,
        redosky = F, 
        box = 1000,
        grid = 1000,
        verbose = F,
        iters = 0
      )
      
      message("Computing delta rot accounting for xy shifts...")
      match_idx = coordmatch(
        coordref = pro_ref$segstats[, c("RAmax", "Decmax")],
        coordcompare = pro_test$segstats[, c("RAmax", "Decmax")],
        rad = 1.0
      )
      
      theta_ref = atan2(
        pro_ref$segstats$ycen[match_idx$bestmatch$refID], pro_ref$segstats$xcen[match_idx$bestmatch$refID]
      )
      
      theta_test = atan2(
        pro_test$segstats$ycen[match_idx$bestmatch$compareID] + tweak_vals_translate[1], pro_test$segstats$xcen[match_idx$bestmatch$compareID] + tweak_vals_translate[2]
      )
      
      delta_phi = -median(theta_ref - theta_test) * 180/pi
      
      tweak_vals = c(tweak_vals_translate, delta_phi)
      
      new_frame = propaneWCSmod(
        warp_frame,
        delta_x = tweak_vals[1], 
        delta_y = tweak_vals[2],
        delta_rot = tweak_vals[3]
      )
      
      message("Making new super reference frame. Super mosaic.")
      new_ref_frame = propaneStackWarpInVar(
        image_list = list(ref_frame, new_frame),
        cores = 1,
        keyvalues_out = wcs,
        magzero_in = 23.9, 
        magzero_out = 23.9
      )
      
      tweak_pars = c(tweak_pars, list(tweak_vals))
      
      ref_frame = new_ref_frame$image
      counter = counter + 1
    }
  }
  return(tweak_pars)
}

deep_stacker = function(input_args){
  
  library(Rfits)
  library(Rwcs)
  library(ProPane)
  library(ProFound)
  library(celestial)
  library(foreach)
  library(doParallel)
  library(data.table)
  library(stringr)
  
  inVar_dir = paste0(input_args$ref_dir, "/InVar_Stacks")
  median_dir = paste0(input_args$ref_dir, "/Median_Stacks")
  patch_dir = paste0(input_args$ref_dir, "/Patch_Stacks")
  dump_dir = paste0(input_args$ref_dir, "/dump")
  
  frames_info = fread(
    paste0(
      input_args$ref_dir, "/FieldFootprints/frames_info.csv"
    )
  )
  
  cal_sky_info = fread(
    paste0(
      input_args$ref_dir, "/Pro1oF/cal_sky_info.csv"
    )
  )
  
  VIDS = paste0(frames_info$VISIT_ID)
  
  message("Searching for frames in ", input_args$RA, " ", input_args$Dec)
  find_frames = propaneFrameFinder(
    dirlist = patch_dir,
    RAcen = input_args$RA,
    Deccen = input_args$Dec,
    rad = input_args$mosaic_size_deg,
    cores = 8,
    plot=F
  )
  
  UNIQUE_VIDS = VIDS[sapply(VIDS, function(x){any(grepl(x,find_frames$file))} )==TRUE]
  
  stack_grid = cal_sky_info[paste0(cal_sky_info$VISIT_ID) %in% UNIQUE_VIDS, c("VISIT_ID", "FILTER", "DETECTOR")]
  
  for(grid_size in c("long")){
    
    file_info = find_frames[grepl(grid_size, find_frames$file), ]
    
    wcs_info = Rfits_key_scan(
      filelist = file_info$full,
      keylist = c("CRVAL1",
                  "CRVAL2",
                  "CD1_1",
                  "CD1_2",
                  "CD2_1",
                  "CD2_2",
                  "NAXIS1",
                  "NAXIS2",
                  "CRPIX1",
                  "CRPIX2"),
      get_pixscale = T,
      get_centre = T,
      cores = 8
    )
    
    arith_mean_RA = mean(wcs_info$centre_RA)
    arith_mean_Dec = mean(wcs_info$centre_Dec)
    
    max_min_RA = (max(wcs_info$centre_RA)+min(wcs_info$centre_RA))/2.0
    max_min_Dec = (max(wcs_info$centre_Dec)+min(wcs_info$centre_Dec))/2.0
    
    RA_val = sum(arith_mean_RA,max_min_RA)/2.0
    Dec_val = sum(arith_mean_Dec,max_min_Dec)/2.0
    
    wcs = Rwcs_keypass(
      CRVAL1 = RA_val,
      CRVAL2 = Dec_val,
      CD1_1 = median(wcs_info$CD1_1),
      CD1_2 = median(wcs_info$CD1_2),
      CD2_1 = median(wcs_info$CD2_1),
      CD2_2 = -median(wcs_info$CD1_1)
    )
    
    wcs$NAXIS1 = ceiling(input_args$mosaic_size_deg / median(wcs_info$pixscale/3600))
    wcs$NAXIS2 = ceiling(input_args$mosaic_size_deg / median(wcs_info$pixscale/3600))
    wcs$CRPIX1 = wcs$NAXIS1/2.0
    wcs$CRPIX2 = wcs$NAXIS2/2.0
  
    vid_filt_stacks = foreach(i = 1:length(UNIQUE_VIDS))%do%{
      message(
        "Stacking VID = ", UNIQUE_VIDS[i], " ", grid_size
      )
      #basically make module A/B mosaic for each VID in 
      # frame_stubs = str_replace(
      #   paste0(patch_dir, "/", grep(UNIQUE_VIDS[i], file_info$file, value = T) ),
      #   "stack_", "stack_patch_"
      # )
      frame_stubs = paste0(patch_dir, "/", grep(UNIQUE_VIDS[i], file_info$file, value = T) )
      
      wcs_info_temp = Rfits_key_scan(
        filelist = frame_stubs,
        keylist = c("CRVAL1",
                    "CRVAL2",
                    "CD1_1",
                    "CD1_2",
                    "CD2_1",
                    "CD2_2",
                    "NAXIS1",
                    "NAXIS2",
                    "CRPIX1",
                    "CRPIX2")
      )
      wcs_temp = Rwcs_keypass(
        CRVAL1 = median(wcs_info_temp$CRVAL1),
        CRVAL2 = median(wcs_info_temp$CRVAL2),
        CD1_1 = median(wcs_info_temp$CD1_1),
        CD1_2 = median(wcs_info_temp$CD1_2),
        CD2_1 = median(wcs_info_temp$CD2_1),
        CD2_2 = -median(wcs_info_temp$CD1_1)
      )
      
      load_images = Rfits_make_list(
        filelist = frame_stubs,
        pointer = T,
        extlist = 1
      )
      
      non_nan_frac = 1 - sum(is.na(load_images[[1]][,]$imDat) / prod(dim(load_images[[1]])))
      
      wcs_temp$NAXIS1 = max(wcs_info_temp$NAXIS1) * non_nan_frac * 3.5
      wcs_temp$NAXIS2 = max(wcs_info_temp$NAXIS2) * non_nan_frac * 1.5
      wcs_temp$CRPIX1 = wcs_temp$NAXIS1/2.0
      wcs_temp$CRPIX2 = wcs_temp$NAXIS2/2.0
      
      load_weights = Rfits_make_list(
        filelist = frame_stubs,
        pointer = T,
        extlist = 4,
        header = F
      )
      
      load_inVar = Rfits_make_list(
        filelist = frame_stubs,
        pointer = T,
        extlist = 3,
        header = F
      )
      
      stack_all_filt = propaneStackWarpInVar(
        image_list = load_images,
        weight_list = load_weights,
        inVar_list = load_inVar,
        keyvalues_out = wcs_temp,
        magzero_in = 23.9,
        magzero_out = 23.9,
        multitype = "cluster",
        cores = 1
      )
      
      stack_all_filt$image
    }
    vid_filt_stacks_stubs = foreach(i = 1:length(UNIQUE_VIDS), .combine = "c") %do%{
      Rfits_write_image(
        vid_filt_stacks[[i]],
        filename = paste0(
          dump_dir, "/", UNIQUE_VIDS[i], "/temp_", i, ".fits"
        )
      )
      paste0(
        dump_dir, "/", UNIQUE_VIDS[i], "/temp_", i, ".fits"
      )
    }
    rm(vid_filt_stacks)
    gc()
    vid_filt_stacks = Rfits_make_list(
      vid_filt_stacks_stubs,
      pointer = T
    )
    
    propane_tweak = propane_tweak_widget(load_frames = vid_filt_stacks, wcs = wcs)

    file.remove(
      vid_filt_stacks_stubs
    )
    
    UNIQUE_FILTERS = unique(stack_grid$FILTER)
    for(i in 1:length(UNIQUE_FILTERS)){
      
      filenames = file_info$full[grepl(paste0(UNIQUE_FILTERS[i]), file_info$full)]
      
      tweaked_frames = foreach(j = 1:length(UNIQUE_VIDS)) %do% {
        
        VID_filenames = filenames[grepl(paste0(UNIQUE_VIDS[j]), filenames)]
        if(length(VID_filenames)==0){
          message("No filters")
          return(NULL)
        }else{
          message("Tweaking: ", VID_filenames)
          
          image_extlist = sapply(VID_filenames, function(x)Rfits_extname_to_ext(x, extname = "image"))
          inVar_extlist = sapply(VID_filenames, function(x)Rfits_extname_to_ext(x, extname = "inVar"))
          weight_extlist = sapply(VID_filenames, function(x)Rfits_extname_to_ext(x, extname = "weight"))

          image_list = Rfits_make_list(
            filelist = VID_filenames,
            extlist = image_extlist,
            pointer = T
          )
          weight_list = Rfits_make_list(
            filelist = VID_filenames,
            extlist = weight_extlist,
            pointer = T,
            header = F
          )
          inVar_list = Rfits_make_list(
            filelist = VID_filenames,
            extlist = inVar_extlist,
            pointer = T,
            header = F
          )

          temp_stack = propaneStackWarpInVar(
            image_list = image_list,
            inVar_list = inVar_list,
            weight_list = weight_list,
            magzero_in = 23.9,
            magzero_out = 23.9,
            keyvalues_out = wcs
          )
          
          tweak_stack = temp_stack
          for(key in c("image","inVar","weight")){
            tweak_stack[[key]] = propaneWCSmod(
              tweak_stack[[key]],
              delta_x = propane_tweak[[j]][1],
              delta_y = propane_tweak[[j]][2],
              delta_rot = propane_tweak[[j]][3]
            )
          }
          tweak_stack
        }
      }
      
      tweaked_frames = tweaked_frames[!sapply(tweaked_frames, is.null)]
      
      loc = which(!is.na(tweaked_frames[[1]]$image$imDat),arr.ind=TRUE)
      xlo = min(loc[,1])-100
      xhi = max(loc[,1])+100
      ylo = min(loc[,2])-100
      yhi = max(loc[,2])+100
      
      wcs_crop_temp = tweaked_frames[[1]]$image[xlo:xhi, ylo:yhi]$keyvalues
      
      wcs_crop = Rwcs_keypass(
        CRVAL1 = wcs_crop_temp$CRVAL1,
        CRVAL2 = wcs_crop_temp$CRVAL2,
        
        CD1_1 = wcs_crop_temp$CD1_1,
        CD1_2 = wcs_crop_temp$CD1_2,
        CD2_1 = wcs_crop_temp$CD2_1,
        CD2_2 = wcs_crop_temp$CD2_2,
      )
      wcs_crop$NAXIS1 = wcs_crop_temp$NAXIS1
      wcs_crop$NAXIS2 = wcs_crop_temp$NAXIS2
      wcs_crop$CRPIX1 = wcs_crop_temp$CRPIX1
      wcs_crop$CRPIX2 = wcs_crop_temp$CRPIX2
      
      deep_frame = propaneStackWarpInVar(
        image_list = lapply(tweaked_frames, function(x)x$image),
        inVar_list = lapply(tweaked_frames, function(x)x$inVar),
        weight_list = lapply(tweaked_frames, function(x)x$weight),
        keyvalues_out = wcs_crop,
        magzero_in = 23.9,
        magzero_out = 23.9,
        direction="backward",
        dump_frames = T,
        dump_dir = paste0(dump_dir, "/", input_args$super_name, "/", grid_size, "/", UNIQUE_FILTERS[i], "/"),
        multitype = "cluster",
        cores = 1,
        cores_warp = 1
      )

      Rfits_write(
        deep_frame,
        paste0(inVar_dir, "/stack_", input_args$super_name, "_", UNIQUE_FILTERS[i], "_", grid_size, ".fits")
      )
      
      med_stack = propaneStackWarpMed(
        dirlist = deep_frame$dump_dir,
        pattern = glob2rx("*image*.fits"),
        keyvalues_out = deep_frame$image$keyvalues
      )

      Rfits_write(
        med_stack,
        paste0(median_dir, "/med_", input_args$super_name, "_", UNIQUE_FILTERS[i], "_", grid_size, ".fits"),
      )
      
    }
    
  }

}

deep_stacker(input_args = input_args)



## spare stuff
propane_tweak_widget_2 = function(load_frames, wcs){
  #chain up tweak solutions
  
  N_frames = length(load_frames)
  tweak_pars = {}
  counter = 0
  while(counter < N_frames){
    message("Counter = ", counter)
    if(counter == 0){
      message("Initialise reference frame.")
      ref_frame = propaneWarp(
        image_in = load_frames[[1]],
        keyvalues_out = wcs
      )
      tweak_pars = c(tweak_pars, list(c(0,0,0)))
      counter = counter + 1
      # store_frames = c(list(ref_frame), store_frames)
    }else{
      message("Tweak test frame.")
      
      message("Warp test to target WCS")
      warp_frame = propaneWarp(
        load_frames[[counter+1]],
        keyvalues_out = wcs
      )
      
      message("Make initial stack for weight selection")
      init_stack = propaneStackWarpInVar(
        image_list = list(ref_frame, warp_frame),
        cores = 1,
        keyvalues_out = wcs,
        magzero_in = 23.9,
        magzero_out = 23.9
      )
      
      array_idx = which(init_stack$weight$imDat >= 2, arr.ind = T)
      
      xmid_cent = quantile(array_idx[,1], 0.5, na.rm = T)
      ymid_cent = quantile(array_idx[,2], 0.5, na.rm = T)
      
      box = c(4000,2000)
      
      if(
        sum(
          is.na(ref_frame[xmid_cent, ymid_cent, box = box]$imDat)
        )/prod(
          dim(
            ref_frame[xmid_cent, ymid_cent, box = box]$imDat
          )
        ) >= 0.8
      ){
        xmid_cent = quantile(array_idx[,1], 0.16, na.rm = T)
        ymid_cent = quantile(array_idx[,2], 0.16, na.rm = T) # should hopefully be robust against fields like the IDF in program 2738 :D
      }
      
      message("Tweaking for integer shifts...")
      tweak_soln = propaneTweak(
        image_pre_fix = warp_frame[xmid_cent, ymid_cent, box = box],
        image_ref = ref_frame[xmid_cent, ymid_cent, box = box],
        WCS_match = T,
        cutcheck = F,
        shift_int = F, 
        delta_max = c(10, 0.5),
        cores = 8
      )
      
      tweak_vals = tweak_soln$optim_out$par
      # tweak_vals = c(tweak_soln$optim_out$par)
      
      new_frame = propaneWCSmod(
        warp_frame,
        delta_x = tweak_vals[1], 
        delta_y = tweak_vals[2],
        delta_rot = tweak_vals[3]
      )
      
      message("Making new super reference frame. Super mosaic.")
      new_ref_frame = propaneStackWarpInVar(
        image_list = list(ref_frame, new_frame),
        cores = 1,
        keyvalues_out = wcs,
        magzero_in = 23.9, 
        magzero_out = 23.9
      )
      
      tweak_pars = c(tweak_pars, list(tweak_vals))
      
      ref_frame = new_ref_frame$image
      counter = counter + 1
    }
  }
  
  
  loc = which(!is.na(init_stack$image$imDat),arr.ind=TRUE)
  xlo = min(loc[,1])-100
  xhi = max(loc[,1])+100
  ylo = min(loc[,2])-100
  yhi = max(loc[,2])+100
  
  png("~/Desktop/foo2_align.png", width = dim(new_ref_frame$image)[1], height = dim(new_ref_frame$image)[2]/2.0, units="px", res = 100)
  par(mar = rep(0,4), oma = rep(0,4))
  plot(new_ref_frame$image[xlo:xhi, ylo:yhi], qdiff=T, sparse = 1)
  dev.off()
  
  return(tweak_pars)
}
profound_tweak_widget = function(load_frames, wcs){
  
  N_frames = length(load_frames)
  
  tweak_pars = {}
  counter = 0
  profound_tweak_pars = while(counter < N_frames){
    message("Counter = ", counter)
    if(counter == 0){
      message("Initialise reference frame.")
      ref_frame = propaneWarp(
        image_in = load_frames[[1]],
        keyvalues_out = wcs
      )
      tweak_pars = c(tweak_pars, list(c(0,0,0)))
      counter = counter + 1
      # store_frames = c(list(ref_frame), store_frames)
    }else{
      message("Tweak test frame.")
      message("ProFound ref.")
      pro_ref = profoundProFound(
        image = ref_frame,
        skycut = 5.0,
        pixcut = 20,
        cliptol = 100,
        tolerance = Inf,
        box = 100,
        magzero = 23.9,
        rem_mask = T
      )
      
      message("Warp test to target WCS")
      warp_frame = propaneWarp(
        load_frames[[counter+1]],
        keyvalues_out = wcs
      )
      
      message("ProFound test.")
      pro = profoundProFound(
        image = warp_frame,
        skycut = 5.0,
        pixcut = 20,
        cliptol = 100,
        tolerance = Inf,
        box = 100,
        magzero = 23.9,
        rem_mask = T
      )
      
      match_idx = coordmatch(
        coordref = pro_ref$segstats[, c("RAmax", "Decmax")],
        coordcompare = pro$segstats[, c("RAmax", "Decmax")],
        rad = 0.5
      )
      
      match_ref = pro_ref$segstats[match_idx$bestmatch$refID, ]
      match_test = pro$segstats[match_idx$bestmatch$compareID, ]
      
      dx = -median(match_test$xcen - match_ref$xcen, na.rm = T)
      dy = -median(match_test$ycen - match_ref$ycen, na.rm = T)
      dphi = 0 ## Ignore rotations.
      
      new_frame_imdat = propaneTran(
        image = warp_frame$imDat, 
        delta_x = dx, 
        delta_y = dy,
        delta_rot = dphi
      )
      new_frame = Rfits_create_image(
        image = new_frame_imdat,
        keyvalues = wcs
      )
      
      new_ref_frame = propaneStackWarpInVar(
        image_list = list(ref_frame, new_frame),
        keyvalues_out = wcs,
        magzero_in = 23.9, 
        magzero_out = 23.9
      )
      
      tweak_pars = c(tweak_pars, list(c(dx,dy,dphi)))
      
      ref_frame = new_ref_frame$image
      counter = counter + 1
    }
  }
  return(tweak_pars)
}



