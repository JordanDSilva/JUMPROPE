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
ref_cat = fread("~/Documents/RefCats/a2744_astro.csv") #example for HST catalogue in EGS field

names(ref_cat) = c("RA", "Dec")

input_args = list(
  ref_dir = "/Volumes/Expansion/totJP/", ## Directory containing the ProPane/ProFound/JUMPROPE stuff
  super_name = "A2744", ## Name of mosaic. Advise against including "_"  

  ref_cat = ref_cat, ## Reference catalogue for alignment
  
  RA = colMeans(ref_cat)[1], ## Coords to search for frames
  Dec = colMeans(ref_cat)[2],     
  mosaic_size = NULL, ## Set the size of the search radius/mosaic size
  grid_size = "short", #default save on memory, 0.06arcsec/pix,
  rotation = "North", ## Rotation argument for ProPaneGenWCS, e.g., "North" for North-Up alignment

  program_id = "2561",  # Maybe we we just want to stack a single program e.g., Primer (1837) and not COSMOS Web (1727)
  cores = 2 ## If NULL then use half the number of cores in the system
) 


deep_stacker = function(input_args){
  
  ## Find aligning frames in search radius, align with ProPaneTweakCat and stack with ProPane
  
  inVar_dir = paste0(input_args$ref_dir, "/InVar_Stacks")
  median_dir = paste0(input_args$ref_dir, "/Median_Stacks")
  patch_dir = paste0(input_args$ref_dir, "/Patch_Stacks")
  dump_dir = paste0(input_args$ref_dir, "/dump")
  
  
  mosaic_dir = paste0(input_args$ref_dir, "/Mosaic_Stacks/")
  mosaic_invar = paste0(mosaic_dir, "/inVar/")
  mosaic_med = paste0(mosaic_dir, "/Median/")
  mosaic_patch = paste0(mosaic_dir, "/Patch/")
  
  dir.create(mosaic_dir, recursive = T)
  dir.create(mosaic_invar, recursive = T)
  dir.create(mosaic_med, recursive = T)
  dir.create(mosaic_patch, recursive = T)
  
  #cal_sky_info = fread(
  #  paste0(
  #    input_args$ref_dir, "/Pro1oF/cal_sky_info.csv"
  #  )
  #)
  
  #VIDS = paste0(cal_sky_info$VISIT_ID)
  
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
    cores = 8,
    plot=F
  )

  
  filters_w2 = str_extract(find_frames$file, pattern = regex("F(\\d+).W2")) ## string match F and W with only digits in between 
  filters_all = str_extract(find_frames$file, pattern = regex("F(\\d+).W|F(\\d+).M"))
  filters_all[which(filters_w2 == "F150W2")] = "F150W2"
  filters_all = na.omit(filters_all)
  
  find_frames$FILT = filters_all
  
  #UNIQUE_VIDS = VIDS[sapply(VIDS, function(x){any(grepl(x,find_frames$file))} )==TRUE]
  # 
#  stack_grid = cal_sky_info[paste0(cal_sky_info$VISIT_ID) %in% UNIQUE_VIDS, c("VISIT_ID", "FILTER", "DETECTOR")]
 
  stack_grid = data.frame("VISIT_ID" = c(rep("2561001002",7), rep("2561001003",7)), 
                          "FILTER" = c("F115W", "F150W", "F200W", "F277W", "F356W", "F410M", "F444W", "F115W", "F150W", "F200W", "F277W", "F356W", "F410M", "F444W"), 
                          "DETECTOR" = rep("NRCB", 14)) 
 
  for(grid_size in c(input_args$grid_size)){
    
    file_info = find_frames[grepl(grid_size, find_frames$file), ]
    
    wcs = propaneGenWCS(filelist = file_info$full, rotation = input_args$rotation)
    
    ref_cat = input_args$ref_cat
    
    UNIQUE_FILTERS = unique(stack_grid$FILTER)
    UNIQUE_FILTERS=UNIQUE_FILTERS[!grepl("CLEAR", UNIQUE_FILTERS) & UNIQUE_FILTERS %in% file_info$FILT]
    for(i in which(UNIQUE_FILTERS == "F115W")){
      
      file_info_filt = file_info[file_info$FILT == paste0(UNIQUE_FILTERS[i]),]
      filenames = file_info_filt$full
 
      print(filenames)
       
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
        pro_test = profoundProFound(
          image = image_list[[k]][,],
          skyRMS = 1/sqrt(inVar_list[[k]][,]$imDat),
          sky = 0,
          redosky = F,
          pixcut = 20,
          skycut = 10.0,
          cliptol = 100,
          tolerance = Inf,
          box = dim(image_list[[k]])/10.0,
          mask = is.infinite(image_list[[k]][,]$imDat),
          rem_mask = T
        )
        
        match_cat = coordmatch(
          coordref = ref_cat,
          coordcompare = pro_test$segstats[, c("RAmax", "Decmax")]
        )
        
        if(any(is.na(match_cat$bestmatch))){
          tweak_cat = list(par = c(0,0,0))
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
      
      cores = ifelse(is.null(input_args$cores), floor(detectCores()/2), input_args$cores)
      deep_frame = propaneStackWarpInVar(
        image_list = image_list,
        inVar_list = inVar_list,
        weight_list = weight_list,
        keyvalues_out = wcs,
        magzero_in = 23.9,
        magzero_out = 23.9,
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


## DEPRECATED
test_wcs = function(input_args){
  
  ## Find optimum search radius to make sure all frames fit inside the target WCS
  patch_dir = paste0(input_args$ref_dir, "/Patch_Stacks/")
  
  mosaic_dir = paste0(input_args$ref_dir, "/mosaic_stacks/")   
  mosaic_invar = paste0(mosaic_dir, "/inVar/")
  mosaic_med = paste0(mosaic_dir, "/Median/")
  mosaic_patch = paste0(mosaic_dir, "/Patch/")
  
  dir.create(mosaic_dir, recursive = T, showWarnings = F)
  dir.create(mosaic_invar, recursive = T, showWarnings = F)   
  dir.create(mosaic_med, recursive = T, showWarnings = F)
  dir.create(mosaic_patch, recursive = T, showWarnings = F)
  
  message("Searching for frames in ", input_args$RA, " ", input_args$Dec)
  find_frames_raw = propaneFrameFinder(
    dirlist = patch_dir,
    RAcen = input_args$RA,
    Deccen = input_args$Dec,
    rad = 2.0,
    cores = 8,
    plot=F
  )
  
  grid_size = input_args$grid_size
  find_frames = find_frames_raw[grepl(grid_size, find_frames_raw$full) & grepl("NRC", find_frames_raw$full), ]
  
  wcs_info = Rfits_key_scan(
    filelist = find_frames$full,
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
    cores = 1
  )
  
  pixscale_deg = mean(find_frames$pixscale)/3600
  
  if(input_args$at_these_coords){
    RA_val = input_args$RA
    Dec_val = input_args$Dec
  }else{
    RA_val = mean(wcs_info$centre_RA)
    Dec_val = mean(wcs_info$centre_Dec)
  }
  
  wcs = Rwcs_keypass(
    CRVAL1 = RA_val,
    CRVAL2 = Dec_val,
    CD1_1 = pixscale_deg,
    CD1_2 = 0,
    CD2_1 = 0,
    CD2_2 = pixscale_deg
  )
  
  foobar = function(mosaic_size_deg){
    message("...")
    foreach(i = seq_along(mosaic_size_deg), .combine = "c")%do%{
      wcs$NAXIS1 = ceiling(mosaic_size_deg[i] / pixscale_deg)
      wcs$NAXIS2 = ceiling(mosaic_size_deg[i] / pixscale_deg)
      wcs$CRPIX1 = wcs$NAXIS1/2.0
      wcs$CRPIX2 = wcs$NAXIS2/2.0
      
      wcs_temp = wcs
      wcs_temp$NAXIS1 = 10.0
      wcs_temp$NAXIS2 = 10.0
      wcs_temp$CRPIX1 = 5.0
      wcs_temp$CRPIX2 = 5.0
      
      cen_coords = Rwcs_s2p(
        RA = RA_val,
        Dec = Dec_val,
        keyvalues = wcs
      )
      
      cen_frames = Rwcs_s2p(
        RA = find_frames$centre_RA[grepl(grid_size, find_frames$stub)], 
        Dec = find_frames$centre_Dec[grepl(grid_size, find_frames$stub)], 
        keyvalues = wcs
      )
      
      boundary_frames = Rwcs_s2p(
        RA = unlist(find_frames[, grep("RA", names(find_frames), value = T)]), 
        Dec = unlist(find_frames[, grep("Dec", names(find_frames), value = T)]), 
        keyvalues = wcs
      )
      
      xcen = as.integer(round(cen_frames[,1],0))
      ycen = as.integer(round(cen_frames[,2],0))
      
      box = as.integer(round(find_frames$dim_1[grepl(grid_size, find_frames$stub)]/2.0,0))
      
      left = (xcen - box)
      right = (xcen + box)
      bottom = (ycen - box)
      up = (ycen + box)
      
      mat = matrix(data = 0, nrow = wcs$NAXIS1, ncol = wcs$NAXIS2)
      
      if(any(
        (left <= 0) | (right >= dim(mat)[1]) | (bottom <= 0) | (up >= dim(mat)[2]) | (boundary_frames <= 0) | (boundary_frames >= dim(mat)[1])
      ) 
      ){
        return(Inf)
      }else{
        for(j in 1:dim(find_frames)[1]){
          mat[((xcen[j]-box[j]):(xcen[j]+box[j])),((ycen[j]-box[j]):(ycen[j]+box[j]))] = 1
        }
        
        fill_fraction = sum(mat==1)/prod(dim(mat))
        
        if((fill_fraction >= 1)){
          return(Inf)
        }else{
          fill_fraction
        }
      }
    }
  }
  
  message("Finding optimum mosaic size")
  suppressWarnings({
    foo_optim <- optim(
      par = c(0.17), fn = foobar, method = "Brent", control = list(fnscale = -1), lower = 0.1, upper = 1.0
    )
  })
  
  message("Rendering footprint plot")
  mosaic_size_deg = foo_optim$par
  # mosaic_size_deg = 0.17
  wcs$NAXIS1 = ceiling(mosaic_size_deg / pixscale_deg)
  wcs$NAXIS2 = ceiling(mosaic_size_deg / pixscale_deg)
  wcs$CRPIX1 = wcs$NAXIS1/2.0
  wcs$CRPIX2 = wcs$NAXIS2/2.0
  
  wcs_temp = wcs
  wcs_temp$NAXIS1 = 10.0
  wcs_temp$NAXIS2 = 10.0
  wcs_temp$CRPIX1 = 5.0
  wcs_temp$CRPIX2 = 5.0
  
  png(paste0(mosaic_dir, "/", input_args$super_name, "_layout.png"), width = 10, height = 8, units="in", res=72) 
  Rwcs_overlap(
    keyvalues_test = wcs_temp, keyvalues_ref = wcs, plot = T
  )
  legend(x="topright", legend = foo_optim$par)
  
  cen_coords = Rwcs_s2p(
    RA = RA_val,
    Dec = Dec_val,
    keyvalues = wcs
  )
  BL = Rwcs_s2p(
    RA = find_frames$corner_BL_RA, Dec = find_frames$corner_BL_Dec, keyvalues = wcs
  )
  BR = Rwcs_s2p(
    RA = find_frames$corner_BR_RA, Dec = find_frames$corner_BR_Dec, keyvalues = wcs
  )
  TL = Rwcs_s2p(
    RA = find_frames$corner_TL_RA, Dec = find_frames$corner_TL_Dec, keyvalues = wcs
  )
  TR = Rwcs_s2p(
    RA = find_frames$corner_TR_RA, Dec = find_frames$corner_TR_Dec, keyvalues = wcs
  )
  
  BL_dist = (BL[,1] - cen_coords[,1])^2 + (BL[,2] - cen_coords[,2])^2
  BL_dist = (BL[,1] - cen_coords[,1])^2 + (BL[,2] - cen_coords[,2])^2
  BL_dist = (BL[,1] - cen_coords[,1])^2 + (BL[,2] - cen_coords[,2])^2
  BL_dist = (BL[,1] - cen_coords[,1])^2 + (BL[,2] - cen_coords[,2])^2
  
  for(i in seq_along(find_frames$full)){
    polygon(c(BL[i,1], TL[i,1], TR[i,1], BR[i,1]),
            c(BL[i,2], TL[i,2], TR[i,2], BR[i,2]),
            # as.numeric(find_frames[i,c('corner_BL_RA','corner_TL_RA','corner_TR_RA','corner_BR_RA')]),
            # as.numeric(find_frames[i,c('corner_BL_Dec','corner_TL_Dec','corner_TR_Dec','corner_BR_Dec')]),
            border = "red",
            col = "red"
    )
  }
  dev.off()
  message("Optimum mosaic size: ", mosaic_size_deg, " deg")
  return(list(radius=mosaic_size_deg, wcs = wcs))
}
