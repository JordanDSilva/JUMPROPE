## Deep stacker + Alignment

# ref_cat = read.table("/Users/22252335/Downloads/hlsp_candels_hst_wfc3_egs_multi_v2_redshift-cat.txt")
# ref_cat = ref_cat[, 3:4]
# names(ref_cat) = c("RA", "Dec")

# farmer_cat = Rfits_read("/Users/22252335/Downloads/COSMOS2020_FARMER_R1_v2.2_p3.fits")
# 
# ref_cat = data.frame(
#   "RA" = farmer_cat[[2]]$ALPHA_J2000,
#   "Dec" = farmer_cat[[2]]$DELTA_J2000
# )

library(Rfits)
library(Rwcs)
library(ProPane)
library(ProFound)
library(celestial)
library(foreach)
library(doParallel)
library(data.table)
library(stringr)

# ref_cat = fread("/Volumes/RAIDY/JWST/test_stack/ref_cats/cosmos_ref_cat.csv")
# ref_cat = fread("/Volumes/RAIDY/JWST/test_stack/ref_cats/mosaic_plckg191_nircam.csv")

ref_cat = read.table("/Users/22252335/Downloads/hlsp_candels_hst_wfc3_egs_multi_v2_redshift-cat.txt")
ref_cat = data.frame(
  RA = ref_cat$V3,
  Dec = ref_cat$V4
)

input_args = list(
  ref_dir = "/Volumes/Expansion/JWST/prof/",
  RA = colMeans(ref_cat)[1], 
  Dec = colMeans(ref_cat)[2],
  super_name = "DDT",
  ref_cat = ref_cat
)

test_wcs = function(input_args){
  
  patch_dir = paste0(input_args$ref_dir, "/Patch_Stacks/")

  message("Searching for frames in ", input_args$RA, " ", input_args$Dec)
  find_frames = propaneFrameFinder(
    dirlist = patch_dir,
    RAcen = input_args$RA,
    Deccen = input_args$Dec,
    rad = 1.0,
    cores = 8,
    plot=F
  )
  
  grid_size = "long"

  wcs_info = Rfits_key_scan(
      filelist = find_frames$full[grepl(grid_size, find_frames$full)],
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
    
    pixscale_deg = mean(find_frames$pixscale)/3600
    
    RA_val = mean(wcs_info$centre_RA)
    Dec_val = mean(wcs_info$centre_Dec)

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
          RA = find_frames$centre_RA, Dec = find_frames$centre_Dec, keyvalues = wcs
        )
        
        xcen = as.integer(round(cen_frames[,1],0))
        ycen = as.integer(round(cen_frames[,2],0))
        
        box = as.integer(round(find_frames$dim_1/2.0,0))
        
        left = xcen-box
        right = xcen+box
        bottom = ycen-box
        up = ycen + box
        
        mat = matrix(data = 0, nrow = wcs$NAXIS1, ncol = wcs$NAXIS2)
        
        if(any(
          (left <= 0) | (right >= dim(mat)[1]) | (bottom <= 0) | (up >= dim(mat)[2]))
        ){
          return(Inf)
        }else{
          for(j in 1:dim(find_frames)[1]){
            mat[((xcen[j]-box[j]):(xcen[j]+box[j])),((ycen[j]-box[j]):(ycen[j]+box[j]))] = 1
          }
          
          fill_fraction = sum(mat==1)/prod(dim(mat))
          
          border_frac = mat
          border_frac[100:(dim(mat)[1]-100), 100:(dim(mat)[2]-100)] = 0
          
          if((fill_fraction >= 1) | any(border_frac==1)){
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
        par = c(0.15), fn = foobar, method = "Brent", control = list(fnscale = -1), lower = 0.1, upper = 0.5
      )
    })
   
    message("Rendering footprint plot")
    mosaic_size_deg = foo_optim$par
    wcs$NAXIS1 = ceiling(mosaic_size_deg / pixscale_deg)
    wcs$NAXIS2 = ceiling(mosaic_size_deg / pixscale_deg)
    wcs$CRPIX1 = wcs$NAXIS1/2.0
    wcs$CRPIX2 = wcs$NAXIS2/2.0
    
    wcs_temp = wcs
    wcs_temp$NAXIS1 = 10.0
    wcs_temp$NAXIS2 = 10.0
    wcs_temp$CRPIX1 = 5.0
    wcs_temp$CRPIX2 = 5.0
    
    par(mar = c(3.5, 3.5, 1.5, 1.5), oma = rep(0,4))
    Rwcs_overlap(
      keyvalues_test = wcs_temp, keyvalues_ref = wcs, plot = T
    )
    
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
    message("Optimum mosaic size: ", mosaic_size_deg, " deg")
    return(list(radius=mosaic_size_deg, wcs = wcs))
}
foo=test_wcs(input_args)

input_args$mosaic_size_deg = foo$radius

deep_stacker = function(input_args){
  
  inVar_dir = paste0(input_args$ref_dir, "/InVar_Stacks")
  median_dir = paste0(input_args$ref_dir, "/Median_Stacks")
  patch_dir = paste0(input_args$ref_dir, "/Patch_Stacks")
  dump_dir = paste0(input_args$ref_dir, "/dump")
  # reftweak_dir = paste0(input_args$ref_dir, "/reftweak")
  # test_dir = paste0(input_args$ref_dir, "/test_stack")
  # dump_dir = "/Users/22252335/Desktop/stack_test/dump"
  # test_dir = "/Users/22252335/Desktop/stack_test/test_stack"
  
  
  # frames_info = fread(
  #   paste0(
  #     input_args$ref_dir, "/FieldFootprints/frames_info.csv"
  #   )
  # )
  # 
  # cal_sky_info = fread(
  #   paste0(
  #     input_args$ref_dir, "/Pro1oF/cal_sky_info.csv"
  #   )
  # )
  
  # VIDS = paste0(frames_info$VISIT_ID)
  
  message("Searching for frames in ", input_args$RA, " ", input_args$Dec)
  find_frames = propaneFrameFinder(
    dirlist = patch_dir,
    RAcen = input_args$RA,
    Deccen = input_args$Dec,
    rad = input_args$mosaic_size_deg,
    cores = 8,
    plot=F
    )
  
  # UNIQUE_VIDS = VIDS[sapply(VIDS, function(x){any(grepl(x,find_frames$file))} )==TRUE]
  # 
  # stack_grid = cal_sky_info[paste0(cal_sky_info$VISIT_ID) %in% UNIQUE_VIDS, c("VISIT_ID", "FILTER", "DETECTOR")]
  
  
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
    
    pixscale_deg = mean(file_info$pixscale)/3600
    
    RA_val = mean(wcs_info$centre_RA)
    Dec_val = mean(wcs_info$centre_Dec)

    wcs = Rwcs_keypass(
      CRVAL1 = RA_val,
      CRVAL2 = Dec_val,
      CD1_1 = pixscale_deg,
      CD1_2 = 0,
      CD2_1 = 0,
      CD2_2 = pixscale_deg
    )
    
    wcs$NAXIS1 = ceiling(input_args$mosaic_size_deg / pixscale_deg)
    wcs$NAXIS2 = ceiling(input_args$mosaic_size_deg / pixscale_deg)
    wcs$CRPIX1 = wcs$NAXIS1/2.0
    wcs$CRPIX2 = wcs$NAXIS2/2.0

    ref_cat = input_args$ref_cat
    
    UNIQUE_FILTERS = c("F115W", "F150W", "F200W", "F277W", "F356W", "F444W")
    # idx_filt = which(UNIQUE_FILTERS == "F277W")
    for(i in 1:length(UNIQUE_FILTERS)){
      
      filenames = file_info$full[grepl(paste0(UNIQUE_FILTERS[i]), file_info$full)]
      
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
      
      for(k in 1:length(image_list)){
        
        message(
          paste0("Tweaking ", filenames[k])
        )
        pro_test = profoundProFound(
          image = image_list[[k]][,],
          skyRMS = 1/sqrt(inVar_list[[k]][,]$imDat),
          pixcut = 20,
          skycut = 10.0,
          cliptol = 100,
          tolerance = Inf,
          box = dim(image_list[[k]])/10.0,
          mask = is.infinite(image_list[[k]][,]$imDat)
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
                                      cat_pre_fix = pro_test$segstats[match_cat$bestmatch$compareID,c("xcen","ycen")],
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
        cores = 1,
        cores_warp = 1
      )
      
      deep_frame$image$imDat[is.infinite(deep_frame$image$imDat)] = NA
      
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
      
      rm(deep_frame)
      rm(med_stack)
      gc()
    }
    
  }

}

deep_stacker(input_args = input_args)


