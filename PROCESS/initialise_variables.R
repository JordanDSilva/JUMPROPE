###########################
## Initialise parameters ##
###########################

initialise_params = function(){
 
  ## Edit this:
  #ref_dir = "~/Desktop/pro/" #Commented out as the scripts will prompt you to create the correct directory structure
  
  # dir_raw =  "Imaging" # Important thing to edit/check. Tells us where the input CAL files are stored. Usually the result of stage 2 of the JWST Calibration Pipeline.
  # 
  # do_NIRISS = F
  # cores_pro = 6
  # cores_stack = 4
  # tasks_stack = 1
  # 
  ## set up the keep_trend data
  keep_trend_data = list(
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
    )
    
  )
  ## Finish editing
  
  return( mget(ls()) )
}
