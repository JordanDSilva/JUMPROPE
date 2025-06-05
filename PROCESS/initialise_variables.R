###########################
## Initialise parameters ##
###########################

initialise_params = function(){
 
  ## Edit this:

  ## set up the keep_trend data
  additional_params = list(
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
   ## I.e., perform "wisp removal" algorithm on all NIRCam short wavelength detectors
   ## and not only wisp affected [A3,A4,B3,B4]
   do_claws = FALSE, 
   
   ## Path to reference astrometric catalogue - to make the ProPane stacks
   tweak_catalogue = NULL, ## NULL will have no internal tweak

   NAXIS_long = NULL, ## Size of the long pixel scales mosaic, keep NULL for default (3000 pixels) 
   NAXIS_short = NULL, ## Size of the short pixel scales mosaic, keep NULL for default (6000 pixels) 
   module_list = NULL ## What modules should we stack, options are ('NRCA_short', 'NRCA_long', 'NRCB_short', 'NRCB_long', 'NIS', 'MIRIMAGE')
  )
  ## Finish editing

  return( mget(ls()) )
}
