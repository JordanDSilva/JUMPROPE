export JUMPROPE_MAST_TOKEN=9999 #Only needed for exclusive access data. 
export JUMPROPE_DOWNLOAD_DIR='/Volumes/Expansion/JWST'
export JUMPROPE_CRDS_PATH='/Volumes/Expansion/crds'
export JUMPROPE_CRDS_CONTEXT='jwst_1140.pmap'


#export JUMPROPE_REF_DIR='/Volumes/Expansion/JWST/prof' ## I recommend you keep this commented. Only add when you are certain the directory structure has been made. The code will auto prompt you to create it on first run.
export JUMPROPE_RAW_DIR='/Volumes/Expansion/JWST/'

export JUMPROPE_do_NIRISS=F
export JUMPROPE_do_MIRI=F

export JUMPROPE_cores_pro=1 #Number of cores to use for 1/f, sky subtr, wisp rem
export JUMPROPE_cores_stack=1 #Number of cores to use internally for propaneWarpStackInVar
export JUMPROPE_tasks_stack=1 #Number of stacking tasks to run


# e.g., If I have 20 cal files to process, I could set cores_pro = 5 and only do 4 lots of computations.
# I could then set cores_stack = 10 and tasks_stack = 2 to stack two things at once while taking advantage of 10 cores for the WCS warping. 
# cores_stack * tasks_stack < Number of Cores. Use detectCores() in R or other to check system resources. 
