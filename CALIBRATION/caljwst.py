###################################################################################################
#JWST Calibration Pipeline, Takes us from Uncal files to Cal files                                #
#Need to have JWST calibration pipeline installed. Activate python environment on Simon's machine.#
###################################################################################################

import numpy as np
import os
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu" #where to download reference files from

from jwst.pipeline import calwebb_detector1 #import detector processor, uncal-->rates
from jwst.pipeline import calwebb_image2 #import image processor, rates-->cal

from pathos import multiprocessing as mp #multiprocessing utility
from pathos.pools import ThreadPool as Pool
import glob
import sys
import argparse #to run from command line, run <python3 calwebb.py -h> to display help
parser = argparse.ArgumentParser()
parser.add_argument("--VISITID",
                    help = "VISITID. e.g., `2736001001` for SMACS")

parser.add_argument("--max_cores", help = "How many cores e.g., 'half'", default = "half")

args = parser.parse_args()

def run1(input_data):
    """Process every uncal to cal"""

    files = input_data
    detector1 = calwebb_detector1.Detector1Pipeline()

    detector1.output_dir = rates_dir
    detector1.ipc.skip = False
    
    detector1.jump.rejection_threshold = 4.0
    detector1.jump.expand_large_events = True

    detector1.sat_required_snowball = False
    detector1.max_jump_to_flag_neighbors = 1
    detector1.min_jump_to_flag_neighbors = 100000

    detector1.jump.maximum_cores = args.max_cores
    detector1.ramp_fit.maximum_cores = args.max_cores

    run_output = detector1.run(files)

    image2 = calwebb_image2.Image2Pipeline()
    image2.output_dir = cal_dir
    image2.save_results = True
    image2.resample.skip = True #we don't need to resample the cals, since they get drizzled later
    image2.run(run_output)

if __name__ == "__main__":

    JWST_QUERY_TOOLS_REF_DIR = os.getenv('JWST_QUERY_TOOLS_REF_DIR')
    CRDS_PATH = os.getenv('JWST_QUERY_TOOLS_CRDS_PATH')
    CRDS_CONTEXT = os.getenv('JWST_QUERY_TOOLS_CRDS_CONTEXT')

    os.environ['CRDS_PATH'] = CRDS_PATH
    os.environ['CRDS_CONTEXT'] = CRDS_CONTEXT

    if None in [JWST_QUERY_TOOLS_REF_DIR, CRDS_PATH, CRDS_CONTEXT]:
        print("Please set ENV variables. ")
        exit()

    print(args)
    N = 6 #N cores/processes 

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    if args.VISITID is not None:
        visitid = args.VISITID
        PID = str(visitid)[0:4]
        ref_dir = JWST_QUERY_TOOLS_REF_DIR
        uncal_dir = os.path.join(ref_dir, PID, 'UNCAL/NIRCAM') # folder must  be present in input directory
        print(uncal_dir)
        files = glob.glob(uncal_dir + "/*fits")

        files = [s for s in files if visitid in s]
        # files = os.listdir(uncal_dir) #list all of the uncals.fits
   
        #dir_frames = "/Volumes/Expansion/JWST/"
        cal_dir = os.path.join(ref_dir, PID, 'CAL/NIRCAM') # cal directory
        rates_dir = os.path.join(ref_dir, PID, 'RATES/NIRCAM') # rates directory

        rates_files = glob.glob(rates_dir + "/*fits")
        cal_files = glob.glob(cal_dir + "/*fits")

        uncal_files = [foo.split("_uncal")[0].split("/")[-1] for foo in files]

        files_redo = []
        for i in range(len(files)):
            idx = [r for r, s in enumerate(cal_files) if uncal_files[i] in s]
            if len(idx)>1:
                print("Something went wrong \n")
                print(cal_files[idx])
                break
            if len(idx)==1:
                cal_size = os.stat(cal_files[idx[0]]).st_size
                if cal_size != 117538560: #size of cal file on disk, hopefully shouldn't change :O
                    files_redo.append(files[i])
            if len(idx)==0:
                files_redo.append(files[i])

        print(files_redo, sep="\n")
        # # Make the cal and rates directories
        os.makedirs(cal_dir, exist_ok=True)
        os.makedirs(rates_dir, exist_ok=True)

        for ff in files_redo:
          run1(ff)

        # #print(cal_files)
        # # Process files!
        #with Pool(processes=N) as p:
        #   p.map(run1, files_redo)
        # files_redo = [files[i] for i in range(len(files)) if not (any(uncal_files[i] in s for s in rates_files)) & (any(uncal_files[i] in s for s in cal_files))]




#I’ve processed all the NIRCam smacs frames with cal-webb. The issue I was having was that two of the dark frames must not have downloaded #properly, so I redownloaded them manually. Process took about ~1 hour with 10 frames running in parallel over 270 frames. These are all with pmap 0953.
#Data summary:
#crds (reference files) = 64.13GB
#uncals = 20.07 GB
#rates = 58.96GB
#cals = 31.74 GB
#There is an option to not save the rate files, but I think Jake does the wisp correction at the rate #level. He also basically runs with #default cal webb parameters.
#I also ran 2 frames of the IDF with identical settings and they worked fine => it doesn’t #automatically redownload reference files. So it #must be that the darks+flats+etc are pmap dependent.
