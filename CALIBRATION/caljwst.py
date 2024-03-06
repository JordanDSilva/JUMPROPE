###################################################################################################
#JWST Calibration Pipeline, Takes us from Uncal files to Cal files                                #
#Need to have JWST calibration pipeline installed. Activate python environment on Simon's machine.#
###################################################################################################

import numpy as np
import os
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu" #where to download reference files from

from jwst.pipeline import calwebb_detector1 #import detector processor, uncal-->rates
from jwst.pipeline import calwebb_image2 #import image processor, rates-->cal
from astroquery.mast import Observations
from astropy.io import ascii

import glob
import sys
import argparse #to run from command line, run <python3 calwebb.py -h> to display help

parser = argparse.ArgumentParser()
parser.add_argument("--VISITID",
                    help = "VISITID. e.g., `2736001001` for SMACS")
parser.add_argument("--redo",
                    type = bool,
                    help = "Should we redo the calibration even if files exist?",
                    default = False)
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

    JUMPROPE_DOWNLOAD_DIR = os.getenv('JUMPROPE_DOWNLOAD_DIR')

    JUMPROPE_CRDS_PATH = os.getenv('JUMPROPE_CRDS_PATH')
    JUMPROPE_CRDS_CONTEXT = os.getenv('JUMPROPE_CRDS_CONTEXT')

    os.environ['CRDS_PATH'] = JUMPROPE_CRDS_PATH
    os.environ['CRDS_CONTEXT'] = JUMPROPE_CRDS_CONTEXT

    if None in [JUMPROPE_DOWNLOAD_DIR, JUMPROPE_CRDS_PATH, JUMPROPE_CRDS_CONTEXT]:
        print("Please set ENV variables. ")
        exit()

    N = 6 #N cores/processes

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    elif args.VISITID is not None:
        visitid = args.VISITID
        PID = str(visitid)[0:4]
        ref_dir = JUMPROPE_DOWNLOAD_DIR
        uncal_dir = os.path.join(ref_dir, PID, 'UNCAL/NIRCAM') # folder must  be present in input directory
        files = glob.glob(uncal_dir + "/*fits")
        files = [s for s in files if visitid in s]
        uncal_files = [foo.split("_uncal")[0].split("/")[-1] for foo in files]

        #dir_frames = "/Volumes/Expansion/JWST/"
        rates_dir = os.path.join(ref_dir, PID, 'RATES/NIRCAM') # rates directory
        rates_files = glob.glob(rates_dir + "/*fits")
        rates_files = [s for s in rates_files if visitid in s]
        rates_files_names = [foo.split("_cal")[0].split("/")[-1] for foo in rates_files]

        cal_dir = os.path.join(ref_dir, PID, 'CAL/NIRCAM') # cal directory
        cal_files = glob.glob(cal_dir + "/*fits")
        cal_files = [s for s in cal_files if visitid in s]
        cal_files_names = [foo.split("_cal")[0].split("/")[-1] for foo in cal_files]

        # # Make the cal and rates directories
        os.makedirs(cal_dir, exist_ok=True)
        os.makedirs(rates_dir, exist_ok=True)

        ## download csv from mast to check expected files and their sizes
        obs_table = Observations.query_criteria(proposal_id=PID,
                                                project='JWST',
                                                dataproduct_type="IMAGE",
                                                instrument_name=["NIRCAM", "NIRCAM/IMAGE"])
        data_products = Observations.get_product_list(obs_table)
        products_stage2 = Observations.filter_products(data_products,
                                                       extension="fits",
                                                       productSubGroupDescription="CAL",
                                                       )
        visit_ids = np.array([obs[3:13] for obs in products_stage2["obs_id"]])
        idx = np.flatnonzero(np.core.defchararray.find(visit_ids, str(visitid)) != -1)
        df = products_stage2[idx].to_pandas().drop_duplicates("productFilename")
        # df = dff[dff["dataRights"] == "PUBLIC"]
        if(len(df) > len(files)):
            print("Some files appear to be missing in VISIT PROGRAM " + visitid + "\n compared to what I can find on MAST.")

            if sum(df["dataRights"] == "EXCLUSIVE_ACCESS") == len(df) - len(uncal_files):
                print("Exclusive access files found.")
                df = df[df["dataRights"] == "PUBLIC"]
            else:
                df = df[df["obs_id"].isin(uncal_files)]

        if args.redo:  ## override previous checks. In case we want to use a new pmap of calibration version for example
            for cal, rate in zip(cal_files, rates_files):
                print("Removing " + str(cal))
                os.remove(cal)
                print("Removing " + str(rate))
                os.remove(rate)
            files_redo = files
        else:
            files_redo = []
            if sum(df["obs_id"].isin(cal_files_names)) != len(df["obs_id"]): ## check that all of the expected cal files have been produced

                idx = list(~df["obs_id"].isin(cal_files_names))

                for i,j in enumerate(idx):
                    if j:
                        files_redo.append(files[i])

        print("Running " + str(len(files_redo)) + " files")
        for ff in files_redo:
          run1(ff)

#I’ve processed all the NIRCam smacs frames with cal-webb. The issue I was having was that two of the dark frames must not have downloaded #properly, so I redownloaded them manually. Process took about ~1 hour with 10 frames running in parallel over 270 frames. These are all with pmap 0953.
#Data summary:
#crds (reference files) = 64.13GB
#uncals = 20.07 GB
#rates = 58.96GB
#cals = 31.74 GB
#There is an option to not save the rate files, but I think Jake does the wisp correction at the rate #level. He also basically runs with #default cal webb parameters.
#I also ran 2 frames of the IDF with identical settings and they worked fine => it doesn’t #automatically redownload reference files. So it #must be that the darks+flats+etc are pmap dependent.
