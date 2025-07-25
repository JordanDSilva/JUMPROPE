#####################################################################
# JWST Calibration Pipeline, Takes us from Uncal files to Cal files.#
# Need to have JWST calibration pipeline installed.                 #
#####################################################################

import numpy as np
import os
import pandas as pd
import glob
import sys
import argparse

os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"  # where to download reference files from
from jwst.pipeline import calwebb_detector1  # import detector processor, uncal-->rates
from jwst.pipeline import calwebb_image2  # import image processor, rates-->cal

parser = argparse.ArgumentParser()
parser.add_argument('--RA', type=float, nargs="?",
                    help='Central right ascension of catalogue (deg)', default=110.8375)
parser.add_argument('--DEC', type=float, nargs="?", help="Central declination of catalogue (deg)", default=-73.4391)
parser.add_argument('--RAD', type=float, nargs="?", help="Search radius (deg)", default=1.0)

parser.add_argument('--VISITID', type=str, help="VISIT ID (e.g., 2736 for SMACS)", default = None)

parser.add_argument("--redo",
                    type=bool,
                    help="Should we redo the calibration even if files exist?",
                    default=False)
parser.add_argument("--max_cores", help="How many cores. Must be written as e.g., 'half', 'quarter' or 'all'", default="half")
args = parser.parse_args()

def run1(input_data, rate_dir, cal_dir, files_bad_dir):
    """ Process every uncal to cal """

    # detector1.clean_flicker_noise.fit_method = 'fft'
    # detector1.jump.rejection_threshold = 4.0
    # detector1.jump.expand_large_events = True
    # detector1.sat_required_snowball = False
    # detector1.max_jump_to_flag_neighbors = 1
    # detector1.min_jump_to_flag_neighbors = 100000

    ## Make the cal and rates directories
    os.makedirs(rate_dir, exist_ok=True)
    os.makedirs(cal_dir, exist_ok=True)

    files = input_data
    detector1 = calwebb_detector1.Detector1Pipeline()

    ## EDIT THE CODE WITH THE JWST CALIBRATION OPTIONS HERE 
    detector1.output_dir = rate_dir
    detector1.jump.maximum_cores = args.max_cores
    detector1.ramp_fit.maximum_cores = args.max_cores
    # detector1.emicorr.skip = False ## Do the EMI banding correction for MIRI 

    try:
        run_output = detector1.run(files)
        image2 = calwebb_image2.Image2Pipeline()
        image2.output_dir = cal_dir
        image2.save_results = True
        image2.resample.skip = True  # we don't need to resample the cals, since they get drizzled later
        image2.run(run_output)
    except Exception as e:
        ## in case the uncal mucks up
        files_bad_stub = (files.split("UNCAL")[1].split("/"))[-1]
        os.rename(files, files_bad_dir + "/" + files_bad_stub)
        print(f"Bad file :( {files_bad_stub}")
        print(f"Error: {e}")
        return None

def cal(ref_dir, ra, dec, rad, visitid):

    if visitid is None:
        print("CALJWST in region")
        uncal_dir = os.path.join(ref_dir, "JWST",
                                 str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                                 str(round(rad, 2)),
                                 "UNCAL") + str("/")
        files_bad_dir = os.path.join(ref_dir, "JWST",
                                 str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                                 str(round(rad, 2)),
                                 "bad") + str("/")
    else:
        print("CALJWST in region")
        uncal_dir = os.path.join(ref_dir, str(visitid), "UNCAL") + str("/")
        files_bad_dir = os.path.join(ref_dir, str(visitid), "bad") + str("/")

    os.makedirs(files_bad_dir, exist_ok=True) ## in case uncal files are corrupted

    files = glob.glob(uncal_dir + "/**/*fits", recursive = True)
    uncal_files = [foo.split("_uncal")[0].split("/")[-1] for foo in files]


    rate_dir = []
    rate_files = []

    cal_dir = []
    cal_files = []
    for i in range(len(files)):
        getUncalPath = os.path.dirname(files[i])

        getRatePath = getUncalPath.replace("UNCAL", "RATE")
        rate_dir.append(getRatePath)
        getRateFile = os.path.join(getRatePath, uncal_files[i] + "_trapsfilled.fits")
        rate_files.append(getRateFile)

        getCalPath = getUncalPath.replace("UNCAL", "CAL")
        cal_dir.append(getCalPath)
        getCalFile = os.path.join(getCalPath, uncal_files[i] + "_cal.fits")
        cal_files.append(getCalFile)

    ## Read info tables from cone search
    info_tables = glob.glob(
        uncal_dir + "*info.csv"
    )
    print(uncal_dir)
    df = pd.concat((pd.read_csv(f) for f in info_tables), ignore_index=True)

    if (len(df) > len(files)):
        print("Some files appear to be missing in this region \n compared to what was downloaded with cone-search.")

    if args.redo:  ## override previous checks. In case we want to use a new pmap of calibration version for example
        for cal in cal_files:
            if os.path.isfile(cal):
                print("Removing " + str(cal))
                os.remove(cal)
        files_redo = files
        rate_dir_redo = rate_dir
        cal_dir_redo = cal_dir
    else:
        files_redo = []
        rate_dir_redo = []
        cal_dir_redo = []
        for i in range(len(files)):
            if not os.path.isfile(cal_files[i]):
                files_redo.append(files[i])
                rate_dir_redo.append(rate_dir[i])
                cal_dir_redo.append(cal_dir[i])

    print("Running " + str(len(files_redo)) + " files")

    if len(files_redo) > 0:
        for i, ff in enumerate(files_redo):
            run1(ff, rate_dir_redo[i], cal_dir_redo[i], files_bad_dir)
    else:
        exit()

if __name__ == "__main__":
    JUMPROPE_DOWNLOAD_DIR = os.getenv('JUMPROPE_DOWNLOAD_DIR')

    JUMPROPE_CRDS_PATH = os.getenv('JUMPROPE_CRDS_PATH')
    JUMPROPE_CRDS_CONTEXT = os.getenv('JUMPROPE_CRDS_CONTEXT')

    os.environ['CRDS_PATH'] = JUMPROPE_CRDS_PATH
    os.environ['CRDS_CONTEXT'] = JUMPROPE_CRDS_CONTEXT

    if None in [JUMPROPE_DOWNLOAD_DIR, JUMPROPE_CRDS_PATH, JUMPROPE_CRDS_CONTEXT]:
        print("Please set ENV variables. ")
        exit()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        cal(JUMPROPE_DOWNLOAD_DIR, args.RA, args.DEC, args.RAD, args.VISITID)
