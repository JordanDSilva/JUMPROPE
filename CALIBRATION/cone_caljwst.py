#####################################################################
# JWST Calibration Pipeline, Takes us from Uncal files to Cal files.#
# Need to have JWST calibration pipeline installed.                 #
#####################################################################

import numpy as np
import os
import pandas as pd

os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"  # where to download reference files from
from jwst.pipeline import calwebb_detector1  # import detector processor, uncal-->rates
from jwst.pipeline import calwebb_image2  # import image processor, rates-->cal

import glob
import sys
import argparse  # to run from command line, run <python3 calwebb.py -h> to display help

parser = argparse.ArgumentParser()
parser.add_argument('--RA', type=float, nargs="?",
                    help='Central right ascension of catalogue (deg)', default=110.8375)
parser.add_argument('--DEC', type=float, nargs="?", help="Central declination of catalogue (deg)", default=-73.4391)
parser.add_argument('--RAD', type=float, nargs="?", help="Search radius (deg)", default=1.0)

parser.add_argument("--redo",
                    type=bool,
                    help="Should we redo the calibration even if files exist?",
                    default=False)
parser.add_argument("--max_cores", help="How many cores. Must be written as e.g., 'half', 'quarter' or 'all'", default="half")
args = parser.parse_args()


def run1(input_data, rates_dir, cal_dir):
    """Process every uncal to cal"""

    files = input_data
    detector1 = calwebb_detector1.Detector1Pipeline()

    detector1.output_dir = rates_dir
    # detector1.ipc.skip = False

    # detector1.clean_flicker_noise.fit_method = 'fft'

    # detector1.jump.rejection_threshold = 4.0
    detector1.jump.expand_large_events = True

    # detector1.sat_required_snowball = False
    # detector1.max_jump_to_flag_neighbors = 1
    # detector1.min_jump_to_flag_neighbors = 100000

    detector1.jump.maximum_cores = args.max_cores
    detector1.ramp_fit.maximum_cores = args.max_cores

    run_output = detector1.run(files)

    image2 = calwebb_image2.Image2Pipeline()
    image2.output_dir = cal_dir
    image2.save_results = True
    image2.resample.skip = True  # we don't need to resample the cals, since they get drizzled later
    image2.run(run_output)

def cal(ref_dir, ra, dec, rad):
    print("CALJWST in region")
    uncal_dir = os.path.join(ref_dir, "JWST",
                             str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                             str(round(rad, 2)),
                             "UNCAL") + str("/")

    files = glob.glob(uncal_dir + "/**/*fits", recursive = True)
    uncal_files = [foo.split("_uncal")[0].split("/")[-1] for foo in files]

    # dir_frames = "/Volumes/Expansion/JWST/"
    rates_dir = os.path.join(ref_dir, "JWST",
                             str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                             str(round(rad, 2)),
                             "RATES") + str("/")
    rates_files = glob.glob(rates_dir + "/**/*fits", recursive = True)
    rates_files_names = [foo.split("_cal")[0].split("/")[-1] for foo in rates_files]

    cal_dir = os.path.join(ref_dir, "JWST",
                           str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                           str(round(rad, 2)),
                           "CAL") + str("/")  # cal directory
    cal_files = glob.glob(cal_dir + "/**/*fits", recursive = True)
    # cal_files = [s for s in cal_files if visitid in s]
    cal_files_names = [foo.split("_cal")[0].split("/")[-1] for foo in cal_files]

    # # Make the cal and rates directories
    os.makedirs(cal_dir, exist_ok=True)
    os.makedirs(rates_dir, exist_ok=True)

    ## Read info tables from cone search
    info_tables = glob.glob(
        uncal_dir + "*info.csv"
    )
    df = pd.concat((pd.read_csv(f) for f in info_tables), ignore_index=True)

    if (len(df) > len(files)):
        print("Some files appear to be missing in this region \n compared to what was downloaded with cone-search.")

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
        if sum(df["obs_id"].isin(cal_files_names)) != len(df["obs_id"]):  ## check that all of the expected cal files have been produced

            idx = list(~df["obs_id"].isin(cal_files_names))

            for i, j in enumerate(idx):
                if j:
                    files_redo.append(files[i])

    files_redo = files
    print("Running " + str(len(files_redo)) + " files")
    for ff in files_redo:
        run1(ff, rates_dir, cal_dir)


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
        cal(JUMPROPE_DOWNLOAD_DIR, args.RA, args.DEC, args.RAD)