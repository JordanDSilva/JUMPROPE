##Query data from MAST and download a .sh script
## Execute shell script to download data

from astroquery.mast import Observations
import os
import numpy as np
import pandas as pd
import argparse
from astropy.io import ascii
import glob
import sys
import shutil
import astropy.units as u
import astropy.table
from astropy.coordinates import SkyCoord
from astropy.table import vstack
import time

parser = argparse.ArgumentParser(description='Download from MAST')

## Set up input args to run from command line
parser.add_argument('--RA', type=float, nargs="?",
                    help='Central right ascension of catalogue (deg)', default=110.8375)
parser.add_argument('--DEC', type=float, nargs="?", help="Central declination of catalogue (deg)", default=-73.4391)
parser.add_argument('--RAD', type=float, nargs="?", help="Search radius (deg)", default=1.0)

parser.add_argument('--TELESCOPE', type=str, nargs="*", help="Telescope to download")
parser.add_argument('--CAMERA', type=str, nargs="*", help="Camera to download", default = ["ALL"])
parser.add_argument('--FILTER', type=str, nargs="*", help="Filter to download", default = ["ALL"])

parser.add_argument('--STAGE', type=str, nargs="*", default = False,
                    help='Product level (e.g., UN/CAL for un/cal files)')

parser.add_argument('--dry_run', type=bool, nargs=1,
                    help='Only download info table.')


args = parser.parse_args()

def query(ref_dir, ra, dec, rad, telescope, camera, filter, stage, dl_products):
    print("Querying observation")

    ## Set up the search cone
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")

    ## convenience stub name for saving the astropy tables
    stub_name = str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" + str(round(rad, 2)) + \
                "_" + telescope + "_" + camera + "_" + filter + "_" + stage

    ## where we will put the downloaded frames
    dl_dir = os.path.join(ref_dir, telescope,
                          str(round(ra, 2)) + "_" + str(round(dec, 2)) + "_" +
                          str(round(rad, 2)),
                          stage) + str("/")

    ## Make sure that we get all HST frames that are part of HAP or HLA products
    telescopes = [telescope]
    if "HST" in telescope:
        telescopes.append("HAP")
        telescopes.append("HLA")

    ## Now the meat of the code. Query MAST.
    if camera == filter == "ALL":
        obs_table = Observations.query_criteria(coordinates=coord,
                                                radius=rad * u.degree,
                                                project=telescopes,
                                                dataproduct_type=["IMAGE", "image"]
                                                )
    elif filter == "ALL" and camera != "ALL":
        obs_table = Observations.query_criteria(coordinates=coord,
                                        radius=rad * u.degree,
                                        project=telescopes,
                                        dataproduct_type=["IMAGE", "image"],
                                        instrument_name=[camera, camera + "/IMAGE"])
    elif camera == "ALL" and filter != "ALL":
        obs_table = Observations.query_criteria(coordinates=coord,
                                                radius=rad * u.degree,
                                                project=telescopes,
                                                dataproduct_type=["IMAGE", "image"],
                                                filters = filter)
    else:
        obs_table = Observations.query_criteria(coordinates=coord,
                                                radius=rad * u.degree,
                                                project=telescopes,
                                                dataproduct_type=["IMAGE", "image"],
                                                instrument_name=[camera, camera + "/IMAGE"],
                                                filters = filter)


    if len(obs_table) == 0:
        pass
    else:
        # Observations.enable_cloud_dataset(provider='AWS')
        #data_products = Observations.get_product_list(obs_table)
        product_list = [Observations.get_product_list(obs) for obs in obs_table]
        data_products = vstack(product_list)

        products_stage2 = Observations.filter_products(data_products,
                                                       extension="fits",
                                                       productSubGroupDescription=stage
                                                       )

        if "HAP-SVM" in set(products_stage2["project"]):
            products_stage2 = products_stage2[products_stage2["project"] == "HAP-SVM"]

        if len(products_stage2) == 0:
            pass
        else:

            os.makedirs(dl_dir, exist_ok=True)

            sh_files = glob.glob(dl_dir + "*.sh")
            for file in sh_files:
                os.remove(file)

            # Observations.download_products(products_stage2, download_dir=dl_dir, curl_flag=True)

            sh_files = glob.glob(dl_dir + "*.csv")
            for file in sh_files:
                os.remove(file)

            ascii.write(obs_table, dl_dir + stub_name +
                        '_obs_table.csv', overwrite=True, format='csv')
            ascii.write(products_stage2, dl_dir + stub_name +
                        '_info.csv', overwrite=True, format='csv')

            if dl_products:
                Observations.download_products(products_stage2,
                                              download_dir=dl_dir,
                                              curl_flag=False,
                                              cache=True)
            return obs_table

def main(args_dict):

    JUMPROPE_MAST_TOKEN = os.getenv('JUMPROPE_MAST_TOKEN')
    JUMPROPE_DOWNLOAD_DIR = os.getenv('JUMPROPE_DOWNLOAD_DIR')

    if None in [JUMPROPE_MAST_TOKEN, JUMPROPE_DOWNLOAD_DIR]:
        print("Please set ENV variables. ")
        exit()

    my_session = Observations.login(token=str(JUMPROPE_MAST_TOKEN))
    ref_dir = str(JUMPROPE_DOWNLOAD_DIR)
    args_grid = [(a,b,c,d) for a in args_dict["TELESCOPE"]
                           for b in args_dict["CAMERA"]
                           for c in args_dict["FILTER"]
                           for d in args_dict["STAGE"]]

    for i in range(len(args_grid)):
        query(
            ref_dir, args_dict["RA"], args_dict["DEC"], args_dict["RAD"],
            args_grid[i][0], args_grid[i][1], args_grid[i][2], args_grid[i][3],
            dl_products= not args_dict["dry_run"]
        )

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:

        args_dict = {
            "RA" : args.RA,
            "DEC" : args.DEC,
            "RAD" : args.RAD,
            "TELESCOPE" : args.TELESCOPE,
            "CAMERA" : args.CAMERA,
            "FILTER" : args.FILTER,
            "STAGE" : args.STAGE,
            "dry_run" : args.dry_run
        }

        main(args_dict)
