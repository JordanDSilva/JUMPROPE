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

parser = argparse.ArgumentParser(description='Download from MAST')
parser.add_argument('--STAGE', type=str, nargs='+',
                    help='Product level (e.g., UN/CAL for un/cal files)')
parser.add_argument('--PID', type=str, nargs="+", help="Proposal ID (e.g., 2736 for SMACS)")
parser.add_argument('--INSTRUMENT', type=str, nargs="+", help="JWST Detector (e.g., MIRI/NIRCAM)")

args = parser.parse_args()

def query_filename(pid, stage, instrument, filenames):
    dl_dir = "/Volumes/GAMA/jwst/MAST-Query/"+pid + "/" + stage + "/" + instrument + "/"
    my_session = Observations.login(token="6ff31bec13ae4e0eb16d90e0103969e9")
    obs_table = Observations.query_criteria(proposal_id=pid,
                                            project='JWST',
                                            dataproduct_type="IMAGE",
                                            instrument_name=[instrument, instrument+"/IMAGE"]
                                            )
    Observations.enable_cloud_dataset(provider='AWS')
    data_products = Observations.get_product_list(obs_table)

    for filename in filenames:
        products_stage2 = Observations.filter_products(data_products,
                                                       extension="fits",
                                                       productFilename=filename,
                                                       productSubGroupDescription=stage)
       # print(filename, products_stage2)
        Observations.download_products(products_stage2,
                                       download_dir=dl_dir,
                                       curl_flag=False,
                                       cache=True)
    file_new = glob.glob(dl_dir + "/mastDownload/JWST/**/*fits")
    #print(file_new)
    for ff in file_new:
        shutil.copy(ff, dl_dir)
    #shutil.copy(*file_new, dl_dir)
    shutil.rmtree(dl_dir + "/mastDownload/")

def query(pid, instrument):

    obs_table = Observations.query_criteria(proposal_id=pid,
                                            project='JWST',
                                            dataproduct_type="IMAGE",
                                            instrument_name=[instrument, instrument+"/IMAGE"]
                                            )
    return obs_table

def check_all_files_downloaded(all_files, csv_file, query_obj, stage, instrument):

    pid = str(*set(query_obj["proposal_id"]))
    files_in_csv = csv_file["productFilename"]
    base_names = [os.path.basename(file) for file in all_files]

    files_to_download = []
    for file in files_in_csv:
        if file not in base_names:
            print(file)
            files_to_download.append(file)
    if len(files_to_download)>0:
        query_filename(pid, stage, instrument, filenames=files_to_download)
    else:
        pass


def check_files(all_files, csv_file, query_obj, stage, instrument):
    files_to_download = []
    pid = str(*set(query_obj["proposal_id"]))
    for file in all_files:
        temp_size = os.stat(file).st_size
        base_name = os.path.basename(file)
        check_temp = csv_file[csv_file["productFilename"] == base_name]
        if check_temp['size'].values == temp_size:
            pass
        else:
            files_to_download.append(base_name)

    if len(files_to_download)>0:
        query_filename(pid, stage, instrument, filenames=files_to_download)
    else:
        pass

def main(pid, stage, instrument):

    JWST_QUERY_TOOLS_MAST_TOKEN = os.getenv('JWST_QUERY_TOOLS_MAST_TOKEN')
    JWST_QUERY_TOOLS_DOWNLOAD_DIR = os.getenv('JWST_QUERY_TOOLS_DOWNLOAD_DIR')
    if None in [JWST_QUERY_TOOLS_MAST_TOKEN, JWST_QUERY_TOOLS_DOWNLOAD_DIR]:
        print("Please set ENV variables. ")
        exit()

    my_session = Observations.login(token=str(JWST_QUERY_TOOLS_MAST_TOKEN))

    ref_dir = str(JWST_QUERY_TOOLS_DOWNLOAD_DIR)
    dl_dir = os.path.join(ref_dir, pid, stage, instrument)
    os.makedirs(dl_dir, exist_ok=True)
    qpid = query(pid, instrument)

    all_files = glob.glob(dl_dir + "*fits")
    check_csv = pd.read_csv(*glob.glob(dl_dir + "*csv"))
    check_all_files_downloaded(all_files=all_files,
                csv_file=check_csv,
                query_obj=qpid,
                stage=stage,
                instrument=instrument)
#    check_files(all_files=all_files,
#                csv_file=check_csv,
#                query_obj=qpid,
#                stage=stage,
#                instrument=instrument)


if __name__ == "__main__":
    if len(sys.argv)==1:
        parser.print_help()
        print("E.g., >> python3 query_jwst.py --PID `2736` --INSTRUMENT `NIRCAM` --STAGE `UNCAL` --check_miri False")
        sys.exit(1)
    else:
        main(*args.PID, *args.STAGE, *args.INSTRUMENT)
