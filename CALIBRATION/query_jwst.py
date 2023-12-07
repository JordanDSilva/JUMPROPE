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
import time

parser = argparse.ArgumentParser(description='Download from MAST')
parser.add_argument('--STAGE', type=str, nargs='+',
                    help='Product level (e.g., UN/CAL for un/cal files)')
parser.add_argument('--VISITID', type=str, nargs="+", help="VISIT ID (e.g., 2736 for SMACS)")
parser.add_argument('--INSTRUMENT', type=str, nargs="+", help="JWST Detector (e.g., MIRI/NIRCAM)")
parser.add_argument('--check_miri', type=bool, nargs="+", help="Look for MIRI frames in NIRCAM pid. Default is False", default=False)

args = parser.parse_args()
def query_miri(nircam_query, miri_dl_dir):
    center_ra = pd.unique(nircam_query['s_ra'])
    center_dec = pd.unique(nircam_query['s_dec'])
    pid = pd.unique(nircam_query["proposal_id"])
    rad = u.Quantity(5.0 / 60.0, unit="deg")
    k = 1
    query_list = []
    query_tables = []
    for ra, dec in zip(center_ra, center_dec):
        coord = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
        miri_table = Observations.query_criteria(coordinates=coord,
                                                 radius=rad,
                                                 project='JWST',
                                                 dataproduct_type="IMAGE",
                                                 instrument_name=["MIRI", "MIRI/IMAGE"],
                                                 )
        data_products = Observations.get_product_list(miri_table)
        miri_stage2 = Observations.filter_products(data_products,
                                                   extension="fits",
                                                   productSubGroupDescription="CAL")
        query_tables.append(miri_stage2)
        Observations.download_products(miri_stage2,
                                       download_dir=miri_dl_dir,
                                       curl_flag=True,
                                       )
        file = glob.glob(miri_dl_dir + '/mastDownload*.sh')
        shutil.move(*file, miri_dl_dir + str(*pid) + '_' + "MIRI" + '_' + "CAL" + "_" + str(k) + '.sh')

        query_list.append(miri_table)
        Observations.download_products(miri_stage2,
                                           download_dir=miri_dl_dir,
                                           curl_flag=False,
                                           cache=True)
        k += 1
    stack_table_temp = astropy.table.vstack(query_tables)
    stack_table = astropy.table.unique(stack_table_temp, keys='obsID')
    ascii.write(stack_table, miri_dl_dir + str(*pid) + '_' + "MIRI" + '_' + "CAL" + '.csv',
                overwrite=True,
                format='csv')
    return query_list

def count_unique(l):
    seen = set()
    for x in l:
        if x not in seen:
            seen.add(x)
    print(len(seen))
def read_table():
    files = glob.glob('./*.csv')
    x = [pd.read_csv(file) for file in files]

    for n, file in enumerate(x):
        print(files[n])
        print(len(file['obs_id']))

        count_unique(file['obs_id'])


def query_filename(dl_dir, visit_id, stage, instrument, filename):
    pid = str(visit_id)[0:4]
    obs_table = Observations.query_criteria(proposal_id=pid,
                                            project='JWST',
                                            dataproduct_type="IMAGE",
                                            instrument_name=[instrument, instrument+"/IMAGE"]
                                            )
    # Observations.enable_cloud_dataset(provider='AWS')
    data_products = Observations.get_product_list(obs_table)
    products_stage2 = Observations.filter_products(data_products,
                                                   extension="fits",
                                                   productFilename=filename,
                                                   productSubGroupDescription=stage)
    Observations.download_products(products_stage2,
                                   download_dir=dl_dir,
                                   curl_flag=False,
                                   cache=True)
    file_new = glob.glob(dl_dir + "/mastDownload/JWST/**/*fits")
    for ff in file_new:
        shutil.move(ff, os.path.join(dl_dir, os.path.basename(ff)))

    shutil.rmtree(dl_dir + "/mastDownload/")


def query_table(dl_dir, table):
    """Get a pandas table and convert to astropy table for easy downloading"""
    products_stage2 = astropy.table.Table.from_pandas(table)
    Observations.download_products(products_stage2,
                                   download_dir=dl_dir,
                                   curl_flag=False,
                                   cache=True)
    file_new = glob.glob(dl_dir + "/mastDownload/JWST/**/*fits")
    for ff in file_new:
        shutil.move(ff, os.path.join(dl_dir, os.path.basename(ff)))

    shutil.rmtree(dl_dir + "/mastDownload/")

def query(visit_id, stage, instrument, dl_dir, dl_products=False):
    pid = str(visit_id)[0:4]
    print("Querying observation")
    obs_table = Observations.query_criteria(proposal_id=pid,
                                            project='JWST',
                                            dataproduct_type="IMAGE",
                                            instrument_name=[instrument, instrument+"/IMAGE"])

    # Observations.enable_cloud_dataset(provider='AWS')
    data_products = Observations.get_product_list(obs_table)
    products_stage2 = Observations.filter_products(data_products,
                                                   extension="fits",
                                                   productSubGroupDescription=stage,
                                                   )

    visit_ids = np.array([obs[3:13] for obs in products_stage2["obs_id"]])
    # idx = np.where(visit_ids == str(visit_id))
    idx = np.flatnonzero(np.core.defchararray.find(visit_ids,str(visit_id))!=-1)

    print("Moving the table")
    Observations.download_products(products_stage2[idx], download_dir=dl_dir, curl_flag=True)
    file = glob.glob(dl_dir + '/mastDownload*.sh')
    for ff in file:
        shutil.move(ff, dl_dir + pid + '_' + instrument + '_' + stage + '.sh')
    ascii.write(products_stage2, dl_dir + pid + '_' + instrument + '_' + stage + '.csv', overwrite=True, format='csv')

    if dl_products:
        Observations.download_products(products_stage2[idx],
                                   download_dir=dl_dir,
                                   curl_flag=False,
                                   cache=True)

    return obs_table


def check_files(dl_dir, csv_file, visit_id, stage, instrument):
    all_files = glob.glob(dl_dir + "*fits")
    all_files_redo = [ii for ii in all_files if visit_id in ii]
    all_files_redo_basenames = [os.path.basename(ii) for ii in all_files_redo]

    csv_file_redo = [ii for ii in csv_file["productFilename"] if visit_id in ii]

    # print(all_files_redo_basenames)
    for file in csv_file_redo:
        if file not in all_files_redo_basenames:
            print("Missing file, downloading: " + file)
            query_filename(dl_dir, visit_id, stage, instrument, filename=file)

    all_files = glob.glob(dl_dir + "*fits")
    all_files_redo = [ii for ii in all_files if visit_id in ii]

    for file in all_files_redo:
        temp_size = os.stat(file).st_size
        base_name = os.path.basename(file)
        check_temp = csv_file[csv_file["productFilename"] == base_name]
        if check_temp['size'].values != temp_size:
            print(check_temp['size'].values, temp_size)
            print("File wrong size, redownloading: " + file)
            query_filename(dl_dir, visit_id, stage, instrument, filename=base_name)

def check_filesV2(dl_dir, csv_file, visit_id):
    all_files = glob.glob(dl_dir + "*fits")
    all_files_redo = [ii for ii in all_files if visit_id in ii]
    all_files_redo_basenames = [os.path.basename(ii) for ii in all_files_redo]

    csv_file_redo = csv_file[csv_file["productFilename"].str.contains(visit_id)]

    df1 = csv_file_redo[~csv_file_redo["productFilename"].isin(all_files_redo_basenames)]
    if len(df1)>0:
        print("Missing " + str(len(df1)) + " files, downloading!")
        query_table(dl_dir, df1)
    # print(df1)

    all_files = glob.glob(dl_dir + "*fits")
    all_files_redo = [ii for ii in all_files if visit_id in ii]
    all_files_redo_basenames = [os.path.basename(ii) for ii in all_files_redo]
    temp_size = [os.stat(file).st_size for file in all_files_redo]

    df1 = csv_file_redo[csv_file_redo["productFilename"].isin(all_files_redo_basenames)]
    df2 = df1.loc[~(df1["size"].values  == temp_size)]
    print(df2)
    print(sum(df1["size"].values  == temp_size))
    if len(df2)>0:
        print(str(len(df2)) + " files wrong sizes, redownloading!")
        query_table(dl_dir, df2)

def main(visit_id, stage, instrument, check_miri):

    JUMPROPE_MAST_TOKEN = os.getenv('JUMPROPE_MAST_TOKEN')
    JUMPROPE_DOWNLOAD_DIR = os.getenv('JUMPROPE_DOWNLOAD_DIR')

    if None in [JUMPROPE_MAST_TOKEN, JUMPROPE_DOWNLOAD_DIR]:
        print("Please set ENV variables. ")
        exit()

    my_session = Observations.login(token=str(JUMPROPE_MAST_TOKEN))
    pid = str(visit_id)[0:4]
    ref_dir = str(JUMPROPE_DOWNLOAD_DIR)

    dl_dir = os.path.join(ref_dir, pid, stage, instrument) + str("/")
    print(dl_dir)
    os.makedirs(dl_dir, exist_ok=True)
    all_files = glob.glob(
        dl_dir + "*fits")
    # print(all_files)
    all_files_redo = [ii for ii in all_files if visit_id in ii]

    qpid = query(visit_id, stage, instrument, dl_dir, dl_products=False) ##produce the csv file with all exposure information
    check_csv = pd.read_csv(*glob.glob(dl_dir + "*csv"))

    if len(all_files_redo)==0:
        ## if the target directory has nothing in it, then download everything in that visitid
        print("Never seen visit " + visit_id + " before. Downloading everything!")
        query(visit_id, stage, instrument, dl_dir, dl_products=True)
        time.sleep(300) #wait 5 minuts and check files
        check_filesV2(dl_dir = dl_dir,
                    csv_file=check_csv,
                    visit_id = visit_id)
    if len(all_files_redo)>0:
        ## if there are already files in the directory
        ## check that we've downloaded everything
        check_filesV2(dl_dir=dl_dir,
                      csv_file=check_csv,
                      visit_id=visit_id)
        mast_files = glob.glob(dl_dir + "/mastDownload/JWST/**/*fits")
        # print(mast_files) ## check if there are any files in the mastDownload folder
        if len(mast_files) == 0:
            pass
        else:
            for file in mast_files:
                shutil.move(file, dl_dir)
            shutil.rmtree(dl_dir + "/mastDownload/")

    if instrument == "MIRI" or not check_miri: #not because the default is false
        pass
    else:
        miri_dl_dir = dl_dir.replace("NIRCAM", "MIRI")
        miri_dl_dir = miri_dl_dir.replace("UNCAL", "CAL") + "/"
        os.makedirs(miri_dl_dir, exist_ok=True)
        qmiri = query_miri(nircam_query=qpid, miri_dl_dir=miri_dl_dir)
        mast_files = glob.glob(miri_dl_dir + "/mastDownload/JWST/**/*fits")
        if len(mast_files) == 0:
            pass
        else:
            for file in mast_files:
                shutil.move(file, miri_dl_dir)
            shutil.rmtree(miri_dl_dir + "/mastDownload/")
        all_files = glob.glob(miri_dl_dir + "*fits")
        check_csv = pd.read_csv(*glob.glob(miri_dl_dir + "*csv"))
        for i in range(len(qmiri)):
            check_files(csv_file=check_csv,
                        query_obj=qmiri[i],
                        stage="CAL",
                        instrument="MIRI")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        print("E.g., >> python3 query_jwst.py --VISITID `2736001001` --INSTRUMENT `NIRCAM` --STAGE `UNCAL` --check_miri False")
        sys.exit(1)
    else:
        main(*args.VISITID, *args.STAGE, *args.INSTRUMENT, args.check_miri)
