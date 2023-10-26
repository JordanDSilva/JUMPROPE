#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:27:19 2022

@author: Jun01ee
"""
from zipfile import ZipFile
import zipfile
import shutil
import pandas as pd
import requests
# import numpy as np
import argparse
import os
import glob
from http.client import IncompleteRead


parser = argparse.ArgumentParser(description='Getting HAP cutout from MAST with HTTP request')
parser.add_argument('--VID', type=str, help='Observation ID (e.g., 2736001001 for SMACS)')
parser.add_argument('--MODULE', type=str, help='Targeted channel of NIRCam (A or B, case insensitive)')
args = parser.parse_args()

ref_dir = os.getenv('JUMPROPE_REF_DIR')

lookup_path = os.path.join(ref_dir, 'ProFound', 'long_warp_info.csv')

zip_dir = os.path.join(ref_dir, 'ProFound', 'Zip')
zip_path = os.path.join(zip_dir, 'astrocut_{0}_{1}.zip'.format(str(args.VID), str(args.MODULE).upper()) )

data_path = os.path.join(ref_dir, 'ProFound', 'HST_cutout', args.VID, args.MODULE )

#if not os.path.exists(data_path):
#    os.mkdir(data_path)
os.makedirs(data_path, exist_ok=True)
os.makedirs(zip_dir, exist_ok=True)

if os.path.isfile(lookup_path):
    df = pd.read_csv(lookup_path)
    cat = df.loc[(df['VISIT_ID'] == args.VID) &
                 (df['MODULE'] == args.MODULE.upper()), :]
    print(cat)
    targ_ra = cat['RA_CEN'].values[0]
    targ_dec = cat['DEC_CEN'].values[0]

    targ_ra_tl = cat['RA_TL'].values[0]
    targ_dec_tl = cat['DEC_TL'].values[0]
    targ_ra_tr = cat['RA_TR'].values[0]
    targ_dec_tr = cat['DEC_TR'].values[0]
    targ_ra_br = cat['RA_BR'].values[0]
    targ_dec_br = cat['DEC_BR'].values[0]
    targ_ra_bl = cat['RA_BL'].values[0]
    targ_dec_bl = cat['DEC_BL'].values[0]

    targ_ra_mbl = cat['RA_MBL'].values[0]
    targ_dec_mbl = cat['DEC_MBL'].values[0]
    targ_ra_mtl = cat['RA_MTL'].values[0]
    targ_dec_mtl = cat['DEC_MTL'].values[0]
    targ_ra_mtr = cat['RA_MTR'].values[0]
    targ_dec_mtr = cat['DEC_MTR'].values[0]
    targ_ra_mbr = cat['RA_MBR'].values[0]
    targ_dec_mbr = cat['DEC_MBR'].values[0]

    targ_ra_cl = cat['RA_CL'].values[0]
    targ_dec_cl = cat['DEC_CL'].values[0]
    targ_ra_ct = cat['RA_CT'].values[0]
    targ_dec_ct = cat['DEC_CT'].values[0]
    targ_ra_cb = cat['RA_CB'].values[0]
    targ_dec_cb = cat['DEC_CB'].values[0]
    targ_ra_cr = cat['RA_CR'].values[0]
    targ_dec_cr = cat['DEC_CR'].values[0]

    print(targ_ra, targ_dec)
else:
    print(lookup_path+' Do not exist!')

link_0 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra,3), round(targ_dec,3))

link_1 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_tl,3), round(targ_dec_tl,3))
link_2 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_tr,3), round(targ_dec_tr,3))
link_3 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_br,3), round(targ_dec_br,3))
link_4 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_bl,3), round(targ_dec_bl,3))

link_5 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_mbl,3), round(targ_dec_mbl,3))
link_6 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_mtl,3), round(targ_dec_mtl,3))
link_7 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_mtr,3), round(targ_dec_mtr,3))
link_8 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_mbr,3), round(targ_dec_mbr,3))

link_9 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_cl,3), round(targ_dec_cl,3))
link_10 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_ct,3), round(targ_dec_ct,3))
link_11 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_cb,3), round(targ_dec_cb,3))
link_12 = 'https://mast.stsci.edu/hapcut/api/v0.1/astrocut?ra={0}&dec={1}&x=5000&y=5000'.format(round(targ_ra_cr,3), round(targ_dec_cr,3))

# from urllib.request import urlretrieve
# urlretrieve(link, zip_path, str(args.MODULE).upper()))

#with requests.get(link, stream=True) as r:
#    print(r.status_code)
#    print(r.json())
#    with open(zip_path, 'wb') as f:
#        shutil.copyfileobj(r.raw, f)
#f.close()

import wget
for link in [link_0, link_1, link_2,
             link_3, link_4, link_5,
             link_6, link_7, link_8,
             link_9, link_10, link_11,
             link_12]:

    for ii in range(5):
        try:
            wget.download(link, out = zip_path)
            break
        except IncompleteRead:
            if ii == 4:
                raise
            else:
                continue


    if zipfile.is_zipfile(zip_path):
       with ZipFile(zip_path, 'r') as zip_ref:
           listOfFileNames = zip_ref.namelist()
           # Iterate over the file names
           for fileName in listOfFileNames:
                # Check filename endswith csv
              if 'coarse' not in fileName:
               # Extract a single file from zip
                  zip_ref.extract(fileName, data_path)
       if os.path.isfile(zip_path):
          os.remove(zip_path)
       else:
          pass
    else:
        if os.path.isfile(zip_path):
           os.remove(zip_path)
        else:
           pass

