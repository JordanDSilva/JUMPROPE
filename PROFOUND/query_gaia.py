from astroquery.gaia import Gaia
import astropy.units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
from astropy.io import ascii
import argparse
import pandas as pd
import numpy as np
import sys
import os
import glob
import shutil

parser = argparse.ArgumentParser(description='Query GAIA')
parser.add_argument('--RA', type=float, nargs='+',
                    help='Central right ascension of catalogue (deg)', default=110.8375)
parser.add_argument('--DEC', type=float, nargs="+", help="Central declination of catalogue (deg)", default=-73.4391)
parser.add_argument('--VID', type=str, help="VISIT ID (e.g., 2736 for SMACS)", default = None)
args = parser.parse_args()

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
Gaia.ROW_LIMIT = -1

def query_gaia(coord_dict, dl_dir = "./"):

    centre_ra = coord_dict["ra"]
    centre_dec = coord_dict["dec"]
    visit_id = coord_dict["visit_id"]
    module_label = coord_dict["module"]

    for ra, dec, visit, module in zip(centre_ra, centre_dec, visit_id, module_label):
        if visit == -999:
           fname_temp = "./GAIA_" + str(ra) + "_" + str(dec)+ ".csv"
        else:
           fname_temp = os.path.join(dl_dir, "GAIA_"+str(visit)+"_"+module+".csv")

        coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")
        radius = u.Quantity(2.0, unit = u.degree)
        jj = Gaia.cone_search(coord, radius)
        gaia_table = jj.get_results()
        # gaia_table_trim = gaia_table
        gaia_table_trim = gaia_table[["ra",
                                      "dec",
                                      "phot_g_mean_mag",
                                      "classprob_dsc_combmod_star",
                                      "phot_bp_rp_excess_factor",
                                      "pm",
                                      "pmra",
                                      "pmdec",
                                      "ref_epoch"
                                      ]]
        ascii.write(gaia_table_trim, fname_temp, overwrite=True,
                    format='csv')

        fname_temp = ""


if __name__ == "__main__":
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    if args.VID is not None:

        ref_dir = os.getenv('JUMPROPE_REF_DIR')
        if None in [ref_dir]:
           print("Please set ENV variables")
           exit()

        df = pd.read_csv( os.path.join(ref_dir, 'ProFound', 'long_warp_info.csv') )
        dl_dir = os.path.join(ref_dir, 'ProFound', 'GAIA_Cats', str(args.VID) )
        os.makedirs(dl_dir, exist_ok=True)

        dff = df.where(df["VISIT_ID"] == args.VID).dropna()

        coord_dict = {"ra" : dff["RA_CEN"], "dec" : dff["DEC_CEN"], "visit_id" : dff["VISIT_ID"], "module" : dff["MODULE"]}
        query_gaia(coord_dict, dl_dir)
    else:
        coord_dict = {"ra" : args.RA, "dec" : args.DEC, "visit_id" : [-999], "module" : [-999]}
        query_gaia(coord_dict)
