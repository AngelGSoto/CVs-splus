'''
This is a simply script to make table with the SDSS format for cross-match.
'''
from astropy.table import Table, vstack
import numpy as np
import argparse
import sys
import os
import pandas as pd

parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("source", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")

cmd_args = parser.parse_args()
file_ = cmd_args.source + ".csv"

tab = pd.read_csv(file_)

n = len(tab)

sep = np.linspace(2.0/60., 2.0/60., num=n)
ra = tab["RA"]
dec = tab["DEC"]
table = Table([ra, dec, sep], names=('ra', 'dec', 'sep'), meta={'name': 'first table'})

#Save the file
asciifile = file_.replace(".csv", "-match-sdss.dat")
table.write(asciifile, format="ascii", delimiter=',', overwrite=True)
    
