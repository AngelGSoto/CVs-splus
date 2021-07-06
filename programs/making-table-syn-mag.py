'''
This script makes table with synthectic fotometry from the json files 
'''
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table, QTable
from astropy.coordinates import SkyCoord 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
import sys
import os
from astropy.visualization import hist
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.crossmatch import crossmatch_angular
from collections import OrderedDict
import glob
import json
import seaborn as sn
sn.set_context("poster")

# Read FITS files
pattern = "*-spectros/*.fits"
file_fits = glob.glob(pattern)

# Read JSON files
pattern1 = "*-spectros/*-SPLUS21-magnitude.json"
file_list = glob.glob(pattern1)

shape1 = (len(file_fits), 4)

inffits = []
for name_fit in file_fits:
    inffits.append(name_fit.split("s/")[-1].split(".fit")[0])
    hdulist = fits.open(name_fit)
    c = SkyCoord(ra=float(hdulist[1].header["RAOBJ"])*u.degree, dec=float(hdulist[1].header["DECOBJ"])*u.degree) 
    inffits.append('SDSSJ{0}{1}'.format(c.ra.to_string(sep='', precision=2, pad=True), c.dec.to_string(sep='', precision=1, alwayssign=True)))
    inffits.append(float(hdulist[1].header["RAOBJ"]))
    inffits.append(float(hdulist[1].header["DECOBJ"]))

shape = (len(file_list), 13)

mag = []
for file_name in file_list:
    mag.append(file_name.split("s/")[-1].split("-catB")[0])
    with open(file_name) as f:
        data = json.load(f)
        data = OrderedDict((k, v) for k, v in sorted(data.items(), key=lambda x: x[0]))
    mag.append(float(data["F0352_uJAVA"]))
    mag.append(float(data["F0378"]))
    mag.append(float(data["F0395"]))
    mag.append(float(data["F0410"]))
    mag.append(float(data["F0430"]))
    mag.append(float(data["F0475_gSDSS"]))
    mag.append(float(data["F0515"]))
    mag.append(float(data["F0626_rSDSS"]))
    mag.append(float(data["F0660"]))
    mag.append(float(data["F0769_iSDSS"]))
    mag.append(float(data["F0861"]))
    mag.append(float(data["F0883_zSDSS"]))


XX_fits = np.array(inffits).reshape(shape1)
print("Data shape:", XX_fits.shape)

XX = np.array(mag).reshape(shape)
print("Data shape:", XX.shape)

# Tables with all information 
new_order = ['ID', 'RA', 'DEC', 'U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I', 'F861', 'Z']

t1 = Table(XX_fits, names=('IDsdss', 'ID', 'RA', 'DEC'), meta={'name': 'first table'}, dtype=('S', 'S', 'f8', 'f8'))
t2 = Table(XX, names=('IDsdss', 'U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I', 'F861', 'Z'),
                                                           meta={'name': 'first table'}, dtype=('S', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))
t1.sort('IDsdss')
t2.sort('IDsdss')

t2['ID'] = t1['ID']
t2['RA'] = t1['RA']
t2['DEC'] = t1['DEC']

# Order of the columns of the tables
table_new = t2[new_order]

table_new.write("CV-synthetic-mag-splus.ecsv", format="ascii.ecsv", overwrite=True)
table_new.write("CV-synthetic-mag-splus.fits", format="fits", overwrite=True)
