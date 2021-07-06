'''
This script makes the SDSS spectra overlapped on SPLUS photometry
'''
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sn
sn.set_context("poster")
import glob
import argparse
import sys
import os
from astropy.visualization import hist
from astroML.datasets import fetch_imaging_sample, fetch_sdss_S82standards
from astroML.crossmatch import crossmatch_angular

# Read the file
parser = argparse.ArgumentParser(
    description="""Make a spectras""")

parser.add_argument("fileSdss", type=str,
                    default="teste-program",
                    help="Name of file, taken the prefix")

parser.add_argument("TableSplus", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix")

parser.add_argument("--ymin", required=False, type=float, default=None,
                    help="""Value y-axis min""")

parser.add_argument("--ymax", required=False, type=float, default=None,
                    help="""Value y-axis max""")

cmd_args = parser.parse_args()
file_spec = cmd_args.fileSdss + ".fits"
file_table = cmd_args.TableSplus + ".ecsv"

hdu = fits.open(file_spec)
datadir = "../"
table = Table.read(os.path.join(datadir, file_table), format="ascii.ecsv")

# name of the file of sdss
idsdss = file_spec.split('.fit')[0]
idsdss1 = table["ID"]
print(idsdss1)
#ID from table
idsplus = table["IDsdss"]

# Data from the SDSS spectra
hdudata = hdu[1].data
wl = hdudata.field("wavelength")
#wl = 10**hdudata.field("loglam")
Flux = 1E-17*hdudata.field("flux")

# Data of the SPLUs list
mag, mag_err = [], []
wl_sp = [3485, 3785, 3950, 4100, 4300, 4803, 5150, 6250, 6600, 7660, 8610, 9110]
color = ["#CC00FF", "#9900FF", "#6600FF", "#0000FF", "#009999", "#006600", "#DD8000", "#FF0000", "#CC0066", "#990033", "#660033", "#330034"]
marker = ["s", "o", "o", "o", "o", "s", "o", "s", "o", "s", "o", "s"] ### tienen todos los filtros

mag.append(table["U"][ind]) 
mag.append(table["F378"][ind])
mag.append(table["F395"][ind])
mag.append(table["F410"][ind])
mag.append(table["F430"][ind])
mag.append(table["G"][ind])
mag.append(table["F515"][ind]) 
mag.append(table["R"][ind]) 
mag.append(table["F660"][ind])
mag.append(table["I"][ind]) 
mag.append(table["F861"][ind]) 
mag.append(table["Z"][ind])

# ff = (10**(-(table["R_PStotal"][ind] + 2.41) / 2.5)) / 6250.0**2
# print(ff)
# for i, ii in zip(wl, Flux):
#     if i> 6000 and i< 6300:
#          print(i, ii)

# Find scale factor
m = wl == 6250.289 
wl_part = wl[m]
flux_part = Flux[m]
Fsp = (10**(-(table["R"][ind] + 2.41) / 2.5)) / 6250.0**2
factor = flux_part / Fsp


# PLOTS
fig, ax = plt.subplots(figsize=(12, 5))
ax.spines["top"].set_visible(False)  
ax.spines["right"].set_visible(False)
ax.set(xlim=[3350,9300])

# set Y-axis range (if applicable)
if cmd_args.ymin is not None and cmd_args.ymax is not None:
    plt.ylim(cmd_args.ymin,cmd_args.ymax)
elif cmd_args.ymin is not None:
    plt.ylim(ymin=cmd_args.ymin)
elif cmd_args.ymax is not None:
    plt.ylim(ymax=cmd_args.ymax)
    
#plt.ylim(ymin=-50.0,ymax=200)
ax.set(xlabel='Wavelength $(\AA)$')
ax.set(ylabel=r'F$(\mathrm{10^{-15} erg\ s^{-1} cm^{-2} \AA^{-1}})$')
Flux /=1e-15
ax.plot(wl, Flux, c = "gray", linewidth=1.3, alpha=0.6, zorder=5)
for wl1, mag, colors, marker_ in zip(wl_sp, mag, color, marker): #
    F = (10**(-(mag + 2.41) / 2.5)) / wl1**2
    F /= 1e-15
    #F *= factor
    ax.scatter(wl1, F, c = colors, marker=marker_, s=80, zorder=4)
#ax.axvline(4686, color='r', linewidth=0.3, linestyle='-', zorder = 6, label="He II")
# plt.text(0.70, 0.19, table["ID"].split("R3.")[-1]).replace(".", "-"),
#              transform=ax.transAxes, fontsize=25, weight='bold')

if cmd_args.ymax is not None:
    ax.annotate(str(table["ID"][ind]), xy=(9000, cmd_args.ymax),  xycoords='data', size=13,
            xytext=(-120, -60), textcoords='offset points', 
            bbox=dict(boxstyle="round4,pad=.5", fc="0.94"),)
    ax.annotate("r=" + format(float(table["R"][ind]), '.2f'), xy=(9000, 0.86*cmd_args.ymax),  xycoords='data', size=13,
            xytext=(-120, -60), textcoords='offset points', 
            bbox=dict(boxstyle="round4,pad=.5", fc="0.94"),)
else:
    None
    
ax.legend()
plt.tight_layout()
#asciifile = "sdss-" + file_spec.replace(".fits", ".jpg")
asciifile = idsdss1 + ".jpg"
plt.savefig(asciifile)
