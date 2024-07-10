import sys
import pylab
import math
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm#, colors
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.io as io
from scipy.io import readsav
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from IPython import display

from PIL import Image

textwidth = 17.
pylab.rcParams['figure.figsize'] = textwidth, textwidth*1   # = x, y
pylab.rcParams['axes.labelsize'] = textwidth*1.5
pylab.rcParams['xtick.labelsize'] = textwidth*1.
pylab.rcParams['ytick.labelsize'] = textwidth*1.
pylab.rcParams['legend.fontsize'] = textwidth*1.5
pylab.rcParams['axes.titlesize'] = textwidth*1.5
pylab.rcParams['axes.linewidth'] = textwidth*0.2
pylab.rcParams['lines.linewidth'] = textwidth*0.2
pylab.rcParams['lines.markersize'] = textwidth*0.2
pylab.rcParams['xtick.direction'] = 'in'
pylab.rcParams['ytick.direction'] = 'in'
pylab.rcParams['ytick.right'] = 'true'
pylab.rcParams['xtick.top'] = 'true'
pylab.rcParams['legend.frameon'] = 'False'
pylab.rcParams['legend.labelspacing'] = 1

pylab.rcParams['xtick.major.pad'] = textwidth*1.
pylab.rcParams['xtick.major.size'] = textwidth*0.7
pylab.rcParams['xtick.major.width'] = textwidth*0.15
pylab.rcParams['xtick.minor.size'] = textwidth*0.4
pylab.rcParams['xtick.minor.width'] = textwidth*0.15
pylab.rcParams['xtick.minor.visible'] = True
#
pylab.rcParams['ytick.major.pad'] = textwidth*1.
pylab.rcParams['ytick.major.size'] = textwidth*0.7
pylab.rcParams['ytick.major.width'] = textwidth*0.15
pylab.rcParams['ytick.minor.size'] = textwidth*0.4
pylab.rcParams['ytick.minor.width'] = textwidth*0.15
pylab.rcParams['ytick.minor.visible'] = True

# wavelengths to be plotted
waves = ["GALEX_NUV", "SDSS_g", "2MASS_Ks"]
# from left to right: masked data, aziumthally averaged data, model maps
version = ["obs_", "av_", "model_"]
distance = 8.58e6   # pc
pixsize = [[1.5,1.5,1.5],
           [1.,1.,1.],
           [2.85,2.85,3.],
           [14, 14, 14]]  # arcsec/pix

# controls the levels in the map plots
# adjust these until the maps are legible.
vmin = [-0.04,-1,-4]
vmax = [0.18,4.5,18]


fig, axs = plt.subplots(3, 3, figsize=(textwidth, textwidth),# constrained_layout=True,
                        sharex=True, sharey=True)
plt.subplots_adjust(left=0.11, bottom=0.07, right=None, top=None, wspace=0., hspace=0.)

for i in range(0,3):
    for j in range(0,3):
        ax = axs[i,j]
        data = fits.open(f"./maps/{version[j]}{waves[i]}.fits")
        image = data[0]
        hdr = data[0].header   # read fits header
        data = image.data
        data = np.flipud(data)

        #ax = plt.subplot2grid((4,3), (i,j), rowspan=1, colspan=1)
        lim = 15
        extent=[-lim,lim,-lim,lim]
        ax.imshow(data, cmap='inferno', vmin=vmin[i], vmax=vmax[i], extent=extent)
        #ax.imshow(data, cmap='inferno', extent=extent)

	# plot ellipses at certain radii as fiducial markers
        #radius= [0.8,2.7,7]
        #colours = ["green","green","white"]
        #ellipses = [Ellipse(
        #xy=[0.1,0.1],
        #width= radius[ii]*2,
        #height= radius[ii]*2*np.cos(20.3*(np.pi/180)),
        #angle=90+12,
        #edgecolor=colours[ii],
        #facecolor="none",
        #alpha=1,
        #linewidth=2)
        #for ii in range(0,len(radius))]

#        for ii, ell in enumerate(ellipses):
#            ax.add_patch(ell)

#        if i != 2:
#            ax.set_xticklabels([])
#        if j != 0:
#            ax.set_yticklabels([])
        if i == 0:
            if j == 0:
                ax.set_title("Data")
            if j == 1:
                ax.set_title("Data - averaged")
            if j == 2:
                ax.set_title("Model")
        if j == 0:
            ax.set_ylabel(waves[i].replace("_"," "))
        if i == 2 and j==1:
            ax.set_xlabel("Radius (kpc)")

        ax.spines['bottom'].set_color('w')
        ax.spines['top'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')
        #ax.xaxis.label.set_color('w')
        ax.tick_params(colors='w', which='both', labelcolor='k')
        #plt.axis([xc[i][j]-lim, xc[i][j]+lim, yc[i][j]-lim, yc[i][j]+lim])

#plt.tight_layout(pad=0.1)
plt.savefig("./maps.pdf")
