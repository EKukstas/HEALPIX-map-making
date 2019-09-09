from __future__ import division
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits # https://pythonhosted.org/pyfits/
from matplotlib import cm
from collections import Counter
import pickle
import h5py 
import sys
from scipy.io import readsav
from astropy import units as u
from astropy.coordinates import SkyCoord
from decimal import Decimal
# Initialise the module with these two lines
import safe_colours
safe_colours = safe_colours.initialise()
cmap = plt.get_cmap('viridis')
cmap.set_over(cmap(1.0))
cmap.set_under('w')
cmap.set_bad('gray')
####################################################################
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

file = readsav('/hpcdata3/ariekuks/2019_02_28/sim_vars/bahamas/power_spectra/BAHAMAS_subfind_y_z_lt_0p15_image_ngp_0p008z0p150_logM11p0_fov25_smooth_20_ntot+nq+fq_rebin+ms+beam+noise.sav')
#print file.keys()

ra_t = np.arange(0.0, 25.0, 25.0/9000)
dec_t = np.arange(0.0, 25.0, 25.0/9000)
ra = []
dec = []
for i in range(0,9000):
    ra = np.append(ra, ra_t)
    dec = np.append(dec, np.ones(9000) * dec_t[i])


sz_y = np.array(file['image_y'])
sz_y = np.reshape(sz_y, np.size(sz_y))
filename = 'tSZ_BAHAMAS_subfind_y_z_lt_0p15_image_ngp'
##########################################################################################
map_to_analyse = sz_y
##########################################################################################
print 'read in map, making HEALPix'
NSIDE = 2048

# set up a blank map
m = hp.ma(np.arange(hp.nside2npix(NSIDE), dtype=np.double))
pixel_theta, pixel_phi = hp.pix2ang(NSIDE, np.arange(hp.nside2npix(NSIDE)))
mask = 1
m.mask = mask

gal = np.empty(shape=(np.size(dec), 2))
c_icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
c_gal = c_icrs.galactic

gal[:,0] = c_gal.l*u.degree
gal[:,1] = c_gal.b*u.degree
gal[:,1] = 90. - gal[:,1]
pix_index = hp.ang2pix(NSIDE, np.radians(gal[:,1]), np.radians(gal[:,0])) # get pixel indices for all input coordinates, many will be repeated


y = Counter(pix_index) # count the number of repeated pix_index, their values, etc.

print 'No. of pixels used, total: ', np.size(pix_index)
print 'No. of pixels assigned galaxies: ', np.size(y.values())
print 'Mean pixel occupation: ', np.mean(y.values())

ind_d = np.where(np.array(y.values()) > 1)[0] # index in y for multiply-occupied pixels
ind_nd = np.where(np.array(y.values()) == 1)[0] # index in y for 

print 'No. of multiply-occupied pixels: ', np.size(np.where(np.array(y.values()) > 1)[0])
print 'No. of singly-occupied pixels: ', np.size(y.values()) - np.size(np.where(np.array(y.values()) > 1)[0])

N = np.array(list(y)).astype(int) # List of unique pixels; repeated or not

ind_sorted = np.argsort(pix_index).astype(int) # sort pix_index in increasing order, repeated indices bunch up
pix_index_sorted = pix_index[ind_sorted] # apply this sorting
map_to_analyse_sorted = map_to_analyse[ind_sorted] # sort their values to preserve order

in_points_nd = np.searchsorted(pix_index_sorted, N[ind_nd]) # find the sub-indices for singly-occupied pixel indices
m[pix_index_sorted[in_points_nd]] = map_to_analyse_sorted[in_points_nd] # locate these singly-occupied pixels in the HP map and assign their value

# now do the same to multiply-occupied pixels
in_points_d = np.searchsorted(pix_index_sorted, N[ind_d]) # find the sub-index for multiply-occupied pixel indices
for i in range(0, len(in_points_d)): # loop over each one, locating its value and adding it to the rest, divinding by the number of constituent square pixels to get the mean.
    #print i
    j = 0
    sum = 0.0
    while (pix_index_sorted[in_points_d[i]] == pix_index_sorted[in_points_d[i] + j]): # while HP pixel index is the same as the next (they are sorted)
        #print in_points[i] + j
        sum = sum + map_to_analyse_sorted[in_points_d[i]] # this code is only aimed at computing total value and the mean. Change this line to do something else, like the median if that's what's needed.
        j = j + 1
        if in_points_d[i] + j >= np.size(pix_index_sorted):
            print in_points_d[i], in_points_d[i] + j
            break
        m[pix_index_sorted[in_points_d[i]]] = sum/j # assign the mean value to the healpix map. sum==the sum of constituent pixel values; j==the number of constituent pixels
                                                    # replace sum/j by a desired operation, such as just sum.
hp.mollview(m, title='BAHAMAS tSZ', cmap=cmap) # display this map
plt.show()

map_max = 2e-5
map_min = 1e-10
hp.mollview(m, cbar=None, rot=[90, -60.0], cmap=cmap, norm='hist', min=map_min, max=map_max, title='')

# Add a colour bar to the mollweide plot
fig = plt.gcf()
ax = plt.gca()
image = ax.get_images()[0]
c_min = map_min
c_max = map_max
cbar = fig.colorbar(image, ax=ax, ticks=[c_min, c_max], orientation='horizontal', fraction=0.036, pad=0.04, shrink=1.0, aspect=35)
#cbar = fig.colorbar(image, ax=ax, ticks=[c_min, 0.00, c_max], orientation='horizontal', fraction=0.036, pad=0.04, shrink=1.0, aspect=35)
cbar.ax.get_xaxis().labelpad = -13.
cbar.ax.set_xlabel('BAHAMAS tSZ')
#cbar.ax.set_xticklabels([np.round(c_min, decimals=4), 0.000, np.round(c_max, decimals=4)])  # horizontal colorbar
cbar.ax.set_xticklabels(["{:.2E}".format(Decimal(c_min)), "{:.2E}".format(Decimal(c_max))])  # horizontal colorbar
# cbar.ax.set_xticklabels(["{:.5f}".format(Decimal(c_min)), "{:.5f}".format(Decimal(c_max))])  # horizontal colorbar

##################################################################################################
print 'Compare totals', np.mean(map_to_analyse), np.mean(m) # compare the input and output map means, should be very close if mean is preserved. Compare sums if that's what's required.
# plt.savefig('/hpcdata3/ariekuks/2019_02_28/figures/BAHAMAS_hp/map_'+str(filename)+'.png', format='png', dpi=300, bbox_inches='tight')

#print 'Finished computing ns 20 maps'
plt.show()

