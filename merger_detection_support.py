# imports
import numpy as np
# import os

# from astropy.coordinates import Angle
# from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
# import astropy.units as u
from scipy import ndimage

# constants
num_bootraps_runs = 1000

def as_params(file_list, scale_list) :
    
    asymms, asymms_l, asymms_u = [], [], []
    clumps, clumps_l, clumps_u = [], [], []
    for i in range(len(file_list)) :
        with fits.open('SDSS_images/' + file_list[i]) as hdu :
            science_image = hdu[0].data # get the science data that will be used
        
        rot_image = np.rot90(science_image, 2)
            # rotate the original image by 180 degrees
        asymm, asymm_l, asymm_u = asymmetry(science_image, rot_image)
        
        smooth_image = ndimage.gaussian_filter(science_image, scale_list[i])
            # smooth the original image by the Gaussian kernel
        clump, clump_l, clump_u = clumpiness(science_image, smooth_image)
        
        asymms.append(asymm)
        asymms_l.append(asymm-asymm_l)
        asymms_u.append(asymm_u-asymm)
        clumps.append(clump)
        clumps_l.append(clump-clump_l)
        clumps_u.append(clump_u-clump)
    
    
    table = Table([asymms, asymms_l, asymms_u,
                   clumps, clumps_l, clumps_u],
                  names=('Asymm', 'Asymm_lo', 'Asymm_hi',
                          'Clump', 'Clump_lo', 'Clump_hi'))
    # table.pprint(max_width=-1)
    
    # plt.plot_with_errorbars(table['Asymm'], [table['Asymm_lo'],table['Asymm_hi']], 'Asymm',
    #                         table['Clump'], [table['Clump_lo'],table['Clump_hi']], 'Clump',
    #                         xmin=0, xmax=1.4, ymin=0, ymax=2)
    
    return

def asymmetry(image, rotated) :
    
    numer = (image - rotated)**2
    denom = image**2
    asymm = 0.5*np.sum(numer)/np.sum(denom)
    
    numer_flat = numer.flatten()
    denom_flat = denom.flatten()
    numer_len = len(numer_flat)
    denom_len = len(denom_flat)
    
    # bootstrap
    asymms = []
    for i in range(num_bootraps_runs) :
        numer_resample = np.random.choice(numer_flat,size=numer_len,replace=True)
        denom_resample = np.random.choice(denom_flat,size=denom_len,replace=True)
        asymm_resample = 0.5*np.sum(numer_resample)/np.sum(denom_resample)
        asymms.append(asymm_resample)
    
    asymms = np.array(asymms)
    # plt.histo(asymms, 'Asymmetry')
    
    # std_dev = np.std(asymms)
    # return asymm, std_dev
    
    asymm = np.median(asymms)
    asymm_l = np.percentile(asymms, 16)
    asymm_u = np.percentile(asymms, 84)
    return asymm, asymm_l, asymm_u

def clumpiness(image, smoothed) :
    
    numer = image - smoothed
    numer = np.ma.masked_where(numer < 0, numer)
    np.ma.set_fill_value(numer, 0)
    numer = numer.filled()
    denom = image
    clumpy = np.sum(numer)/np.sum(denom)
    
    numer_flat = numer.flatten()
    denom_flat = denom.flatten()
    numer_len = len(numer_flat)
    denom_len = len(denom_flat)
    
    # bootstrap
    clumpys = []
    for i in range(num_bootraps_runs) :
        numer_resample = np.random.choice(numer_flat,size=numer_len,replace=True)
        denom_resample = np.random.choice(denom_flat,size=denom_len,replace=True)
        clumpy_resample = np.sum(numer_resample)/np.sum(denom_resample)
        clumpys.append(clumpy_resample)
    
    clumpys = np.array(clumpys)
    # plt.histo(clumpys, 'Clumpiness')
    
    # std_dev = np.std(clumpys)
    # return clumpy, std_dev
    
    clumpy = np.median(clumpys)
    clumpy_l = np.percentile(clumpys, 16)
    clumpy_u = np.percentile(clumpys, 84)
    return clumpy, clumpy_l, clumpy_u