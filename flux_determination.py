
# imports
import numpy as np

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS

import plots as plt

# constants
image_dir = 'muse_r-band_images/'
eline_dir = 'eline_tables/'
targets = ['HE0040-1105', 'HE0114-0015', 'HE0119-0118', 'HE0203-0031',
           'HE0212-0059', 'HE0227-0913', 'HE0232-0900', 'HE0345+0056',
           'HE0853+0102', 'HE0934+0119', 'HE2222-0026', 'HE2302-0857']
centers = [ [192,152], [168,163], [114,81], [83,79],
            [161,159], [116,88],  [68,52],  [80,81],
            [163,117], [156,163], [93,91],  [140,161] ]

def circular_mask(height, width, center, radius) :
    
    YY, XX = np.ogrid[:height, :width]
    dist_squared = ((XX - center[0])**2 + (YY - center[1])**2)
    
    return (dist_squared <= radius**2)

def line_flux_map(name, line, center) :
    
    infile = eline_dir + 'full/' + name + '.eline_table.fits' # unbinned table
    
    catalog = Table.read(infile)
    
    eline_x_cor = catalog['x_cor']
    eline_y_cor = catalog['y_cor']
    
    if line == 'HA' :
        flux = catalog['Halpha_flux']
        mask = flux / catalog['Halpha_flux_err'] >= 3
        label = r'$\rm H\alpha$ Flux'
    if line == 'HB' :
        flux = catalog['Hbeta_flux']
        mask = flux / catalog['Hbeta_flux_err'] >= 3
        label = r'$\rm H\beta$ Flux'
    if line == 'NII' :
        flux = catalog['NII6583_flux']
        mask = flux / catalog['NII6583_flux_err'] >= 3
        label = r'[N II] $\lambda 6583$ Flux'
    if line == 'OIII' :
        flux = catalog['OIII5007_flux']
        mask = flux / catalog['OIII5007_flux_err'] >= 3
        label = r'[O III] $\lambda 5007$ Flux'
    
    with fits.open(image_dir+name+'_FOV_SDSSr.fits') as hdu :
        wcs = WCS(hdu[0].header)
        dim = hdu[0].data.shape
    
    image = np.full( (dim[0], dim[1]), np.nan )
    image[eline_y_cor[mask], eline_x_cor[mask]] = flux[mask]
    # plt.emission_line_map(image, wcs, label)
    
    radius = 1.5*u.arcsec/(0.2*u.arcsec/u.pixel) # SDSS has 3" diameter fibers,
     # and MUSE has a spatial resolution of 0.2"/pixel
    fiber_mask = circular_mask(dim[0], dim[1], center, radius.value)
    fiber_image = image[fiber_mask]
    
    return np.nansum(fiber_image)

def determine_BPT_fluxes() :
    
    log_OIII_HB = []
    log_NII_HA = []
    
    for i in range(len(targets)) :
        OIII = line_flux_map(targets[i], 'OIII', centers[i])
        HB = line_flux_map(targets[i], 'HB', centers[i])
        NII = line_flux_map(targets[i], 'NII', centers[i])
        HA = line_flux_map(targets[i], 'HA', centers[i])
        log_OIII_HB.append( np.log10(OIII/HB) )
        log_NII_HA.append( np.log10(NII/HA) )
    
    # print(log_OIII_HB)
    # print(log_NII_HA)
    
    return
