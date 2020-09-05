
# imports
import numpy as np
import os

from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from scipy import ndimage

import requests
import statmorph
import photutils

# import plots as plt

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
# SDSS_path = 'catalogs/CARS_SDSS/SDSS_allObjects_logMass10_complete_0-52113.fits'
GAMA_path = 'catalogs/CARS_GAMA/GAMA_allObjects_logMass10_complete_0-1364.fits'
pixel_scale = (0.396*u.arcsec)/(1*u.pix) # 1 pixel = 0.396" for SDSS
scale_of_interest = 15*u.kpc # bubbles and cavities are often found on such scales

def main(calculate_morphologies=False, download=False) :
    
    import matplotlib.pyplot as plt
    
    GAMA = Table.read(GAMA_path)
    
    RAs = GAMA['RA']
    Decs = GAMA['DEC']
    redshifts = GAMA['Z']
    
    # scale = cosmo.arcsec_per_kpc_proper(redshifts)*scale_of_interest/pixel_scale
        # determine number of pixels that correspond to 15 kpc in projected
        # size, to highlight AGN driven features like bubbles and cavities
    # kpc_per_pixel = scale_of_interest/scale
    
    files = download_images(RAs, Decs, download=download)
    
    if calculate_morphologies == True :
        asymmetries = []
        smoothnesses = []
        ginis = []
        M20s = []
        for i in range(len(files)) :
            try :
                with fits.open('SDSS_images/' + files[i]) as hdu :
                    image = hdu[0].data
                    image[image < 0] = 0        
                    sigma = np.sqrt(image)
                    # plt.imshow(image)
                
                    threshold = photutils.detect_threshold(image, 1.5)
                    npixels = 5
                    segm = photutils.detect_sources(image, threshold, npixels)
                    label = np.argmax(segm.areas) + 1
                    segmap = segm.data == label
                    segmap_float = ndimage.uniform_filter(np.float64(segmap), size=10)
                    segmap = segmap_float > 0.5
                    # plt.imshow(segmap, origin='lower', cmap='gray')
                    
                    source_morphs = statmorph.source_morphology(image, segmap, weightmap=sigma)
                    morph = source_morphs[0]
                    # params = [morph.gini, morph.m20, morph.asymmetry, morph.smoothness]
                    asymmetries.append(morph.asymmetry)
                    smoothnesses.append(morph.smoothness)
                    ginis.append(morph.gini)
                    M20s.append(morph.m20)
                    
                    from statmorph.utils.image_diagnostics import make_figure
                    make_figure(morph, 'SDSS_images/PNGs/' + files[i][:-5] + '.png')
            except :
                asymmetries.append(-99)
                smoothnesses.append(-99)
                ginis.append(-99)
                M20s.append(-99)
    
    GAMA['Asymm'] = asymmetries
    GAMA['Smooth'] = smoothnesses
    GAMA['Gini'] = ginis
    GAMA['M20'] = M20s
    
    GAMA.write('catalogs/CARS_GAMA/GAMA_allObjects_logMass10_complete_0-1364_morphologies.fits')
    
    # files = []
    # for i in range(len(small)) :
    #     file = ('cutout_%.4f_%.4f.fits') % (RAs[i], Decs[i])
    #     files.append(file)
    
    # as_params(files, scale.value)
    
    return

def image_urls(ra_list, dec_list) :
    
    filenames = []
    url_list = []
    for i in range(len(ra_list)) :
        filename = 'cutout_%.4f_%.4f.fits' % (ra_list[i], dec_list[i])
        filenames.append(filename)
        
        url = ('http://legacysurvey.org/viewer/cutout.fits?ra=%s&dec=%s'
               % (str(ra_list[i]), str(dec_list[i])) )
        url += '&layer=sdss&pixscale=0.396&size=120&bands=r'
        url_list.append(url)
    
    return filenames, url_list

def download_images(ra_list, dec_list, download=False) :
    
    files, urls = image_urls(ra_list, dec_list)
    
    if download == True :
        for i in range(len(ra_list)) :
            r = requests.get(urls[i], stream=True)        
            with open('SDSS_images/' + files[i], 'wb') as f :
                for chunk in r.iter_content(chunk_size=128):
                    f.write(chunk)
    
    return files

main(calculate_morphologies=False)
