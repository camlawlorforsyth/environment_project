
# imports
import numpy as np

from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u
import scipy.stats as sp

import plots as plt

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
currentFig = 1

def mass_comparison(cat_name) :
    
    if (cat_name == 'SDSS') :
        cat_path = 'catalogs/joined_cats/SDSS_gal-info_totlgm_masked-besides-z-0075.fits'
    if (cat_name == 'GAMA') :
        cat_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_masked-besides-z-0075.fits'
    
    catalog = Table.read(cat_path)
    
    if (cat_name == 'SDSS') :
        mask = ( (catalog['KCOR_MAG'][:,0] > 0) &
                 (catalog['KCOR_MAG'][:,2] > 0) &
                 (catalog['MEDIAN'] >= 6) )
        catalog = catalog[mask]
        
        D_A = cosmo.angular_diameter_distance(catalog['Z'])
        D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)/u.pc
        M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        logMass = (1.15 + 0.7*(catalog['KCOR_MAG'][:,0] -
                               catalog['KCOR_MAG'][:,2]) - 0.4*M_i)*u.solMass
        catalog['colourMass'] = logMass
        
        lastmask = (catalog['colourMass'] >= 6*u.solMass)
        catalog = catalog[lastmask]
        
        catmass = catalog['MEDIAN']*u.solMass
    
    if (cat_name == 'GAMA') :
        D_A = cosmo.angular_diameter_distance(catalog['Z_1'])
        fluxscale = catalog['fluxscale']
        condition = ((fluxscale >= 0.1) & (fluxscale <= 10))
        M_i_corrected = catalog['absmag_i'] - 2.5*np.log10(catalog['fluxscale'])
        M_i = np.where(condition, M_i_corrected, catalog['absmag_i'])
            # see required aperture correction for stellar mass or absolute magnitude
            # http://www.gama-survey.org/dr3/data/cat/StellarMasses/v20/StellarMasses.notes
        logMass = (1.15 + 0.7*catalog['gminusi'] - 0.4*M_i)/u.mag*u.solMass
                    # based on relation from Taylor+ 2011, MNRAS, 418, 1587
        catalog['colourMass'] = logMass
        
        catmass = catalog['logmstar']
    
    logmass = catalog['colourMass']
    residual = logmass - catmass
    logmass = logmass/u.solMass
    residual = residual/u.solMass
    
    nbins = 30
    meds, bin_edge, bin_num = sp.binned_statistic(logmass, residual, 'median', nbins, [(8.45,12)])
    lo, lo_edge, lo_bn = sp.binned_statistic(logmass, residual, lower, nbins, [(8.45,12)])
    hi, hi_edge, hi_bn = sp.binned_statistic(logmass, residual, upper, nbins, [(8.45,12)])
    
    plt.plot_with_errors(bin_edge[:-1], r'$\log{(M_*/M_\odot)} = 1.15 + 0.70(g-i) - 0.4M_i$',
                          meds, r'$\Delta \log{(M_*/M_\odot)}$', lo, hi,
                          xmin=8.3, xmax=12, ymin=-1, ymax=1)
    
    return

def upper(array) :
    return np.percentile(array, 84)

def lower(array) :
    return np.percentile(array, 16)

# mass_comparison('SDSS')
# mass_comparison('GAMA')
