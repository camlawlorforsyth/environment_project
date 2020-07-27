
# imports
import numpy as np

from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u
import scipy.stats as sp

import functions as funcs
import plots as plt

import warnings
warnings.filterwarnings('ignore')

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
currentFig = 1

def mass_comparison(cat_name, intercept_corr=False, linear_corr=False, n_bins=25) :
    
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
        
        # D_A = cosmo.angular_diameter_distance(catalog['Z'])
        D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)/u.pc
        
        if intercept_corr :
            g_int_corr = intercept_corr[0]
            i_int_corr = intercept_corr[1]
            M_g = catalog['KCOR_MAG'][:,0] - 5*np.log10(D_L) + 5 + g_int_corr
            M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 + i_int_corr
        elif linear_corr :
            g_slope_corr = linear_corr[0]
            g_int_corr = linear_corr[1]
            i_slope_corr = linear_corr[2]
            i_int_corr = linear_corr[3]
            M_g = g_slope_corr*(catalog['KCOR_MAG'][:,0] - 5*np.log10(D_L) + 5) + g_int_corr
            M_i = i_slope_corr*(catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5) + i_int_corr
        else :
            M_g = catalog['KCOR_MAG'][:,0] - 5*np.log10(D_L) + 5 # absolute g mag.
            M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        
        logMass = (1.15 + 0.7*(M_g - M_i) - 0.4*M_i)*u.solMass
        catalog['colourMass'] = logMass
        
        lastmask = (catalog['colourMass'] >= 6*u.solMass)
        catalog = catalog[lastmask]
        
        catmass = catalog['MEDIAN']*u.solMass
    
    if (cat_name == 'GAMA') :
        # D_A = cosmo.angular_diameter_distance(catalog['Z_1'])
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
    
    nbins = n_bins
    meds, bin_edge, bin_num = sp.binned_statistic(logmass, residual, 'median', nbins, [(9,12)])
    lo, lo_edge, lo_bn = sp.binned_statistic(logmass, residual, lower, nbins, [(9,12)])
    hi, hi_edge, hi_bn = sp.binned_statistic(logmass, residual, upper, nbins, [(9,12)])
    
    plt.plot_with_errors(bin_edge[:-1], r'$\log{(M_*/M_\odot)} = 1.15 + 0.70(g-i) - 0.4M_i$',
                          meds, r'$\Delta \log{(M_*/M_\odot)}$', lo, hi,
                          xmin=8.9, xmax=12, ymin=-1, ymax=1)
    
    return

def mass_correction() :
    
    cat_path = 'catalogs/joined_cats/SDSS_gal-info_totlgm_GAMA_GaussFitSimple_StellarMasses_masked-besides-z-0075.fits'
    catalog = Table.read(cat_path)
    
    mask = ( (catalog['KCOR_MAG'][:,0] > 0) &
             (catalog['KCOR_MAG'][:,2] > 0) &
             (catalog['MEDIAN'] >= 6) )
    catalog = catalog[mask]
    
    D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)/u.pc # redshift from SDSS
    
    # SDSS photometry
    M_g_sdss = catalog['KCOR_MAG'][:,0] - 5*np.log10(D_L) + 5 # absolute g mag.
    catalog['M_g_sdss'] = M_g_sdss
    
    M_i_sdss = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
    catalog['M_i_sdss'] = M_i_sdss
    
    # GAMA photometry
    fluxscale = catalog['fluxscale']
    condition = ((fluxscale >= 0.1) & (fluxscale <= 10))
    
    M_g_corrected = catalog['absmag_g'] - 2.5*np.log10(catalog['fluxscale'])
    M_g_gama = np.where(condition, M_g_corrected, catalog['absmag_g'])
    catalog['M_g_gama'] = M_g_gama
    
    M_i_corrected = catalog['absmag_i'] - 2.5*np.log10(catalog['fluxscale'])
    M_i_gama = np.where(condition, M_i_corrected, catalog['absmag_i'])
    catalog['M_i_gama'] = M_i_gama
    
    mask = ( (catalog['M_g_gama'] > -70) & (catalog['M_i_gama'] > -70) )
    catalog = catalog[mask]
    
    from scipy.optimize import curve_fit
    # from scipy.stats import chisquare
    
    # fit1 to g magnitude data
    print('\nUsing linear function with slope of unity')
    popt_int, pcov_int = curve_fit(funcs.intercept, catalog['M_g_sdss'], catalog['M_g_gama'])
    print('b=%s' % str(popt_int[0]))
    
    print('\nUsing linear function with free slope')
    popt_lin, pcov_lin = curve_fit(funcs.line, catalog['M_g_sdss'], catalog['M_g_gama'])
    print('m=%s  b=%s' % (str(popt_lin[0]), str(popt_lin[1]) ))
    
    plt.plot(catalog['M_g_sdss'], r'Absolute $g$ Magnitude from SDSS',
              catalog['M_g_gama'], r'Absolute $g$ Magnitude from GAMA',
              hist2d=True, fit1=True, xmin=-26, xmax=-11, ymin=-26, ymax=-11)
    
    # fit2 to i magnitude data
    print('\nUsing linear function with slope of unity')
    popt_int, pcov_int = curve_fit(funcs.intercept, catalog['M_i_sdss'], catalog['M_i_gama'])
    print('b=%s' % str(popt_int[0]))
    
    print('\nUsing linear function with free slope')
    popt_lin, pcov_lin = curve_fit(funcs.line, catalog['M_i_sdss'], catalog['M_i_gama'])
    print('m=%s  b=%s' % (str(popt_lin[0]), str(popt_lin[1]) ))
    
    plt.plot(catalog['M_i_sdss'], r'Absolute $i$ Magnitude from SDSS',
              catalog['M_i_gama'], r'Absolute $i$ Magnitude from GAMA',
              hist2d=True, fit2=True, xmin=-26, xmax=-11, ymin=-26, ymax=-11)
    
    return

def upper(array) :
    return np.percentile(array, 84)

def lower(array) :
    return np.percentile(array, 16)

# mass_comparison('SDSS', n_bins=25)
# mass_comparison('SDSS', intercept_corr=[-1.2722659777062275, -1.2179114940900249])
# mass_comparison('SDSS', linear_corr=[0.8125759463914866, -4.961460348702913,
#                                       0.8492931852404171, -4.343099872187622], n_bins=23)
# mass_comparison('GAMA', n_bins=23)

# mass_correction()
