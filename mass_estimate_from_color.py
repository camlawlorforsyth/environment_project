
# imports
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u

import functions as funcs
import plots as plt

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
currentFig = 1

def mass_comparison(cat_name) :
    
    if (cat_name == 'SDSS') :
        cat_path = 'catalogs/SDSS_gal_info_Mstar_SFR.fits'
    if (cat_name == 'GAMA') :
        cat_path = 'catalogs/GAMA_GFS_Mstar.fits'
    
    catalog = Table.read(cat_path)
    
    if (cat_name == 'SDSS') :
        test = catalog[(catalog['RA'] >= 0) & (catalog['Z'] > 0) &
                       (catalog['KCOR_MAG'][:,0] > 0) &
                       (catalog['KCOR_MAG'][:,2] > 0) &
                       (catalog['MEDIAN_2'] > 0)]
        
        D_L = cosmo.luminosity_distance(test['Z']).to(u.pc)/u.pc
        M_i = test['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        g_i_Mstar = (1.15 + 0.7*(test['KCOR_MAG'][:,0] -
                                 test['KCOR_MAG'][:,2]) - 0.4*M_i)
        # based on relation from Taylor+ 2011, MNRAS, 418, 1587
        
        funcs.plot(g_i_Mstar, r'$\log(M_*/M_\odot)$ = 1.15 + 0.70($g-i$) - 0.4$M_i$',
             test['MEDIAN_2'], r'$M_*$ from MPA/JHU [$\log(M_\odot)$]',
             cat_name, hist2d=True, xmin=4, xmax=12, ymin=6, ymax=12)
    
    if (cat_name == 'GAMA') :
        test = catalog[(catalog['logmstar'] > 0) & (catalog['absmag_i'] > -30)]
        
        g_i_Mstar = 1.15 + 0.7*test['gminusi'] - 0.4*test['absmag_i']
        # based on relation from Taylor+ 2011, MNRAS, 418, 1587
        
        plt.plot(g_i_Mstar, r'$\log(M_*/M_\odot)$ = 1.15 + 0.70($g-i$) - 0.4$M_i$',
             test['logmstar'], r'$M_*$ from GAMA [$\log(M_\odot)]$',
             cat_name, hist2d=True, nbins=50, xmin=6, xmax=13, ymin=6, ymax=13)
    
    return

# mass_comparison('SDSS')
# mass_comparison('GAMA')
