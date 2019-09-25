
import numpy as np
from astropy.io import ascii
#import astropy.units as u
import matplotlib.pyplot as plt
import scipy.stats as sp

dat = ascii.read('cars_catalog.csv')

name, alt_name, RA, Dec = dat['name'], dat['alt_name'], dat['RA'], dat['Dec']
zz, large_radii = dat['z'], dat['2Mpc_radius']
H2_mass = np.power(10, dat['log_M(H2)'])
HI_mass = np.power(10, dat['log_M(HI)'])

# would like masses of CARS galaxies themselves, and star formation rates

currentFig = 1

LABELS = {
          'RA':'Right Ascension (deg)',
          'Dec':'Declination (deg)',
          'zz':'Redshift',
          'large_radii':'Radius corresponding to 2 Mpc projected (arcmin)',
          'H2_mass':'Mass of Molecular Hydrogen ($M_\odot$)',
          'HI_mass':'Mass of Atomic Hydrogen ($M_\odot$)'
          }

def plot(x_param, x_label, y_param, y_label, linear=False) :
    
    global currentFig
    spear = sp.spearmanr(x_param, y_param, nan_policy='omit')
    print("Figure %2.1d   %13s vs %-13s   Spearman: %8.3g   pvalue: %8.2g" % 
        (currentFig, y_label, x_label, spear[0], spear[1]) )
    
    fig = plt.figure(currentFig)  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111) # set axes, figure location
    
    if linear == True :
        ax.plot(x_param, y_param, 'o')
    else :
        ax.loglog(x_param, y_param, 'o')
    
    ax.set_xlabel("%s" % LABELS[x_label], fontsize = 15 )
    ax.set_ylabel("%s" % LABELS[y_label], fontsize = 15 )
    
    plt.tight_layout()
    plt.show()
    
    return

#plot(zz, 'zz', RA, 'RA', linear=True)
#plot(zz, 'zz', Dec, 'Dec', linear=True)
#plot(H2_mass, 'H2_mass', HI_mass, 'HI_mass')
#plot(zz, 'zz', H2_mass, 'H2_mass')
#plot(zz, 'zz', HI_mass, 'HI_mass')


