# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 4
    ASSIGNMENT: 
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2019-Mar-10
    
    PURPOSE: Plot the parameters from the fitted SEDs to look for correlations.
"""

# imports
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii
from astropy.table import Table

# constants
currentFig = 1

LABELS = {'age':'Age (yr)',
          'irlum':r'Luminosity ($L_\odot)$',
          'optsfr':r'Optical SFR ($M_\odot$ yr$^{-1}$)',
          'mstar':r'Mass of the Stellar Population ($M_\odot$)',
          'lir':r'$L_{IR}$ (erg s$^{-1}$)',
          'irsfr':r'IR SFR ($M_\odot$ yr$^{-1}$)'
          }

#......................................................................plotting
def plotting(xvals, xlab, yvals, ylab) :
#             yarrays, labels, xlab, ylab, location='upper left') :
    """
    This function plots multiple sets of data on the same labelled plot,
    effectively supplied as a list of lists (yarrays). These data sets
    can be either continuous distributions, or discrete values.
    """
    
    global currentFig
    
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111)
    
#    i = 0
#    while i < len(yarrays) : # add each dataset in a list to the plot
#        ax.plot(xvals, yarrays[i], label = labels[i])
#        i += 1
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.errorbar(xvals, yvals, xerr=UNCERTS[xlab],
                yerr=UNCERTS[ylab], fmt='ko', elinewidth=0.3,
                capsize=1.5, errorevery=1)
    
    ax.set_xlabel("%s" % LABELS[xlab], fontsize = 15 )
    ax.set_ylabel("%s" % LABELS[ylab], fontsize = 15 )
    
#    ax.set_xlim(10**39, 10**45)
#    ax.set_ylim(10**39, 10**45)
    
#    plt.legend(loc = location)
    
    plt.tight_layout()
    plt.show()
    
    return
#.....................................................................writedata
def writedata(headerlist, paramlist) :
    
    table = {}
    
    if len(headerlist) == len(paramlist) :
        for i in range(len(paramlist)) :
            table.update({headerlist[i]:paramlist[i]})
    else :
        print("Lists aren't same length.")
        
    newdata = Table( table )
    
    ascii.write(newdata, 'params_to_plot.txt', overwrite=True)
    
    return
#..............................................................end of functions

age, age_lo, age_hi = [], [], []
irlum, irlum_lo, irlum_hi = [], [], []
mstar, mstar_lo, mstar_hi = [], [], []
optsfr, optsfr_lo, optsfr_hi = [], [], []
lir, lir_lo, lir_hi = [], [], []
irsfr, irsfr_lo, irsfr_hi = [], [], []

for i in range(1, 6) :
    filename = 'parameter_outvalues_'
    identifier = str(i)
    filename += identifier
    filename += '.txt'
    dat = ascii.read(filename) # requires columns to have unique names
    
    age.append(np.power(10,dat['age'][2]))
    age_lo.append(np.power(10,dat['age'][2])-np.power(10,dat['age'][1]))
    age_hi.append(np.power(10,dat['age'][3])-np.power(10,dat['age'][2]))
    
    irlum.append(np.power(10,dat['irlum'][2]))
    irlum_lo.append(np.power(10,dat['irlum'][2])-np.power(10,dat['irlum'][1]))
    irlum_hi.append(np.power(10,dat['irlum'][3])-np.power(10,dat['irlum'][2]))
    
    mstar.append(np.power(10,dat['log_Mstar'][2]))
    mstar_lo.append(np.power(10,dat['log_Mstar'][2])-np.power(10,dat['log_Mstar'][1]))
    mstar_hi.append(np.power(10,dat['log_Mstar'][3])-np.power(10,dat['log_Mstar'][2]))
    
    optsfr.append(dat['SFR_opt'][2])
    optsfr_lo.append(dat['SFR_opt'][2]-dat['SFR_opt'][1])
    optsfr_hi.append(dat['SFR_opt'][3]-dat['SFR_opt'][2])
    
    lir.append(np.power(10,dat['LIR(8-1000)'][2]))
    lir_lo.append(np.power(10,dat['LIR(8-1000)'][2])-np.power(10,dat['LIR(8-1000)'][1]))
    lir_hi.append(np.power(10,dat['LIR(8-1000)'][3])-np.power(10,dat['LIR(8-1000)'][2]))
    
    irsfr.append(dat['SFR_IR'][2])
    irsfr_lo.append(dat['SFR_IR'][2]-dat['SFR_IR'][1])
    irsfr_hi.append(dat['SFR_IR'][3]-dat['SFR_IR'][2])

UNCERTS = {'age':[age_lo, age_hi],
           'irlum':[irlum_lo, irlum_hi],
           'optsfr':[optsfr_lo, optsfr_hi],
           'mstar':[mstar_lo, mstar_hi],
           'lir':[lir_lo, lir_hi],
           'irsfr':[irsfr_lo, irsfr_hi]
           }

# plots
#plotting(age, 'age', irlum, 'irlum')
#plotting(age, 'age', mstar, 'mstar')
#plotting(age, 'age', optsfr, 'optsfr')
#plotting(age, 'age', lir, 'lir')
#plotting(age, 'age', irsfr, 'irsfr')

#plotting(irlum, 'irlum', mstar, 'mstar')
#plotting(irlum, 'irlum', optsfr, 'optsfr')
#plotting(irlum, 'irlum', lir, 'lir')
#plotting(irlum, 'irlum', irsfr, 'irsfr')

#plotting(mstar, 'mstar', optsfr, 'optsfr')
#plotting(mstar, 'mstar', lir, 'lir')
#plotting(mstar, 'mstar', irsfr, 'irsfr')

#plotting(optsfr, 'optsfr', lir, 'lir')
#plotting(optsfr, 'optsfr', irsfr, 'irsfr')

#plotting(lir, 'lir', irsfr, 'irsfr') # irsfr = 3.88e-44*lir, see Calistro-Rivera paper

# write the data to a file
headers = ['age','age_lo','age_hi','irlum','irlum_lo','irlum_hi',
           'optsfr','optsfr_lo','optsfr_hi','mstar','mstar_lo','mstar_hi',
           'lir','lir_lo','lir_hi','irsfr','irsfr_lo','irsfr_hi']
params = [age, age_lo, age_hi, irlum, irlum_lo, irlum_hi,
          optsfr, optsfr_lo, optsfr_hi, mstar, mstar_lo, mstar_hi,
          lir, lir_lo, lir_hi, irsfr, irsfr_lo, irsfr_hi]
writedata(headers, params)

