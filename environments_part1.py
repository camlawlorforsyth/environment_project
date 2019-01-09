# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 1
    ASSIGNMENT: Search for physically close companions to CARS host galaxies
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2018-Dec-20
    
    PURPOSE: Search for physically close companion objects to CARS host
             galaxies, within 2 Mpc projected, and +/-1500 km/s along the LOS.
"""

# imports
import numpy as np

import astropy.constants as const
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore") # ignore warnings about division by 0 when
# taking z/z_err >= 3, and ignore integration warning for the cosmology

# the following are all in SDSS as per Bernd
galaxies = {
           'HE0040-1105':{'RA':'10.65358542','dec':'-10.82278903','z':0.04199},
           'RBS175':{'RA':'19.26495342','dec':'7.63966751E-3','z':0.04563},
           'Mrk1503':{'RA':'20.49924553','dec':'-1.04009986','z':0.05435},
           'Mrk1018':{'RA':'31.56660166','dec':'-0.29145178','z':0.04298},
           'Mrk1044':{'RA':'37.52302517','dec':'-8.99810955','z':0.016451},
               # in GAMA?, no DR7,8 spectrum
           'Mrk1048':{'RA':'38.65766625','dec':'-8.78779752','z':0.043143},
               # in GAMA?, no DR7,8 spectrum
           'HE0345+0056':{'RA':'56.9174638','dec':'1.0872126','z':0.031000},
               # no DR7,8 spectrum
           'HE0853+0102':{'RA':'133.97614006','dec':'0.85306112','z':0.05247,
                          'GAMA_CATAID':278841}, # in GAMA
           'Mrk707':{'RA':'144.25436933','dec':'1.0954803','z':0.05025},
           'HE2222-0026':{'RA':'336.14705483','dec':'-0.18442192','z':0.05873},
           'Mrk926':{'RA':'346.18116195','dec':'-8.68573563','z':0.04702},
           'Mrk590':{'RA':'33.63983722','dec':'-0.76672072','z':0.02609}
               # no DR7 spectrum
        } # coordinates from SDSS DR8 catalog, z from SDSS DR8, otherwise NED
#galaxies = {'HE0040-1105':{'RA':'00:42:36.860','dec':'-10:49:22.03','z':0.041962},
#           'RBS175':{'RA':'01:17:03.587','dec':'+00:00:27.41','z':0.045605},
#           'Mrk1503':{'RA':'01:21:59.827','dec':'-01:02:24.08','z':0.054341},
#           'Mrk1018':{'RA':'02:06:15.990','dec':'-00:17:29.20','z':0.042436},
#           'Mrk1044':{'RA':'02:30:05.525','dec':'-08:59:53.29','z':0.016451},
                # in GAMA?
#           'Mrk1048':{'RA':'02:34:37.769','dec':'-08:47:15.44','z':0.043143},
                # in GAMA?
#           'HE0345+0056':{'RA':'03:47:40.188','dec':'+01:05:14.02','z':0.031000},
#           'HE0853+0102':{'RA':'08:55:54.268','dec':'+00:51:10.60','z':0.052000},
                # in GAMA
#           'Mrk707':{'RA':'09:37:01.030','dec':'+01:05:43.48','z':0.050338},
#           'HE2222-0026':{'RA':'22:24:35.292','dec':'-00:11:03.89','z':0.059114},
#           'Mrk926':{'RA':'23:04:43.478','dec':'-08:41:08.62','z':0.046860},
#           'Mrk590':{'RA':'02:14:33.562','dec':'-00:46:00.09','z':0.026385}
#        } # coordinates/z from NED, note that 'ras' must be changed to u.hour

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
currentFig = 1 # first figure will be numbered as 'Figure 1'

#..........................................................................main
def main(catalog, index, sample_table=False, printtable=False,
         return_vals=False) :
    
    IDs = list( galaxies.keys() ) # the object name/identifier
    ras = Angle([ galaxy['RA'] for galaxy in galaxies.values() ], u.deg)
    decs = Angle([ galaxy['dec'] for galaxy in galaxies.values() ], u.deg)
    zs = np.array([ galaxy['z'] for galaxy in galaxies.values() ])
    
    lower_z = np.array(zs) - 1500*(u.km/u.s)/const.c.to('km/s')
    upper_z = np.array(zs) + 1500*(u.km/u.s)/const.c.to('km/s')
    
    dists = cosmo.angular_diameter_distance(zs) # compute D_A
    radii = (2*u.Mpc/dists)*(180*60*u.arcmin/np.pi) # find the radii = 2 Mpc
    
    if (sample_table == True) :
        
        print('\nSample Properties\n')
    
        print('Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 ' +
              'Mpc\u207B\u00B9, \u03A9\u2098 = 0.3\n')
        
        table = Table([IDs, ras.to_string(unit=u.hour),
                       decs.to_string(unit=u.degree), zs,
                       zs*const.c.to('km/s'),
                       dists, radii],
                      names=('Object Name','RA','Dec','z','Velocity','D_A',
                             '2 Mpc Radius') )
        table['Object Name'].format = '14s'
        table['RA'].format = '13s'
        table['z'].format = '10.6f'
        table['Velocity'].format = '8.0f'
        table['D_A'].format = '8.2f'
        table['2 Mpc Radius'].format = '12.3f'
#        table.sort('z')
        table.pprint(max_lines=-1, max_width=-1)
    
    if (catalog == 'SDSS') :
        catname = 'SDSS'
        
        # open the SDSS catalog and populate relevant information
#        SDSS_catalog = fits.open('gal_info_dr7_v5_2.fit.gz') # SDSS DR7
        SDSS_catalog = fits.open('galSpecInfo-dr8.fits') # SDSS DR8
#        info = SDSS_catalog.info()
#        header = SDSS_catalog[1].header
#        print(header) # to see what is in the actual data table
        data = SDSS_catalog[1].data
        RAs = data.field('RA')*u.deg # get all the RAs, Decs, redshifts, etc.
        Decs = data.field('DEC')*u.deg
        redshifts = data.field('Z')
        z_errs = data.field('Z_ERR') # we only want z/z_err >= 3
        SDSS_catalog.close() # close the SDSS catalog fits file
        
        badIndex = np.where(Decs==-9999.0*u.deg) # the one bad RA/Dec value
        RAs = np.delete(RAs, badIndex) # remove bad value for all parameters
        Decs = np.delete(Decs, badIndex)
        redshifts = np.delete(redshifts, badIndex)
        z_errs = np.delete(z_errs, badIndex)
        
        mask = ( (lower_z[index] <= redshifts) & (redshifts <= upper_z[index] )
                & ((redshifts / z_errs) >= 3) ) # mask the data based on
                # velocity cuts, redshift quality
        
        if return_vals == True :
            catlen = cat_search(IDs[index], ras[index],decs[index], zs[index],
                                dists[index], lower_z[index], upper_z[index],
                                radii[index], RAs, Decs, redshifts, catname,
                                mask, printtable=printtable,
                                return_vals=return_vals)
        else :
            cat_search(IDs[index], ras[index],decs[index], zs[index],
                       dists[index], lower_z[index], upper_z[index],
                       radii[index], RAs, Decs, redshifts, catname, mask,
                       printtable=printtable)
        
    else :
        catname = 'GAMA'
        
        # open the GAMA catalog and populate relevant information
        GAMA_catalog = fits.open('GaussFitSimple.fits') # SpecLineSFRv05
#        info = GAMA_catalog.info()
#        header = GAMA_catalog[1].header
#        print(header)
        data = GAMA_catalog[1].data
        RAs = data.field('RA')*u.deg
        Decs = data.field('DEC')*u.deg
        redshifts = data.field('Z')
        redshift_quality = data.field('NQ')
        GAMA_catalog.close() # close the GAMA catalog fits file
        
        mask = ( (lower_z[index] <= redshifts) & (redshifts <= upper_z[index] )
                & (redshift_quality >= 3) ) # mask the data based on
                # velocity cuts, redshift quality
        
        if return_vals == True :
            catlen = cat_search(IDs[index], ras[index],decs[index], zs[index],
                                dists[index], lower_z[index], upper_z[index],
                                radii[index], RAs, Decs, redshifts, catname,
                                mask, printtable=printtable,
                                return_vals=return_vals)
        else :
            cat_search(IDs[index], ras[index],decs[index], zs[index],
                       dists[index], lower_z[index], upper_z[index],
                       radii[index], RAs, Decs, redshifts, catname, mask,
                       printtable=printtable)
    
    if return_vals == True :
        close = cosmo.angular_diameter_distance(lower_z[index])
        far = cosmo.angular_diameter_distance(upper_z[index])
        volume = np.pi*(far-close)*(2*u.Mpc)**2
        return zs[index], catlen, volume
    else :
        return

#....................................................................cat_search
def cat_search(galaxy_ID, RA_c, Dec_c, zs_c, dists_c, low_z, high_z, radius,
               RAs, Decs, redshifts, catname, mask, printtable=False,
               return_vals=False) :
    
    center = SkyCoord(ra=RA_c, dec=Dec_c, distance=dists_c) #galaxy of interest
    
    distances = cosmo.angular_diameter_distance(redshifts) # compute D_A
    catalog = SkyCoord(ra=RAs, dec=Decs, distance=distances,
                       unit=(u.deg, u.deg, u.Mpc) ) # create the catalog
    
    d2d = center.separation(catalog) # find projected separations on the sky
    mask = mask & ( d2d <= radius ) # mask the data based on 2D separation
    catalog = catalog[mask]
    
    if printtable == True :
        print('\n-----ANALYSIS FOR {0:s}-----\n'.format(galaxy_ID) )
        
        print('Reference Galaxy: {}\n'.format(galaxy_ID) )
        
        print('Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 ' +
              'Mpc\u207B\u00B9, \u03A9\u2098 = 0.3\n')
        
        print('Velocity between: {0:.0f} and {1:.0f} km/s\n'.format(
                low_z*const.c.to('km/s'), high_z*const.c.to('km/s') ) )
        
        print('{0:g} objects found in {1:s} within {2:.3f} of {3:s}\n'.format(
                len(catalog), catname, radius, galaxy_ID) )
    
    table = Table([Angle(RAs[mask], u.deg), #.to_string(unit=u.hour),
                   Angle(Decs[mask], u.deg), #.to_string(unit=u.degree),
                   redshifts[mask],
                   redshifts[mask]*const.c.to('km/s'),distances[mask],
                   d2d[mask]*60*u.arcmin/u.deg],
                  names=('RA','Dec','z','Velocity','D_A','Separation') )
#    table.add_row( [RA_c.to_string(unit=u.hour)+'*',
#            Dec_c.to_string(unit=u.degree), zs_c, zs_c*const.c.to('km/s'),
#            dists_c/u.Mpc, 0.0] ) # add the galaxy of interest to the table
#    table['RA'].format = '15s' # only for h:m:s
    table['z'].format = '10.6f'
    table['Velocity'].format = '8.0f'
    table['D_A'].format = '8.2f'
    table['Separation'].format = '10.3f'
    table.sort('Separation') # sort the table by the separation
    
    if printtable == True :
        table.pprint(max_lines=-1) # print full table, use print(table) for reg
#    print("\nNote: '*' in the first line denotes the object of interest.")
    
    if return_vals == True :
        return len(catalog)-1
    else :
        return

#..................................................................density_plot
def density_plot(catalog) :
    
    if catalog == 'both' :
        redshifts = []
        SDSS_companions = []
        SDSS_volumes = []
        GAMA_companions = []
        GAMA_volumes = []
        for i in range(0, 12) :
            redshift, number, volume = main('SDSS', i, return_vals=True)
            redshifts.append(redshift)
            SDSS_companions.append(number)
            SDSS_volumes.append(volume.value)
            if i in [4,5,7] :
                redshift, number, volume = main('GAMA', i, return_vals=True)
                GAMA_companions.append(number)
                GAMA_volumes.append(volume.value)
            else :
                GAMA_companions.append(np.nan)
                GAMA_volumes.append(np.nan)
        multi(redshifts, 'Redshift',
              np.array(SDSS_companions)/np.array(SDSS_volumes),
              np.array(GAMA_companions)/np.array(GAMA_volumes),
              'Companions Density (Mpc$^{-3}$)')
    
    redshifts = []
    number_of_companions = []
    volumes = []
    
    if catalog == 'SDSS' :
        for i in range(0, 12) :
            redshift, number, volume = main('SDSS', i, return_vals=True)
            redshifts.append(redshift)
            number_of_companions.append(number)
            volumes.append(volume.value)
        plot(redshifts, 'Redshift',
             np.array(number_of_companions)/np.array(volumes),
             'Companions Density (Mpc$^{-3}$)')
    
    if catalog == 'GAMA' :
        for i in [4,5,7] :
            redshift, number, volume = main('GAMA', i, return_vals=True)
            redshifts.append(redshift)
            number_of_companions.append(number)
            volumes.append(volume.value)
        plot(redshifts, 'Redshift',
             np.array(number_of_companions)/np.array(volumes),
             'Companions Density (Mpc$^{-3}$)')
    
    return

#.........................................................................multi
def multi(xvals, xlab, yvals1, yvals2, ylab, xmin=None, xmax=None, ymin=None,
          ymax=None, location='upper left') :
    
    global currentFig        
    fig = plt.figure(currentFig)  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111)
    
    ax.plot(xvals, yvals1, 'ko', label = "%s" % 'SDSS', zorder=2)
#    ax.axhline(np.nanmean(yvals1), color='k', linestyle='-', label='mean')
    ax.axhline(np.nanmedian(yvals1), color='k', linestyle='--', label='median',
               zorder=1)
    
    ax.plot(xvals, yvals2, 'ro', label = "%s" % 'GAMA', zorder=2)
#    ax.axhline(np.nanmean(yvals2), color='r', linestyle='-', label='mean')
    ax.axhline(np.nanmedian(yvals2), color='r', linestyle='--', label='median',
               zorder=1)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel("%s" % xlab, fontsize = 15 )
    ax.set_ylabel("%s" % ylab, fontsize = 15 )
    
    plt.legend(loc = location)
    
    plt.tight_layout()
    plt.show()
    
    return

#..........................................................................plot
def plot(xvals, xlab, yvals, ylab, xmin=None, xmax=None, ymin=None,
         ymax=None) : #, logx=False, logy=False, linear=False) :
    
    global currentFig
    fig = plt.figure(currentFig)  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111) # set axes, figure location
    
#    if (logx == True) and (logy == False) and (linear == False) :
#        ax.semilogx(xvals, yvals, 'ko')
#    elif (logx == False) and (logy == True) and (linear == False) :
#        ax.semilogy(xvals, yvals, 'ko')
#    elif (logx == False) and (logy == False) and (linear == True) :
    ax.plot(xvals, yvals, 'ko')
#    elif (logx == True) and (logy == True) and (linear == False) :
#        ax.loglog(xvals, yvals, 'ko')
#    else :
#        ax.loglog(xvals, yvals, 'ko')
            
    ax.set_xlabel("%s" % xlab, fontsize = 15 )
    ax.set_ylabel("%s" % ylab, fontsize = 15 )
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.tight_layout()
    plt.show()
    
    return
#..............................................................end of functions

#main('SDSS', 0, printtable=True, sample_table=True)

#for i in range(1, 12) :
#    main('SDSS', i, printtable=True)

#for i in [4,5,7] :
#    main('GAMA', i)
