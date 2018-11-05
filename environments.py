# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 1
    ASSIGNMENT: Search for physically close companions to CARS host galaxies
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2018-Nov-4
    
    PURPOSE: Search for physically close companion objects to CARS host
             galaxies, within 2 Mpc projected, and +/-1500 km/s along the LOS.
"""

# imports
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np

# the following are all in SDSS as per Bernd
galaxies = {'HE0040-1105':{'RA':'10.65358672','dec':'-10.82278912','z':0.041962},
           'RBS175':{'RA':'19.26494998','dec':'7.61312E-3','z':0.045605},
           'Mrk1503':{'RA':'20.49921832','dec':'-1.04010828','z':0.054341},
           'Mrk1018':{'RA':'31.56661419','dec':'-0.29144322','z':0.042436},
           'Mrk1044':{'RA':'37.52302021','dec':'-8.99813544','z':0.016451}, # in GAMA?
           'Mrk1048':{'RA':'38.6576586','dec':'-8.78777681','z':0.043143}, # in GAMA?
           'HE0345+0056':{'RA':'56.91745592','dec':'1.08722631','z':0.031000},
           'HE0853+0102':{'RA':'133.97612964','dec':'0.85305195','z':0.052000}, # in GAMA
           'Mrk707':{'RA':'144.25436927','dec':'1.09547822','z':0.050338},
           'HE2222-0026':{'RA':'336.14705331','dec':'-0.18441524','z':0.059114},
           'Mrk926':{'RA':'346.18116153','dec':'-8.68572479','z':0.046860},
           'Mrk590':{'RA':'33.63982968','dec':'-0.76672567','z':0.026385}
        } # coordinates from SDSS DR7 catalog
#galaxies = {'HE0040-1105':{'RA':'00:42:36.860','dec':'-10:49:22.03','z':0.041962},
#           'RBS175':{'RA':'01:17:03.587','dec':'+00:00:27.41','z':0.045605},
#           'Mrk1503':{'RA':'01:21:59.827','dec':'-01:02:24.08','z':0.054341},
#           'Mrk1018':{'RA':'02:06:15.990','dec':'-00:17:29.20','z':0.042436},
#           'Mrk1044':{'RA':'02:30:05.525','dec':'-08:59:53.29','z':0.016451}, # in GAMA?
#           'Mrk1048':{'RA':'02:34:37.769','dec':'-08:47:15.44','z':0.043143}, # in GAMA?
#           'HE0345+0056':{'RA':'03:47:40.188','dec':'+01:05:14.02','z':0.031000},
#           'HE0853+0102':{'RA':'08:55:54.268','dec':'+00:51:10.60','z':0.052000}, # in GAMA
#           'Mrk707':{'RA':'09:37:01.030','dec':'+01:05:43.48','z':0.050338},
#           'HE2222-0026':{'RA':'22:24:35.292','dec':'-00:11:03.89','z':0.059114},
#           'Mrk926':{'RA':'23:04:43.478','dec':'-08:41:08.62','z':0.046860},
#           'Mrk590':{'RA':'02:14:33.562','dec':'-00:46:00.09','z':0.026385}
#        } # coordinates from NED, note that the ras must be changed to u.hour

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
c = 299792.458 # km/s

#...........................................................................main
def main(catalog, index, sample_table=False) :
    
    IDs = list( galaxies.keys() ) # the object name/identifier
    ras = Angle([ galaxy['RA'] for galaxy in galaxies.values() ], u.deg)
    decs = Angle([ galaxy['dec'] for galaxy in galaxies.values() ], u.deg)
    zs = np.array([ galaxy['z'] for galaxy in galaxies.values() ])
    
    lower_z = np.array(zs) - 1500/c
    upper_z = np.array(zs) + 1500/c
    
    dists = cosmo.angular_diameter_distance(zs) # compute D_A
    radii = (2*u.Mpc/dists)*(180*60*u.arcmin/np.pi) # find the radii = 2 Mpc
    
    if (sample_table == True) :
        
        print('\nSample Properties\n')
    
        print('Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 ' +
              'Mpc\u207B\u00B9, \u03A9\u2098 = 0.3\n')
        
        table = Table([IDs, ras.to_string(unit=u.hour),
                       decs.to_string(unit=u.degree), zs, zs*c*u.km/u.s,
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
        SDSS_catalog = fits.open('gal_info_dr7_v5_2.fit.gz') # SDSS DR7
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
        RAs = np.delete(RAs, badIndex) # remove the bad value for all parameters
        Decs = np.delete(Decs, badIndex)
        redshifts = np.delete(redshifts, badIndex)
        z_errs = np.delete(z_errs, badIndex)
        
        mask = ( (lower_z[index] <= redshifts) & (redshifts <= upper_z[index] )
                & ((redshifts / z_errs) >= 3) ) # mask the data based on
                # velocity cuts, redshift quality
        
        cat_search(IDs[index], ras[index], decs[index], zs[index], dists[index],
                   lower_z[index], upper_z[index], radii[index],
                   RAs, Decs, redshifts, catname, mask)
        
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
        
        cat_search(IDs[index], ras[index], decs[index], zs[index], dists[index],
                   lower_z[index], upper_z[index], radii[index],
                   RAs, Decs, redshifts, catname, mask)
        
    return

#.....................................................................cat_search
def cat_search(galaxy_ID, RA_c, Dec_c, zs_c, dists_c, low_z, high_z, radius,
               RAs, Decs, redshifts, catname, mask) :
    
    print('\n-----ANALYSIS FOR {0:s}-----\n'.format(galaxy_ID) )
    center = SkyCoord(ra=RA_c, dec=Dec_c, distance=dists_c) # galaxy of interest
    
    distances = cosmo.angular_diameter_distance(redshifts) # compute D_A
    catalog = SkyCoord(ra=RAs, dec=Decs, distance=distances,
                       unit=(u.deg, u.deg, u.Mpc) ) # create the catalog
    
    d2d = center.separation(catalog) # find the projected separations on the sky
    mask = mask & ( d2d <= radius ) # mask the data based on 2D separation
    catalog = catalog[mask]
    
    print('Reference Galaxy: {}\n'.format(galaxy_ID) )
    
    print('Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 ' +
          'Mpc\u207B\u00B9, \u03A9\u2098 = 0.3\n')
    
    print('Velocity between: {0:.0f} and {1:.0f} km/s\n'.format(
            low_z*c, high_z*c) )
    
    print('{0:g} objects found in {1:s} within {2:.3f} of {3:s}\n'.format(
            len(catalog), catname, radius, galaxy_ID) )
    
    table = Table([Angle(RAs[mask], u.deg).to_string(unit=u.hour),
                   Angle(Decs[mask], u.deg).to_string(unit=u.degree),
                   redshifts[mask], redshifts[mask]*c*u.km/u.s, distances[mask],
                   d2d[mask]*60*u.arcmin/u.deg],
                  names=('RA','Dec','z','Velocity','D_A','Separation') )
#    table.add_row( [RA_c.to_string(unit=u.hour)+'*',
#            Dec_c.to_string(unit=u.degree), zs_c, zs_c*c,
#            dists_c/u.Mpc, 0.0] ) # add the galaxy of interest to the table
    table['RA'].format = '15s'
    table['z'].format = '10.6f'
    table['Velocity'].format = '8.0f'
    table['D_A'].format = '8.2f'
    table['Separation'].format = '10.3f'
    table.sort('Separation') # sort the table by the separation
    table.pprint(max_lines=-1) # print full table, use print(table) for regular
    print("\nNote: '*' in the first line denotes the object of interest.")
    
    return
#...............................................................end of functions

main('SDSS', 0)
