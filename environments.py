
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import Table

c = 299792.458 # km/s

# the following are all in SDSS as per Bernd
galaxies = {'HE0040-1105':{'RA':'00:42:36.860','dec':'-10:49:22.03','z':0.041962},
           'RBS175':{'RA':'01:17:03.587','dec':'+00:00:27.41','z':0.045605},
           'Mrk1503':{'RA':'01:21:59.827','dec':'-01:02:24.08','z':0.054341},
           'Mrk1018':{'RA':'02:06:15.990','dec':'-00:17:29.20','z':0.042436},
           'Mrk1044':{'RA':'02:30:05.525','dec':'-08:59:53.29','z':0.016451},
           'Mrk1048':{'RA':'02:34:37.769','dec':'-08:47:15.44','z':0.043143},
           'HE0345+0056':{'RA':'03:47:40.188','dec':'+01:05:14.02','z':0.031000},
           'HE0853+0102':{'RA':'08:55:54.268','dec':'+00:51:10.60','z':0.052000}, # also in GAMA
           'Mrk707':{'RA':'09:37:01.030','dec':'+01:05:43.48','z':0.050338},
           'HE2222-0026':{'RA':'22:24:35.292','dec':'-00:11:03.89','z':0.059114},
           'Mrk926':{'RA':'23:04:43.478','dec':'-08:41:08.62','z':0.046860},
           'Mrk590':{'RA':'02:14:33.562','dec':'-00:46:00.09','z':0.026385}
        }

IDs = list( galaxies.keys() )
ras = Angle([ galaxy['RA'] for galaxy in galaxies.values() ], u.hour)
decs = Angle([ galaxy['dec'] for galaxy in galaxies.values() ], u.deg)
zs = [ galaxy['z'] for galaxy in galaxies.values() ]

lower_z = np.array(zs) - 1500/c
upper_z = np.array(zs) + 1500/c

cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
dists = cosmo.angular_diameter_distance(zs) # compute the angular size distances
radii = (2*u.Mpc/dists)*(180*60*u.arcmin/np.pi) # find the radii corresponding to 2 Mpc

#...........................................................................SDSS
def SDSS(i) :
    
    # open the SDSS catalog and populate relevant information (ie. RA+Dec, z)
    SDSS_catalog = fits.open('SDSS_-_gal_info_dr7_v5_2.fit.gz') # SDSS DR7, per Yjan
    #info = SDSS_catalog.info()
    header = SDSS_catalog[1].header
    #print(header) # to see what is in the actual data table
    data = SDSS_catalog[1].data
    RAs = data.field('RA') # get all the RAs, Decs, redshifts, etc.
    Decs = data.field('DEC')
    redshifts = data.field('Z')
    z_errs = data.field('Z_ERR') # we only want z/z_err >= 3 for reliability
    SDSS_catalog.close() # close the SDSS catalog fits file
    
    badIndex = np.where(Decs==-9999) # have to remove the one RA and Dec value of -9999.0
    RAs = np.delete(RAs, badIndex)*u.deg # remove the bad value for all parameters
    Decs = np.delete(Decs, badIndex)*u.deg
    redshifts = np.delete(redshifts, badIndex)
    z_errs = np.delete(z_errs, badIndex)
    
    # select which galaxy to investigate
    catname = 'SDSS'
    
    print('\n-----ANALYSIS FOR {0:s}-----\n'.format( IDs[i] ) )
    center = SkyCoord(ra=ras[i], dec=decs[i], distance=dists[i]) # galaxy of interest
    
    distances = cosmo.angular_diameter_distance(redshifts) # compute the angular size distances
    catalog = SkyCoord(ra=RAs, dec=Decs, distance=distances,
                       unit=(u.deg, u.deg, u.Mpc) ) # create the catalog
    
    d2d = center.separation(catalog) # find the projected separations on the sky (ie. 2D)
    mask = ( (lower_z[i] <= redshifts) & (redshifts <= upper_z[i] )
                & ( d2d <= radii[i] ) & ((redshifts / z_errs) >= 3) ) # mask the data
            # based on redshift reliability, velocity cuts, 2D separation
    catalog = catalog[mask]
    
    print('Reference Galaxy: {}\n'.format(IDs[i]) )
    
    print('Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 ' +
          'Mpc\u207B\u00B9, \u03A9\u2098 = 0.3\n')
    
    print('Velocity between: {0:.0f} and {1:.0f} km/s\n'.format(
            lower_z[i]*c, upper_z[i]*c) )
    
    print('{0:g} objects found in {1:s} within {2:.3f} of {3:s}\n'.format(
            len(catalog), catname, radii[i], IDs[i]) )
    
    table = Table([Angle(RAs[mask]).to_string(unit=u.hour),
                   Angle(Decs[mask]).to_string(unit=u.degree),
                   redshifts[mask],
                   redshifts[mask]*c*u.km/u.s,
                   distances[mask],
                   d2d[mask]*60*u.arcmin/u.deg],
                  names=('RA','Dec','z','Velocity','D_A','Separation') )
    table.add_row( [ras[i].to_string(unit=u.hour)+'*',
            decs[i].to_string(unit=u.degree), zs[i], zs[i]*c,
            dists[i]/u.Mpc, 0.0] ) # add the galaxy of interest to the table
    table['RA'].format = '15s'
    table['z'].format = '10.6f'
    table['Velocity'].format = '8.0f'
    table['D_A'].format = '8.2f'
    table['Separation'].format = '10.3f'
    table.sort('Separation') # sort the table by the separation
    table.pprint(max_lines=-1) # to print full table, use print(table) for regular
    print("\nNote: '*' in the first line denotes the object of interest.")
    
    return

#...........................................................................GAMA
def GAMA(i) :
    
    
    
    
    
    return
#...............................................................end of functions

SDSS(8)

'''
#centers = SkyCoord(ra=ras, dec=decs, distance=dists,
#                   unit=(u.hour, u.deg, u.Mpc) ) # for full sample
'''
"""
#..........................................................................table
def coords(galaxy) :
    print('%-13s %14s %14s %10s' %
          (galaxy, galaxies[galaxy]['RA'], galaxies[galaxy]['dec'],
           galaxies[galaxy]['z']) )
    return

#..........................................................................table
def coordtable() :
    print('Galaxy               RA            Dec           z\n')
    gals = list( galaxies.keys() )
    newtable = np.vectorize(coords)
    newtable(gals)
    return

#...............................................................................
coordtable()
"""
