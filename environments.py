
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import Table

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

ras = [ galaxy['RA'] for galaxy in galaxies.values() ]
decs = [ galaxy['dec'] for galaxy in galaxies.values() ]
zs = [ galaxy['z'] for galaxy in galaxies.values() ]

cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
dists = cosmo.angular_diameter_distance(zs)

SDSS_catalog = fits.open('SDSS_-_gal_info_dr7_v5_2.fit.gz')
#info = SDSS_catalog.info()
header = SDSS_catalog[1].header
#print(header) # to see what is in the actual data table
data = SDSS_catalog[1].data
RAs_all = data.field('RA')
Decs_all = data.field('DEC')
redshifts_all = data.field('Z')
z_errs_all = data.field('Z_ERR') # we only want z/z_err >= 3 for reliability
SDSS_catalog.close() # close the SDSS catalog fits file

mask = (redshifts_all / z_errs_all) >= 3 # use this mask for reliability
RAs = RAs_all[mask] # and apply it to all the data of interest
Decs = Decs_all[mask]
redshifts = redshifts_all[mask]
z_errs = z_errs_all[mask]




#centers = SkyCoord(ra=ras, dec=decs, distance=dists,
#                   unit=(u.hour, u.deg, u.Mpc) ) # for full sample

#print('-----ANALYSIS FOR {}-----'.format( list(galaxies.keys())[0] ) )
center = SkyCoord(ra=ras[0], dec=decs[0], distance=dists[0],
                  unit=(u.hour, u.deg, u.Mpc) ) # for first
#print(center)




badIndex = np.where(Decs==-9999)
RAs = np.delete(RAs, badIndex)
Decs = np.delete(Decs, badIndex)
redshifts = np.delete(redshifts, badIndex)
distances = cosmo.angular_diameter_distance(redshifts) # compute the angular size distances

catalog = SkyCoord(ra=RAs, dec=Decs, distance=distances,
                   unit=(u.deg, u.deg, u.Mpc) )

#idx, idxcatalog, d2d, d3d = 



'''
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
'''
