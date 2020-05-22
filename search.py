
# imports
import numpy as np

import astropy.constants as const
from astropy.coordinates import Angle, SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)

def catalog_search(catalog_path, cat_name, CARS_host, CARS_sky_coords,
                   redshift, delta_v, radius, self_in_search) :
    
    lower_z = redshift - delta_v*(u.km/u.s)/const.c.to('km/s')
    upper_z = redshift + delta_v*(u.km/u.s)/const.c.to('km/s')
    
    D_A = cosmo.angular_diameter_distance(redshift)
    radius = Angle( ((radius*u.Mpc)/D_A).value, u.radian ) # find the angular radius
    radius = radius.to('arcmin')
    
    catalog = Table.read(catalog_path)
    
    if (cat_name == 'SDSS') :
        mask = (
                (lower_z <= catalog['Z']) & (catalog['Z'] <= upper_z) &
                ((catalog['Z'] / catalog['Z_ERR']) >= 3)
                ) # mask based on velocity cuts, redshift quality
    if (cat_name == 'GAMA') or (cat_name == 'GAMA_alt') :
        catalog.rename_column('Z_1', 'Z')
        # catalog.rename_column('RA_1', 'RA')
        # catalog.rename_column('DEC_1', 'DEC')
        catalog.rename_column('NQ_1', 'NQ')
        mask = (
                (lower_z <= catalog['Z']) & (catalog['Z'] <= upper_z) &
                (catalog['NQ'] >= 3)
                ) # mask based on velocity cuts, redshift quality
    catalog = catalog[mask] # mask the catalog once to speed things up
    
    D_A = cosmo.angular_diameter_distance(catalog['Z'])
    positions = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], distance=D_A,
                         unit=(u.deg, u.deg, u.Mpc) ) # get sky positions
    
    d2d = CARS_sky_coords.separation(positions).to(u.arcmin) # find sky sep.s
    new_mask = (d2d <= radius) # mask based on 2D sky separation
    catalog = catalog[new_mask]
    
    d3d = CARS_sky_coords.separation_3d(positions)[new_mask] # find physical sep.s
    
    if (cat_name == 'SDSS') :
        D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)/u.pc
        M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        g_i_Mstar = (1.15 + 0.7*(catalog['KCOR_MAG'][:,0] -
                                 catalog['KCOR_MAG'][:,2]) - 0.4*M_i)*u.solMass
    if (cat_name == 'GAMA') or (cat_name == 'GAMA_alt') :
        g_i_Mstar = (1.15 + 0.7*catalog['gminusi'] -
                     0.4*catalog['absmag_i'])/u.mag*u.solMass
            # based on relation from Taylor+ 2011, MNRAS, 418, 1587
    
    # see required aperture correction for stellar mass or absolute magnitude
    # http://www.gama-survey.org/dr3/data/cat/StellarMasses/v20/StellarMasses.notes
    
    host_list = []
    for i in range(len(catalog)) :
        host_list.append(str(CARS_host))
    
    # add the additional columns to the catalog
    catalog['D_A'] = D_A[new_mask]
    catalog['d2d'] = d2d[new_mask]
    catalog['d3d'] = d3d
    catalog['g_i_Mstar'] = g_i_Mstar
    catalog['CARS_host'] = host_list
    
    if self_in_search == True :
        if len(catalog) > 0 :
            catalog.sort('d2d')
            catalog.remove_row(0)
    
#    if (cat_name == 'SDSS') :
#        catalog['log_mass'] = catalog['MEDIAN_2']
#    if (cat_name == 'GAMA') :
#        catalog['log_mass'] = catalog['logmstar']
    
    return catalog
