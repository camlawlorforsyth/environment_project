
# imports
import astropy.constants as const
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table
import astropy.units as u

# constants
mass_limit = 8.452021*u.solMass

def catalog_search(catalog, redshift, D_A, delta_v, radius) :
    
    radius = Angle( ((radius*u.Mpc)/D_A).value, u.radian) # find the angular radius
    radius = radius.to('arcmin')
    
    lower_z = redshift - delta_v*(u.km/u.s)/const.c.to('km/s')
    upper_z = redshift + delta_v*(u.km/u.s)/const.c.to('km/s')
    
    mask = ( (lower_z <= catalog['Z']) & (catalog['Z'] <= upper_z) &
             (catalog['d2d'] <= radius) ) # mask based on velocity cuts and
        # projected sky separation
    catalog = catalog[mask]
    
    return catalog

def search_prep(catalog_path, host, CARS_sky_coords, self_in_search) :
    
    catalog = Table.read(catalog_path)
    
    positions = SkyCoord(ra=catalog['RA_1'], dec=catalog['DEC_1'],
                         distance=catalog['D_A'],
                         unit=(u.deg, u.deg, u.Mpc) ) # get sky positions
    
    d2d = CARS_sky_coords.separation(positions).to(u.arcmin) # find sky sep.s
    d3d = CARS_sky_coords.separation_3d(positions) # find physical sep.s
    
    # add the additional columns to the catalog
    catalog['d2d'] = d2d
    catalog['d3d'] = d3d
    catalog['target_host'] = [host] * len(catalog)
    
    if self_in_search == True :
        if len(catalog) > 0 :
            catalog.sort('d2d')
            if host == catalog['Name'][0] :
                catalog.remove_row(0)
    
    return catalog
