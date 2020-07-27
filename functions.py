
# imports
import numpy as np

from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)

def intercept(xx, b) :
    return xx + b

def line(xx, m, b) :
    return m*xx + b

def parabola(xx, amp, translation, offset) :
    return amp*(xx - translation)**2 + offset

def gaussian(xx, mu, sigma_sqr) :
    exponent = -0.5*np.square(xx - mu)/sigma_sqr
    norm = 1/np.sqrt(2*np.pi*sigma_sqr)
    return norm*np.exp(exponent)

def image_scale() :
    
    galaxies = Table.read('catalogs/raw_cats/CARS_thesis_sample.fits')
    redshifts = galaxies['Z']
    Dists = cosmo.angular_diameter_distance(redshifts)
    arcsec_scale = cosmo.arcsec_per_kpc_proper(redshifts)
    
    radius = 2000*u.kpc
    image_side_length = 3000*u.pix
    size_of_interest = image_side_length/2
    
    pixel_scale = radius*arcsec_scale/(size_of_interest)
    # print(pixel_scale)
    # print( Angle( ((2*u.Mpc)/Dists).value, u.radian ).to('arcmin') )
    
    used = np.array([1.23731709, 1.47521605, 1.24310033, 1.58221001,
                     0.85964875, 1.19281806, 1.56419326, 1.9152773,
                     1.19039626, 1.27362488, 1.05597032, 1.2068972])*u.arcsec/u.pix
    
    # print((used*size_of_interest/arcsec_scale*2).to(u.Mpc))
    
    return
