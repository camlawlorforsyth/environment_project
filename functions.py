
# imports
from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)

def intercept(xx, b) :
    return xx + b

def line(xx, m, b) :
    return m*xx + b

def parabola(xx, amp, translation, offset) :
    return amp*(xx - translation)**2 + offset

def image_scale() :
    
    radius = 2000*u.kpc
    image_side_length = 1000*u.pix
    size_of_interest = image_side_length/2
    # arcsec_scale = cosmo.arcsec_per_kpc_proper(redshifts)
    # pixel_scale = radius*arcsec_scale/(image_side_length/2)
    # print(pixel_scale)
    # print( Angle( ((2*u.Mpc)/Dists).value, u.radian ).to('arcmin') )
    
    return