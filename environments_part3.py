# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 3
    ASSIGNMENT: Fit the SED to obtain stellar masses.
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2019-Feb-25
    
    PURPOSE: Using the obtained multi-wavelength photometry from the relevant
             catalogs (GALEX, SDSS/Pan-STARRS/GAMA, 2MASS, and WISE), fit
             the spectral energy distribution (SED) to obtain stellar masses
             and the star formation rate (SFR) for individual companions.
"""

# imports
import numpy as np

from astropy.io import ascii
from astropy.table import Table

# constants
pogson = np.power(100, 0.2)

#.......................................................................convert
def convert(mag) :
    
    flux = (3631000)*np.power(10, -mag/pogson) # result in mJy
    
    return flux
#..............................................................................

dat = ascii.read('photometry.csv') # requires columns to have unique names

# GALEX
FUV = convert( dat['fuv_mag']-dat['e_bv'] )
FUV_err = np.log(10)/pogson*FUV*dat['fuv_magerr']
FUV.fill_value = -99
FUV = FUV.filled()
FUV_err.fill_value = -99
FUV_err = FUV_err.filled()

NUV = convert( dat['nuv_mag']-dat['e_bv'] )
NUV_err = np.log(10)/pogson*NUV*dat['nuv_magerr']
NUV.fill_value = -99
NUV = NUV.filled()
NUV_err.fill_value = -99
NUV_err = NUV_err.filled()

# SDSS
u = convert( dat['petroMag_u']-dat['extinction_u']-0.04 )
u_err = np.log(10)/pogson*u*dat['petroMagErr_u']
u.fill_value = -99
u = u.filled()
u_err.fill_value = -99
u_err = u_err.filled() 

g = convert( dat['petroMag_g']-dat['extinction_g'] )
g_err = np.log(10)/pogson*g*dat['petroMagErr_g']
g.fill_value = -99
g = g.filled()
g_err.fill_value = -99
g_err = g_err.filled()

r = convert( dat['petroMag_r']-dat['extinction_r'] )
r_err = np.log(10)/pogson*r*dat['petroMagErr_r']
r.fill_value = -99
r = r.filled()
r_err.fill_value = -99
r_err = r_err.filled()

i = convert( dat['petroMag_i']-dat['extinction_i'] )
i_err = np.log(10)/pogson*i*dat['petroMagErr_i']
i.fill_value = -99
i = i.filled()
i_err.fill_value = -99
i_err = i_err.filled()

z = convert( dat['petroMag_z']-dat['extinction_z']+0.02 )
z_err = np.log(10)/pogson*z*dat['petroMagErr_z']
z.fill_value = -99
z = z.filled()
z_err.fill_value = -99
z_err = z_err.filled()

# 2MASS
J = convert( dat['j_m']+0.898 )
J_err = np.log(10)/pogson*J*dat['j_cmsig']
J.fill_value = -99
J = J.filled()
J_err.fill_value = -99
J_err = J_err.filled()

H = convert( dat['h_m']+1.381 )
H_err = np.log(10)/pogson*H*dat['h_cmsig']
H.fill_value = -99
H = H.filled()
H_err.fill_value = -99
H_err = H_err.filled()

K = convert( dat['k_m']+1.849 )
K_err = np.log(10)/pogson*K*dat['k_cmsig']
K.fill_value = -99
K = K.filled()
K_err.fill_value = -99
K_err = K_err.filled()

# WISE
W1 = convert( dat['w1mpro']+2.699 )
W1_err = np.log(10)/pogson*W1*dat['w1sigmpro']
W1.fill_value = -99
W1 = W1.filled()
W1_err.fill_value = -99
W1_err = W1_err.filled()

W2 = convert( dat['w2mpro']+3.339 )
W2_err = np.log(10)/pogson*W2*dat['w2sigmpro']
W2.fill_value = -99
W2 = W2.filled()
W2_err.fill_value = -99
W2_err = W2_err.filled()

W3 = convert( dat['w3mpro']+5.174 )
W3_err = np.log(10)/pogson*W3*dat['w3sigmpro']
W3.fill_value = -99
W3 = W3.filled()
W3_err.fill_value = -99
W3_err = W3_err.filled()

W4 = convert( dat['w4mpro']+6.620 )
W4_err = np.log(10)/pogson*W4*dat['w4sigmpro']
W4.fill_value = -99
W4 = W4.filled()
W4_err.fill_value = -99
W4_err = W4_err.filled()

# wavelengths and redshifts
ids = np.arange(1, 276)
FUV_wl = np.full(len(dat), 1538.6)
NUV_wl = np.full(len(dat), 2315.7)
u_wl = np.full(len(dat), 3551)
g_wl = np.full(len(dat), 4686)
r_wl = np.full(len(dat), 6166)
i_wl = np.full(len(dat), 7480)
z_wl = np.full(len(dat), 8932)
J_wl = np.full(len(dat), 12350)
H_wl = np.full(len(dat), 16620)
K_wl = np.full(len(dat), 21590)
W1_wl = np.full(len(dat), 34000)
W2_wl = np.full(len(dat), 46000)
W3_wl = np.full(len(dat), 120000)
W4_wl = np.full(len(dat), 220000)

# create the table to be used in SED fitting
newdata = Table( {
        '#id':ids,
        'redshift':dat['Spec_z'],
        'FUV_wl':FUV_wl, 'FUV':FUV, 'FUV_err':FUV_err,
        'NUV_wl':NUV_wl, 'NUV':NUV, 'NUV_err':NUV_err,
        'u_wl':u_wl, 'u':u, 'u_err':u_err,
        'g_wl':g_wl, 'g':g, 'g_err':g_err,
        'r_wl':r_wl, 'r':r, 'r_err':r_err,
        'i_wl':i_wl, 'i':i, 'i_err':i_err,
        'z_wl':z_wl, 'z':z, 'z_err':z_err,
        'J_wl':J_wl, 'J':J, 'J_err':J_err,
        'H_wl':H_wl, 'H':H, 'H_err':H_err,
        'K_wl':K_wl, 'K':K, 'K_err':K_err,
        'W1_wl':W1_wl, 'W1':W1, 'W1_err':W1_err,
        'W2_wl':W2_wl, 'W2':W2, 'W2_err':W2_err,
        'W3_wl':W3_wl, 'W3':W3, 'W3_err':W3_err,
        'W4_wl':W4_wl, 'W4':W4, 'W4_err':W4_err
        } )


ascii.write(newdata, 'catalog.txt', overwrite=True)


#tbl = ascii.read('catalog.txt')
#print(tbl.colnames)



