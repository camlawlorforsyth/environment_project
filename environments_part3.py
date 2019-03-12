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
import matplotlib.pyplot as plt
import numpy.ma as ma

# constants
currentFig = 1
pogson = np.power(100, 0.2)

#.......................................................................convert
def convert(mag) :
    
    flux = (3631000)*np.power(10, -mag/pogson) # result in mJy
    
    return flux

#.........................................................................histo
def histo(param, label, num_bins) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    
    ax = fig.add_subplot(111)
    ax.hist(param, bins=num_bins, density=True, color='k')
    plt.xlabel("%s" % label, fontsize = 15)
    
    plt.tight_layout()
    plt.show()
    
    return

#..............................................................................

dat = ascii.read('photometry_v1.csv') # requires columns to have unique names
#histo(2*dat['petroRad_r'], r'$2r_p$', 100) # show petroRad_r distribution

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

i_mag = convert( dat['petroMag_i']-dat['extinction_i'] )
i_err = np.log(10)/pogson*i_mag*dat['petroMagErr_i']
i_mag.fill_value = -99
i_mag = i_mag.filled()
i_err.fill_value = -99
i_err = i_err.filled()

z = convert( dat['petroMag_z']-dat['extinction_z']+0.02 )
z_err = np.log(10)/pogson*z*dat['petroMagErr_z']
z.fill_value = -99
z = z.filled()
z_err.fill_value = -99
z_err = z_err.filled()

# 2MASS
J, H, K = [], [], []
J_err, H_err, K_err = [], [], []
for i in range(len(dat)) :
    if 2*(dat['petroRad_r'][i]) < 5 :
        J.append(convert(dat['j_m_5'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_5'][i])
        H.append(convert(dat['h_m_5'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_5'][i])
        K.append(convert(dat['k_m_5'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_5'][i])
    elif (2*dat['petroRad_r'][i]) < 7 :
        J.append(convert(dat['j_m_7'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_7'][i])
        H.append(convert(dat['h_m_7'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_7'][i])
        K.append(convert(dat['k_m_7'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_7'][i])
    elif (2*dat['petroRad_r'][i]) < 10 :
        J.append(convert(dat['j_m_10'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_10'][i])
        H.append(convert(dat['h_m_10'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_10'][i])
        K.append(convert(dat['k_m_10'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_10'][i])
    elif (2*dat['petroRad_r'][i]) < 15 :
        J.append(convert(dat['j_m_15'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_15'][i])
        H.append(convert(dat['h_m_15'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_15'][i])
        K.append(convert(dat['k_m_15'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_15'][i])
    elif (2*dat['petroRad_r'][i]) < 20 :
        J.append(convert(dat['j_m_20'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_20'][i])
        H.append(convert(dat['h_m_20'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_20'][i])
        K.append(convert(dat['k_m_20'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_20'][i])
    elif (2*dat['petroRad_r'][i]) < 25 :
        J.append(convert(dat['j_m_25'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_25'][i])
        H.append(convert(dat['h_m_25'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_25'][i])
        K.append(convert(dat['k_m_25'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_25'][i])
    elif (2*dat['petroRad_r'][i]) < 30 :
        J.append(convert(dat['j_m_30'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_30'][i])
        H.append(convert(dat['h_m_30'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_30'][i])
        K.append(convert(dat['k_m_30'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_30'][i])
    elif (2*dat['petroRad_r'][i]) < 40 :
        J.append(convert(dat['j_m_40'][i]+0.898))
        J_err.append(np.log(10)/pogson*J[i]*dat['j_msig_40'][i])
        H.append(convert(dat['h_m_40'][i]+1.381))
        H_err.append(np.log(10)/pogson*H[i]*dat['h_msig_40'][i])
        K.append(convert(dat['k_m_40'][i]+1.849))
        K_err.append(np.log(10)/pogson*K[i]*dat['k_msig_40'][i])
    else :
        J.append(ma.masked)
        J_err.append(ma.masked)
        H.append(ma.masked)
        H_err.append(ma.masked)
        K.append(ma.masked)
        K_err.append(ma.masked)
J = ma.array(J)
J_err= ma.array(J_err)
J.fill_value = -99
J = J.filled()
J_err.fill_value = -99
J_err = J_err.filled()
H = ma.array(H)
H_err= ma.array(H_err)
H.fill_value = -99
H = H.filled()
H_err.fill_value = -99
H_err = H_err.filled()        
K = ma.array(K)
K_err= ma.array(K_err)
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
ids = np.arange(0, len(dat))
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
        'i_wl':i_wl, 'i':i_mag, 'i_err':i_err,
        'z_wl':z_wl, 'z':z, 'z_err':z_err,
        'J_wl':J_wl, 'J':J, 'J_err':J_err,
        'H_wl':H_wl, 'H':H, 'H_err':H_err,
        'K_wl':K_wl, 'K':K, 'K_err':K_err,
        'W1_wl':W1_wl, 'W1':W1, 'W1_err':W1_err,
        'W2_wl':W2_wl, 'W2':W2, 'W2_err':W2_err,
        'W3_wl':W3_wl, 'W3':W3, 'W3_err':W3_err,
        'W4_wl':W4_wl, 'W4':W4, 'W4_err':W4_err
        } )


#ascii.write(newdata, 'catalog.txt', overwrite=True)


#tbl = ascii.read('catalog.txt')
#print(tbl.colnames)



