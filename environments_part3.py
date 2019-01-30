# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 3
    ASSIGNMENT: Fit the SED to obtain stellar masses.
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2019-Jan-29
    
    PURPOSE: Using the obtained multi-wavelength photometry from the relevant
             catalogs (GALEX, SDSS/Pan-STARRS/GAMA, 2MASS, and WISE), fit
             the spectral energy distribution (SED) to obtain stellar masses
             and the star formation rate (SFR) for individual companions.
"""

# imports
import numpy as np

index = 0 # control which CARS galaxy's companions are being studied

filelist = {0:'HE0040-1105.csv'}
#filelist = {0:'HE0040-1105.csv', 1:'RBS175.csv', 3:'Mrk1503.csv',
#            4:'Mrk1018.csv', 5:'Mrk1044.csv', 6:'Mrk1048.csv',
#            7:'HE0345+0056.csv', 8:'HE0853+0102.csv', 9:'Mrk707.csv',
#            10:'HE2222-0026.csv', 11:'Mrk926.csv', 12:'Mrk590.csv'}

# create numpy arrays from *.txt data tables
(spectro_ra, spectro_dec, spectro_z, # from initial search
 spectro_vel, spectro_D_A, spectro_sep,
 specobjid, sdss_ra, sdss_dec, sdss_z, sdss_objid, sdss_type, # SDSS values
 modelmag_u, modelmag_u_err, modelmag_g, modelmag_g_err,
 modelmag_r, modelmag_r_err, modelmag_i, modelmag_i_err,
 modelmag_z, modelmag_z_err, sdss_sci_prim,
 twomass_ra, twomass_dec, j_mag, j_mag_err, h_mag, h_mag_err, # 2MASS values
 k_mag, k_mag_err, twomass_photo_qual, twomass_read_flag,
 twomass_cont_flag, twomass_gal_contam ) = np.genfromtxt(
    filelist[index], delimiter = ',', unpack = True)
