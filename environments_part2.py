# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 2
    ASSIGNMENT: Get the multi-wavelength photometric data for companions
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2018-Dec-20
    
    PURPOSE: Given the list of physically close companions to CARS host
             galaxies, now retrieve the multi-wavelength photometry from
             relevant catalogs (ex. UV: Galex; optical: SDSS/Pan-STARRS/GAMA;
             NIR: 2MASS; MIR: WISE) for those companions.
"""

# imports
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np
import warnings
warnings.filterwarnings("ignore") # ignore warnings about division by 0 when
# taking z/z_err >= 3, and ignore integration warning for the cosmology



#xsc_aaa.gz # 2MASS extended source catalog (declination < 0.0°)
#xsc_baa.gz # 2MASS extended source catalog (declination > 0.0°)
#more info here:
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/explsup.html
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec1_3.html
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_3.html
#    https://irsa.ipac.caltech.edu/2MASS/download/allsky/format_xsc.html