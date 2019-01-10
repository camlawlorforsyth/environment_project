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


# GALEX
# see here for bandpass info:
#   https://asd.gsfc.nasa.gov/archive/galex/Documents/instrument_summary.html
#   http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=GALEX/GALEX.FUV
#   https://galex.stsci.edu/GR6/?page=userfaq
# data access here:
#   http://www.gama-survey.org/dr3/data/cat/GalexPhotometry/v02/
# data access available through online form, here:
#   http://galex.stsci.edu/GR6/?page=mastform
#   https://galex.stsci.edu/GR6/?page=sqlform
# schema browser: https://galex.stsci.edu/GR6/?page=dbinfo

# GAMA
# see here for bandpass info:
#
# more info here:
#   http://www.gama-survey.org/dr3/
#
# data access here:
#   need info/help from Yjan - where is this photometry?

# SDSS # ugriz
# see here for bandpass info:
#   http://www.sdss3.org/instruments/camera.php#Filters
#   http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
#   http://www.sdss3.org/dr8/algorithms/magnitudes.php
# data access available through online form, here:
#   http://skyserver.sdss.org/dr8/en/tools/crossid/crossid.asp
# schema browser: http://skyserver.sdss.org/dr8/en/help/browser/browser.asp

# Pan-STARRS1 # uses grizy
# more info here:
#   https://outerspace.stsci.edu/display/PANSTARRS/PS1+FAQ+-+Frequently+asked+questions
# see here for bandpass info:
#   https://outerspace.stsci.edu/display/PANSTARRS/PS1+Filter+properties
#   http://adsabs.harvard.edu/abs/2012ApJ...750...99T
# with updated zero-points here:
#   http://adsabs.harvard.edu/abs/2012ApJ...756..158S
# data access available through online form, here:
#   https://outerspace.stsci.edu/display/PANSTARRS/How+to+retrieve+and+use+PS1+data
#   http://archive.stsci.edu/panstarrs/search.php

# 2MASS # uses 2MASS JHK_s
#xsc_aaa.gz # 2MASS extended source catalog (declination < 0.0°)
#xsc_baa.gz # 2MASS extended source catalog (declination >= 0.0°)
# see here for bandpass info:
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
# more info here:
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/explsup.html
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec1_3.html
#    https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_3.html
#    https://irsa.ipac.caltech.edu/2MASS/download/allsky/format_xsc.html
# data access available through bulk download, here:
#    https://irsa.ipac.caltech.edu/2MASS/download/allsky/

# WISE # uses AB system
#WISECat.fits # WISE photometric catalog
# see here for bandpass info:
#   http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
# more info here:
#    http://wise2.ipac.caltech.edu/docs/release/allwise/
#    http://wise.ssl.berkeley.edu/astronomers.html
# data access here:
#   http://www.gama-survey.org/dr3/data/cat/WISEPhotometry/v02/WISECat.fits
