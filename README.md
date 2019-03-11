# environment_project #

CARS environment project, looking at optical photometry for objects that are physically close to CARS host galaxies.

## Required Catalog Files ##

The Python script 'environments_part1.py' requires the SDSS and GAMA catalog files to be in the `catalog/` directory. The files used in the preliminary analysis can be found here: http://data.sdss3.org/sas/dr12/sdss/spectro/redux/specObj-dr12.fits for the SDSS DR12 catalog (4.53 GB), and http://www.gama-survey.org/dr3/data/cat/SpecCat/v27/SpecObj.fits for the GAMA catalog (55.8 MB).

## SDSS Photometry ##

As SDSS can be a bit tricky to navigate (at the best of times), instead of getting a whole photometric catalog and using Astropy to search through it based on the spectroscopic results from [environments_part1.py](environments_part1.py), the user is encouraged to use the CrossID matching tool at http://skyserver.sdss.org/dr12/en/tools/crossid/crossid.aspx, using the printed output from [environments_part1.py](environments_part1.py), along with the SQL query available in [SDSS_SQL_query.txt](SQL_queries/SDSS_SQL_query.txt). Note that this search queries the SpecPhotoAll table, which combines the spectroscopic and photometric parameters of an object. It is a precomputed join between the 'PhotoObjAll' and 'SpecObjAll' tables.

This should greatly simplify the SDSS photometric search.

## 2MASS and WISE Photometry ##

In a similar fashion to SDSS, the user is encouraged to use the NASA/IPAC Infrared Science Archive (IRSA) catalog search at https://irsa.ipac.caltech.edu/frontpage/ to obtain both NIR and MIR photometry from 2MASS and WISE, respectively. This should again simplify the IR photometric search, as many photometries can be selected and compared, with useful information regarding if a measurement is contaminated or associated with a very extended object.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
