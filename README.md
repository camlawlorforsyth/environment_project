# environment_project #

CARS environment project, looking at optical photometry for objects that are physically close to CARS host galaxies.

## Required Catalog Files ##

The Python script 'environments.py' requires the SDSS and GAMA catalog files to be in the same directory as the script itself. The files used in the preliminary analysis can be found here: http://data.sdss3.org/sas/dr8/common/sdss-spectro/redux/galSpecInfo-dr8.fits for the SDSS DR8 catalog (363 MB), and http://www.gama-survey.org/dr3/data/cat/SpecLineSFR/v05/GaussFitSimple.fits for the GAMA catalog (165 MB).

## SDSS Photometry ##

As SDSS can be a bit tricky to navigate (at the best of times), instead of getting a whole photometric catalog and using Astropy to search through it based on the spectroscopic results from [environments_part1.py](environments_part1.py), the user is encouraged to use the CrossID matching tool at http://skyserver.sdss.org/dr8/en/tools/crossid/crossid.asp, using the printed output from [environments_part1.py](environments_part1.py), along with the SQL query available in [SDSS_SQL_query_for_spec_and_phot_matching.txt](SDSS_SQL_query_for_spec_and_phot_matching.txt).

This should greatly simplify the SDSS photometric search.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
