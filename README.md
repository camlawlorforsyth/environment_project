# environment_project #

CARS environment project, looking at optical photometry for objects that are physically close to CARS host galaxies.

## Required Catalog Files ##

The Python script 'environments_part1.py' requires the SDSS and GAMA catalog files to be in the `catalog/` directory. The files used in the preliminary analysis can be found here: http://data.sdss3.org/sas/dr12/sdss/spectro/redux/specObj-dr12.fits for the SDSS DR12 spectro catalog (4.53 GB), and http://www.gama-survey.org/dr3/data/cat/SpecCat/v27/SpecObj.fits for the GAMA spectro catalog (55.8 MB). Information regarding the contents of these files are here: http://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=description+SpecObj+V and http://www.gama-survey.org/dr3/schema/table.php?id=31, respectively. Note that this analysis uses SDSS DR12, in order to match to GALEX.

## SDSS Photometry ##

As SDSS can be a bit tricky to navigate (at the best of times), instead of getting a whole photometric catalog and using Astropy to search through it based on the spectroscopic results from [environments_part1.py](environments_part1.py), the user is encouraged to use the CrossID matching tool at http://skyserver.sdss.org/dr12/en/tools/crossid/crossid.aspx, using the printed output from [environments_part1.py](environments_part1.py), along with the SQL query available in [SDSS_SQL_query.txt](SQL_queries/SDSS_SQL_query.txt). Note that this search queries the SpecPhotoAll table, which combines the spectroscopic and photometric parameters of an object. It is a precomputed join between the 'PhotoObjAll' and 'SpecObjAll' tables. Schema regarding this file is available at http://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=description+SpecPhoto+V.

This should greatly simplify the SDSS photometric search.

## 2MASS and WISE Photometry ##

In a similar fashion to SDSS, the user is encouraged to use the NASA/IPAC Infrared Science Archive (IRSA) catalog search at https://irsa.ipac.caltech.edu/frontpage/ to obtain both NIR and MIR photometry from 2MASS and WISE, respectively. This should again simplify the IR photometric search, as many photometries can be selected and compared, with useful information regarding if a measurement is contaminated.

The 2MASS XSC photometric catalog is available at https://irsa.ipac.caltech.edu/Missions/2mass.html, with its schema available at https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2a.html.

The WISE (AllWISE) photometric catalog is available at https://irsa.ipac.caltech.edu/Missions/wise.html, with its schema available at http://wise2.ipac.caltech.edu/docs/release/allwise/expsup/sec2_1a.html

## GALEX Photometry ##

The GALEX GR6/GR7 photometry can be obtained from the catalog at https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html, whose schema is available at https://galex.stsci.edu/GR6/?page=dbinfo, with additional information at http://www.galex.caltech.edu/researcher/files/mcat_columns_long.txt.

## GAMA Photometry ##

The GAMA aperture matched photometry (matched to SDSS and VIKING) for HE 0853+0102 (located in the equatorial region) is available in the file 'GAMA_ApMatchedCat.fits' from http://www.gama-survey.org/dr3/data/cat/ApMatchedPhotom/v06/ApMatchedCat.fits, with schema available at http://www.gama-survey.org/dr3/schema/dmu.php?id=10, and http://www.gama-survey.org/dr3/data/cat/ApMatchedPhotom/v06/ApMatchedCat.par, with further information regarding the matching process at http://www.gama-survey.org/dr3/data/cat/ApMatchedPhotom/v06/ApMatchedCat.notes. This should be used as a consistency check for this object, once the SDSS photometry has been collected.

# GAMA Environmental Parameters #

Some environmental parameters of interest regarding GAMA are available in the file 'GAMA_EnvironmentMeasures.fits' from http://www.gama-survey.org/dr3/data/cat/EnvironmentMeasures/v05/EnvironmentMeasures.fits, with schema available at http://www.gama-survey.org/dr3/schema/dmu.php?id=16, and http://www.gama-survey.org/dr3/data/cat/EnvironmentMeasures/v05/EnvironmentMeasures.par, with further information regarding the envrionmental measures at http://www.gama-survey.org/dr3/data/cat/EnvironmentMeasures/v05/EnvironmentMeasures.notes.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
