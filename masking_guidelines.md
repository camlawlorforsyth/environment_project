
# Masking guidelines for SDSS and GAMA catalogs #

## SDSS ##

### Preliminary Masking ###

We being with the 'gal_info_dr7_v5_2.fit.gz' file and the 'gal_line_dr7_v5_2.fit.gz' file, both available from the MPA/JHU collaboration here: https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/raw_data.html

'gal_info_dr7_v5_2.fit.gz' has 927552 rows and 25 columns, while 'gal_line_dr7_v5_2.fit.gz' has 927552 rows and 239 columns.

Next we use TOPCAT to join the two tables by employing a pair match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in each table, we write '$index' which will match the two tables side by side (ie. like using an 'hstack' in Python with Astropy).

The resulting table will have 927552 rows and 264 columns. We will now mask the resulting merged table based on several criteria (either in the table itself, but more conveniently using the Row Subsets interface):
1) Redshift (column 'Z')
   -We only want objects with a redshift between 0 and 0.075 (note that we will further mask based on redshift using various Python functions).
   --This removes 690631 rows, leaving us with 236921 good rows.
2) Redshift quality (column 'Z_ERR')
   -We only want objects with Z / Z_ERR >=3.
   --This removes 313 rows, leaving us with 236608 good rows.
3) Target type (column 'TARGETTYPE')
   -We only want targets with TARGETTYPE==GALAXY.
   --This removes 10551 rows, leaving us with 226057 good rows.
4) Spectro type (column 'SPECTROTYPE')
   -We only want targets with SPECTROTYPE==GALAXY.
   --This removes 158 rows, leaving us with 225889 good rows.
5) Primary target (column 'PRIMTARGET')
   -We only want targets with PRIMTARGET==32,64,96,128,192,224,256,288,320 (more information available here: https://users.obs.carnegiescience.edu/yshen/SDSS/Princeton_MIT%20SDSS%20Spectroscopy%20Home%20Page.htm#targetflags)
   --This removies 1084 rows, leaving us with 224815 good rows.
6) Redshift warning (column 'Z_WARNING')
   -We only want targets with Z_WARNING==0.
   --This removes 810 rows, leaving us with 224005 good rows.

Save the resulting final table as 'SDSS_gal_info_gal_line_masked.fits'. This file should take up ~244 MB.

### Subsequent file size reduction (optimization) ###

We will now mask out columns we do not need before creating a final useable table that is smaller in size and thus faster to load in Python.

In TOPCAT, use the Display column metadata option to view the column metadata. In this interface, we can also deselect columns that are shown. For our final table, we require only the following columns:
-PLATEID_1
-MJD
-FIBERID_1
-RA
-DEC
-PRIMTARGET
-Z
-Z_ERR
-NQ
-V_DISP
-V_DISP_ERR
-E_BV_SFD
-SPECTRO_MAG
-KCOR_MAG
-RELEASE
-SIGMA_BALMER
-SIGMA_BALMER_ERR
-SIGMA_FORBIDDEN
-SIGMA_FORBIDDEN_ERR
-H_BETA_CONT
-H_BETA_CONT_ERR
-H_BETA_REQW
-H_BETA_REQW_ERR
-H_BETA_EQW
-H_BETA_EQW_ERR
-H_BETA_FLUX
-H_BETA_FLUX_ERR
-OIII_5007_CONT
-OIII_5007_CONT_ERR
-OIII_5007_REQW
-OIII_5007_REQW_ERR
-OIII_5007_EQW
-OIII_5007_EQW_ERR
-OIII_5007_FLUX
-OIII_5007_FLUX_ERR
-H_ALPHA_CONT
-H_ALPHA_CONT_ERR
-H_ALPHA_REQW
-H_ALPHA_REQW_ERR
-H_ALPHA_EQW
-H_ALPHA_EQW_ERR
-H_ALPHA_FLUX
-H_ALPHA_FLUX_ERR
-NII_6584_CONT
-NII_6584_CONT_ERR
-NII_6584_REQW
-NII_6584_REQW_ERR
-NII_6584_EQW
-NII_6584_EQW_ERR
-NII_6584_FLUX
-NII_6584_FLUX_ERR
-SPECTOFIBER

The resulting final table will thus have 224005 rows with 52 columns. Save the resulting table as 'SDSS_gal_info_gal_line_vCam.fits'. This file should take up ~48.3 MB.

## GAMA ##

### Preliminary Masking ###

We begin with the 'GaussFitSimple.fits' file, the 'StellarMasses.fits' file, both available from GAMA here: http://www.gama-survey.org/dr3/data/cat/, along with the 'SpecClass_Gordon2019_v1.fits' file, from Y. Gordon (private communication).

'GaussFitSimple.fits' has 200228 rows and 214 columns, while 'StellarMasses.fits' has 120739 rows and 92 columns, while 'SpecClass_Gordon2019_v1.fits' has 305529 rows and 57 columns.

Next we use TOPCAT to join the three tables by employing a triple match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in each table, we select 'CATAID' which will match the three tables based on the unique GAMA identifier.

The resulting table will have 119301 rows and 363 columns. We will now mask the resulting merged table based on several criteria (either in the table itself, but more conveniently using the Row Subsets interface):
1) Redshift (column 'Z_1')
   -We only want objects with a redshift between 0 and 0.075 (note that we will further mask based on redshift using various Python functions).
   --This removes 110316 rows, leaving us with 8985 good rows.
2) Redshift quality (column 'NQ_1')
   -We only want objects with NQ_1 >= 3.
   --This removes 416 rows, leaving us with 8569 good rows.
3) Survey class (column 'SURVEY_CLASS')
   -We only want objects with SURVEY_CLASS >= 3.
   --This removes 30 rows, leaving us with 8539 good rows.
4) Fringing in AATSpecAllv27 (column 'BAD_FLAG')
   -We only want objects with BAD_FLAG==0.
   --This removes 66 rows, leaving us with 8473 good rows.

Save the resulting table as 'GAMA_GaussFitSimple_StellarMasses_SpecClassGordon_masked.fits'. This file should take up ~12.7 MB.

### Subsequent file size reduction (optimization) ###

We will now mask out columns we do not need before creating a final useable table that is smaller in size and thus faster to load in Python.

In TOPCAT, use the Display column metadata option to view the column metadata. In this interface, we can also deselect columns that are shown. For our final table, we require only the following columns:
-CATAID_1
-RA_1
-DEC_1
-Z_1
-NQ_1
-SURVEY_1
-SURVEY_CODE_1
-FSCALE
-HB_CONT
-HB_CONT_ERR
-HB_GRAD
-HB_GRAD_ERR
-POS_HB
-POS_HB_ERR
-SIG_HB
-SIG_HB_ERR
-HB_FLUX_1
-HB_FLUX_ERR_1
-HB_EW_1
-HB_EW_ERR
-POS_OIIIR
-POS_OIIIR_ERR
-SIG_OIIIR
-SIG_OIIIR_ERR
-OIIIR_FLUX_1
-OIIIR_FLUX_ERR_1
-OIIIR_EW_1
-OIIIR_EW_ERR
-HA_CONT
-HA_CONT_ERR
-HA_GRAD
-HA_GRAD_ERR
-POS_HA
-POS_HA_ERR
-SIG_HA
-SIG_HA_ERR
-HA_FLUX_1
-HA_FLUX_ERR_1
-HA_EW_1
-HA_EW_ERR
-POS_NIIR
-POS_NIIR_ERR
-SIG_NIIR
-SIG_NIIR_ERR
-NIIR_FLUX_1
-NIIR_FLUX_ERR_1
-NII_EW_1
-NII_EW_ERR
-fluxscale
-logmstar
-dellogmstar
-logmoverl_i
-dellogmoverl_i
-logmintsfh
-dellogmintsfh
-logmremnants
-dellogmremnants
-extBV
-delextBV
-gminusi
-delgminusi
-gminusi_stars
-absmag_g
-delabsmag_g
-absmag_g_stars
-absmag_i
-delabsmag_i
-absmag_i_stars
-BAD_FLAG
-broadHB
-broadHA
-log_OIII_HB
-log_NII_HA
-EmLineType
-EmLineMethod

The resulting final table will thus have 8473 rows with 75 columns. Save the resulting table as 'GAMA_GaussFitSimple_StellarMasses_SpecClassGordon_vCam.fits'. This file should take up ~2.58 MB.
