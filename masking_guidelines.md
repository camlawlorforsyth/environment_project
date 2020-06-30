
# Masking guidelines for SDSS and GAMA catalogs #

## SDSS ##

### Preliminary Masking ###

We begin with the 'gal_info_dr7_v5_2.fit.gz' file and the 'gal_line_dr7_v5_2.fit.gz' file, both available from the MPA/JHU collaboration here: https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/raw_data.html, as well as the 'gal_totsfr_dr7_v5_2.fits.gz' file available from the MPA/JHU collaboration here: https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/sfrs.html, and lastly the 'totlgm_dr7_v5_2b.fit' file available from the MPA/JHU collaboration here: https://home.strw.leidenuniv.nl/~jarle/SDSS/

'gal_info_dr7_v5_2.fit.gz' has 927552 rows and 25 columns, while 'gal_line_dr7_v5_2.fit.gz' has 927552 rows and 239 columns, while 'totlgm_dr7_v5_2b.fit' has 927552 rows and 7 columns, while 'gal_totsfr_dr7_v5_2.fits.gz' has 927552 rows and 9 columns.

Next we use TOPCAT to join the four tables by employing a quadruple match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in each table, we write '$index' which will match the two tables side by side (ie. like using an 'hstack' in Python with Astropy).

The resulting table will have 927552 rows and 280 columns. We will now mask the resulting merged table based on several criteria (either in the table itself, but more conveniently using the Row Subsets interface):
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

Save the resulting final table as 'SDSS_gal-info_gal-line_totlgm_totsfr_masked.fits'. This file should take up ~258 MB.

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
-MEDIAN_logmass
-P16_logmass
-P84_logmass
-MEDIAN_sfr
-P16_sfr
-P84_sfr
-FLAG

The resulting final table will thus have 224005 rows with 59 columns. Save the resulting table as 'SDSS_gal-info_gal-line_totlgm_totsfr_vCam.fits'. This file should take up ~54.2 MB.

### Secondary Masking ###

#### Cross-match table preparation ####

We fist must prepare the catalogue files that we will cross-match with our 'SDSS_gal-info_gal-line_totlgm_totsfr_vCam.fits' file. This will encompass 4 steps.

We begin with the 'OSSY_emission_flux_errors.fits' file, and the 'OSSY_SDSS_parameters.fits' file, and the 'Oh+_2015_newType1AGNcatalog.fits' file, all available from the OSSY database here: http://gem.yonsei.ac.kr/~ksoh/wordpress/

In addition we need the 'ALPAKA_v1.fits.gz' file from the ALPAKA website here: https://sites.google.com/site/sdssalpaka/, and lastly the 'Liu+_J_ApJS_243_21_table1.dat.gz.fits' and 'Liu+_J_ApJS_243_21_table2.dat.gz.fits' files from here: https://ui.adsabs.harvard.edu/link_gateway/2019ApJS..243...21L/CDS

'OSSY_emission_flux_errors.fits' has 664187 rows and 75 columns, while 'OSSY_SDSS_parameters.fits' has 664187 rows and 7 columns, while 'Oh+_2015_newType1AGNcatalog.fits' has 5553 rows and 14 columns, while 'ALPAKA_v1.fits.gz' has 25670 rows and 220 columns, while 'Liu+_J_ApJS_243_21_table1.dat.gz.fits' has 14584 rows and 9 columns, while 'Liu+_J_ApJS_243_21_table2.dat.gz.fits' has 14584 rows and 65 columns.

1)
Next we use TOPCAT to join the 'OSSY_SDSS_parameters.fits' and 'OSSY_emission_flux_errors.fits' tables by employing a pair match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in 'OSSY_SDSS_parameters.fits' we write 'SDSS_ID' and for the Matched Value column in 'OSSY_emission_flux_errors.fits' we write 'SDSS_ID', which will cross-match the two tables based on the unique SDSS object identifier.

The resulting table will have 664187 rows and 82 columns. We will now hide various columns that we do not care about in our analysis. In TOPCAT, use the Display column metadata option to view the column metadata. In this interface, we can also deselect columns that are shown. For our final table, we require only the following columns:
-RA
-DEC
-REDSHIFT
-PLATE
-MJD
-FIBERID
-SIG_BALMER
-SIG_BALMER_E
-EBV_STAR
-EBV_GAS
-E_EBV_STAR
-E_EBV_GAS
-HB_4861
-OIII_5006
-HA_6562
-NII_6583
-HB_4861_ERROR
-OIII_5006_ERROR
-HA_6562_ERROR
-NII_6583_ERROR

The resulting final table will thus have 664187 rows with 20 columns. Save the resulting table as 'OSSY_EmissionLine_final.fits'. This file should take up ~93.7 MB.

2)
Next we use TOPCAT to join the 'OSSY_SDSS_parameters.fits' and 'Oh+_2015_newType1AGNcatalog.fits' tables by employing a pair match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in 'OSSY_SDSS_parameters.fits' we write 'SDSS_ID' and for the Matched Value column in 'Oh+_2015_newType1AGNcatalog.fits' we write 'SDSS_OBJID', which will cross-match the two tables based on the unique SDSS object identifier.

The resulting table will have 5553 rows and 21 columns. We will now hide various columns that we do not care about in our analysis. In TOPCAT, use the Display column metadata option to view the column metadata. In this interface, we can also deselect columns that are shown. For our final table, we require only the following columns:
-RA_1
-DEC_1
-REDSHIFT_1
-PLATE
-MJD
-FIBERID
-FWHM
-FWHM_ERR
-EW_BHA

The resulting final table will thus have 5553 rows with 6 columns. Save the resulting table as 'Oh+_Type1AGN_final.fits'. This file should take up ~334 KB.

3)
Next we use TOPCAT to mask the 'ALPAKA_v1.fits.gz' table. We will now hide various columns that we do not care about in our analysis. In TOPCAT, use the Display column metadata option to view the column metadata. In this interface, we can also deselect columns that are shown. For our final table, we require only the following columns:
-RA
-DEC
-Z_1
-AGN_TYPE
-SDSS_HA_EW
-SDSS_HA_EW_ERR
-SDSS_HA_CONT
-SDSS_HA_SIGMA
-SDSS_HA_SIGMA_ERR
-SDSS_HB_EW
-SDSS_HB_EWERR
-SDSS_HB_CONT
-SDSS_HB_SIGMA
-SDSS_HB_SIGMAERR
-SDSS_OIII_EW
-SDSS_OIII_EWERR
-SDSS_OIII_CONT
-SDSS_OIII_SIGMA
-SDSS_OIII_SIGMAERR
-SDSS_NII_EW
-SDSS_NII_EWERR
-SDSS_NII_CONT
-SDSS_NII_SIGMA
-SDSS_NII_SIGMAERR
-HB_WAV
-HB_FWHM
-HB_FWHM_ERR
-HB_FLUX
-HB_FLUX_ERR
-HBB_WAV
-HBB_FWHM
-HBB_FWHM_ERR
-HBB_FLUX
-HBB_FLUX_ERR
-HBBB_WAV
-HBBB_FWHM
-HBBB_FWHM_ERR
-HBBB_FLUX
-HBBB_FLUX_ERR
-OIII_5007_WAV
-OIII_5007_FWHM
-OIII_5007_FWHM_ERR
-OIII_5007_FLUX
-OIII_5007_FLUX_ERR
-OIII_5007B_WAV
-OIII_5007B_FWHM
-OIII_5007B_FWHM_ERR
-OIII_5007B_FLUX
-OIII_5007B_FLUX_ERR
-HA_WAV
-HA_FWHM
-HA_FWHM_ERR
-HA_FLUX
-HA_FLUX_ERR
-HAB_WAV
-HAB_FWHM
-HAB_FWHM_ERR
-HAB_FLUX
-HAB_FLUX_ERR
-HABB_WAV
-HABB_FWHM
-HABB_FWHM_ERR
-HABB_FLUX
-HABB_FLUX_ERR
-NII_6584_WAV
-NII_6584_FWHM
-NII_6584_FWHM_ERR
-NII_6584_FLUX
-NII_6584_FLUX_ERR
-NII_6584B_WAV
-NII_6584B_FWHM
-NII_6584B_FWHM_ERR
-NII_6584B_FLUX
-NII_6584B_FLUX_ERR
-R_EDD

The resulting final table will thus have 25670 rows with 75 columns. Save the resulting table as 'ALPAKA_final.fits'. This file should take up ~7.51 MB.

4)
Next we use TOPCAT to join the 'Liu+_J_ApJS_243_21_table1.dat.gz.fits' and 'Liu+_J_ApJS_243_21_table2.dat.gz.fits' tables by employing a pair match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in 'Liu+_J_ApJS_243_21_table1.dat.gz.fits' we write 'Seq' and for the Matched Value column in 'Liu+_J_ApJS_243_21_table2.dat.gz.fits' we write 'Seq', which will cross-match the two tables based on the unique sequence object identifier.

The resulting table will have 14584 rows and 74 columns. We will now hide various columns that we do not care about in our analysis. In TOPCAT, use the Display column metadata option to view the column metadata. In this interface, we can also deselect columns that are shown. For our final table, we require only the following columns:
-RAdeg
-DECdeg
-z
-Plate
-MJD
-Fiber
-lamBHa
-e_lamBHa
-FWHMBHa
-e_FWHMBHa
-EWBHa
-e_EWBHa
-lam5007
-e_lam5007
-FWHM5007
-e_FWHM5007
-EW5007
-e_EW5007

The resulting final table will thus have 14584 rows with 18 columns. Save the resulting table as 'Liu+_BLAGN_final.fits'. This file should take up ~1.01 MB.

#### Cross-match table execution ####

We are now able to finally include the relevant information from the above 4 tables into the 'SDSS_gal-info_gal-line_totlgm_totsfr_vCam.fits' file.

We begin with the 'SDSS_gal-info_gal-line_totlgm_totsfr_vCam.fits' file, the 'OSSY_EmissionLine_final.fits' file, the 'Oh+_Type1AGN_final.fits', the 'ALPAKA_final.fits' file, and the 'Liu+_BLAGN_final.fits' file, as previously saved above.

Next we use TOPCAT to join the various tables by employing successive pair matches. The algorithm we will be using is 'Exact Value'. For the Matched Value column in 'SDSS_gal-info_gal-line_totlgm_totsfr_vCam.fits' we write 'concat(PLATEID_1, "_", MJD, "_", FIBERID_1)' and for 'OSSY_EmissionLine_final.fits' we write 'concat(PLATE, "_", MJD, "_", FIBERID)'. We ensure that the Join Type is 'All from 1' so that we still have our 224005 galaxies in SDSS.

The resulting table will have 224005 rows and 79 columns. Of the columns from 'OSSY_EmissionLine_final.fits', we only require the following:
-SIG_BALMER
-SIG_BALMER_E

Thus the first matched table will have 224005 rows and 61 columns. We now employ another pair match. For the Matched Value column in 'match(1,3)' we write 'concat(PLATEID_1, "_", MJD_1, "_", FIBERID_1)' and for 'Oh+_Type1AGN_final.fits' we write 'concat(PLATE, "_", MJD, "_", FIBERID)'. We ensure that the Join Type is 'All from 1' so that we still have our 224005 galaxies in SDSS.

The resulting table will have 224005 rows and 70 columns. Of the columns from 'Oh+_Type1AGN_final.fits', we only require the following:
-FWHM
-FWHM_ERR

Thus the second matched table will have 224005 rows and 63 columns. We now employ another pair match, but this time using the 'Sky' algorithm. For the RA column in 'match(6,4)' we write 'RA_1_1' and for the Dec column we write 'DEC_1_1' and for the RA column in 'ALPAKA_final.fits' we write 'RA' and for the Dec column we write 'DEC'. We ensure that the Join Type is 'All from 1' so that we still have our 224005 galaxies in SDSS.

The resulting table will have 224005 rows and 139 columns. Of the columns from 'ALPAKA_final.fits', we only require the following:
-SDSS_HA_SIGMA
-SDSS_HA_SIGMA_ERR
-SDSS_OIII_SIGMA
-SDSS_OIII_SIGMAERR
-OIII_5007_FWHM
-OIII_5007_FWHM_ERR
-HABB_FWHM
-HABB_FWHM_ERR

Thus the third matched table will have 224005 rows and 71 columns. We now employ another pair match. For the Matched Value column in 'match(7,2)' we write 'concat(PLATEID_1, "_", MJD_1, "_", FIBERID_1)' and for 'Liu+_BLAGN_final.fits' we write 'concat(Plate, "_", MJD, "_", Fiber)'. We ensure that the Join Type is 'All from 1' so that we still have our 224005 galaxies in SDSS.

The resulting table will have 224005 rows and 89 columns. Of the columns from 'Liu+_BLAGN_final.fits', we only require the following:
-FWHMBHa
-e_FWHMBHa
-FWHM5007
-e_FWHM5007

Thus the fourth and final matched table will have 224005 rows and 75 columns. We now rename the 'RA_1_1' column to 'RA_1', the 'DEC_1_1' column to 'DEC_1', and the 'Z_1' column to 'Z'. We finally include some information in the column descriptions for the various columns we've added above. For those columns from 'OSSY_EmissionLine_final.fits' we write 'from Oh et al. 2011, ApJS, 195, 13', for those columns from 'Oh+_Type1AGN_final.fits' we write 'from Oh et al. 2015, ApJS, 219, 1', for those columns from 'ALPAKA_final.fits' we write 'from Mullaney et al. 2013, MNRAS, 433, 622', and lastly for those columns from 'Liu+_BLAGN_final.fits' we write 'from Liu et al. 2019, ApJS, 243, 21'. We next define a new synthetic column 'BL' with the expression as follows: "((SIG_BALMER/SIG_BALMER_E>=3) && (2.3548200450309493*SIG_BALMER>=1200)) || ((FWHM/FWHM_ERR>=3) && (FWHM>=1200)) || ((SDSS_HA_SIGMA/SDSS_HA_SIGMA_ERR>=3) && (2.3548200450309493*2.99792458e5*SDSS_HA_SIGMA/6562.8>=1200)) || ((HABB_FWHM/HABB_FWHM_ERR>=3) && (HABB_FWHM>=1200)) || ((FWHMBHa/e_FWHMBHa>=3) && (FWHMBHa>=1200))". Save the resulting table as 'SDSS_gal-info_gal-line_totlgm_totsfr_SpecClassInitial_vCam.fits'. This file should take up ~71.6 MB.

## GAMA ##

### Preliminary Masking ###

We begin with the 'GaussFitSimple.fits' file, the 'StellarMasses.fits' file, and the 'MagPhys.fits' file, all available from GAMA here: http://www.gama-survey.org/dr3/data/cat/, along with the 'SpecClass_Gordon2019_v1.fits' file, from Y. Gordon (private communication).

'GaussFitSimple.fits' has 200228 rows and 214 columns, while 'StellarMasses.fits' has 120739 rows and 92 columns, while 'MagPhys.fits' has 120114 rows and 181 columns, while 'SpecClass_Gordon2019_v1.fits' has 305529 rows and 57 columns.

Next we use TOPCAT to join the four tables by employing a quadruple match. The algorithm we will be using is 'Exact Value'. For the Matched Value column in each table, we select 'CATAID' which will match the four tables based on the unique GAMA identifier.

The resulting table will have 119190 rows and 544 columns. We will now mask the resulting merged table based on several criteria (either in the table itself, but more conveniently using the Row Subsets interface):
1) Redshift (column 'Z_1')
   -We only want objects with a redshift between 0 and 0.075 (note that we will further mask based on redshift using various Python functions).
   --This removes 110244 rows, leaving us with 8946 good rows.
2) Redshift quality (column 'NQ_1')
   -We only want objects with NQ_1 >= 3.
   --This removes 414 rows, leaving us with 8532 good rows.
3) Survey class (column 'SURVEY_CLASS')
   -We only want objects with SURVEY_CLASS >= 3.
   --This removes 10 rows, leaving us with 8522 good rows.
4) Fringing in AATSpecAllv27 (column 'BAD_FLAG')
   -We only want objects with BAD_FLAG==0.
   --This removes 66 rows, leaving us with 8456 good rows.

Save the resulting table as 'GAMA_GaussFitSimple_StellarMasses_MagPhys_SpecClassGordon_masked.fits'. This file should take up ~20.2 MB.

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
-HB_FITFAIL
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
-OIIIR_NPEG
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
-NIIR_EW_1
-NIIR_EW_ERR
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
-uminusr
-deluminusr
-gminusi_stars
-uminusr_stars
-absmag_u
-delabsmag_u
-absmag_u_stars
-absmag_g
-delabsmag_g
-absmag_g_stars
-absmag_r
-delabsmag_r
-absmag_r_stars
-absmag_i
-delabsmag_i
-absmag_i_stars
-absmag_z
-delabsmag_z
-absmag_z_stars
-SFR_0_1Gyr_best_fit
-SFR_0_1Gyr_percentile16
-SFR_0_1Gyr_percentile50
-SFR_0_1Gyr_percentile84
-BAD_FLAG
-broadHB
-broadHA
-log_OIII_HB
-log_NII_HA
-EmLineType
-EmLineMethod

The resulting final table will thus have 8456 rows with 93 columns. Save the resulting table as 'GAMA_GaussFitSimple_StellarMasses_MagPhys_SpecClassGordon_vCam.fits'. This file should take up ~3.10 MB.
