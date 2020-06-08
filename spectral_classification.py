
# imports
import numpy as np

from astropy.table import Table

def spectral_classification(sample) :
    
    if sample == 'SDSS' :
        in_path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassPrelim_vCam.fits'
        out_path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassCam_vCam.fits'
        catalog = Table.read(in_path)
        
        log_OIII_HB = np.log10( catalog['OIII_5007_FLUX'] / catalog['H_BETA_FLUX'] )
        log_NII_HA = np.log10( catalog['NII_6584_FLUX'] / catalog['H_ALPHA_FLUX'] )
        HA_EW = np.log10( catalog['H_ALPHA_FLUX'] / catalog['H_ALPHA_CONT'] )
    
    if sample == 'GAMA' :
        in_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassPrelim_vCam.fits'
        out_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassCam_vCam.fits'
        catalog = Table.read(in_path)
        
        HB_flux_corr = (1 + 2.5/catalog['HB_EW_1'])*catalog['HB_FLUX_1']
        log_OIII_HB = np.log10( catalog['OIIIR_FLUX_1'] / HB_flux_corr )
        HA_flux_corr = (1 + 2.5/catalog['HA_EW_1'])*catalog['HA_FLUX_1']
        log_NII_HA = np.log10( catalog['NIIR_FLUX_1'] / HA_flux_corr )
        HA_EW = np.log10(catalog['HA_EW_1'])
    
    BPT_SFG = ( (log_OIII_HB < 0.61/(log_NII_HA - 0.05) + 1.3) &
                (log_NII_HA < 0.05) )
    BPT_Comp = ( (log_OIII_HB > 0.61/(log_NII_HA - 0.05) + 1.3) &
                 (log_OIII_HB < 0.61/(log_NII_HA - 0.47) + 1.19) &
                 (log_NII_HA < 0.47) )
    BPT_Sey = ( (log_OIII_HB > 0.61/(log_NII_HA - 0.47) + 1.19) &
                (log_OIII_HB > 1.05*log_NII_HA + 0.45) )
    BPT_LIN = ( (log_OIII_HB > 0.61/(log_NII_HA - 0.47) + 1.19) &
                (log_OIII_HB < 1.05*log_NII_HA + 0.45) )
    
    WHAN_passives = ( (HA_EW < -log_NII_HA + np.log10(0.5)) |
                      (HA_EW < np.log10(0.5)) )
    WHAN_SFG = (HA_EW > -log_NII_HA + np.log10(0.5)) & (log_NII_HA < -0.4)
    WHAN_Comp = ( (HA_EW > -log_NII_HA + np.log10(0.5)) &
                  (log_NII_HA > -0.4) & (log_NII_HA < -0.32) )
    WHAN_Sey = (HA_EW > np.log10(6)) & (log_NII_HA > -0.32)
    WHAN_LIN = ( (HA_EW > -log_NII_HA + np.log10(0.5)) &
                 (HA_EW > np.log10(0.5)) & (log_NII_HA > -0.32) & (HA_EW < np.log10(6)) )
    
    EmLineType = []
    EmLineMethod = []
    bad_x = []
    bad_y = []
    for i in range(len(catalog)) :
        if catalog['BL'][i] :
            EmLineType.append('BLAGN')
            EmLineMethod.append('BL')
        elif catalog['BPT'][i] :
            if BPT_SFG[i] :
                EmLineType.append('SFG')
            elif BPT_Comp[i] :
                EmLineType.append('Comp')
            elif BPT_Sey[i] :
                EmLineType.append('Seyfert')
            elif BPT_LIN[i] :
                EmLineType.append('LINER')
            else :
                bad_x.append(log_NII_HA[i])
                bad_y.append(log_OIII_HB[i])
                EmLineType.append('not_ELG')
            EmLineMethod.append('BPT')
        elif catalog['WHAN'][i] :
            if WHAN_passives[i] :
                EmLineType.append('Passive')
            elif WHAN_SFG[i] :
                EmLineType.append('SFG')
            elif WHAN_Comp[i] :
                EmLineType.append('Comp')
            elif WHAN_Sey[i] :
                EmLineType.append('Seyfert')
            elif WHAN_LIN[i] :
                EmLineType.append('LINER')
            else :
                EmLineType.append('not_ELG')
            EmLineMethod.append('WHAN')
        else :
            EmLineType.append('not_ELG')
            EmLineMethod.append('not_ELG')
    
    catalog['EmLineType'] = EmLineType
    catalog['EmLineMethod'] = EmLineMethod
    
    print(bad_x)
    print(bad_y)
    
    catalog.write(out_path, overwrite=True)
    
    return
