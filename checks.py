
# imports
import numpy as np

from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u

import environmental_parameters as env_params
import plots as plt

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
GAMA_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassGordon_vCam.fits'

def comparison() :
    
    catalog = Table.read('catalogs/joined_cats/GAMA_SpecClassGordon_EnvMeas.fits')
    
    mask = (
            (catalog['EmLineType'] == 'Seyfert') | (catalog['EmLineType'] == 'BLAGN  ')
            )
    catalog = catalog[mask]
    
    # dist = catalog['DistanceTo5nn']
    surf_dens = catalog['SurfaceDensity']
    counts = catalog['CountInCyl']
    ageden_par = catalog['AGEDenPar']
    
    # plt.plot(np.log10(dist), 'log10(Distance to 5NN)',
    #            np.log10(surf_dens), r'log10(Surface Density (Mpc$^{-2}$))',
    #            cat_name='GAMA')
    
    # plt.histo2d(np.log10(dist), 'log10(Distance to 5NN)',
    #               np.log10(surf_dens), r'log10(Surface Density (Mpc$^{-2}$))')
    
    return surf_dens, counts, ageden_par

def consistency_check(cat_name) :
    
    length = 3
    
    cataid_catalog = Table.read('catalogs/GAMA_EnvironmentMeasures.fits')
    cataids = cataid_catalog['CATAID']
    
    np.random.seed(0)
    random_10 = np.random.choice(cataids, length)
#    print(random_10)
    
    base_cat = Table.read('catalogs/GAMA_GaussFitSimple.fits')
    
    indexes = []
    for i in range(length) :
        index = np.where( base_cat['CATAID'] == random_10[i] )[0][0]
        indexes.append(index)
    
    ids = base_cat['CATAID'][indexes]
    ras = Angle(base_cat['RA'][indexes], u.deg)
    decs = Angle(base_cat['DEC'][indexes], u.deg)
    zs = base_cat['Z'][indexes]
    dists = cosmo.angular_diameter_distance(zs)
    
    list_of_sub_cats = []
    rows = []
    for index in range(length) :
        sub_catalog, row = env_params.gama_params(cat_name, GAMA_path, ras[index],
                                                  decs[index], zs[index], dists[index],
                                                  ids[index])
        list_of_sub_cats.append(sub_catalog)
        rows.append(row)
    
    RAND = Table(rows = rows, names=('CATAID', '# Comp.', 'SurfaceDensity',
                                     'SurfaceDensityErr', 'CountInCyl',
                                     'CountInCylErr', 'Overdens.',
                                     'AGEPar', 'AGEParErr', 'AGEScale'))
    RAND['SurfaceDensity'] = RAND['SurfaceDensity']*u.Mpc**(-2)
    RAND['SurfaceDensityErr'] = RAND['SurfaceDensityErr']*u.Mpc**(-2)
    RAND['SurfaceDensity'].format = '.6f'
    RAND['SurfaceDensityErr'].format = '.6f'
    RAND['CountInCylErr'].format = '.3g'
    RAND['Overdens.'].format = '.3f'
    RAND['AGEPar'] = RAND['AGEPar']*u.Mpc**(-1)
    RAND['AGEPar'].format = '.3f'
    RAND['AGEParErr'].format = '.3g'
    RAND['AGEScale'].format = '.2g'
    RAND.sort('CATAID')
    print("")
    RAND.pprint(max_lines = -1, max_width=-1)
    
    return

def magnitude_check() :
    
    colour_array = np.arange(-2, 5, 0.1)
    D_Ls = cosmo.luminosity_distance( np.arange(0.005,0.08,0.005) ).to(u.pc)/u.pc
    absmag_array = 17.8 - 5*np.log10(D_Ls) + 5
    
    XX, YY = np.meshgrid(colour_array, absmag_array)
    logmass = 1.15 + 0.7*XX - 0.4*YY
    
    plt.contour3d(XX, YY, logmass, 'g-i', r'$M_i$', r'$\log$ (Mass)')
    
    # limit = 1.15 + 0.7*0 - 0.4*absmag_array
    
    return

def specClass_check_BPT(sample) :
    
    if sample == 'SDSS' :
        in_path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassCam_vCam.fits'
        catalog = Table.read(in_path)
        
        log_OIII_HB = np.log10( catalog['OIII_5007_FLUX'] / catalog['H_BETA_FLUX'] )
        log_NII_HA = np.log10( catalog['NII_6584_FLUX'] / catalog['H_ALPHA_FLUX'] )
        
        SFG = (catalog['EmLineType'] == 'SFG') & (catalog['EmLineMethod'] == 'BPT')
        Comp = (catalog['EmLineType'] == 'Comp') & (catalog['EmLineMethod'] == 'BPT')
        Sey = (catalog['EmLineType'] == 'Seyfert') & (catalog['EmLineMethod'] == 'BPT')
        LINER = (catalog['EmLineType'] == 'LINER') & (catalog['EmLineMethod'] == 'BPT')
    
    if sample == 'GAMA' :
        in_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassCam_vCam.fits'
        catalog = Table.read(in_path)
        
        HB_flux_corr = (1 + 2.5/catalog['HB_EW_1'])*catalog['HB_FLUX_1']
        log_OIII_HB = np.log10( catalog['OIIIR_FLUX_1'] / HB_flux_corr )
        HA_flux_corr = (1 + 2.5/catalog['HA_EW_1'])*catalog['HA_FLUX_1']
        log_NII_HA = np.log10( catalog['NIIR_FLUX_1'] / HA_flux_corr )
        
        SFG = (catalog['EmLineType'] == 'SFG') & (catalog['EmLineMethod'] == 'BPT')
        Comp = (catalog['EmLineType'] == 'Comp') & (catalog['EmLineMethod'] == 'BPT')
        Sey = (catalog['EmLineType'] == 'Seyfert') & (catalog['EmLineMethod'] == 'BPT')
        LINER = (catalog['EmLineType'] == 'LINER') & (catalog['EmLineMethod'] == 'BPT')
    
    # totBPT = np.sum(catalog['BPT'])
    # print(totBPT)
    # totSFG, totComp, totSey, totLIN = np.sum(SFG), np.sum(Comp), np.sum(Sey), np.sum(LINER)
    # print(totSFG, totComp, totSey, totLIN)
    # print(totBPT - (totSFG + totComp + totSey + totLIN) )
    # not_ELG = (catalog['EmLineType'] == 'not_ELG') & (catalog['BPT'] == True)
    # print(np.sum(not_ELG))
    
    plt.diagram_BPT(log_NII_HA[SFG], log_OIII_HB[SFG],
                    log_NII_HA[Comp], log_OIII_HB[Comp],
                    log_NII_HA[Sey], log_OIII_HB[Sey],
                    log_NII_HA[LINER], log_OIII_HB[LINER])
    
    return

def specClass_check_WHAN(sample) :
    
    if sample == 'SDSS' :
        in_path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassCam_vCam.fits'
        catalog = Table.read(in_path)
        
        HA_width = np.log10( catalog['H_ALPHA_FLUX'] / catalog['H_ALPHA_CONT'] )
        log_NII_HA = np.log10( catalog['NII_6584_FLUX'] / catalog['H_ALPHA_FLUX'] )
        
        Pass = (catalog['EmLineType'] == 'Passive') & (catalog['EmLineMethod'] == 'WHAN')
        SFG = (catalog['EmLineType'] == 'SFG') & (catalog['EmLineMethod'] == 'WHAN')
        Comp = (catalog['EmLineType'] == 'Comp') & (catalog['EmLineMethod'] == 'WHAN')
        Sey = (catalog['EmLineType'] == 'Seyfert') & (catalog['EmLineMethod'] == 'WHAN')
        LINER = (catalog['EmLineType'] == 'LINER') & (catalog['EmLineMethod'] == 'WHAN')
    
    if sample == 'GAMA' :
        in_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassCam_vCam.fits'
        catalog = Table.read(in_path)
        
        HA_width = np.log10(catalog['HA_EW_1'])
        HA_flux_corr = (1 + 2.5/catalog['HA_EW_1'])*catalog['HA_FLUX_1']
        log_NII_HA = np.log10( catalog['NIIR_FLUX_1'] / HA_flux_corr )
        
        Pass = (catalog['EmLineType'] == 'Passive') & (catalog['EmLineMethod'] == 'WHAN')
        SFG = (catalog['EmLineType'] == 'SFG') & (catalog['EmLineMethod'] == 'WHAN')
        Comp = (catalog['EmLineType'] == 'Comp') & (catalog['EmLineMethod'] == 'WHAN')
        Sey = (catalog['EmLineType'] == 'Seyfert') & (catalog['EmLineMethod'] == 'WHAN')
        LINER = (catalog['EmLineType'] == 'LINER') & (catalog['EmLineMethod'] == 'WHAN')
    
    # totWHAN = np.sum(catalog['WHAN'])
    # print(totWHAN)
    # totPass, totSFG, totComp, totSey, totLIN = np.sum(Pass), np.sum(SFG), np.sum(Comp), np.sum(Sey), np.sum(LINER)
    # print(totPass, totSFG, totComp, totSey, totLIN)
    # print(totWHAN - (totPass + totSFG + totComp + totSey + totLIN) )
    # not_ELG = (catalog['EmLineType'] == 'not_ELG') & (catalog['WHAN'] == True)
    # print(np.sum(not_ELG))
    
    plt.diagram_WHAN(log_NII_HA[Pass], HA_width[Pass],
                     log_NII_HA[SFG], HA_width[SFG],
                     log_NII_HA[Comp], HA_width[Comp],
                     log_NII_HA[Sey], HA_width[Sey],
                     log_NII_HA[LINER], HA_width[LINER])
    
    return

#consistency_check('GAMA_alt') # this doesn't work

# specClass_check_BPT('SDSS')
# specClass_check_WHAN('SDSS')
# specClass_check_BPT('GAMA')
# specClass_check_WHAN('GAMA')
