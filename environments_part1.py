# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 1
    ASSIGNMENT: Search for physically close companions to CARS host galaxies
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2020-Apr-11
    
    PURPOSE: Search for physically close companion objects to CARS host
             galaxies, within 2 Mpc projected, and +/-1500 km/s along the LOS.
"""

# imports
import numpy as np

# import astropy.constants as const
from astropy.coordinates import Angle
# from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
# from astropy.io import ascii
from astropy.table import Table, vstack, hstack#, Column
import astropy.units as u
import warnings
warnings.filterwarnings("ignore") # ignore warnings about division by 0 when
# taking z/z_err >= 3, and ignore integration warning for the cosmology

import plots as plt
import functions as funcs
import checks as chk
import environmental_parameters as env_params
import search as srch

# the following are all in SDSS as per Bernd
galaxies = { # coordinates from SDSS DR7 catalog, z from SDSS DR7, otherwise NED
            'HE 0040-1105':{'RA':'10.65358672','dec':'-10.82278912','z':0.0419,
                            'g':15.873565, 'i':15.084779},
            'HE 0114-0015':{'RA':'19.26494998','dec':'7.61312E-3','z':0.0456,
                            'g':15.659997, 'i':14.723861}, # RBS 175
            'HE 0119-0118':{'RA':'20.49921832','dec':'-1.04010828','z':0.0542,
                            'g':14.495553, 'i':13.845582}, # Mrk 1503
            'HE 0203-0031':{'RA':'31.56661419','dec':'-0.29144322','z':0.0426,
                            'g':14.032228, 'i':12.812394}, # Mrk 1018
            'HE 0212-0059':{'RA':'33.63982968','dec':'-0.76672567','z':0.0261,
                            'g':13.756615, 'i':12.667005}, # Mrk 590
            'HE 0227-0913':{'RA':'37.52302021','dec':'-8.99813544','z':0.01645,
                            'g':14.035069, 'i':13.183728}, # Mrk 1044, in GAMA? no DR7/8 spectra
            'HE 0232-0900':{'RA':'38.6576586','dec':'-8.78777681','z':0.04314,
                            'g':13.537529, 'i':12.712051}, # Mrk 1048, in GAMA?, no DR7/8 spectra
            'HE 0345+0056':{'RA':'56.91745592','dec':'1.08722631','z':0.03100,
                            'g':15.105369, 'i':14.4509535}, # no DR7/8 spectra
            'HE 0853+0102':{'RA':'133.97612964','dec':'0.85305195','z':0.0524,
                            'g':16.010803, 'i':15.319082,'GAMA_CATAID':278841},
            'HE 0934+0119':{'RA':'144.25436927','dec':'1.09547822','z':0.0505,
                            'g':15.794184, 'i':15.187816}, # Mrk 707
            'HE 2222-0026':{'RA':'336.14705331','dec':'-0.18441524','z':0.0581,
                            'g':17.061754, 'i':16.162113},
            'HE 2302-0857':{'RA':'346.18116153','dec':'-8.68572479','z':0.0469,
                            'g':14.508565, 'i':13.438229}, # Mrk 926
            
#            'HE 0853-0126':{'RA':'134.07430518','dec':'-1.63535577','z':0.05981}, # no DR7 spectra
#            'HE 0949-0122':{'RA':'148.07959493','dec':'-1.61209535','z':0.01993}, # Mrk 1239, no DR7 spectra
#            'HE 2128-0221':{'RA':'322.69512987','dec':'-2.14690352','z':0.05248} # no DR7 spectra
           }

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
currentFig = 1 # first figure will be numbered as 'Figure 1'

IDs = list( galaxies.keys() ) # the object name/identifier
RAs = Angle([ galaxy['RA'] for galaxy in galaxies.values() ], u.deg)
Decs = Angle([ galaxy['dec'] for galaxy in galaxies.values() ], u.deg)
redshifts = np.array([ galaxy['z'] for galaxy in galaxies.values() ])
Dists = cosmo.angular_diameter_distance(redshifts) # compute D_A
g_colors = np.array([ galaxy['g'] for galaxy in galaxies.values() ])
i_colors = np.array([ galaxy['i'] for galaxy in galaxies.values() ])
lum_dists = cosmo.luminosity_distance(redshifts).to(u.pc)/u.pc
abs_i_colors = i_colors - 5*np.log10(lum_dists) + 5
masses = 1.15 + 0.7*(g_colors - i_colors) - 0.4*abs_i_colors

CARS_base = Table([IDs, RAs, Decs, redshifts, Dists, (lum_dists*u.pc).to(u.Mpc),
                   Angle( ((2*u.Mpc)/Dists).value, u.radian ).to('arcmin'),
                   g_colors*u.mag, i_colors*u.mag, abs_i_colors*u.mag,
                   masses*u.solMass],
                  names=('CARS_Host', 'RA', 'DEC', 'Z', 'D_A', 'D_L',
                         '2_Mpc_Radius', 'g', 'i', 'M_i', 'log_mass'))
CARS_base.meta['comments'] = ['Flat \u039BCDM cosmology: H\u2080 = 70 km ' +
                              's\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3']
#CARS_base.write('catalogs/CARS_SDSS/CARS_SDSS_base.fits', overwrite=True)
#CARS_base.write('catalogs/CARS_SDSS/CARS_SDSS_base.tex', format='ascii.latex', overwrite=True)

CARS_GAMA_base = Table(rows=CARS_base[8],
                       names=('CARS_Host', 'RA', 'DEC', 'Z', 'D_A', 'D_L',
                              '2_Mpc_Radius', 'g', 'i', 'M_i', 'log_mass'))
CARS_GAMA_base.meta['comments'] = ['Flat \u039BCDM cosmology: H\u2080 = 70 km ' +
                                   's\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3']
#CARS_GAMA_base.write('catalogs/CARS_GAMA/CARS_GAMA_base.fits', overwrite=True)
#CARS_GAMA_base.write('catalogs/CARS_GAMA/CARS_GAMA_base.tex', format='ascii.latex', overwrite=True)

# mass-limit using color-based mass estimate
mass_limit = 8.452021 # 104 km/s away from HE 2222-0026
if (mass_limit > 0) :
    limited = '_masslimited'
else :
    limited = ''

SDSS_path = 'catalogs/joined_cats/SDSS_gal_info_gal_line_vCam.fits'
#SDSS_path = 'catalogs/joined_cats/SDSS_gal_info_Mstar_SFR.fits' # v2 removes -9999 value
# GAMA_path = 'catalogs/joined_cats/GAMA_GFS_StelMass.fits'
#GAMA_alt = 'catalogs/joined_cats/GAMA_GFS_StelMass_EnvMeas.fits'
#GAMA_path = 'catalogs/raw_cats/GAMA_GaussFitSimple.fits'
#GAMA_alt = 'catalogs/raw_cats/GAMA_EnvironmentMeasures.fits'
GAMA_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassGordon_vCam.fits'

# further information regarding GAMA environmental parameters
# www.gama-survey.org/dr3/data/cat/EnvironmentMeasures/v05/EnvironmentMeasures.notes

#..........................................................................main
def main(cat_name, mass_check=False) :
    
    # complete the search and determine values for all CARS host galaxies
    
    if (cat_name == 'SDSS') :
        indexes = range(len(galaxies))
        path = SDSS_path
        base = CARS_base
    if (cat_name == 'GAMA') :
        indexes = [8]
        path = GAMA_path
        base = CARS_GAMA_base
    
    list_of_sub_cats = []
    tables = []
    for index in indexes :
        sub_catalog, table = env_params.gama_params(cat_name, path, RAs[index],
                                         Decs[index], redshifts[index],
                                         Dists[index], IDs[index],
                                         np.power(10, masses[index])*u.solMass,
                                         self_in_search=False)
        list_of_sub_cats.append(sub_catalog)
#        bins = int(np.ceil(np.sqrt( len(sub_catalog) )))
#        plt.histo(sub_catalog['g_i_Mstar'], r'$\log_{10}(\rm M_{*}/M_{\odot})$',
#                    IDs[index], masses[index])        
        tables.append(table)
    
#    if (mass_check==True) :
#        for index in [10] :
#            sub_catalog, row = gama_params(cat_name, path, RAs[index], Decs[index],
#                                           redshifts[index],Dists[index],IDs[index])
#            list_of_sub_cats.append(sub_catalog)
#    #        bins = int(np.ceil(np.sqrt( len(sub_catalog) )))
##            plt.histo(sub_catalog['g_i_Mstar'], r'$\log_{10}(\rm M_{*}/M_{\odot})$',
##                        IDs[index], masses[index])
#            
#            print(redshifts[10]*const.c.to('km/s'))
#            print()
#            v_diff = (redshifts[10]-sub_catalog['Z'])*const.c.to('km/s')
#            mass = sub_catalog['g_i_Mstar']
#            t = Table([v_diff, mass], names=('v_diff', 'mass'))
#            t.sort('v_diff')
#            print(t)
#            
#            rows.append(row)
    
    # table that contains all information from the catalogs for the companions 
#    master_table = vstack(list_of_sub_cats)
#    master_table.pprint(max_lines=-1, max_width=-1) # print full table
#    print(master_table) # print small table
    
    # table that contains the determined parameters for the CARS galaxies
    CARS = vstack(tables)
    CARS['Companions'].description = 'number of companions in 2Mpc, +/-1500km/s cylinder'
    CARS['Most_Massive_Mass'].description='mass of most massive companion'
    CARS['Most_Massive_Distance'].description='distance to most massive companion'
    CARS['Closest_Mass'].description='mass of closest companion'
    CARS['Closest_Distance'].description='distance to closest companion'
    CARS['SurfaceDensity'].description='surface density based on the distance to the 5th nearest neighbour'
    CARS['SurfaceDensityErr'].description='surface density uncertainty'
    CARS['CountInCyl'].description='number of (other) galaxies within cylinder of radius 1 co-moving Mpc'
    CARS['CountInCylErr'].description='Poisson error on number of galaxies in cylinder'
    CARS['Overdensity'].description='ratio of CountInCyl over the average number of galaxies within the considered volume'
    CARS['Excess'].description='difference between CountInCyl and average number of galaxies within the considered volume'
    CARS['AGEPar'].description='adaptive Gaussian environment parameter'
    CARS['AGEParErr'].description='Poisson error on the number of galaxies used to calculate AGEPar'
    CARS['AGEScale'].description='adaptive scaling factor used for the adaptive Gaussian ellipsoid'
    outpath = ('catalogs/CARS_' + cat_name + '/CARS_' + cat_name +
               '_environments' + limited + '.fits')
    CARS.write(outpath, overwrite=True)
    
    # table that contains the basic information and the environment parameters
    combined = hstack([base, CARS])
    combined_outpath = ('catalogs/CARS_' + cat_name + '/CARS_' + cat_name +
                        '_complete' + limited + '.fits')
    combined.write(combined_outpath, overwrite=True)
    
#    plt.plot(CARS['AGEPar'], 'AGEPar', CARS['AGEParErr'], 'AGEParErr')
    # use linear scaling to illustrate relationship
    
    # see Ned Taylor's 2011 paper for mass comparison between SDSS and GAMA
#    g_i_Mstar = master_table['g_i_Mstar']/u.solMass
#    mass = master_table['log_mass']
#    plt.plot(g_i_Mstar, 'mass from colour', mass, 'mass from catalog')
    
#    from scipy.optimize import curve_fit
#    from scipy.stats import chisquare
    
#    print("\nUsing linear function with free slope")
#    popt_lin, pcov_lin = curve_fit(funcs.line, g_i_Mstar, mass)
#    print("m=%.4g  b=%.4g" % tuple(popt_lin))
#    expected = funcs.line(g_i_Mstar, popt_lin[0], popt_lin[1])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
#    print("R^2=%.4g" % R_squared)
#    plt.plot(g_i_Mstar, 'mass from colour', (mass - popt_lin[1])/popt_lin[0],
#               'corrected mass from catalog')
    
#    print("\nUsing linear function with slope of unity")
#    popt_int, pcov_int = curve_fit(funcs.intercept, g_i_Mstar, mass)
#    print("b=%.4g" % popt_int)
#    expected = funcs.intercept(g_i_Mstar, popt_int[0])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
#    print("R^2=%.4g" % R_squared)
#    plt.plot(g_i_Mstar, 'mass from colour', mass - popt_int[0],
#               'corrected mass from catalog ')
    
#    print("\nUsing parabola with free parameters")
#    popt_parab, pcov_parab = curve_fit(funcs.parabola, g_i_Mstar, mass)
#    print("A=%.4g  t=%.4g  b=%.4g" % tuple(popt_parab))
#    expected = funcs.parabola(g_i_Mstar, popt_parab[0], popt_parab[1], popt_parab[2])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
#    print("R^2=%.4g" % R_squared)
#    plt.plot(g_i_Mstar, 'mass from colour',
#               np.sqrt( (mass - popt_parab[2])/popt_parab[0] ) + popt_parab[1],
#               'corrected mass from catalog')
    
#    cat_surf_dens, cat_ageden, cat_counts = comparison()
#    plt.multi2(np.log10(cat_surf_dens), cat_ageden, 'GAMA',
#                 np.log10(CARS['SurfaceDensity']),
#                 CARS['AGEPar'], 'CARS', r'$\log_{10} (\frac{Surface Density}{Mpc^{-2}})$',
#                 r'AGE Density (Mpc$^{-1}$)')
    
    return

def spectral_classification() :
    
    in_path = 'catalogs/joined_cats/SDSS_gal_info_gal_line_SpecClassPrelim_vCam.fits'
    
    catalog = Table.read(in_path)
    
    log_OIII_HB = np.log10( catalog['OIII_5007_FLUX'] / catalog['H_BETA_FLUX'] )
    log_NII_HA = np.log10( catalog['NII_6584_FLUX'] / catalog['H_ALPHA_FLUX'] )
    BPT_SFG = ( (log_OIII_HB < 0.61/(log_NII_HA - 0.05) + 1.3) &
                (log_NII_HA < 0.05) )
    BPT_Comp = ( (log_OIII_HB > 0.61/(log_NII_HA - 0.05) + 1.3) &
                 (log_OIII_HB < 0.61/(log_NII_HA - 0.47) + 1.19) &
                 (log_NII_HA < 0.47) )
    BPT_Sey = ( (log_OIII_HB > 0.61/(log_NII_HA - 0.47) + 1.19) &
                (log_OIII_HB > 1.05*log_NII_HA + 0.45) )
    BPT_LIN = ( (log_OIII_HB > 0.61/(log_NII_HA - 0.47) + 1.19) &
                (log_OIII_HB < 1.05*log_NII_HA + 0.45) )
    
    HA_EW = np.log10(np.absolute(catalog['H_ALPHA_EQW']))
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
                print(log_NII_HA[i], log_OIII_HB[i])
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
    
    catalog.write('catalogs/joined_cats/SDSS_gal_info_gal_line_SpecClassCam_vCam.fits',
                   overwrite=True)
    
    return

#..............................................................whole_cat_params
def whole_cat_params(cat_name) :
    
    if (cat_name == 'SDSS') :
        path = 'catalogs/joined_cats/SDSS_gal_info_Mstar_SFR.fits'
    
    catalog = Table.read(path)
    
    if (cat_name == 'SDSS') :
        
        mask1 = ( (catalog['RA'] > -9999) & (catalog['DEC'] > -9999) &
                  (catalog['Z'] > 2e-07 ) )
        catalog = catalog[mask1]
        
        mask2 = ( (catalog['KCOR_MAG'][:,0] > 0) &
                  (catalog['KCOR_MAG'][:,0] < 30) &
                  (catalog['KCOR_MAG'][:,2] > 0) &
                  (catalog['KCOR_MAG'][:,2] < 30) )
        catalog = catalog[mask2]
        
        mass_mask = (catalog['MEDIAN_2'] > 0) # solar masses
        SFR_mask = (catalog['MEDIAN_3'] > -99) # star formation rate
        
        stellar_mass = catalog['MEDIAN_2']
        stellar_mass_label = r'$\log$ Stellar Mass from MPA/JHU ($M_{\odot}$)'
        
        SFR = catalog['MEDIAN_3']
        SFR_label = r'SFR from MPA/JHU ($M_{\odot}$ year$^{-1}$)'
        
        D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)
        D_L_cm = D_L.to('cm')
        M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L/u.pc) + 5 # absolute i mag.
        g_i_Mstar = (1.15 + 0.7*(catalog['KCOR_MAG'][:,0] -
                                 catalog['KCOR_MAG'][:,2]) - 0.4*M_i)*u.solMass
        g_i_Mstar_label = r'$\log$ Stellar Mass from Colour ($M_{\odot}$)'
        
        cat2 = Table.read('catalogs/raw_cats/SDSS_gal_line_dr7_v5_2.fit.gz')
        cat2 = cat2[mask1]
        cat2 = cat2[mask2]
        Ha_flux = cat2['H_ALPHA_FLUX']*1e-17*u.erg/u.s/(u.cm**2)
        SFR_alt_per_D_L = 7.9e-42*(u.solMass/u.yr*u.s/u.erg)*4*np.pi*D_L_cm*Ha_flux
        SFR_alt = SFR_alt_per_D_L*D_L_cm
        SFR_alt_label = r'$\log$ SFR from H$\mathrm{\alpha}$ flux ($M_{\odot}$ year$^{-1}$)'
        SFR_alt_mask = (SFR_alt > 0)
        
        m_Ha = 0.76*catalog['KCOR_MAG'][:,1] + 0.24*catalog['KCOR_MAG'][:,1]
        M_Ha = m_Ha - 5*np.log10(D_L/u.pc) + 5
        L_Ha = (100**((4.62 - M_Ha)/5))*u.solLum
        SFR_other = ((7.9e-42*L_Ha).to('erg/s'))*(u.solMass/u.yr*u.s/u.erg)
        SFR_other_label = r'$\log$ SFR from H$\mathrm{\alpha}$ magnitude ($M_{\odot}$ year$^{-1}$)'
        SFR_other_mask = ( (m_Ha > 10) & (m_Ha < 30) )
        
    if (cat_name == 'GAMA') or (cat_name == 'GAMA_alt') :
        g_i_Mstar = (1.15 + 0.7*catalog['gminusi'] -
                      0.4*catalog['absmag_i'])/u.mag*u.solMass
            # based on relation from Taylor+ 2011, MNRAS, 418, 1587
    
    # using plot - probably not the way to go
    # plt.plot(g_i_Mstar.value, 'Color-Based Mass',
                # catalog['MEDIAN_2'], 'Mass from Catalog', hist2d=True)
    # plt.plot(g_i_Mstar, 'Color-Based Mass',
                # catalog['MEDIAN_3'], 'SFR', hist2d=True)
    
    # using histo2d - probably the better option
    plt.histo2d( (g_i_Mstar[mass_mask]).value, g_i_Mstar_label,
                  stellar_mass[mass_mask], stellar_mass_label,
                  xmin=4, xmax=12, ymin=6, ymax=13)
    
    plt.histo2d(stellar_mass[mass_mask&SFR_mask], stellar_mass_label,
                  SFR[mass_mask&SFR_mask], SFR_label, xmin=5)
    
    plt.histo2d( (g_i_Mstar[SFR_mask]).value, g_i_Mstar_label,
                  SFR[SFR_mask], SFR_label)
    
    # sSFR = (SFR[mass_mask&SFR_mask])/(stellar_mass[mass_mask&SFR_mask])
    
    # plt.histo2d(stellar_mass[mass_mask&SFR_mask], stellar_mass_label,
                  # sSFR, 'sSFR', nbins=300, xmin=5, ymin=-1)
    
    plt.histo2d(SFR[SFR_mask&SFR_alt_mask], SFR_label,
                  np.log10(SFR_alt[SFR_mask&SFR_alt_mask].value), SFR_alt_label)
    
    plt.histo2d(SFR[SFR_mask&SFR_other_mask], SFR_label,
                  np.log10(SFR_other[SFR_mask&SFR_other_mask].value),
                  SFR_other_label)
    
    plt.histo2d(np.log10(SFR_alt[SFR_alt_mask&SFR_other_mask].value),
                  SFR_alt_label,
                  np.log10(SFR_other[SFR_alt_mask&SFR_other_mask].value),
                  SFR_other_label)
    
    return
#..............................................................end of functions

# main('SDSS')
# main('GAMA')

# whole_cat_params('SDSS')
