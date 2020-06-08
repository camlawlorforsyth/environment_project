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

from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
from astropy.table import hstack, Table, vstack
import astropy.units as u
import warnings
warnings.filterwarnings('ignore') # ignore warnings about division by 0 when
# taking z/z_err >= 3, and ignore integration warning for the cosmology

import checks as chk
import environmental_parameters as env_params
import functions as funcs
import plots as plt
import search as srch

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
currentFig = 1 # first figure will be numbered as 'Figure 1'
mass_limit = 8.452021 # [dex(M_*)] ie. log(M_*) = 8.452021, 104 km/s away from HE 2222-0026
SDSS_path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassCam_logMass_vCam.fits'
GAMA_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassCam_logMass_vCam.fits'

# mass-limit using color-based mass estimate
if (mass_limit > 0) :
    limited = '_masslimited'
else :
    limited = ''

# open the CARS catalog and populate values
galaxies = Table.read('catalogs/raw_cats/CARS_thesis_sample.fits')
IDs = galaxies['Name']
RAs = Angle(galaxies['RA'], u.deg)
Decs = Angle(galaxies['DEC'], u.deg)
redshifts = galaxies['Z']
g_colors = galaxies['g_mag']
i_colors = galaxies['i_mag']

Dists = cosmo.angular_diameter_distance(redshifts) # compute D_A
lum_dists = cosmo.luminosity_distance(redshifts).to(u.pc)
M_g = g_colors - 5*np.log10(lum_dists.value) + 5
M_i = i_colors - 5*np.log10(lum_dists.value) + 5
masses = 1.15 + 0.7*(g_colors - i_colors) - 0.4*M_i

def main(cat_name, mass_check=False) :
    
    # complete the search and determine values for all CARS host galaxies
    
    base = save_base_table(cat_name)
    
    if (cat_name == 'SDSS') :
        indexes = range(len(galaxies))
        path = SDSS_path
    if (cat_name == 'GAMA') :
        indexes = [8]
        path = GAMA_path
    
    list_of_sub_cats = []
    tables = []
    for index in indexes :
        sub_catalog, table = env_params.gama_params(cat_name, path, RAs[index],
                                                    Decs[index], redshifts[index],
                                                    Dists[index], IDs[index],
                                                    masses[index], self_in_search=False)
        list_of_sub_cats.append(sub_catalog)
#        bins = int(np.ceil(np.sqrt( len(sub_catalog) )))
#        plt.histo(sub_catalog['logMass'], r'$\log_{10}(\rm M_{*}/M_{\odot})$',
#                    IDs[index], masses[index])        
        tables.append(table)
    
#    if (mass_check==True) :
#        for index in [10] :
#            sub_catalog, row = gama_params(cat_name, path, RAs[index], Decs[index],
#                                           redshifts[index],Dists[index],IDs[index])
#            list_of_sub_cats.append(sub_catalog)
#    #        bins = int(np.ceil(np.sqrt( len(sub_catalog) )))
##            plt.histo(sub_catalog['logMass'], r'$\log_{10}(\rm M_{*}/M_{\odot})$',
##                        IDs[index], masses[index])
#            
#            print(redshifts[10]*const.c.to('km/s'))
#            print()
#            v_diff = (redshifts[10]-sub_catalog['Z'])*const.c.to('km/s')
#            mass = sub_catalog['logMass']
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
    CARS['sigma_2D'].description='number of companions in default aperture over aperture area'
    CARS['rho_3D'].description='number of companions in default aperture over volume of aperture'
    CARS['d_CoM'].description='projected distance to the center of mass of the system from the CARS target'
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
#    logMass = master_table['logMass']/u.solMass
#    mass = master_table['log_mass']
#    plt.plot(logMass, 'mass from colour', mass, 'mass from catalog')
    
#    from scipy.optimize import curve_fit
#    from scipy.stats import chisquare
    
#    print("\nUsing linear function with free slope")
#    popt_lin, pcov_lin = curve_fit(funcs.line, logMass, mass)
#    print("m=%.4g  b=%.4g" % tuple(popt_lin))
#    expected = funcs.line(logMass, popt_lin[0], popt_lin[1])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(logMass))**2)
#    print("R^2=%.4g" % R_squared)
#    plt.plot(logMass, 'mass from colour', (mass - popt_lin[1])/popt_lin[0],
#               'corrected mass from catalog')
    
#    print("\nUsing linear function with slope of unity")
#    popt_int, pcov_int = curve_fit(funcs.intercept, logMass, mass)
#    print("b=%.4g" % popt_int)
#    expected = funcs.intercept(logMass, popt_int[0])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(logMass))**2)
#    print("R^2=%.4g" % R_squared)
#    plt.plot(logMass, 'mass from colour', mass - popt_int[0],
#               'corrected mass from catalog ')
    
#    print("\nUsing parabola with free parameters")
#    popt_parab, pcov_parab = curve_fit(funcs.parabola, logMass, mass)
#    print("A=%.4g  t=%.4g  b=%.4g" % tuple(popt_parab))
#    expected = funcs.parabola(logMass, popt_parab[0], popt_parab[1], popt_parab[2])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(logMass))**2)
#    print("R^2=%.4g" % R_squared)
#    plt.plot(logMass, 'mass from colour',
#               np.sqrt( (mass - popt_parab[2])/popt_parab[0] ) + popt_parab[1],
#               'corrected mass from catalog')
    
#    cat_surf_dens, cat_ageden, cat_counts = comparison()
#    plt.multi2(np.log10(cat_surf_dens), cat_ageden, 'GAMA',
#                 np.log10(CARS['SurfaceDensity']),
#                 CARS['AGEPar'], 'CARS', r'$\log_{10} (\frac{Surface Density}{Mpc^{-2}})$',
#                 r'AGE Density (Mpc$^{-1}$)')
    
    return

def final_catalog_operations(cat_name) :
    
    if cat_name == 'SDSS' :
        path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassCam_vCam.fits'
    if cat_name == 'GAMA' :
        path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassCam_vCam.fits'
    
    catalog = Table.read(path)
    
    if (cat_name == 'SDSS') :
        D_A = cosmo.angular_diameter_distance(catalog['Z'])
        D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)/u.pc
        M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        logMass = (1.15 + 0.7*(catalog['KCOR_MAG'][:,0] -
                               catalog['KCOR_MAG'][:,2]) - 0.4*M_i)*u.solMass
    
    if (cat_name == 'GAMA') :
        D_A = cosmo.angular_diameter_distance(catalog['Z_1'])
        catalog.rename_column('Z_1', 'Z')
        fluxscale = catalog['fluxscale']
        condition = ((fluxscale >= 0.1) & (fluxscale <= 10))
        M_i_corrected = catalog['absmag_i'] - 2.5*np.log10(catalog['fluxscale'])
        M_i = np.where(condition, M_i_corrected, catalog['absmag_i'])
            # see required aperture correction for stellar mass or absolute magnitude
            # http://www.gama-survey.org/dr3/data/cat/StellarMasses/v20/StellarMasses.notes
        logMass = (1.15 + 0.7*catalog['gminusi'] - 0.4*M_i)/u.mag*u.solMass
                    # based on relation from Taylor+ 2011, MNRAS, 418, 1587
    
    RAs = Angle(catalog['RA_1'], u.deg)
    decs = Angle(catalog['DEC_1'], u.deg)
    julianIDs = []
    for i in range(len(catalog)) :
        RA_str = RAs[i].to_string(unit=u.hour, sep='', precision=0, pad=True)
        dec_str = decs[i].to_string(unit=u.deg, sep='', precision=0,
                                    alwayssign=True, pad=True)
        ID = 'J' + RA_str + dec_str
        julianIDs.append(ID)
    
    # add the additional columns to the catalog
    catalog['D_A'] = D_A
    catalog['logMass'] = logMass
    catalog['Name'] = julianIDs
    
    mask = (logMass >= mass_limit*u.solMass) & (logMass < 13.5*u.solMass) # mask
        # based on unphysical stellar masses and the defined stellar mass cut
    catalog = catalog[mask]
    
    if cat_name == 'SDSS' :
        outpath = 'catalogs/joined_cats/SDSS_gal-info_gal-line_SpecClassCam_logMass_vCam.fits'
    if cat_name == 'GAMA' :
        outpath = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_SpecClassCam_logMass_vCam.fits'
    
    catalog.write(outpath, overwrite=False)
    
    return

def save_base_table(sample) :
    
    columns = [IDs, RAs, Decs, redshifts, Dists, lum_dists.to(u.Mpc),
               Angle( ((2*u.Mpc)/Dists).value, u.radian ).to('arcmin'),
               M_g*u.mag, M_i*u.mag, masses*u.solMass]
    
    column_headings = ('CARS_Host', 'RA', 'DEC', 'Z', 'D_A', 'D_L', 
                       '2_Mpc_Radius', 'M_g', 'M_i', 'logMass')
    
    comment = 'Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3'
    
    if sample == 'SDSS' :
        CARS_base = Table(columns, names=column_headings)
        CARS_base.meta['comments'] = comment
        CARS_base.write('catalogs/CARS_SDSS/CARS_SDSS_base.fits',
                        overwrite=True)
        # CARS_base.write('catalogs/CARS_SDSS/CARS_SDSS_base.tex',
                        # format='ascii.latex', overwrite=True)
        base_table = CARS_base
    
    if sample == 'GAMA' :
        CARS_GAMA_base = Table.read('catalogs/CARS_SDSS/CARS_SDSS_base.fits')
        CARS_GAMA_base = Table(CARS_GAMA_base[8], names=column_headings)
        CARS_GAMA_base.write('catalogs/CARS_GAMA/CARS_GAMA_base.fits',
                             overwrite=True)
        # CARS_GAMA_base.write('catalogs/CARS_GAMA/CARS_GAMA_base.tex',
                             # format='ascii.latex', overwrite=True)
        base_table = CARS_GAMA_base
    
    return base_table

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

# main('SDSS')
# main('GAMA')

# whole_cat_params('SDSS')
