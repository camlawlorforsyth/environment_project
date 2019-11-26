# -*- coding: utf-8 -*-
"""
    MASTERS PROJECT - PART 1
    ASSIGNMENT: Search for physically close companions to CARS host galaxies
    AUTHOR:     Cam Lawlor-Forsyth (lawlorfc@myumanitoba.ca)
    SUPERVISOR: Chris O'Dea
    VERSION:    2019-Sept-23
    
    PURPOSE: Search for physically close companion objects to CARS host
             galaxies, within 2 Mpc projected, and +/-1500 km/s along the LOS.
"""

# imports
import numpy as np

import astropy.constants as const
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import ascii
from astropy.table import Table, vstack, hstack, Column
import astropy.units as u
import warnings
warnings.filterwarnings("ignore") # ignore warnings about division by 0 when
# taking z/z_err >= 3, and ignore integration warning for the cosmology

import functions as funcs

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
                   (2*u.Mpc/Dists)*(180*60*u.arcmin/np.pi),
                   g_colors*u.mag, i_colors*u.mag, abs_i_colors*u.mag,
                   masses*u.solMass],
                  names=('CARS_Host', 'RA', 'DEC', 'Z', 'D_A', 'D_L',
                         '2_Mpc_Radius', 'g', 'i', 'M_i', 'log_mass'))
CARS_base.meta['comments'] = ['Flat \u039BCDM cosmology: H\u2080 = 70 km ' +
                              's\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3']
#CARS_base.write('catalogs/CARS_base.fits', overwrite=True)
#CARS_base.write('CARS_base.tex', format='ascii.latex', overwrite=True)

CARS_GAMA_base = Table(rows=CARS_base[8],
                       names=('CARS_Host', 'RA', 'DEC', 'Z', 'D_A', 'D_L',
                              '2_Mpc_Radius', 'g', 'i', 'M_i', 'log_mass'))
CARS_GAMA_base.meta['comments'] = ['Flat \u039BCDM cosmology: H\u2080 = 70 km ' +
                                   's\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3']
#CARS_GAMA_base.write('catalogs/CARS_GAMA_base.fits', overwrite=True)
#CARS_GAMA_base.write('CARS_GAMA_base.tex', format='ascii.latex', overwrite=True)

# using color-based mass estimate
#mass_limit = 9.942327 # 17 km/s away from CARS target z
#mass_limit = 9.253037 # 65 km/s away from CARS target z
#mass_limit = 9.726674 # 89 km/s away from CARS target z
#mass_limit = 8.585107 # 91 km/s away from CARS target z
mass_limit = 8.452021 # 104 km/s away from CARS target z

SDSS_path = 'catalogs/gal_info_dr7_v5_2_vCam.fits' # vCam removes -9999 and large z
#SDSS_path = 'catalogs/SDSS_gal_info_Mstar_SFR.fits' # v2 removes -9999 value
GAMA_path = 'catalogs/GAMA_GFS_StelMass.fits'
#GAMA_alt = 'catalogs/GAMA_GFS_StelMass_EnvMeas.fits'
#GAMA_path = 'catalogs/GAMA_GaussFitSimple.fits'
GAMA_alt = 'catalogs/GAMA_EnvironmentMeasures.fits'

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
        sub_catalog, table = gama_params(cat_name, path, RAs[index],
                                         Decs[index], redshifts[index],
                                         Dists[index], IDs[index])
        list_of_sub_cats.append(sub_catalog)
#        bins = int(np.ceil(np.sqrt( len(sub_catalog) )))
#        funcs.histo(sub_catalog['g_i_Mstar'], r'$\log_{10}(\rm M_{*}/M_{\odot})$',
#                    IDs[index], masses[index])        
        tables.append(table)
    
#    if (mass_check==True) :
#        for index in [10] :
#            sub_catalog, row = gama_params(cat_name, path, RAs[index], Decs[index],
#                                           redshifts[index],Dists[index],IDs[index])
#            list_of_sub_cats.append(sub_catalog)
#    #        bins = int(np.ceil(np.sqrt( len(sub_catalog) )))
##            funcs.histo(sub_catalog['g_i_Mstar'], r'$\log_{10}(\rm M_{*}/M_{\odot})$',
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
    CARS['Most_Massive'].description='mass of most massive companion'
    CARS['Dist_to_Massive'].description='distance to most massive companion'
    CARS['Closest_Mass'].description='mass of closest companion'
    CARS['Dist_to_Closest'].description='distance to closest companion'
    CARS['SurfaceDensity'].description='surface density based on the distance to the 5th nearest neighbour'
    CARS['SurfaceDensityErr'].description='surface density uncertainty'
    CARS['CountInCyl'].description='number of (other) galaxies within cylinder of radius 1 co-moving Mpc'
    CARS['CountInCylErr'].description='Poisson error on number of galaxies in cylinder'
    CARS['Overdensity'].description='ratio of CountInCyl over the average number of galaxies within the considered volume'
    CARS['Excess'].description='difference between CountInCyl and average number of galaxies within the considered volume'
    CARS['AGEPar'].description='adaptive Gaussian environment parameter'
    CARS['AGEParErr'].description='Poisson error on the number of galaxies used to calculate AGEPar'
    CARS['AGEScale'].description='adaptive scaling factor used for the adaptive Gaussian ellipsoid'
    CARS.write('catalogs/CARS_environments_masslimited_new.fits', overwrite=True)
    
    # table that contains the basic information and the environment parameters
    combined = hstack([base, CARS])
    combined.write('catalogs/CARS_complete_masslimited_new.fits', overwrite=True)
    
#    funcs.plot(CARS['AGEPar'], 'AGEPar', CARS['AGEParErr'], 'AGEParErr')
    # use linear scaling to illustrate relationship
    
    # see Ned Taylor's 2011 paper for mass comparison between SDSS and GAMA
#    g_i_Mstar = master_table['g_i_Mstar']/u.solMass
#    mass = master_table['log_mass']
#    funcs.plot(g_i_Mstar, 'mass from colour', mass, 'mass from catalog')
    
#    from scipy.optimize import curve_fit
#    from scipy.stats import chisquare
    
#    print("\nUsing linear function with free slope")
#    popt_lin, pcov_lin = curve_fit(line, g_i_Mstar, mass)
#    print("m=%.4g  b=%.4g" % tuple(popt_lin))
#    expected = line(g_i_Mstar, popt_lin[0], popt_lin[1])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
#    print("R^2=%.4g" % R_squared)
#    funcs.plot(g_i_Mstar, 'mass from colour', (mass - popt_lin[1])/popt_lin[0],
#               'corrected mass from catalog')
    
#    print("\nUsing linear function with slope of unity")
#    popt_int, pcov_int = curve_fit(intercept, g_i_Mstar, mass)
#    print("b=%.4g" % popt_int)
#    expected = intercept(g_i_Mstar, popt_int[0])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
#    print("R^2=%.4g" % R_squared)
#    funcs.plot(g_i_Mstar, 'mass from colour', mass - popt_int[0],
#               'corrected mass from catalog ')
    
#    print("\nUsing parabola with free parameters")
#    popt_parab, pcov_parab = curve_fit(parabola, g_i_Mstar, mass)
#    print("A=%.4g  t=%.4g  b=%.4g" % tuple(popt_parab))
#    expected = parabola(g_i_Mstar, popt_parab[0], popt_parab[1], popt_parab[2])
#    chisq, pval = chisquare(mass, expected)
#    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
#    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
#    print("R^2=%.4g" % R_squared)
#    funcs.plot(g_i_Mstar, 'mass from colour',
#               np.sqrt( (mass - popt_parab[2])/popt_parab[0] ) + popt_parab[1],
#               'corrected mass from catalog')
    
#    cat_surf_dens, cat_ageden, cat_counts = comparison()
#    funcs.multi2(np.log10(cat_surf_dens), cat_ageden, 'GAMA',
#                 np.log10(CARS['SurfaceDensity']),
#                 CARS['AGEPar'], 'CARS', r'$\log_{10} (\frac{Surface Density}{Mpc^{-2}})$',
#                 r'AGE Density (Mpc$^{-1}$)')
    
    return

#.............................................................adaptive_gaussian
def adaptive_gaussian(CARS_coords, CARS_dist, catalog) :
    
    sigma = 2*u.Mpc
    
    distances = cosmo.angular_diameter_distance(catalog['Z'])
    r_zs = abs(distances - CARS_dist)
    scale = cosmo.kpc_proper_per_arcmin(catalog['Z'])
    r_as = (scale*catalog['d2d']).to(u.Mpc)
    
    nn = sum(catalog['d3d'].to(u.Mpc) < sigma)
    if nn > 10 :
        nn = 10
    AGEScale = 1 + 0.2*nn
    
    ellipse = (r_as/(3*sigma))**2 + (r_zs/(AGEScale*3*sigma))**2
    mask = (ellipse <= 1)
    catalog = catalog[mask]
    AGEErr = np.sqrt(len(catalog))
    
    distances = cosmo.angular_diameter_distance(catalog['Z'])
    new_r_zs = abs(distances - CARS_dist)
    scale = cosmo.kpc_proper_per_arcmin(catalog['Z'])
    new_r_as = (scale*catalog['d2d']).to(u.Mpc)
    
    indv = np.exp( -0.5*( (new_r_as/sigma)**2 + (new_r_zs/(AGEScale*sigma))**2 ) )
    AGEDenPar = 1/np.sqrt(2*np.pi)/sigma * np.sum(indv)
    
    return AGEDenPar, AGEErr, AGEScale

#................................................................catalog_search
def catalog_search(catalog_path, cat_name, CARS_host, CARS_sky_coords,
                   redshift, delta_v, radius) :
    
    lower_z = redshift - delta_v*(u.km/u.s)/const.c.to('km/s')
    upper_z = redshift + delta_v*(u.km/u.s)/const.c.to('km/s')
    
    D_A = cosmo.angular_diameter_distance(redshift)
    radius = (radius*u.Mpc/D_A)*(180*60*u.arcmin/np.pi) # find the angular radius
    
    catalog = Table.read(catalog_path)
    
    if (cat_name == 'SDSS') :
        mask = (
                (lower_z <= catalog['Z']) & (catalog['Z'] <= upper_z) &
                ((catalog['Z'] / catalog['Z_ERR']) >= 3)
                ) # mask based on velocity cuts, redshift quality
    if (cat_name == 'GAMA') or (cat_name == 'GAMA_alt') :
        catalog.rename_column('Z_1', 'Z')
        catalog.rename_column('RA_1', 'RA')
        catalog.rename_column('DEC_1', 'DEC')
        catalog.rename_column('NQ_1', 'NQ')
        mask = (
                (lower_z <= catalog['Z']) & (catalog['Z'] <= upper_z) &
                (catalog['NQ'] >= 3)
                ) # mask based on velocity cuts, redshift quality
    catalog = catalog[mask] # mask the catalog once to speed things up
    
    D_A = cosmo.angular_diameter_distance(catalog['Z'])
    positions = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], distance=D_A,
                         unit=(u.deg, u.deg, u.Mpc) ) # get sky positions
    
    d2d = CARS_sky_coords.separation(positions).to(u.arcmin) # find sky sep.s
    new_mask = (d2d <= radius) # mask based on 2D sky separation
    catalog = catalog[new_mask]
    
    d3d = CARS_sky_coords.separation_3d(positions)[new_mask] # find physical sep.s
    
    if (cat_name == 'SDSS') :
        D_L = cosmo.luminosity_distance(catalog['Z']).to(u.pc)/u.pc
        M_i = catalog['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        g_i_Mstar = (1.15 + 0.7*(catalog['KCOR_MAG'][:,0] -
                                 catalog['KCOR_MAG'][:,2]) - 0.4*M_i)*u.solMass
    if (cat_name == 'GAMA') or (cat_name == 'GAMA_alt') :
        g_i_Mstar = (1.15 + 0.7*catalog['gminusi'] -
                     0.4*catalog['absmag_i'])/u.mag*u.solMass
            # based on relation from Taylor+ 2011, MNRAS, 418, 1587
    
    # see required aperture correction for stellar mass or absolute magnitude
    # http://www.gama-survey.org/dr3/data/cat/StellarMasses/v20/StellarMasses.notes
    
    host_list = []
    for i in range(len(catalog)) :
        host_list.append(str(CARS_host))
    
    # add the additional columns to the catalog
    catalog['D_A'] = D_A[new_mask]
    catalog['d2d'] = d2d[new_mask]
    catalog['d3d'] = d3d
    catalog['g_i_Mstar'] = g_i_Mstar
    catalog['CARS_host'] = host_list
    
#    if (cat_name == 'SDSS') :
#        catalog['log_mass'] = catalog['MEDIAN_2']
    if (cat_name == 'GAMA') :
        catalog['log_mass'] = catalog['logmstar']
    
    return catalog

#....................................................................comparison
def comparison() :
    
    catalog = Table.read('catalogs/GAMA_SpecClassGordon_EnvMeas.fits')
    
    mask = (
            (catalog['EmLineType'] == 'Seyfert') | (catalog['EmLineType'] == 'BLAGN  ')
            )
    catalog = catalog[mask]
    
#    dist = catalog['DistanceTo5nn']
    surf_dens = catalog['SurfaceDensity']
    counts = catalog['CountInCyl']
    ageden_par = catalog['AGEDenPar']
    
#    funcs.plot(np.log10(dist), 'log10(Distance to 5NN)',
#               np.log10(surf_dens), r'log10(Surface Density (Mpc$^{-2}$))',
#               cat_name='GAMA')
    
#    funcs.histo2d(np.log10(dist), 'log10(Distance to 5NN)',
#                  np.log10(surf_dens), r'log10(Surface Density (Mpc$^{-2}$))')
    
    return surf_dens, counts, ageden_par

#.............................................................consistency_check
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
        sub_catalog, row = gama_params(cat_name, GAMA_path, ras[index],
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

#............................................................counts_in_cylinder
def counts_in_cylinder(redshift, catalog) :
    
    # GAMA Counts in Cylinder Environmental Measure
    # requires 1 Mpc radius, delta(velocity) = 1000 km/s
    
    CountInCyl = len(catalog)
    CountInCylErr = np.sqrt(CountInCyl)
    
    lo = redshift - 1000*(u.km/u.s)/const.c.to('km/s')
    hi = redshift + 1000*(u.km/u.s)/const.c.to('km/s')
    D_As = cosmo.angular_diameter_distance([lo, hi])
    volume = (np.pi)*(1*u.Mpc)*(1*u.Mpc)*abs(D_As[1] - D_As[0])
    nbar_ref = 65/(972.1576816002101*u.Mpc**3) # total counts in cylinders over
        # total (1Mpc**2)*(z-1000km/s)*(z+1000km/s) volume for each target
#    nbar_ref = 6.916794832380327e-05*u.Mpc**(-3) # from MPA/JHU and D_A of
        # highest redshift object as radius of sphere
#    nbar_ref = 0.00911*u.Mpc**(-3) # from GAMA catalog
    avg_count = nbar_ref*volume
    OverDensity = CountInCyl/avg_count
    Excess = CountInCyl - avg_count
    
    return CountInCyl, CountInCylErr, OverDensity, Excess

#...................................................................helper_fxns
def intercept(xx, b) :
    return xx + b

def line(xx, m, b) :
    return m*xx + b

def parabola(xx, amp, translation, offset) :
    return amp*(xx - translation)**2 + offset

#...................................................................gama_params
def gama_params(cat_name, path, alpha, delta, zz, D_A, ID) :
    
    center = SkyCoord(ra=alpha, dec=delta, distance=D_A)
        # galaxy of interest
    
    default_catalog = catalog_search(path, cat_name, ID, center, zz, 1500, 2)
    last_mask = (default_catalog['g_i_Mstar'] >= mass_limit)
    default_catalog = default_catalog[last_mask]
    default_catalog.sort('g_i_Mstar')
    massive = default_catalog['g_i_Mstar'][-1]
    massive_dist = (default_catalog['d3d'].to(u.kpc))[-1]
    default_catalog.sort('d3d')
    closest_mass = default_catalog['g_i_Mstar'][0]
    closest = (default_catalog['d3d'].to(u.kpc))[0]
    
    age_catalog = catalog_search(path, cat_name, ID, center, zz, 4931.58, 7)
    age_par, age_par_err, age_scale = adaptive_gaussian(center, D_A, age_catalog)
    
    surface_catalog = catalog_search(path, cat_name, ID, center, zz, 1000, 7)
    surface_dens, surface_dens_err = surface_density(surface_catalog)
    
    cylinder_catalog = catalog_search(path, cat_name, ID, center, zz, 1000, 1)
    counts, counts_err, overdensity, excess = counts_in_cylinder(zz, cylinder_catalog)
    
#    funcs.plot(companion_catalog['g_i_Mstar']/u.solMass,
#               r'$\log(M_*/M_\odot)$ = 1.15 + 0.70($g-i$) - 0.4$M_i$',
#               companion_catalog['MEDIAN_2'],
#               r'$M_*$ from catalog [$\log(M_\odot)$]', cat_name)
#    funcs.histo(companion_catalog['MEDIAN_2'], r'$M_*$ from catalog [$\log(M_\odot)$]')
#    funcs.histo(companion_catalog['d3d']/u.Mpc, 'Physical Separation (Mpc)')
    
    # create small table for CARS host that contains relevent information
#    data = (len(default_catalog), massive, closest.value, closest_mass,
#            surface_dens.value,
#            surface_dens_err.value, counts, counts_err, overdensity.value,
#            excess, age_par.value, age_par_err, age_scale)
    
    envs = Table()
    envs['Companions'] = Column([len(default_catalog)])
    envs['Most_Massive'] = Column([massive], unit=u.solMass)
    envs['Dist_to_Massive'] = Column([massive_dist.value], unit=u.kpc)
    envs['Closest_Mass'] = Column([closest_mass], unit=u.solMass)
    envs['Dist_to_Closest'] = Column([closest.value], unit=u.kpc)
    envs['SurfaceDensity'] = Column([surface_dens.value], unit=u.Mpc**(-2))
    envs['SurfaceDensityErr'] = Column([surface_dens_err.value], unit=u.Mpc**(-2))
    envs['CountInCyl'] = Column([counts])
    envs['CountInCylErr'] = Column([counts_err])
    envs['Overdensity'] = Column([overdensity.value])
    envs['Excess'] = Column([excess])
    envs['AGEPar'] = Column([age_par.value], unit=u.Mpc**(-1))
    envs['AGEParErr'] = Column([age_par_err])
    envs['AGEScale'] = Column([age_scale])
    
    return default_catalog, envs

#...................................................................print_table
def print_table(table_to_print) :
    
    table = Table([host_list, Angle(catalog['RA'],u.deg),
                   Angle(catalog['DEC'],u.deg),
            catalog['Z'], catalog['Z']*const.c.to('km/s'), D_A[new_mask],
            d2d[new_mask], d3d
                   ],
            names=('CARS Host','RA','Dec','z','Velocity','D_A','Proj. Sep.',
                   'Phys. Sep.')
                  )
    table['z'].format = '10.6f'
    table['Velocity'].format = '8.0f'
    table['D_A'].format = '8.2f'
    table['Proj. Sep.'].format = '.3f'
    table['Phys. Sep.'].format = '.3f'
    
    if (cat_name == 'SDSS') :
        table['Type'] = catalog['TARGETTYPE']
        table['log(g-i M*)'] = g_i_Mstar
        table['log(M*)'] = catalog['MEDIAN_2']*u.solMass
        table['Type'].format = '.7s'        
    if (cat_name == 'GAMA') :
        table['log(g-i M*)'] = g_i_Mstar
        table['log(M*)'] = catalog['logmstar']
        table['log(M*)'].unit = u.solMass
    
    table.sort('Proj. Sep.') # sort the table by projected separation
    
    table.pprint(max_lines=-1, max_width=-1) # print full table
    print(table)
    
    return

#...............................................................surface_density
def surface_density(catalog) :
    
    # see Sarah Brough's 2013 paper for the definition of the GAMA surface density
    
    # GAMA Surface Density Environmental Measure
    # requires delta(velocity) = 1000 km/s
    
    catalog.sort('d2d')
    if len(catalog) >= 6 :
        DistanceTo4nn = catalog['d2d'][3]*u.Mpc # Python is 0-indexed
        DistanceTo5nn = catalog['d2d'][4]*u.Mpc
        DistanceTo6nn = catalog['d2d'][5]*u.Mpc
        SurfaceDensity4nn = 4 / (np.pi * DistanceTo4nn**2)
        SurfaceDensity5nn = 5 / (np.pi * DistanceTo5nn**2)
        SurfaceDensity6nn = 6 / (np.pi * DistanceTo6nn**2)
        SurfaceDensityErr = max(abs(SurfaceDensity5nn - SurfaceDensity4nn),
                                abs(SurfaceDensity6nn - SurfaceDensity5nn))
        return SurfaceDensity5nn, SurfaceDensityErr
    else :
        return 0/(u.Mpc**2), -999/(u.Mpc**2)
#..............................................................end of functions

#cat, CARS_data = gama_params('SDSS', SDSS_path, 0)

#cat, CARS_data = gama_params('GAMA', GAMA_path, 8)
#cat.sort('d2d')
#print(cat['d3d'])

#main('SDSS')
#main('GAMA')

#consistency_check('GAMA_alt') # this doesn't work
