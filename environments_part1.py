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
#from astropy.io import fits
from astropy.table import Table, vstack
import astropy.units as u
import warnings
warnings.filterwarnings("ignore") # ignore warnings about division by 0 when
# taking z/z_err >= 3, and ignore integration warning for the cosmology

import functions as funcs

# the following are all in SDSS as per Bernd
galaxies = { # coordinates from SDSS DR7 catalog, z from SDSS DR7, otherwise NED
            'HE0040-1105':{'RA':'10.65358672','dec':'-10.82278912','z':0.0419},
            'RBS175':{'RA':'19.26494998','dec':'7.61312E-3','z':0.0456},
            'Mrk1503':{'RA':'20.49921832','dec':'-1.04010828','z':0.0542},
            'Mrk1018':{'RA':'31.56661419','dec':'-0.29144322','z':0.0426},
            'Mrk590':{'RA':'33.63982968','dec':'-0.76672567','z':0.0261},
            'Mrk1044':{'RA':'37.52302021','dec':'-8.99813544','z':0.01645}, # in GAMA? no DR7/8 spectra
            'Mrk1048':{'RA':'38.6576586','dec':'-8.78777681','z':0.04314}, # in GAMA?, no DR7/8 spectrum
            'HE0345+0056':{'RA':'56.91745592','dec':'1.08722631','z':0.03100}, # no DR7/8 spectrum
            'HE0853+0102':{'RA':'133.97612964','dec':'0.85305195','z':0.0524,
                          'GAMA_CATAID':278841}, # in GAMA
            'Mrk707':{'RA':'144.25436927','dec':'1.09547822','z':0.0505},
            'HE2222-0026':{'RA':'336.14705331','dec':'-0.18441524','z':0.0581},
            'Mrk926':{'RA':'346.18116153','dec':'-8.68572479','z':0.0469}
           }

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # specify the cosmology being used
currentFig = 1 # first figure will be numbered as 'Figure 1'

IDs = list( galaxies.keys() ) # the object name/identifier
ras = Angle([ galaxy['RA'] for galaxy in galaxies.values() ], u.deg)
decs = Angle([ galaxy['dec'] for galaxy in galaxies.values() ], u.deg)
zs = np.array([ galaxy['z'] for galaxy in galaxies.values() ])

dists = cosmo.angular_diameter_distance(zs) # compute D_A

SDSS_path = 'catalogs/SDSS_gal_info_Mstar_SFR.fits'
GAMA_path = 'catalogs/GAMA_GFS_Mstar.fits'

# further information regarding GAMA environmental parameters
# www.gama-survey.org/dr3/data/cat/EnvironmentMeasures/v05/EnvironmentMeasures.notes

#..........................................................................main
def main(cat_name, path, index) :
    
    center = SkyCoord(ra=ras[index], dec=decs[index], distance=dists[index])
        # galaxy of interest
    
    default_catalog = catalog_search(path, cat_name, IDs[index], center,
                                     index, zs[index], 1500, 2)
    table = Table([default_catalog['d2d'], default_catalog['d3d']],
                  names=('2D Separation', '3D Separation'))
    table.sort('2D Separation')
    table.pprint(max_lines=-1, max_width=-1)
    print('')
    
    age_par, age_par_err, age_scale = adaptive_gaussian(center, dists[index],
                                                        default_catalog)
#    print("%.3g +/- %.3g" % (age_par.value, age_par_err))
    
    surface_catalog = catalog_search(path, cat_name, IDs[index], center,
                                     index, zs[index], 1000, 7)
    surface_dens, surface_dens_err = surface_density(surface_catalog)
#    print("%.3g +/- %.3g" % (surface_dens.value, surface_dens_err.value))
    
    cylinder_catalog = catalog_search(path, cat_name, IDs[index], center,
                                      index, zs[index], 1000, 1)
    counts, counts_err, overdensity = counts_in_cylinder(zs[index],
                                                         cylinder_catalog)
#    print("%.3g +/- %.3g, overdensity = %.3g" % (counts, counts_err, overdensity))
    
#    funcs.plot(companion_catalog['g_i_Mstar']/u.solMass,
#               r'$\log(M_*/M_\odot)$ = 1.15 + 0.70($g-i$) - 0.4$M_i$',
#               companion_catalog['MEDIAN_2'],
#               r'$M_*$ from catalog [$\log(M_\odot)$]', cat_name)
#    funcs.histo(companion_catalog['MEDIAN_2'], r'$M_*$ from catalog [$\log(M_\odot)$]')
#    funcs.histo(companion_catalog['d3d']/u.Mpc, 'Physical Separation (Mpc)')
    
    # create small table for CARS host that contains relevent information
    data = (IDs[index], len(default_catalog), surface_dens.value,
            surface_dens_err.value, counts, counts_err, overdensity,
            age_par.value, age_par_err, age_scale)
    
    return default_catalog, data

#.............................................................adaptive_gaussian
def adaptive_gaussian(CARS_coords, CARS_dist, catalog) :
    
    sigma = 2*u.Mpc
    
    distances = cosmo.angular_diameter_distance(catalog['Z'])
    r_zs = abs(distances - CARS_dist)
    
    scale = cosmo.kpc_proper_per_arcmin(catalog['Z'])    
    r_as = (scale*catalog['d2d']).to(u.Mpc)
    
    mask = ( (r_zs < sigma) & (r_as < sigma) )
    
    catalog = catalog[mask]
    AGEErr = np.sqrt(len(catalog)) # this calculation is currently incorrect
    # compare with entire GAMA catalog for AGEErr and AGEPar for clarification
    
    nn = sum(catalog['d3d'].to(u.Mpc) < sigma)
    if nn > 10 :
        nn = 10
    AGEScale = 1 + 0.2*nn
    
    indv = np.exp( -0.5*( (r_as/sigma)**2 + (r_zs/(AGEScale*sigma))**2 ) )
    AGEDenPar = 1/np.sqrt(2*np.pi)/sigma * np.sum(indv)
    
    return AGEDenPar, AGEErr, AGEScale

#................................................................catalog_search
def catalog_search(catalog_path, cat_name, CARS_host, CARS_sky_coords, index,
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
    if (cat_name == 'GAMA') :
        catalog.rename_column('Z_1', 'Z')
        mask = (
                (lower_z <= catalog['Z']) & (catalog['Z'] <= upper_z) &
                (catalog['NQ_1'] >= 3)
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
    if (cat_name == 'GAMA') :
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
    
    return catalog

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
    nbar_ref = 0.00911*u.Mpc**(-3)
    OverDensity = CountInCyl/(nbar_ref * volume)
    
    return CountInCyl, CountInCylErr, OverDensity

#...................................................................full_search
def full_search(cat_name) :
    
    # complete the search and determine values for all CARS host galaxies
    
    if (cat_name == 'SDSS') :
        indexes = range(0, 12)
        path = SDSS_path
    if (cat_name == 'GAMA') :
        indexes = [8]
        path = GAMA_path
    
    list_of_sub_cats = []
    rows = []
    for index in indexes :
        sub_catalog, row = main(cat_name, path, index)
        list_of_sub_cats.append(sub_catalog)
        rows.append(row)
    
    # table that contains all information from the catalogs for the companions 
#    master_table = vstack(list_of_sub_cats)
#    master_table.pprint(max_lines=-1, max_width=-1) # print full table
#    print(master_table) # print small table
    
    # table that contains the determined parameters for the CARS galaxies
    CARS = Table(rows = rows, names=('CARS Host', '# Comp.', 'SurfaceDensity',
                                     'SurfaceDensityErr', 'CountInCyl',
                                     'CountInCylErr', 'Overdens.',
                                     'AGEPar', 'AGEParErr', 'AGEScale'))
    CARS['SurfaceDensity'] = CARS['SurfaceDensity']*u.Mpc**(-2)
    CARS['SurfaceDensityErr'] = CARS['SurfaceDensityErr']*u.Mpc**(-2)
    CARS['SurfaceDensity'].format = '.6f'
    CARS['SurfaceDensityErr'].format = '.6f'
    CARS['CountInCylErr'].format = '.3g'
    CARS['Overdens.'].format = '.3f'
    CARS['AGEPar'] = CARS['AGEPar']*u.Mpc**(-1)
    CARS['AGEPar'].format = '.3f'
    CARS['AGEParErr'].format = '.3g'
    CARS['AGEScale'].format = '.2g'
#    print("")
#    CARS.pprint(max_lines = -1, max_width=-1)
    
#    funcs.plot(CARS['AGEPar'], 'AGEPar', CARS['AGEParErr'], 'AGEParErr')
    
    # see Ned Taylor's 2011 paper for mass comparison between SDSS and GAMA
    '''
    g_i_Mstar = master_table['g_i_Mstar']/u.solMass
    mass = master_table['MEDIAN_2']
    funcs.plot(g_i_Mstar, 'mass from colour', mass, 'mass from catalog')
    
    from scipy.optimize import curve_fit
    from scipy.stats import chisquare
    
    print("\nUsing linear function with free slope")
    popt_lin, pcov_lin = curve_fit(line, g_i_Mstar, mass)
    print("m=%.4g  b=%.4g" % tuple(popt_lin))
    expected = line(g_i_Mstar, popt_lin[0], popt_lin[1])
    chisq, pval = chisquare(mass, expected)
    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
    print("R^2=%.4g" % R_squared)
    funcs.plot(g_i_Mstar, 'mass from colour', (mass - popt_lin[1])/popt_lin[0],
               'corrected mass from catalog')
    
    print("\nUsing linear function with slope of unity")
    popt_int, pcov_int = curve_fit(intercept, g_i_Mstar, mass)
    print("b=%.4g" % popt_int)
    expected = intercept(g_i_Mstar, popt_int[0])
    chisq, pval = chisquare(mass, expected)
    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
    print("R^2=%.4g" % R_squared)
    funcs.plot(g_i_Mstar, 'mass from colour', mass - popt_int[0],
               'corrected mass from catalog ')
    
    print("\nUsing parabola with free parameters")
    popt_parab, pcov_parab = curve_fit(parabola, g_i_Mstar, mass)
    print("A=%.4g  t=%.4g  b=%.4g" % tuple(popt_parab))
    expected = parabola(g_i_Mstar, popt_parab[0], popt_parab[1], popt_parab[2])
    chisq, pval = chisquare(mass, expected)
    print("chisq=%.4g  pval=%.4g" % (chisq, pval))
    R_squared = 1 - np.sum((mass - expected)**2)/np.sum((mass - np.mean(g_i_Mstar))**2)
    print("R^2=%.4g" % R_squared)
    funcs.plot(g_i_Mstar, 'mass from colour',
               np.sqrt( (mass - popt_parab[2])/popt_parab[0] ) + popt_parab[1],
               'corrected mass from catalog')
    '''
    
#    cat_surf_dens, cat_ageden, cat_counts = comparison()
#    funcs.multi2(np.log10(cat_surf_dens), cat_ageden, 'GAMA',
#                 np.log10(CARS['SurfaceDensity']),
#                 CARS['AGEPar'], 'CARS', r'log10(Surface Density (Mpc$^{-2}$))',
#                 r'AGE Density (Mpc$^{-1}$)')
    
    return

#...................................................................helper_fxns
def intercept(xx, b) :
    return xx + b

def line(xx, m, b) :
    return m*xx + b

def parabola(xx, amp, translation, offset) :
    return amp*(xx - translation)**2 + offset

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
        return 0, -999

#..................................................................sample_table
def table_of_sample() :
    
    param_list = [IDs, ras, decs, zs, zs*const.c.to('km/s'), dists, radii]
    headings = ('Object Name','RA','Dec','z','Velocity','D_A', '2 Mpc Radius')
    
    print('\nSample Properties\n')
    
    print('Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 ' +
          'Mpc\u207B\u00B9, \u03A9\u2098 = 0.3\n')
    
    table = Table(param_list, names=headings)
    table['Object Name'].format = '14s'
    table['RA'].format = '12.8f'
    table['Dec'].format = '12.8f'
    table['z'].format = '10.6f'
    table['Velocity'].format = '8.0f'
    table['D_A'].format = '8.2f'
    table['2 Mpc Radius'].format = '12.3f'
    table.sort('RA')
    table.pprint(max_lines=-1, max_width=-1)
    
    print('')
    
    return

#....................................................................comparison
def comparison() :
    
    catalog = Table.read('catalogs/GAMA_SpecClassGordon_EnvMeas.fits')
    
    mask = (
            (catalog['EmLineType'] == 'Seyfert') | (catalog['EmLineType'] == 'BLAGN')
            )
    catalog = catalog[mask]
    
    surf_dens = catalog['SurfaceDensity']
    counts = catalog['CountInCyl']
    ageden_par = catalog['AGEDenPar']
    
    return surf_dens, counts, ageden_par
#..............................................................end of functions

#table_of_sample()

#cat, CARS_data = main('SDSS', SDSS_path, 0)

#cat, CARS_data = main('GAMA', GAMA_path, 8)
#cat.sort('d2d')
#print(cat['d3d'])

full_search('SDSS')
