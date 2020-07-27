
# imports
import numpy as np

import astropy.constants as const
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table, Column
import astropy.units as u

import search as srch

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
mass_limit = 9.071297

def adaptive_gaussian(CARS_coords, CARS_dist, catalog) :
    
    # GAMA 'Adaptive Gaussian' Environmental Measure
    
    sigma = 2*u.Mpc
    
    r_zs = np.absolute(catalog['D_A'] - CARS_dist).to(u.Mpc)
    scale = cosmo.kpc_proper_per_arcmin(catalog['Z'])
    r_as = (scale*catalog['d2d']).to(u.Mpc)
    
    nn = np.sum(catalog['d3d'].to(u.Mpc) < sigma)
    if nn > 10 :
        nn = 10
    AGEScale = 1 + 0.2*nn
    
    ellipse = ( r_as/(3*sigma) )**2 + ( r_zs/(AGEScale*3*sigma) )**2
    mask = (ellipse <= 1)
    catalog = catalog[mask]
    
    indv = np.exp( -0.5*( (r_as/sigma)**2 + (r_zs/(AGEScale*sigma))**2 ) )
    AGEDenPar = 1/np.sqrt(2*np.pi)/sigma * np.sum(indv)
    AGEErr = np.sqrt(len(catalog))
    
    return AGEDenPar, AGEErr, AGEScale

def basic_params(catalog) :
    
    if len(catalog) > 0 :
        mass_catalog = catalog
        mass_catalog.sort('logMass')
        most_massive_mass = mass_catalog['logMass'][-1]
        most_massive_dist = (mass_catalog['d3d'].to(u.kpc))[-1]
        
        closeness_catalog = catalog
        closeness_catalog.sort('d3d')
        closest_mass = closeness_catalog['logMass'][0]
        closest_dist = (closeness_catalog['d3d'].to(u.kpc))[0]
        
        massive_1Mpc_catalog = catalog
        mask = (catalog['d3d'].to(u.kpc) <= 1000*u.kpc)
        massive_1Mpc_catalog = massive_1Mpc_catalog[mask]
        if len(massive_1Mpc_catalog) > 0 :
            massive_1Mpc_catalog.sort('logMass')
            massive_1Mpc_mass = massive_1Mpc_catalog['logMass'][-1]
            massive_1Mpc_dist = (massive_1Mpc_catalog['d3d'].to(u.kpc))[-1]
        else :
            massive_1Mpc_mass = np.nan
            massive_1Mpc_dist = np.nan*u.kpc
    else :
        most_massive_mass = np.nan
        most_massive_dist = np.nan*u.kpc
        closest_mass = np.nan
        closest_dist = np.nan*u.kpc
        massive_1Mpc_mass = np.nan
        massive_1Mpc_dist = np.nan*u.kpc
    
    return most_massive_mass, most_massive_dist, closest_mass, closest_dist, massive_1Mpc_mass, massive_1Mpc_dist

def center_of_mass_calc(companions_RA, companions_DEC, companions_D_As,
                        companions_log_masses, CARS_position, CARS_mass) :
    
    positions = SkyCoord(ra = companions_RA, dec = companions_DEC,
                         distance = companions_D_As,
                         unit=(u.deg, u.deg, u.Mpc) )
    
    CARS_mass = np.power(10, CARS_mass)*u.solMass
    masses = np.power(10, companions_log_masses)
    total_mass = np.sum(masses)*u.solMass + CARS_mass
    
    x_cm = (np.sum(positions.cartesian.x*masses) +
            CARS_position.cartesian.x*CARS_mass)/total_mass
    y_cm = (np.sum(positions.cartesian.y*masses) +
            CARS_position.cartesian.y*CARS_mass)/total_mass
    z_cm = (np.sum(positions.cartesian.z*masses) +
            CARS_position.cartesian.z*CARS_mass)/total_mass
    
    center_of_mass = SkyCoord(x=x_cm.value, y=y_cm.value, z=z_cm.value,
                              unit='Mpc', representation_type='cartesian')
    sep = CARS_position.separation(center_of_mass).to(u.arcmin)
    sep_3d = CARS_position.separation_3d(center_of_mass).to(u.kpc)
    
    return sep

def center_of_mass_map(path, cat_name, ID, center, zz, D_A, CARS_mass) :
    
    # compute the 2D distance from the CARS host galaxy to the center of mass
    # as a function of the radius of the aperture and the velocity difference
    
    catalog = srch.search_prep(path, cat_name, ID, center)
    
    min_val = 0
    max_val = 10
    XX = 250*np.arange(min_val, max_val, 1)
    YY = 0.25*np.arange(min_val, max_val, 1)
    XX, YY = np.meshgrid(XX, YY)
    ZZ = []
    companions = []
    
    for i in range(min_val, max_val) :
        row = []
        companions_row = []
        for j in range(min_val, max_val) :
            subcat = srch.catalog_search(catalog, zz, D_A, 250*i, 0.25*j)
            num_companions = len(subcat)
            companions_row.append(num_companions)
            sub_sep = center_of_mass_calc(subcat['RA_1'], subcat['DEC_1'],
                                          subcat['D_A'], subcat['logMass'],
                                          center, CARS_mass)
            row.append(sub_sep.value)
        companions.append(companions_row)
        ZZ.append(row)
    companions = np.array(companions)
    ZZ = np.array(ZZ)
    
    outfile = ('center_of_mass_arrays/' + cat_name + '_' + ID.replace(' ', '') +
               '_' + str(len(XX)) + '_' + str(len(YY)) )
    np.savez(outfile, x=XX, y=YY, z=ZZ, a=companions) # save the arrays to a file
    
    return

def counts_in_cylinder(redshift, catalog) :
    
    # GAMA 'Counts in Cylinder' Environmental Measure
    # requires 1 Mpc radius, delta(velocity) = 1000 km/s
    
    CountInCyl = len(catalog)
    CountInCylErr = np.sqrt(CountInCyl)
    
    lo = redshift - 1000*(u.km/u.s)/const.c.to('km/s')
    hi = redshift + 1000*(u.km/u.s)/const.c.to('km/s')
    D_As = cosmo.angular_diameter_distance([lo, hi])
    volume = (np.pi)*(1*u.Mpc)*(1*u.Mpc)*abs(D_As[1] - D_As[0])
    nbar_ref = (65+12)/(972.1576816002101*(u.Mpc**3)) # total counts in cylinders over
        # total (1Mpc**2)*(z-1000km/s)*(z+1000km/s) volume for each target,
        # including the target
#    nbar_ref = 6.916794832380327e-05*(u.Mpc**(-3)) # from MPA/JHU and D_A of
        # highest redshift object as radius of sphere
#    nbar_ref = 0.00911*(u.Mpc**(-3)) # from GAMA catalog
    expected = nbar_ref*volume # expected number of galaxies given average
        # density and specific volume
    OverDensity = (CountInCyl + 1)/expected # add the CARS host into the calc.
    Excess = (CountInCyl + 1) - expected # add the CARS host into the calc.
    
    return CountInCyl, CountInCylErr, OverDensity, Excess

def gama_params(cat_name, path, alpha, delta, zz, D_A, ID, CARS_mass,
                com_map=False, self_in_search=False) :
    
    # further information regarding GAMA environmental parameters
    # www.gama-survey.org/dr3/data/cat/EnvironmentMeasures/v05/EnvironmentMeasures.notes
    
    center = SkyCoord(ra=alpha, dec=delta, distance=D_A) # galaxy of interest
    
    catalog = srch.search_prep(path, ID, center, self_in_search)
    
    default_catalog = srch.catalog_search(catalog, zz, D_A, 1500, 2)
    
    (most_massive_mass, most_massive_dist,
     closest_mass, closest_dist,
     massive_1Mpc_mass, massive_1Mpc_dist) = basic_params(default_catalog)
    
    sigma = sigma_2d(default_catalog, 2)
    rho = rho_3d(default_catalog, zz, 1500, 2)
    
    sep = center_of_mass_calc(default_catalog['RA_1'], default_catalog['DEC_1'],
                              default_catalog['D_A'], default_catalog['logMass'],
                              center, CARS_mass)
    if com_map == True :
        center_of_mass_map(path, cat_name, ID, center, zz, D_A, CARS_mass)
    
    age_catalog = srch.catalog_search(catalog, zz, D_A, 2995, 10)
    age_par, age_par_err, age_scale = adaptive_gaussian(center, D_A, age_catalog)
    
    surface_catalog = srch.catalog_search(catalog, zz, D_A, 1000, 7)
    surface_dens, surface_dens_err = surface_density(surface_catalog)
    
    cylinder_catalog = srch.catalog_search(catalog, zz, D_A, 1000, 1)
    counts, counts_err, overdensity, excess = counts_in_cylinder(zz, cylinder_catalog)
    
#    plt.plot(companion_catalog['logMass']/u.solMass,
#               r'$\log(M_*/M_\odot)$ = 1.15 + 0.70($g-i$) - 0.4$M_i$',
#               companion_catalog['MEDIAN_2'],
#               r'$M_*$ from catalog [$\log(M_\odot)$]', cat_name)
#    plt.histo(companion_catalog['MEDIAN_2'], r'$M_*$ from catalog [$\log(M_\odot)$]')
#    plt.histo(companion_catalog['d3d']/u.Mpc, 'Physical Separation (Mpc)')
    
    # create small table for CARS host that contains relevent information
#    data = (len(default_catalog), massive, closest.value, closest_mass,
#            surface_dens.value,
#            surface_dens_err.value, counts, counts_err, overdensity.value,
#            excess, age_par.value, age_par_err, age_scale)
    
    envs = Table()
    envs['Companions'] = Column([len(default_catalog)])
    envs['Most_Massive_Mass'] = Column([most_massive_mass], unit=u.solMass)
    envs['Most_Massive_Distance'] = Column([most_massive_dist.value], unit=u.kpc)
    envs['Closest_Mass'] = Column([closest_mass], unit=u.solMass)
    envs['Closest_Distance'] = Column([closest_dist.value], unit=u.kpc)
    envs['Massive_1Mpc_Mass'] = Column([massive_1Mpc_mass], unit=u.solMass)
    envs['Massive_1Mpc_Distance'] = Column([massive_1Mpc_dist.value], unit=u.kpc)
    envs['sigma_2D'] = Column([sigma.value], unit=u.Mpc**(-2))
    envs['rho_3D'] = Column([rho.value], unit=u.Mpc**(-3))
    envs['d_CoM'] = Column([sep.value], unit=u.arcmin)
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

def sigma_2d(catalog, radius) :
    
    radius = radius*u.Mpc
    sigma = len(catalog) / (np.pi * radius**2)
    
    return sigma

def rho_3d(catalog, redshift, delta_v, radius) :
    
    radius = radius*u.Mpc
    lo = redshift - delta_v*(u.km/u.s)/const.c.to('km/s')
    hi = redshift + delta_v*(u.km/u.s)/const.c.to('km/s')
    D_As = cosmo.angular_diameter_distance([lo, hi])
    volume = (np.pi * radius**2)*abs(D_As[1] - D_As[0])
    
    rho = len(catalog) / volume
    
    return rho

def surface_density(catalog) :
    
    # see Sarah Brough's 2013 paper for the definition of the GAMA surface density
    
    # GAMA 'Surface Density' Environmental Measure
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
        SD = SurfaceDensity5nn
        SD_err = SurfaceDensityErr
    elif len(catalog) >= 5 :
        DistanceTo5nn = catalog['d2d'][4]*u.Mpc
        SurfaceDensity5nn = 5 / (np.pi * DistanceTo5nn**2)
        SD = SurfaceDensity5nn
        SD_err = -999/(u.Mpc**2)
    else :
        SD = np.nan/(u.Mpc**2)
        SD_err = -999/(u.Mpc**2)
    
    return SD, SD_err
