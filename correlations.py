
# imports
from astropy.table import Table

import plots as plt

def correlations() :
    
    catalog = Table.read('catalogs/CARS_SDSS/CARS_SDSS_complete_masslimited.fits')
    
    masses = Table.read('catalogs/raw_cats/CARS_thesis_sample_atomicMolecularMasses.fits')
    atomic = masses['HI_mass']
    molecular = masses['H2_mass']
    atomic_label = r'Mass of Atomic Hydrogren $(\log_{10}[M_\odot])$'
    molecular_label = r'Mass of Molecular Hydrogren $(\log_{10}[M_\odot])$'
    
    redshift = catalog['Z']
    g_mag = catalog['M_g']
    i_mag = catalog['M_i']
    mass = catalog['logMass']
    companions = catalog['Companions']
    massive_mass = catalog['Most_Massive_Mass']
    massive_dist = catalog['Most_Massive_Distance']
    closest_mass = catalog['Closest_Mass']
    closest_dist = catalog['Closest_Distance']
    sigma = catalog['sigma_2D']
    rho = catalog['rho_3D']
    d_CoM = catalog['d_CoM']
    SD = catalog['SurfaceDensity']
    # SD_err = catalog['SurfaceDensityErr']
    GAMA_CiC = catalog['CountInCyl']
    # GAMA_CiC_err = catalog['CountInCylErr']
    overdensity = catalog['Overdensity']
    excess = catalog['Excess']
    agep = catalog['AGEPar']
    
    z_label = 'Redshift'
    g_mag_label = 'Absolute g Magnitude'
    i_mag_label = 'Absolute i Magnitude'
    mass_label = r'CARS Host Mass $(\log_{10}[M_\odot])$'
    companions_label = 'Number of Companions'
    massive_mass_label = r'Most Massive Companion Mass $(\log_{10}[M_\odot])$'
    massive_dist_label = 'Distance to Most Massive Companion (kpc)'
    closest_mass_label = r'Closest Companion Mass $(\log_{10}[M_\odot])$'
    closest_dist_label = 'Distance to Closest Companion (kpc)'
    sigma_label = r'Surface Area Density (Mpc$^{-2}$)'
    rho_label = r'Volume Density (Mpc$^{-3}$)'
    d_CoM_label = 'Projected Distance to the Center of Stellar Mass (arcmin)'
    SD_label = r'Surface Density (Mpc$^{-2}$)'
    GAMA_CiC_label = 'Companions in the GAMA Cyl.'
    overdensity_label = 'Overdensity'
    excess_label = 'Excess'
    agep_label = r'AGE Parameter (Mpc$^{-1}$)'
    
    param_list = [redshift, g_mag, i_mag, mass, companions,
                  massive_mass, massive_dist, closest_mass,
                  closest_dist, sigma, rho, d_CoM,
                  SD, GAMA_CiC, overdensity, agep,
                  atomic, molecular]
    param_labels = [z_label, g_mag_label, i_mag_label, mass_label, companions_label,
                    massive_mass_label, massive_dist_label, closest_mass_label,
                    closest_dist_label, sigma_label, rho_label, d_CoM_label,
                    SD_label, GAMA_CiC_label, overdensity_label, agep_label,
                    atomic_label, molecular_label]
    
    for i in [0] :#range(len(param_list)) :
        plt.plot(param_list[i], param_labels[i], redshift, z_label)
        plt.plot(param_list[i], param_labels[i], g_mag, g_mag_label)
        plt.plot(param_list[i], param_labels[i], i_mag, i_mag_label)
        plt.plot(param_list[i], param_labels[i], mass, mass_label)
        plt.plot(param_list[i], param_labels[i], companions, companions_label)
        plt.plot(param_list[i], param_labels[i], massive_mass, massive_mass_label)
        plt.plot(param_list[i], param_labels[i], massive_dist, massive_dist_label)
        plt.plot(param_list[i], param_labels[i], closest_mass, closest_mass_label)
        plt.plot(param_list[i], param_labels[i], closest_dist, closest_dist_label)
        plt.plot(param_list[i], param_labels[i], sigma, sigma_label)
        plt.plot(param_list[i], param_labels[i], rho, rho_label)
        plt.plot(param_list[i], param_labels[i], d_CoM, d_CoM_label)
        plt.plot(param_list[i], param_labels[i], SD, SD_label)
        plt.plot(param_list[i], param_labels[i], GAMA_CiC, GAMA_CiC_label)
        plt.plot(param_list[i], param_labels[i], overdensity, overdensity_label)
        plt.plot(param_list[i], param_labels[i], agep, agep_label)
        plt.plot(param_list[i], param_labels[i], atomic, atomic_label)
        plt.plot(param_list[i], param_labels[i], molecular, molecular_label)
    
    # possible correlations?
    # plt.plot(mass, mass_label, closest_mass, closest_mass_label)
    # plt.plot(mass, mass_label, closest_dist, closest_dist_label)
    # plt.plot(closest_mass, closest_mass_label, closest_dist, closest_dist_label)
    
    # new possible correlations?
    # plt.plot(redshift, z_label, atomic, atomic_label)
    # plt.plot(g_mag, g_mag_label, molecular, molecular_label)
    # plt.plot(i_mag, i_mag_label, molecular, molecular_label)
    # plt.plot(mass, mass_label, molecular, molecular_label)
    # plt.plot(d_CoM, d_CoM_label, companions, companions_label)
    
    return

# correlations()
