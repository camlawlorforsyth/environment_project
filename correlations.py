
def correlations() :
    
    catalog = Table.read('catalogs/CARS_SDSS/CARS_SDSS_complete_masslimited.fits')
    
    mass = catalog['log_mass']
    companions = catalog['Companions']
    massive_mass = catalog['Most_Massive_Mass']
    massive_dist = catalog['Most_Massive_Distance']
    closest_mass = catalog['Closest_Mass']
    closest_dist = catalog['Closest_Distance']
    SD = catalog['SurfaceDensity']
    SD_err = catalog['SurfaceDensityErr']
    GAMA_CiC = catalog['CountInCyl']
    GAMA_CiC_err = catalog['CountInCylErr']
    overdensity = catalog['Overdensity']
    excess = catalog['Excess']
    agep = catalog['AGEPar']
    
    mass_label = r'CARS Host Mass $(\log_{10}[M_\odot])$'
    companions_label = 'Number of Companions'
    massive_mass_label = r'Most Massive Companion Mass $(\log_{10}[M_\odot])$'
    massive_dist_label = 'Distance to Most Massive Companion (kpc)'
    closest_mass_label = r'Closest Companion Mass $(\log_{10}[M_\odot])$'
    closest_dist_label = 'Distance to Closest Companion (kpc)'
    SD_label = r'Surface Density (Mpc$^{-2}$)'
    GAMA_CiC_label = 'Companions in the GAMA Cyl.'
    overdensity_label = 'Overdensity'
    excess_label = 'Excess'
    agep_label = r'AGE Parameter (Mpc$^{-1}$)'
    
    param_list = [mass, companions, massive_mass, massive_dist,
                  closest_mass, closest_dist, SD, GAMA_CiC, overdensity, agep]
    param_labels = [mass_label, companions_label, massive_mass_label,
                    massive_dist_label, closest_mass_label, closest_dist_label,
                    SD_label, GAMA_CiC_label, overdensity_label, agep_label]
    
    # for i in [0] :#range(len(param_list)) :
    #     plt.plot(param_list[i], param_labels[i], mass, mass_label)
    #     plt.plot(param_list[i], param_labels[i], companions, companions_label)
    #     plt.plot(param_list[i], param_labels[i], massive_mass, massive_mass_label)
    #     plt.plot(param_list[i], param_labels[i], massive_dist, massive_dist_label)
    #     plt.plot(param_list[i], param_labels[i], closest_mass, closest_mass_label)
    #     plt.plot(param_list[i], param_labels[i], closest_dist, closest_dist_label)
    #     plt.plot(param_list[i], param_labels[i], SD, SD_label)
    #     plt.plot(param_list[i], param_labels[i], GAMA_CiC, GAMA_CiC_label)
    #     plt.plot(param_list[i], param_labels[i], overdensity, overdensity_label)
    #     plt.plot(param_list[i], param_labels[i], agep, agep_label)
    
    # clear correlations with understandable explanations
    plt.plot(mass, mass_label, closest_mass, closest_mass_label)
    plt.plot(mass, mass_label, closest_dist, closest_dist_label)
    plt.plot(closest_mass, closest_mass_label, closest_dist, closest_dist_label)
    
    return

# correlations()
