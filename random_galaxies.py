
# imports
import numpy as np

import astropy.constants as const
from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
from astropy.table import hstack, Table, vstack
import astropy.units as u

import environmental_parameters as env_params

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
mass_limit_default = 8.452021
fwhm_conv = 2*np.sqrt( 2*np.log(2) )
speed_of_light = const.c.to('km/s').value
SDSS_path = 'catalogs/joined_cats/SDSS_gal-info_gal-line_totlgm_totsfr_SpecClassCam_logMass_vCam.fits'
GAMA_path = 'catalogs/joined_cats/GAMA_GaussFitSimple_StellarMasses_MagPhys_SpecClassCam_logMass_vCam.fits'

def random_galaxy_candidates(GAMA_or_SDSS, outfilename, targettype, mass_limit) :
    
    if GAMA_or_SDSS == 'GAMA' :
        in_path = GAMA_path
    if GAMA_or_SDSS == 'SDSS' :
        in_path = SDSS_path
    
    base_cat = Table.read(in_path)
    
    # EmLineType==BLAGN, LINER, Seyfert, Comp, SFG, Passive, not_ELG
    # mass >= 8.452021, 0.01 <= z <= 0.06
    mask = ( (base_cat['EmLineType'] == targettype) &
             (0.01 <= base_cat['Z']) & (base_cat['Z'] <= 0.06) &
             (base_cat['logMass'] >= mass_limit) )
    
    # create the resulting masked table of good objects and save it
    comparison_objects = base_cat[mask]
    print(len(comparison_objects))
    
    # add the index columns to the catalog
    comparison_objects['index'] = np.arange(len(comparison_objects))
    
    if mass_limit == mass_limit_default :
        mass_string = ''
    if mass_limit == 10 :
        mass_string = '_logMass10'
    
    out_path = 'catalogs/candidate_cats/' + GAMA_or_SDSS + '_' + outfilename + mass_string + '.fits'
    comparison_objects.write(out_path, overwrite=True)
    
    return

def random_galaxy_comparison(GAMA_or_SDSS, infile, outfilename, random=False,
                             num_random=10) :
    
    in_dir_path = 'catalogs/candidate_cats/'
    
    if GAMA_or_SDSS == 'GAMA' :
        out_dir_path_name = 'catalogs/CARS_GAMA/' + outfilename
        path = GAMA_path
    if GAMA_or_SDSS == 'SDSS' :
        out_dir_path_name = 'catalogs/CARS_SDSS/' + outfilename
        path = SDSS_path
    
    inpath = in_dir_path + infile + '.fits'
    catalog = Table.read(inpath)
    length = len(catalog)
    
    if random == True : # calculate the env. params for a random number of galaxies
        rng = np.random.default_rng() # '0' refers to the seed
        randoms = rng.uniform(0, 1, num_random)
        out_end = 'random' + str(num_random)
        index_mask = (np.array(randoms*length)).astype(int)
    else :
        start, stop, step = 0, length, 1 # calculate the env. params for all in file
        out_end = str(start) + '-' + str(stop-1)
        index_mask = np.arange(start, stop, step)
    
    cat = catalog[index_mask]
    
    RAs = Angle(cat['RA_1'], u.deg)
    decs = Angle(cat['DEC_1'], u.deg)
    D_As = cat['D_A'].to('Mpc')
    D_Ls = cosmo.luminosity_distance(cat['Z'])
    julianIDs = cat['Name']
    colorMass = cat['logMass']
    
    if GAMA_or_SDSS == 'GAMA' :
        absmag_g = cat['absmag_g']
        absmag_i = cat['absmag_i']
        catMass = cat['logmstar']
        sfr = np.power(10, cat['SFR_0_1Gyr_percentile50'])
        condition = ((cat['SIG_OIIIR']/cat['SIG_OIIIR_ERR'] >= 3) &
                     (cat['OIIIR_NPEG'] == 0) & (cat['HB_FITFAIL'] == 0))
        OIII_width = (fwhm_conv*speed_of_light*cat['SIG_OIIIR']/5007)
        OIII_FWHM = np.where(condition, OIII_width, np.nan)
    if GAMA_or_SDSS == 'SDSS' :
        absmag_g = (cat['KCOR_MAG'][:,0] - 5*np.log10(D_Ls.to('pc')/u.pc) + 5)*u.mag
        absmag_i = (cat['KCOR_MAG'][:,2] - 5*np.log10(D_Ls.to('pc')/u.pc) + 5)*u.mag
        catMass = cat['MEDIAN_logmass']
        sfr = np.power(10, cat['MEDIAN_sfr'])
        condition = (cat['OIII_5007_FWHM']/cat['OIII_5007_FWHM_ERR'] >= 3)
        OIII_FWHM = np.where(condition, cat['OIII_5007_FWHM'], np.nan)
    
    columns = [julianIDs, RAs, decs, cat['Z'], D_As, D_Ls,
               Angle( ((2*u.Mpc)/D_As).value, u.radian ).to('arcmin'),
               absmag_g, absmag_i, colorMass,
               catMass, sfr, OIII_FWHM]
    
    column_headings = ('Name', 'RA', 'DEC', 'Z', 'D_A', 'D_L',
                       '2_Mpc_Radius', 'M_g', 'M_i', 'logMass',
                       'catlogMass', 'SFR', 'OIII_FWHM')
    
    comment = 'Flat \u039BCDM cosmology: H\u2080 = 70 km s\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3'
    
    base = Table(columns, names=column_headings)
    base.meta['comments'] = comment
    
    base_outpath = out_dir_path_name + '_base_' + out_end + '.fits'
    base.write(base_outpath, overwrite=True)
    
    tables = []
    for index in range(len(cat)) :
        sub_cat, table = env_params.gama_params(GAMA_or_SDSS, path, RAs[index],
                                                decs[index], cat['Z'][index],
                                                D_As[index], julianIDs[index],
                                                colorMass[index], self_in_search=True)
        tables.append(table)
        print(str(index) + ' finished.')
    
    compliments = vstack(tables)
    outpath = out_dir_path_name + '_envs_' + out_end + '.fits'
    compliments.write(outpath, overwrite=True)
    
    complete = hstack([base, compliments])
    complete_outpath = out_dir_path_name + '_complete_' + out_end + '.fits'
    complete.write(complete_outpath, overwrite=True)
    
    return

def random_galaxy_join_tables(ext) :
    
    # ext = 'base', 'envs', or 'complete'
    
    dir_path = 'catalogs/CARS_GAMA/'
    
    file_list = [
                'GAMA_SFG_1st_' + ext + '_random_500.fits',
                'GAMA_SFG_2nd_' + ext + '_random_500.fits',
                'GAMA_SFG_3rd_' + ext + '_random_500.fits',
                'GAMA_SFG_4th_' + ext + '_random_500.fits',
                'GAMA_SFG_5th_' + ext + '_random_500.fits'
                ]
    
    tables = []
    for file in file_list :
        tables.append(Table.read(dir_path + file))
    
    master_table = vstack(tables)
    master_table.write(dir_path + 'GAMA_SFG_' + ext + '_random_2500.fits')
    
    return

# subset the different objects types in GAMA
# random_galaxy_candidates('GAMA', 'comparison_BLAGN', 'BLAGN  ', mass_limit_default)
# random_galaxy_candidates('GAMA', 'comparison_Seyferts', 'Seyfert', mass_limit_default)
# random_galaxy_candidates('GAMA', 'comparison_LINERs', 'LINER  ', mass_limit_default)
# random_galaxy_candidates('GAMA', 'comparison_Comps', 'Comp   ', mass_limit_default)
# random_galaxy_candidates('GAMA', 'comparison_SFGs', 'SFG    ', mass_limit_default)
# random_galaxy_candidates('GAMA', 'comparison_Passives', 'Passive', mass_limit_default)
# random_galaxy_candidates('GAMA', 'comparison_not-ELGs', 'not_ELG', mass_limit_default)

# subset the different objects types in GAMA with logMass >= 10
# random_galaxy_candidates('GAMA', 'comparison_BLAGN', 'BLAGN  ', 10)
# random_galaxy_candidates('GAMA', 'comparison_Seyferts', 'Seyfert', 10)
# random_galaxy_candidates('GAMA', 'comparison_LINERs', 'LINER  ', 10)
# random_galaxy_candidates('GAMA', 'comparison_Comps', 'Comp   ', 10)
# random_galaxy_candidates('GAMA', 'comparison_SFGs', 'SFG    ', 10)
# random_galaxy_candidates('GAMA', 'comparison_Passives', 'Passive', 10)
# random_galaxy_candidates('GAMA', 'comparison_not-ELGs', 'not_ELG', 10)

# subset the different objects types in SDSS
# random_galaxy_candidates('SDSS', 'comparison_BLAGN', 'BLAGN  ', mass_limit_default)
# random_galaxy_candidates('SDSS', 'comparison_Seyferts', 'Seyfert', mass_limit_default)
# random_galaxy_candidates('SDSS', 'comparison_LINERs', 'LINER  ', mass_limit_default)
# random_galaxy_candidates('SDSS', 'comparison_Comps', 'Comp   ', mass_limit_default)
# random_galaxy_candidates('SDSS', 'comparison_SFGs', 'SFG    ', mass_limit_default)
# random_galaxy_candidates('SDSS', 'comparison_Passives', 'Passive', mass_limit_default)
# random_galaxy_candidates('SDSS', 'comparison_not-ELGs', 'not_ELG', mass_limit_default)

# subset the different objects types in SDSS with logMass >= 10
# random_galaxy_candidates('SDSS', 'comparison_BLAGN', 'BLAGN  ', 10)
# random_galaxy_candidates('SDSS', 'comparison_Seyferts', 'Seyfert', 10)
# random_galaxy_candidates('SDSS', 'comparison_LINERs', 'LINER  ', 10)
# random_galaxy_candidates('SDSS', 'comparison_Comps', 'Comp   ', 10)
# random_galaxy_candidates('SDSS', 'comparison_SFGs', 'SFG    ', 10)
# random_galaxy_candidates('SDSS', 'comparison_Passives', 'Passive', 10)
# random_galaxy_candidates('SDSS', 'comparison_not-ELGs', 'not_ELG', 10)

# determine env. params for the different objects types in GAMA with logMass >= 10
# random_galaxy_comparison('GAMA', 'GAMA_comparison_BLAGN_logMass10', 'GAMA_comparison_BLAGN_logMass10')
# random_galaxy_comparison('GAMA', 'GAMA_comparison_Seyferts_logMass10', 'GAMA_comparison_Seyferts_logMass10')
# random_galaxy_comparison('GAMA', 'GAMA_comparison_LINERs_logMass10', 'GAMA_comparison_LINERs_logMass10')
# random_galaxy_comparison('GAMA', 'GAMA_comparison_Comps_logMass10', 'GAMA_comparison_Comps_logMass10')
# random_galaxy_comparison('GAMA', 'GAMA_comparison_SFGs_logMass10', 'GAMA_comparison_SFGs_logMass10')
# random_galaxy_comparison('GAMA', 'GAMA_comparison_Passives_logMass10', 'GAMA_comparison_Passives_logMass10')
# random_galaxy_comparison('GAMA', 'GAMA_comparison_not-ELGs_logMass10', 'GAMA_comparison_not-ELGs_logMass10')

# determine env. params for the different objects types in SDSS with logMass >= 10
# random_galaxy_comparison('SDSS', 'SDSS_comparison_BLAGN_logMass10', 'SDSS_comparison_BLAGN_logMass10')    
# random_galaxy_comparison('SDSS', 'SDSS_comparison_Seyferts_logMass10', 'SDSS_comparison_Seyferts_logMass10')
# random_galaxy_comparison('SDSS', 'SDSS_comparison_Passives_logMass10', 'SDSS_comparison_Passives_logMass10')
# random_galaxy_comparison('SDSS', 'SDSS_comparison_SFGs_logMass10', 'SDSS_comparison_SFGs_logMass10')
# random_galaxy_comparison('SDSS', 'SDSS_comparison_Comps_logMass10', 'SDSS_comparison_Comps_logMass10')



# random_galaxy_comparison('SDSS', 'SDSS_comparison_not-ELGs_logMass10', 'SDSS_comparison_not-ELGs_logMass10')

# to complete:
# random_galaxy_comparison('SDSS', 'SDSS_comparison_LINERs_logMass10', 'SDSS_comparison_LINERs_logMass10')
