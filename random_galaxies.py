
# imports
import numpy as np

from astropy.coordinates import Angle
from astropy.cosmology import FlatLambdaCDM
from astropy.table import hstack, Table, vstack
import astropy.units as u

import environmental_parameters as env_params

# constants
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
mass_limit = 8.452021
SDSS_path = 'catalogs/joined_cats/SDSS_gal_info_gal_line_vCam.fits'

def random_galaxy_candidates(GAMA_or_SDSS, outfilename, targettype) :
    
    out_dir_path = 'catalogs/candidate_cats/'
    
    if GAMA_or_SDSS == 'GAMA' :
        dir_path = 'catalogs/joined_cats/'
        in_path = dir_path + 'GAMA_GaussFitSimple_StellarMasses_SpecClassGordon_vCam.fits'
    if GAMA_or_SDSS == 'SDSS' :
        dir_path = 'catalogs/joined_cats/'
        in_path = dir_path + 'SDSS_gal_info_gal_line_vCam.fits'
    
    base_cat = Table.read(in_path)
    
    # calculate D_L, M_i, colour masses
    if GAMA_or_SDSS == 'SDSS' :
        D_L = cosmo.luminosity_distance(base_cat['Z']).to(u.pc)/u.pc
        M_i = base_cat['KCOR_MAG'][:,2] - 5*np.log10(D_L) + 5 # absolute i mag.
        colour_mass = (1.15 + 0.7*(base_cat['KCOR_MAG'][:,0] -
                                   base_cat['KCOR_MAG'][:,2]) - 0.4*M_i)
    if GAMA_or_SDSS == 'GAMA' :
        colour_mass = 1.15 + 0.7*base_cat['gminusi'] - 0.4*base_cat['absmag_i']
    
    base_cat['colour_mass'] = colour_mass # add colour mass to catalog
    
    # GAMA mask: EmLineType==BLAGN,Comp,LINER,Passive,SFG,Seyfert,not_ELG
    # mass >= 8.452021, 0.01 <= z <= 0.06, z quality >= 3
    if GAMA_or_SDSS == 'GAMA' :
        mask = (
                (base_cat['EmLineType'] == targettype) &
                (0.01 <= base_cat['Z_1']) & (base_cat['Z_1'] <= 0.06) &
                (base_cat['NQ_1'] >= 3) &
                (base_cat['colour_mass'] >= mass_limit)
                )
    
    # SDSS mask: TARGETTYPE==GALAXY
    # mass >= 8.452021, 0.01 <= z <= 0.06, z quality >= 3
    if GAMA_or_SDSS == 'SDSS' :
        mask = (
                (base_cat['TARGETTYPE'] == targettype) &
                (0.01 <= base_cat['Z']) & (base_cat['Z'] <= 0.06) &
                ((base_cat['Z'] / base_cat['Z_ERR']) >= 3) &
                (base_cat['colour_mass'] >= mass_limit)
                )
    
    # create the resulting masked table of good objects and save it
    compliments = base_cat[mask]
    compliments['index'] = np.arange(len(compliments))
    out_path = out_dir_path + outfilename + '.fits'
    compliments.write(out_path, overwrite=False)
    
    return


def random_galaxy_comparison(GAMA_or_SDSS, infile, outfilename, random=True,
                             num_random=10) :
    
    in_dir_path = 'catalogs/candidate_cats/'
    
    if GAMA_or_SDSS == 'GAMA' :
        out_dir_path_name = 'catalogs/CARS_GAMA/' + outfilename
    if GAMA_or_SDSS == 'SDSS' :
        out_dir_path_name = 'catalogs/CARS_SDSS/' + outfilename
    
    inpath = in_dir_path + infile + '.fits'
    catalog = Table.read(inpath)
    length = len(catalog)
    
    if random == True : # calculate the env. params for a random number of galaxies
        rng = np.random.default_rng() # '0' refers to the seed
        randoms = rng.uniform(0, 1, num_random)
        out_end = 'random_' + str(num_random)
        index_mask = (np.array(randoms*length)).astype(int)
    else :
        start, stop, step = 0, length, 1 # calculate the env. params for all in file
        out_end = str(start) + '-' + str(stop-1)
        index_mask = np.arange(start, stop, step)
    
    cat = catalog[index_mask]
    
    if GAMA_or_SDSS == 'GAMA' :
        RAs = Angle(cat['RA_1'], u.deg)
        decs = Angle(cat['DEC_1'], u.deg)
        D_As = cosmo.angular_diameter_distance(cat['Z_1'])
        D_Ls = cosmo.luminosity_distance(cat['Z_1'])
    if GAMA_or_SDSS == 'SDSS' :
        RAs = Angle(cat['RA'], u.deg)
        decs = Angle(cat['DEC'], u.deg)
        D_As = cosmo.angular_diameter_distance(cat['Z'])
        D_Ls = cosmo.luminosity_distance(cat['Z'])
    
    julianIDs = []
    for i in range(len(cat)) :
        RA_str = RAs[i].to_string(unit=u.hour, sep='', precision=0, pad=True)
        dec_str = decs[i].to_string(unit=u.deg, sep='', precision=0,
                                    alwayssign=True, pad=True)
        ID = 'J' + RA_str + dec_str
        julianIDs.append(ID)
    
    if GAMA_or_SDSS == 'GAMA' :
        table_columns = [julianIDs, RAs, decs, cat['Z_1'], D_As, D_Ls,
                         Angle( ((2*u.Mpc)/D_As).value, u.radian ).to('arcmin'),
                         cat['gminusi'], cat['absmag_i'],
                         cat['colour_mass']/u.mag*u.solMass]
        column_hdrs = ('Name', 'RA', 'DEC', 'Z', 'D_A', 'D_L', '2_Mpc_Radius',
                       'g-i', 'M_i', 'log_mass')
    if GAMA_or_SDSS == 'SDSS' :
        table_columns = [julianIDs, RAs, decs, cat['Z'], D_As, D_Ls,
                         Angle( ((2*u.Mpc)/D_As).value, u.radian ).to('arcmin'),
                         cat['KCOR_MAG'][:,0]*u.mag, cat['KCOR_MAG'][:,2]*u.mag,
                         cat['KCOR_MAG'][:,2] - 5*np.log10(D_Ls.to('pc')/u.pc) + 5,
                         cat['colour_mass']*u.solMass]
        column_hdrs = ('Name', 'RA', 'DEC', 'Z', 'D_A', 'D_L', '2_Mpc_Radius',
                       'g', 'i', 'M_i', 'log_mass')
    
    base = Table(table_columns, names=column_hdrs)
    base.meta['comments'] = ['Flat \u039BCDM cosmology: H\u2080 = 70 km ' +
                              's\u207B\u00B9 Mpc\u207B\u00B9, \u03A9\u2098 = 0.3']
    base_outpath = out_dir_path_name + '_base_' + out_end + '.fits'
    base.write(base_outpath, overwrite=False)
    
    tables = []
    for index in range(len(cat)) :
        sub_cat, table = env_params.gama_params('SDSS', SDSS_path, RAs[index],
                                                decs[index], cat['Z'][index],
                                                D_As[index], julianIDs[index],
                                                np.power(10, cat['colour_mass'][index])*u.solMass,
                                                self_in_search=True)
        tables.append(table)
        print(str(index) + ' finished.')
    
    compliments = vstack(tables)
    outpath = out_dir_path_name + '_envs_' + out_end + '.fits'
    compliments.write(outpath, overwrite=False)
    
    complete = hstack([base, compliments])
    complete_outpath = out_dir_path_name + '_complete_' + out_end + '.fits'
    complete.write(complete_outpath, overwrite=False)
    
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

# for other 116901 galaxies in SDSS
# random_galaxy_candidates('SDSS', 'gal_info_dr7_v5_2_vCam_SDSS-galaxies_test',
                            # 'GALAXY             ') # for SDSS general galaxies
# random_galaxy_comparison('gal_info_dr7_v5_2_vCam_SDSS-galaxies',
                          # 'SDSS-galaxies', random=True, num_random=10)

# for other 53 BLAGN in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_BLAGN',
                           # 'BLAGN  ')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_BLAGN',
                         # 'GAMA_BLAGN', random=False)

# for other 160 Comp in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_Comp',
                          # 'Comp   ')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_Comp',
                         # 'GAMA_Comp', random=False)

# for other 53 LINER in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_LINER',
                          # 'LINER  ')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_LINER',
                         # 'GAMA_LINER', random=False)

# for other 195 Passive in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_Passive',
                          # 'Passive')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_Passive',
                          # 'GAMA_Passive', random=False)

# for other 2727 SFG in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_SFG',
                          # 'SFG    ')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_SFG',
                            # 'GAMA_SFG', random=False)

# for other 40 Seyfert in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_Seyfert',
                          # 'Seyfert')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_Seyfert',
                          # 'GAMA_Seyfert', random=False)

# for other 2195 not_ELG in GAMA
# random_galaxy_candidates('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_not-ELG',
                          # 'not_ELG')
# random_galaxy_comparison('GAMA', 'GAMA_GFS_SM_SCG_vCam_GAMA_not-ELG',
                          # 'GAMA_not-ELG', random=False)
