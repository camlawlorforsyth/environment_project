
# imports
import numpy as np

from astropy.table import Table

import functions as funcs
import plots as plt

def error_calc(cc, kk, nn) :
    # as adapted from https://ui.adsabs.harvard.edu/abs/2011PASA...28..128C/abstract
    # cc is the confidence level
    # kk is the success count
    # nn is the sample size
    
    import scipy.stats.distributions as dist
    
    p_lower = dist.beta.ppf( (1 - cc)/2, kk + 1, nn - kk + 1)
    p_upper = dist.beta.ppf( 1 - (1 - cc)/2, kk + 1, nn - kk + 1)
    
    sigmas = [p_upper - p_lower, abs(p_lower - p_upper)]
    
    return max(sigmas)

def population_fracs() :
    
    from scipy.stats import norm
    
    GAMA = Table.read('catalogs/CARS_GAMA/GAMA_allObjects_logMass10_complete_0-1364.fits')
    SDSS = Table.read('catalogs/CARS_SDSS/SDSS_allObjects_logMass10_complete_0-52113.fits')
    
    logCD_gama = np.log10(GAMA['Closest_Distance'])
    logCD_gama = logCD_gama[logCD_gama > 0]
    logCD_gama_len = len(logCD_gama)
    logCD_gama_weight = np.ones(logCD_gama_len)/logCD_gama_len
    
    mu_CD, std_CD = norm.fit(logCD_gama)
    # print(mu_CD, std_CD)
    CD_xs = np.linspace(0, 5, 10000)
    CD_ys = funcs.gaussian_fit(CD_xs, mu_CD, std_CD)
    
    logSD_gama = np.log10(GAMA['SurfaceDensity'])
    logSD_gama = logSD_gama[logSD_gama < 10]
    logSD_gama_len = len(logSD_gama)
    logSD_gama_weight = np.ones(logSD_gama_len)/logSD_gama_len
    
    mu_SD, std_SD = norm.fit(logSD_gama)
    # print(mu_SD, std_SD)
    SD_xs = np.linspace(-6, 2, 10000)
    SD_ys = funcs.gaussian_fit(SD_xs, mu_SD, std_SD)
    
    logCD_sdss = np.log10(SDSS['Closest_Distance'])
    logCD_sdss = logCD_sdss[logCD_sdss > 0]    
    logCD_sdss_len = len(logCD_sdss)
    logCD_sdss_weight = np.ones(logCD_sdss_len)/logCD_sdss_len
    
    logSD_sdss = np.log10(SDSS['SurfaceDensity'])
    logSD_sdss = logSD_sdss[logSD_sdss < 10]
    logSD_sdss_len = len(logSD_sdss)
    logSD_sdss_weight = np.ones(logSD_sdss_len)/logSD_sdss_len
    
    # plt.histo_with_fxn(logCD_gama, r'$\log($Distance to Closest Companion/kpc)',
    #                     logCD_gama_weight, CD_xs, 0.2*CD_ys, 'Fit to GAMA:\n$\mu=2.93$, $\sigma=0.36$',
    #                     loc=2, num_bins=17, outfile='CD_histo_with_fit_GAMA.pdf',
    #                     xmin=0, xmax=4.5, ymin=0, ymax=0.25)
    
    # plt.histo_with_fxn(logCD_sdss, r'$\log($Distance to Closest Companion/kpc)',
    #                     logCD_sdss_weight, CD_xs, 0.2*CD_ys, 'Fit to GAMA:\n$\mu=2.93$, $\sigma=0.36$',
    #                     loc=2, num_bins=25, outfile='CD_histo_with_fit_SDSS.pdf',
    #                     xmin=0, xmax=4.5, ymin=0, ymax=0.25)
    
    # plt.histo_with_fxn(logSD_gama, r'$\log($Surface Density/Mpc$^{-2}$)',
    #                     logSD_gama_weight, SD_xs, 0.14*SD_ys, 'Fit to GAMA:\n$\mu=-2.22$, $\sigma=0.86$',
    #                     loc=1, num_bins=17, outfile='SD_histo_with_fit_GAMA.pdf',
    #                     xmin=-6, xmax=2, ymin=0, ymax=0.15)
    
    # plt.histo_with_fxn(logSD_sdss, r'$\log($Surface Density/Mpc$^{-2}$)',
    #                     logSD_sdss_weight, SD_xs, 0.14*SD_ys, 'Fit to GAMA:\n$\mu=-2.22$, $\sigma=0.86$',
    #                     loc=1, num_bins=45, outfile='SD_histo_with_fit_SDSS.pdf',
    #                     xmin=-6, xmax=2, ymin=0, ymax=0.15)
    
    # based on distance to closest companion
    # print('\nClosest Distance')
    CD_table = Table()
    CD_table['Population'] = ['AGN', 'SFG', 'Passive']
    
    # low-density environments
    CD_gama_LD_mask = np.log10(GAMA['Closest_Distance']) > mu_CD + std_CD
    CD_gama_LD  = GAMA[CD_gama_LD_mask]
    CD_gama_LD_len = len(CD_gama_LD)
    CD_gama_LD_agn = (np.sum(CD_gama_LD['Type'] == 'BLAGN  ') +
                      np.sum(CD_gama_LD['Type'] == 'Seyfert') +
                      np.sum(CD_gama_LD['Type'] == 'LINER  ') )
    CD_gama_LD_sfg = np.sum(CD_gama_LD['Type'] == 'SFG    ')
    CD_gama_LD_pass = (np.sum(CD_gama_LD['Type'] == 'Passive') +
                       np.sum(CD_gama_LD['Type'] == 'not_ELG') )
    
    CD_sdss_LD_mask = np.log10(SDSS['Closest_Distance']) > mu_CD + std_CD
    CD_sdss_LD  = SDSS[CD_sdss_LD_mask]
    CD_sdss_LD_len = len(CD_sdss_LD)
    CD_sdss_LD_agn = (np.sum(CD_sdss_LD['Type'] == 'BLAGN  ') +
                      np.sum(CD_sdss_LD['Type'] == 'Seyfert') +
                      np.sum(CD_sdss_LD['Type'] == 'LINER  ') )
    CD_sdss_LD_sfg = np.sum(CD_sdss_LD['Type'] == 'SFG    ')
    CD_sdss_LD_pass = (np.sum(CD_sdss_LD['Type'] == 'Passive') +
                       np.sum(CD_sdss_LD['Type'] == 'not_ELG') )
    
    # low-density to medium-density environments
    CD_gama_LDMD_mask = ( (np.log10(GAMA['Closest_Distance']) > mu_CD) &
                          (np.log10(GAMA['Closest_Distance']) < mu_CD + std_CD) )
    CD_gama_LDMD  = GAMA[CD_gama_LDMD_mask]
    CD_gama_LDMD_len = len(CD_gama_LDMD)
    CD_gama_LDMD_agn = (np.sum(CD_gama_LDMD['Type'] == 'BLAGN  ') +
                        np.sum(CD_gama_LDMD['Type'] == 'Seyfert') +
                        np.sum(CD_gama_LDMD['Type'] == 'LINER  ') )
    CD_gama_LDMD_sfg = np.sum(CD_gama_LDMD['Type'] == 'SFG    ')
    CD_gama_LDMD_pass = (np.sum(CD_gama_LDMD['Type'] == 'Passive') +
                         np.sum(CD_gama_LDMD['Type'] == 'not_ELG') )
    
    CD_sdss_LDMD_mask = ( (np.log10(SDSS['Closest_Distance']) > mu_CD) &
                          (np.log10(SDSS['Closest_Distance']) < mu_CD + std_CD) )
    CD_sdss_LDMD  = SDSS[CD_sdss_LDMD_mask]
    CD_sdss_LDMD_len = len(CD_sdss_LDMD)
    CD_sdss_LDMD_agn = (np.sum(CD_sdss_LDMD['Type'] == 'BLAGN  ') +
                        np.sum(CD_sdss_LDMD['Type'] == 'Seyfert') +
                        np.sum(CD_sdss_LDMD['Type'] == 'LINER  ') )
    CD_sdss_LDMD_sfg = np.sum(CD_sdss_LDMD['Type'] == 'SFG    ')
    CD_sdss_LDMD_pass = (np.sum(CD_sdss_LDMD['Type'] == 'Passive') +
                         np.sum(CD_sdss_LDMD['Type'] == 'not_ELG') )
    
    # medium-density to high-density environments
    CD_gama_MDHD_mask = ( (np.log10(GAMA['Closest_Distance']) < mu_CD) &
                          (np.log10(GAMA['Closest_Distance']) > mu_CD - std_CD) )
    CD_gama_MDHD  = GAMA[CD_gama_MDHD_mask]
    CD_gama_MDHD_len = len(CD_gama_MDHD)
    CD_gama_MDHD_agn = (np.sum(CD_gama_MDHD['Type'] == 'BLAGN  ') +
                        np.sum(CD_gama_MDHD['Type'] == 'Seyfert') +
                        np.sum(CD_gama_MDHD['Type'] == 'LINER  ') )
    CD_gama_MDHD_sfg = np.sum(CD_gama_MDHD['Type'] == 'SFG    ')
    CD_gama_MDHD_pass = (np.sum(CD_gama_MDHD['Type'] == 'Passive') +
                         np.sum(CD_gama_MDHD['Type'] == 'not_ELG') )
    
    CD_sdss_MDHD_mask = ( (np.log10(SDSS['Closest_Distance']) < mu_CD) &
                          (np.log10(SDSS['Closest_Distance']) > mu_CD - std_CD) )
    CD_sdss_MDHD  = SDSS[CD_sdss_MDHD_mask]
    CD_sdss_MDHD_len = len(CD_sdss_MDHD)
    CD_sdss_MDHD_agn = (np.sum(CD_sdss_MDHD['Type'] == 'BLAGN  ') +
                        np.sum(CD_sdss_MDHD['Type'] == 'Seyfert') +
                        np.sum(CD_sdss_MDHD['Type'] == 'LINER  ') )
    CD_sdss_MDHD_sfg = np.sum(CD_sdss_MDHD['Type'] == 'SFG    ')
    CD_sdss_MDHD_pass = (np.sum(CD_sdss_MDHD['Type'] == 'Passive') +
                         np.sum(CD_sdss_MDHD['Type'] == 'not_ELG') )
    
    # high-density environments
    CD_gama_HD_mask = np.log10(GAMA['Closest_Distance']) < mu_CD - std_CD
    CD_gama_HD  = GAMA[CD_gama_HD_mask]
    CD_gama_HD_len = len(CD_gama_HD)
    CD_gama_HD_agn = (np.sum(CD_gama_HD['Type'] == 'BLAGN  ') +
                      np.sum(CD_gama_HD['Type'] == 'Seyfert') +
                      np.sum(CD_gama_HD['Type'] == 'LINER  ') )
    CD_gama_HD_sfg = np.sum(CD_gama_HD['Type'] == 'SFG    ')
    CD_gama_HD_pass = (np.sum(CD_gama_HD['Type'] == 'Passive') +
                       np.sum(CD_gama_HD['Type'] == 'not_ELG') )
    
    CD_sdss_HD_mask = np.log10(SDSS['Closest_Distance']) < mu_CD - std_CD
    CD_sdss_HD  = SDSS[CD_sdss_HD_mask]
    CD_sdss_HD_len = len(CD_sdss_HD)
    CD_sdss_HD_agn = (np.sum(CD_sdss_HD['Type'] == 'BLAGN  ') +
                      np.sum(CD_sdss_HD['Type'] == 'Seyfert') +
                      np.sum(CD_sdss_HD['Type'] == 'LINER  ') )
    CD_sdss_HD_sfg = np.sum(CD_sdss_HD['Type'] == 'SFG    ')
    CD_sdss_HD_pass = (np.sum(CD_sdss_HD['Type'] == 'Passive') +
                       np.sum(CD_sdss_HD['Type'] == 'not_ELG') )
    
    # SDSS calculations and errors
    CD_sdss_LD = np.array([CD_sdss_LD_agn, CD_sdss_LD_sfg, CD_sdss_LD_pass])/CD_sdss_LD_len
    CD_sdss_LD_agn_err = error_calc(0.683, CD_sdss_LD_agn, CD_sdss_LD_len)
    CD_sdss_LD_sfg_err = error_calc(0.683, CD_sdss_LD_sfg, CD_sdss_LD_len)
    CD_sdss_LD_pass_err = error_calc(0.683, CD_sdss_LD_pass, CD_sdss_LD_len)
    
    CD_sdss_LDMD = np.array([CD_sdss_LDMD_agn, CD_sdss_LDMD_sfg, CD_sdss_LDMD_pass])/CD_sdss_LDMD_len
    CD_sdss_LDMD_agn_err = error_calc(0.683, CD_sdss_LDMD_agn, CD_sdss_LDMD_len)
    CD_sdss_LDMD_sfg_err = error_calc(0.683, CD_sdss_LDMD_sfg, CD_sdss_LDMD_len)
    CD_sdss_LDMD_pass_err = error_calc(0.683, CD_sdss_LDMD_pass, CD_sdss_LDMD_len)
    
    CD_sdss_MDHD = np.array([CD_sdss_MDHD_agn, CD_sdss_MDHD_sfg, CD_sdss_MDHD_pass])/CD_sdss_MDHD_len
    CD_sdss_MDHD_agn_err = error_calc(0.683, CD_sdss_MDHD_agn, CD_sdss_MDHD_len)
    CD_sdss_MDHD_sfg_err = error_calc(0.683, CD_sdss_MDHD_sfg, CD_sdss_MDHD_len)
    CD_sdss_MDHD_pass_err = error_calc(0.683, CD_sdss_MDHD_pass, CD_sdss_MDHD_len)
    
    CD_sdss_HD = np.array([CD_sdss_HD_agn, CD_sdss_HD_sfg, CD_sdss_HD_pass])/CD_sdss_HD_len
    CD_sdss_HD_agn_err = error_calc(0.683, CD_sdss_HD_agn, CD_sdss_HD_len)
    CD_sdss_HD_sfg_err = error_calc(0.683, CD_sdss_HD_sfg, CD_sdss_HD_len)
    CD_sdss_HD_pass_err = error_calc(0.683, CD_sdss_HD_pass, CD_sdss_HD_len)
    
    # GAMA calculations and errors
    CD_gama_LD = np.array([CD_gama_LD_agn, CD_gama_LD_sfg, CD_gama_LD_pass])/CD_gama_LD_len
    CD_gama_LD_agn_err = error_calc(0.683, CD_gama_LD_agn, CD_gama_LD_len)
    CD_gama_LD_sfg_err = error_calc(0.683, CD_gama_LD_sfg, CD_gama_LD_len)
    CD_gama_LD_pass_err = error_calc(0.683, CD_gama_LD_pass, CD_gama_LD_len)
    
    CD_gama_LDMD = np.array([CD_gama_LDMD_agn, CD_gama_LDMD_sfg, CD_gama_LDMD_pass])/CD_gama_LDMD_len
    CD_gama_LDMD_agn_err = error_calc(0.683, CD_gama_LDMD_agn, CD_gama_LDMD_len)
    CD_gama_LDMD_sfg_err = error_calc(0.683, CD_gama_LDMD_sfg, CD_gama_LDMD_len)
    CD_gama_LDMD_pass_err = error_calc(0.683, CD_gama_LDMD_pass, CD_gama_LDMD_len)
    
    CD_gama_MDHD = np.array([CD_gama_MDHD_agn, CD_gama_MDHD_sfg, CD_gama_MDHD_pass])/CD_gama_MDHD_len
    CD_gama_MDHD_agn_err = error_calc(0.683, CD_gama_MDHD_agn, CD_gama_MDHD_len)
    CD_gama_MDHD_sfg_err = error_calc(0.683, CD_gama_MDHD_sfg, CD_gama_MDHD_len)
    CD_gama_MDHD_pass_err = error_calc(0.683, CD_gama_MDHD_pass, CD_gama_MDHD_len)
    
    CD_gama_HD = np.array([CD_gama_HD_agn, CD_gama_HD_sfg, CD_gama_HD_pass])/CD_gama_HD_len
    CD_gama_HD_agn_err = error_calc(0.683, CD_gama_HD_agn, CD_gama_HD_len)
    CD_gama_HD_sfg_err = error_calc(0.683, CD_gama_HD_sfg, CD_gama_HD_len)
    CD_gama_HD_pass_err = error_calc(0.683, CD_gama_HD_pass, CD_gama_HD_len)
    
    CD_table['SDSS LD'] = CD_sdss_LD
    # CD_table['SDSS LD Err'] = [CD_sdss_LD_agn_err, CD_sdss_LD_sfg_err, CD_sdss_LD_pass_err]
    CD_table['SDSS LDMD'] = CD_sdss_LDMD
    # CD_table['SDSS LDMD Err'] = [CD_sdss_LDMD_agn_err, CD_sdss_LDMD_sfg_err, CD_sdss_LDMD_pass_err]
    CD_table['SDSS MDHD'] = CD_sdss_MDHD
    # CD_table['SDSS MDHD Err'] = [CD_sdss_MDHD_agn_err, CD_sdss_MDHD_sfg_err, CD_sdss_MDHD_pass_err]
    CD_table['SDSS HD'] = CD_sdss_HD
    # CD_table['SDSS HD Err'] = [CD_sdss_HD_agn_err, CD_sdss_HD_sfg_err, CD_sdss_HD_pass_err]
    CD_table['GAMA LD'] = CD_gama_LD
    # CD_table['GAMA LD Err'] = [CD_gama_LD_agn_err, CD_gama_LD_sfg_err, CD_gama_LD_pass_err]
    CD_table['GAMA LDMD'] = CD_gama_LDMD
    # CD_table['GAMA LDMD Err'] = [CD_gama_LDMD_agn_err, CD_gama_LDMD_sfg_err, CD_gama_LDMD_pass_err]
    CD_table['GAMA MDHD'] = CD_gama_MDHD
    # CD_table['GAMA MDHD Err'] = [CD_gama_MDHD_agn_err, CD_gama_MDHD_sfg_err, CD_gama_MDHD_pass_err]
    CD_table['GAMA HD'] = CD_gama_HD
    # CD_table['GAMA HD Err'] = [CD_gama_HD_agn_err, CD_gama_HD_sfg_err, CD_gama_HD_pass_err]
    
    CD_table['SDSS LD'].info.format = '.2f'
    # CD_table['SDSS LD Err'].info.format = '.2f'
    CD_table['SDSS LDMD'].info.format = '.2f'
    # CD_table['SDSS LDMD Err'].info.format = '.2f'
    CD_table['SDSS MDHD'].info.format = '.2f'
    # CD_table['SDSS MDHD Err'].info.format = '.2f'
    CD_table['SDSS HD'].info.format = '.2f'
    # CD_table['SDSS HD Err'].info.format = '.2f'
    CD_table['GAMA LD'].info.format = '.2f'
    # CD_table['GAMA LD Err'].info.format = '.2f'
    CD_table['GAMA LDMD'].info.format = '.2f'
    # CD_table['GAMA LDMD Err'].info.format = '.2f'
    CD_table['GAMA MDHD'].info.format = '.2f'
    # CD_table['GAMA MDHD Err'].info.format = '.2f'
    CD_table['GAMA HD'].info.format = '.2f'
    # CD_table['GAMA HD Err'].info.format = '.2f'
    # CD_table.pprint(max_width=-1)
    print('')
    
    gama_agn = [CD_gama_LD[0], CD_gama_LDMD[0], CD_gama_MDHD[0], CD_gama_HD[0]]
    gama_agn_err = [CD_gama_LD_agn_err, CD_gama_LDMD_agn_err, CD_gama_MDHD_agn_err, CD_gama_HD_agn_err]
    gama_sfg = [CD_gama_LD[1], CD_gama_LDMD[1], CD_gama_MDHD[1], CD_gama_HD[1]]
    gama_sfg_err = [CD_gama_LD_sfg_err, CD_gama_LDMD_sfg_err, CD_gama_MDHD_sfg_err, CD_gama_HD_sfg_err]
    gama_pass = [CD_gama_LD[2], CD_gama_LDMD[2], CD_gama_MDHD[2], CD_gama_HD[2]]
    gama_pass_err = [CD_gama_LD_pass_err, CD_gama_LDMD_pass_err, CD_gama_MDHD_pass_err, CD_gama_HD_pass_err]
    plt.population_plot(gama_agn, gama_agn_err, gama_sfg, gama_sfg_err,
                        gama_pass, gama_pass_err,
                        ['GAMA AGN', 'GAMA SFG', 'GAMA Passive'], ymin=0, ymax=0.8)

    sdss_agn = [CD_sdss_LD[0], CD_sdss_LDMD[0], CD_sdss_MDHD[0], CD_sdss_HD[0]]
    sdss_agn_err = [CD_sdss_LD_agn_err, CD_sdss_LDMD_agn_err, CD_sdss_MDHD_agn_err, CD_sdss_HD_agn_err]
    sdss_sfg = [CD_sdss_LD[1], CD_sdss_LDMD[1], CD_sdss_MDHD[1], CD_sdss_HD[1]]
    sdss_sfg_err = [CD_sdss_LD_sfg_err, CD_sdss_LDMD_sfg_err, CD_sdss_MDHD_sfg_err, CD_sdss_HD_sfg_err]
    sdss_pass = [CD_sdss_LD[2], CD_sdss_LDMD[2], CD_sdss_MDHD[2], CD_sdss_HD[2]]
    sdss_pass_err = [CD_sdss_LD_pass_err, CD_sdss_LDMD_pass_err, CD_sdss_MDHD_pass_err, CD_sdss_HD_pass_err]
    plt.population_plot(sdss_agn, sdss_agn_err, sdss_sfg, sdss_sfg_err,
                        sdss_pass, sdss_pass_err,
                        ['SDSS AGN', 'SDSS SFG', 'SDSS Passive'], ymin=0, ymax=0.52)
    
    # based on surface density
    # print('Surface Density')
    SD_table = Table()
    SD_table['Population'] = ['AGN', 'SFG', 'Passive']
    
    # low-density environments
    SD_gama_LD_mask = np.log10(GAMA['SurfaceDensity']) < mu_SD - std_SD
    SD_gama_LD  = GAMA[SD_gama_LD_mask]
    SD_gama_LD_len = len(SD_gama_LD)
    SD_gama_LD_agn = (np.sum(SD_gama_LD['Type'] == 'BLAGN  ') +
                      np.sum(SD_gama_LD['Type'] == 'Seyfert') +
                      np.sum(SD_gama_LD['Type'] == 'LINER  ') )
    SD_gama_LD_sfg = np.sum(SD_gama_LD['Type'] == 'SFG    ')
    SD_gama_LD_pass = (np.sum(SD_gama_LD['Type'] == 'Passive') +
                       np.sum(SD_gama_LD['Type'] == 'not_ELG') )
    
    SD_sdss_LD_mask = np.log10(SDSS['SurfaceDensity']) < mu_SD - std_SD
    SD_sdss_LD  = SDSS[SD_sdss_LD_mask]
    SD_sdss_LD_len = len(SD_sdss_LD)
    SD_sdss_LD_agn = (np.sum(SD_sdss_LD['Type'] == 'BLAGN  ') +
                      np.sum(SD_sdss_LD['Type'] == 'Seyfert') +
                      np.sum(SD_sdss_LD['Type'] == 'LINER  ') )
    SD_sdss_LD_sfg = np.sum(SD_sdss_LD['Type'] == 'SFG    ')
    SD_sdss_LD_pass = (np.sum(SD_sdss_LD['Type'] == 'Passive') +
                       np.sum(SD_sdss_LD['Type'] == 'not_ELG') )
    
    # low-density to medium-density environments
    SD_gama_LDMD_mask = ( (np.log10(GAMA['SurfaceDensity']) < mu_SD) &
                          (np.log10(GAMA['SurfaceDensity']) > mu_SD - std_SD) )
    SD_gama_LDMD  = GAMA[SD_gama_LDMD_mask]
    SD_gama_LDMD_len = len(SD_gama_LDMD)
    SD_gama_LDMD_agn = (np.sum(SD_gama_LDMD['Type'] == 'BLAGN  ') +
                        np.sum(SD_gama_LDMD['Type'] == 'Seyfert') +
                        np.sum(SD_gama_LDMD['Type'] == 'LINER  ') )
    SD_gama_LDMD_sfg = np.sum(SD_gama_LDMD['Type'] == 'SFG    ')
    SD_gama_LDMD_pass = (np.sum(SD_gama_LDMD['Type'] == 'Passive') +
                         np.sum(SD_gama_LDMD['Type'] == 'not_ELG') )
    
    SD_sdss_LDMD_mask = ( (np.log10(SDSS['SurfaceDensity']) < mu_SD) &
                          (np.log10(SDSS['SurfaceDensity']) > mu_SD - std_SD) )
    SD_sdss_LDMD  = SDSS[SD_sdss_LDMD_mask]
    SD_sdss_LDMD_len = len(SD_sdss_LDMD)
    SD_sdss_LDMD_agn = (np.sum(SD_sdss_LDMD['Type'] == 'BLAGN  ') +
                        np.sum(SD_sdss_LDMD['Type'] == 'Seyfert') +
                        np.sum(SD_sdss_LDMD['Type'] == 'LINER  ') )
    SD_sdss_LDMD_sfg = np.sum(SD_sdss_LDMD['Type'] == 'SFG    ')
    SD_sdss_LDMD_pass = (np.sum(SD_sdss_LDMD['Type'] == 'Passive') +
                         np.sum(SD_sdss_LDMD['Type'] == 'not_ELG') )
    
    # medium-density to high-density environments
    SD_gama_MDHD_mask = ( (np.log10(GAMA['SurfaceDensity']) > mu_SD) &
                          (np.log10(GAMA['SurfaceDensity']) < mu_SD + std_SD) )
    SD_gama_MDHD  = GAMA[SD_gama_MDHD_mask]
    SD_gama_MDHD_len = len(SD_gama_MDHD)
    SD_gama_MDHD_agn = (np.sum(SD_gama_MDHD['Type'] == 'BLAGN  ') +
                        np.sum(SD_gama_MDHD['Type'] == 'Seyfert') +
                        np.sum(SD_gama_MDHD['Type'] == 'LINER  ') )
    SD_gama_MDHD_sfg = np.sum(SD_gama_MDHD['Type'] == 'SFG    ')
    SD_gama_MDHD_pass = (np.sum(SD_gama_MDHD['Type'] == 'Passive') +
                         np.sum(SD_gama_MDHD['Type'] == 'not_ELG') )
    
    SD_sdss_MDHD_mask = ( (np.log10(SDSS['SurfaceDensity']) > mu_SD) &
                          (np.log10(SDSS['SurfaceDensity']) < mu_SD + std_SD) )
    SD_sdss_MDHD  = SDSS[SD_sdss_MDHD_mask]
    SD_sdss_MDHD_len = len(SD_sdss_MDHD)
    SD_sdss_MDHD_agn = (np.sum(SD_sdss_MDHD['Type'] == 'BLAGN  ') +
                        np.sum(SD_sdss_MDHD['Type'] == 'Seyfert') +
                        np.sum(SD_sdss_MDHD['Type'] == 'LINER  ') )
    SD_sdss_MDHD_sfg = np.sum(SD_sdss_MDHD['Type'] == 'SFG    ')
    SD_sdss_MDHD_pass = (np.sum(SD_sdss_MDHD['Type'] == 'Passive') +
                         np.sum(SD_sdss_MDHD['Type'] == 'not_ELG') )
    
    # high-density environments
    SD_gama_HD_mask = np.log10(GAMA['SurfaceDensity']) > mu_SD + std_SD
    SD_gama_HD  = GAMA[SD_gama_HD_mask]
    SD_gama_HD_len = len(SD_gama_HD)
    SD_gama_HD_agn = (np.sum(SD_gama_HD['Type'] == 'BLAGN  ') +
                      np.sum(SD_gama_HD['Type'] == 'Seyfert') +
                      np.sum(SD_gama_HD['Type'] == 'LINER  ') )
    SD_gama_HD_sfg = np.sum(SD_gama_HD['Type'] == 'SFG    ')
    SD_gama_HD_pass = (np.sum(SD_gama_HD['Type'] == 'Passive') +
                       np.sum(SD_gama_HD['Type'] == 'not_ELG') )
    
    SD_sdss_HD_mask = np.log10(SDSS['SurfaceDensity']) > mu_SD + std_SD
    SD_sdss_HD  = SDSS[SD_sdss_HD_mask]
    SD_sdss_HD_len = len(SD_sdss_HD)
    SD_sdss_HD_agn = (np.sum(SD_sdss_HD['Type'] == 'BLAGN  ') +
                      np.sum(SD_sdss_HD['Type'] == 'Seyfert') +
                      np.sum(SD_sdss_HD['Type'] == 'LINER  ') )
    SD_sdss_HD_sfg = np.sum(SD_sdss_HD['Type'] == 'SFG    ')
    SD_sdss_HD_pass = (np.sum(SD_sdss_HD['Type'] == 'Passive') +
                       np.sum(SD_sdss_HD['Type'] == 'not_ELG') )
    
    # SDSS calculations and errors
    SD_sdss_LD = np.array([SD_sdss_LD_agn, SD_sdss_LD_sfg, SD_sdss_LD_pass])/SD_sdss_LD_len
    SD_sdss_LD_agn_err = error_calc(0.683, SD_sdss_LD_agn, SD_sdss_LD_len)
    SD_sdss_LD_sfg_err = error_calc(0.683, SD_sdss_LD_sfg, SD_sdss_LD_len)
    SD_sdss_LD_pass_err = error_calc(0.683, SD_sdss_LD_pass, SD_sdss_LD_len)
    
    SD_sdss_LDMD = np.array([SD_sdss_LDMD_agn, SD_sdss_LDMD_sfg, SD_sdss_LDMD_pass])/SD_sdss_LDMD_len
    SD_sdss_LDMD_agn_err = error_calc(0.683, SD_sdss_LDMD_agn, SD_sdss_LDMD_len)
    SD_sdss_LDMD_sfg_err = error_calc(0.683, SD_sdss_LDMD_sfg, SD_sdss_LDMD_len)
    SD_sdss_LDMD_pass_err = error_calc(0.683, SD_sdss_LDMD_pass, SD_sdss_LDMD_len)
    
    SD_sdss_MDHD = np.array([SD_sdss_MDHD_agn, SD_sdss_MDHD_sfg, SD_sdss_MDHD_pass])/SD_sdss_MDHD_len
    SD_sdss_MDHD_agn_err = error_calc(0.683, SD_sdss_MDHD_agn, SD_sdss_MDHD_len)
    SD_sdss_MDHD_sfg_err = error_calc(0.683, SD_sdss_MDHD_sfg, SD_sdss_MDHD_len)
    SD_sdss_MDHD_pass_err = error_calc(0.683, SD_sdss_MDHD_pass, SD_sdss_MDHD_len)
    
    SD_sdss_HD = np.array([SD_sdss_HD_agn, SD_sdss_HD_sfg, SD_sdss_HD_pass])/SD_sdss_HD_len
    SD_sdss_HD_agn_err = error_calc(0.683, SD_sdss_HD_agn, SD_sdss_HD_len)
    SD_sdss_HD_sfg_err = error_calc(0.683, SD_sdss_HD_sfg, SD_sdss_HD_len)
    SD_sdss_HD_pass_err = error_calc(0.683, SD_sdss_HD_pass, SD_sdss_HD_len)
    
    # GAMA calculations and errors
    SD_gama_LD = np.array([SD_gama_LD_agn, SD_gama_LD_sfg, SD_gama_LD_pass])/SD_gama_LD_len
    SD_gama_LD_agn_err = error_calc(0.683, SD_gama_LD_agn, SD_gama_LD_len)
    SD_gama_LD_sfg_err = error_calc(0.683, SD_gama_LD_sfg, SD_gama_LD_len)
    SD_gama_LD_pass_err = error_calc(0.683, SD_gama_LD_pass, SD_gama_LD_len)
    
    SD_gama_LDMD = np.array([SD_gama_LDMD_agn, SD_gama_LDMD_sfg, SD_gama_LDMD_pass])/SD_gama_LDMD_len
    SD_gama_LDMD_agn_err = error_calc(0.683, SD_gama_LDMD_agn, SD_gama_LDMD_len)
    SD_gama_LDMD_sfg_err = error_calc(0.683, SD_gama_LDMD_sfg, SD_gama_LDMD_len)
    SD_gama_LDMD_pass_err = error_calc(0.683, SD_gama_LDMD_pass, SD_gama_LDMD_len)
    
    SD_gama_MDHD = np.array([SD_gama_MDHD_agn, SD_gama_MDHD_sfg, SD_gama_MDHD_pass])/SD_gama_MDHD_len
    SD_gama_MDHD_agn_err = error_calc(0.683, SD_gama_MDHD_agn, SD_gama_MDHD_len)
    SD_gama_MDHD_sfg_err = error_calc(0.683, SD_gama_MDHD_sfg, SD_gama_MDHD_len)
    SD_gama_MDHD_pass_err = error_calc(0.683, SD_gama_MDHD_pass, SD_gama_MDHD_len)
    
    SD_gama_HD = np.array([SD_gama_HD_agn, SD_gama_HD_sfg, SD_gama_HD_pass])/SD_gama_HD_len
    SD_gama_HD_agn_err = error_calc(0.683, SD_gama_HD_agn, SD_gama_HD_len)
    SD_gama_HD_sfg_err = error_calc(0.683, SD_gama_HD_sfg, SD_gama_HD_len)
    SD_gama_HD_pass_err = error_calc(0.683, SD_gama_HD_pass, SD_gama_HD_len)
    
    SD_table['SDSS LD'] = SD_sdss_LD
    # SD_table['SDSS LD Err'] = [SD_sdss_LD_agn_err, SD_sdss_LD_sfg_err, SD_sdss_LD_pass_err]
    SD_table['SDSS LDMD'] = SD_sdss_LDMD
    # SD_table['SDSS LDMD Err'] = [SD_sdss_LDMD_agn_err, SD_sdss_LDMD_sfg_err, SD_sdss_LDMD_pass_err]
    SD_table['SDSS MDHD'] = SD_sdss_MDHD
    # SD_table['SDSS MDHD Err'] = [SD_sdss_MDHD_agn_err, SD_sdss_MDHD_sfg_err, SD_sdss_MDHD_pass_err]
    SD_table['SDSS HD'] = SD_sdss_HD
    # SD_table['SDSS HD Err'] = [SD_sdss_HD_agn_err, SD_sdss_HD_sfg_err, SD_sdss_HD_pass_err]
    SD_table['GAMA LD'] = SD_gama_LD
    # SD_table['GAMA LD Err'] = [SD_gama_LD_agn_err, SD_gama_LD_sfg_err, SD_gama_LD_pass_err]
    SD_table['GAMA LDMD'] = SD_gama_LDMD
    # SD_table['GAMA LDMD Err'] = [SD_gama_LDMD_agn_err, SD_gama_LDMD_sfg_err, SD_gama_LDMD_pass_err]
    SD_table['GAMA MDHD'] = SD_gama_MDHD
    # SD_table['GAMA MDHD Err'] = [SD_gama_MDHD_agn_err, SD_gama_MDHD_sfg_err, SD_gama_MDHD_pass_err]
    SD_table['GAMA HD'] = SD_gama_HD
    # SD_table['GAMA HD Err'] = [SD_gama_HD_agn_err, SD_gama_HD_sfg_err, SD_gama_HD_pass_err]
    
    SD_table['SDSS LD'].info.format = '.2f'
    # SD_table['SDSS LD Err'].info.format = '.2f'
    SD_table['SDSS LDMD'].info.format = '.2f'
    # SD_table['SDSS LDMD Err'].info.format = '.2f'
    SD_table['SDSS MDHD'].info.format = '.2f'
    # SD_table['SDSS MDHD Err'].info.format = '.2f'
    SD_table['SDSS HD'].info.format = '.2f'
    # SD_table['SDSS HD Err'].info.format = '.2f'
    SD_table['GAMA LD'].info.format = '.2f'
    # SD_table['GAMA LD Err'].info.format = '.2f'
    SD_table['GAMA LDMD'].info.format = '.2f'
    # SD_table['GAMA LDMD Err'].info.format = '.2f'
    SD_table['GAMA MDHD'].info.format = '.2f'
    # SD_table['GAMA MDHD Err'].info.format = '.2f'
    SD_table['GAMA HD'].info.format = '.2f'
    # SD_table['GAMA HD Err'].info.format = '.2f'
    # SD_table.pprint(max_width=-1)
    
    gama_agn = [SD_gama_LD[0], SD_gama_LDMD[0], SD_gama_MDHD[0], SD_gama_HD[0]]
    gama_agn_err = [SD_gama_LD_agn_err, SD_gama_LDMD_agn_err, SD_gama_MDHD_agn_err, SD_gama_HD_agn_err]
    gama_sfg = [SD_gama_LD[1], SD_gama_LDMD[1], SD_gama_MDHD[1], SD_gama_HD[1]]
    gama_sfg_err = [SD_gama_LD_sfg_err, SD_gama_LDMD_sfg_err, SD_gama_MDHD_sfg_err, SD_gama_HD_sfg_err]
    gama_pass = [SD_gama_LD[2], SD_gama_LDMD[2], SD_gama_MDHD[2], SD_gama_HD[2]]
    gama_pass_err = [SD_gama_LD_pass_err, SD_gama_LDMD_pass_err, SD_gama_MDHD_pass_err, SD_gama_HD_pass_err]
    plt.population_plot(gama_agn, gama_agn_err, gama_sfg, gama_sfg_err,
                        gama_pass, gama_pass_err,
                        ['GAMA AGN', 'GAMA SFG', 'GAMA Passive'], ymin=0, ymax=0.8)

    sdss_agn = [SD_sdss_LD[0], SD_sdss_LDMD[0], SD_sdss_MDHD[0], SD_sdss_HD[0]]
    sdss_agn_err = [SD_sdss_LD_agn_err, SD_sdss_LDMD_agn_err, SD_sdss_MDHD_agn_err, SD_sdss_HD_agn_err]
    sdss_sfg = [SD_sdss_LD[1], SD_sdss_LDMD[1], SD_sdss_MDHD[1], SD_sdss_HD[1]]
    sdss_sfg_err = [SD_sdss_LD_sfg_err, SD_sdss_LDMD_sfg_err, SD_sdss_MDHD_sfg_err, SD_sdss_HD_sfg_err]
    sdss_pass = [SD_sdss_LD[2], SD_sdss_LDMD[2], SD_sdss_MDHD[2], SD_sdss_HD[2]]
    sdss_pass_err = [SD_sdss_LD_pass_err, SD_sdss_LDMD_pass_err, SD_sdss_MDHD_pass_err, SD_sdss_HD_pass_err]
    plt.population_plot(sdss_agn, sdss_agn_err, sdss_sfg, sdss_sfg_err,
                        sdss_pass, sdss_pass_err,
                        ['SDSS AGN', 'SDSS SFG', 'SDSS Passive'], ymin=0, ymax=0.52)
    
    return

population_fracs()