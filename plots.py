
# imports
import numpy as np

from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from astropy.coordinates import Angle
from astropy.table import Table
import astropy.units as u

# plt.style.use('ggplot') # or 'default'
# plt.rcParams['axes.grid'] = True
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.linewidth'] = 1.5

# constants
currentFig = 1

def all_histograms(param, xlabel, mass_limit=10, loc=0, log=False, sdss_bins=None,
                   xmin=None, xmax=None, ymin=None, ymax=None) :
    
    if mass_limit==10 :
        mass_string = 'logMass10_'
    else :
        mass_string = ''
    
    import scipy.stats as sp
    import warnings
    warnings.filterwarnings('ignore')
    
    # CARS
    main_dir = 'catalogs/CARS_SDSS/'
    SDSS_CARS = Table.read(main_dir + 'CARS_SDSS_complete_masslimited.fits')
    sdss_cars_len = len(SDSS_CARS)
    CARS_weight = np.ones(sdss_cars_len)/sdss_cars_len
    CARS_label = 'CARS in SDSS (%s)' % sdss_cars_len
    
    # GAMA
    gama_dir = 'catalogs/CARS_GAMA/'
    gama_blagn = Table.read(gama_dir + 'GAMA_comparison_BLAGN_logMass10_complete_0-38.fits')
    gama_sey = Table.read(gama_dir + 'GAMA_comparison_Seyferts_logMass10_complete_0-39.fits')
    gama_lin = Table.read(gama_dir + 'GAMA_comparison_LINERs_logMass10_complete_0-109.fits')
    gama_comp = Table.read(gama_dir + 'GAMA_comparison_Comps_logMass10_complete_0-125.fits')
    gama_sfg = Table.read(gama_dir + 'GAMA_comparison_SFGs_logMass10_complete_0-354.fits')
    gama_pass = Table.read(gama_dir + 'GAMA_comparison_Passives_logMass10_complete_0-52.fits')
    gama_elg = Table.read(gama_dir + 'GAMA_comparison_not-ELGs_logMass10_complete_0-641.fits')
    
    gama_blagn_len = len(gama_blagn)
    gama_sey_len = len(gama_sey)
    gama_lin_len = len(gama_lin)
    gama_comp_len = len(gama_comp)
    gama_sfg_len = len(gama_sfg)
    gama_pass_len = len(gama_pass)
    gama_elg_len = len(gama_elg)
    
    gama_blagn_weight = np.ones(gama_blagn_len)/gama_blagn_len
    gama_sey_weight = np.ones(gama_sey_len)/gama_sey_len
    gama_lin_weight = np.ones(gama_lin_len)/gama_lin_len
    gama_comp_weight = np.ones(gama_comp_len)/gama_comp_len
    gama_sfg_weight = np.ones(gama_sfg_len)/gama_sfg_len
    gama_pass_weight = np.ones(gama_pass_len)/gama_pass_len
    gama_elg_weight = np.ones(gama_elg_len)/gama_elg_len
    
    # Compute the Kolmogorov-Smirnov statistic on 2 samples
    ks_gama_blagn, pval_gama_blagn = sp.ks_2samp(SDSS_CARS[param], gama_blagn[param])
    ks_gama_sey, pval_gama_sey = sp.ks_2samp(SDSS_CARS[param], gama_sey[param])
    ks_gama_lin, pval_gama_lin = sp.ks_2samp(SDSS_CARS[param], gama_lin[param])
    ks_gama_comp, pval_gama_comp = sp.ks_2samp(SDSS_CARS[param], gama_comp[param])
    ks_gama_sfg, pval_gama_sfg = sp.ks_2samp(SDSS_CARS[param], gama_sfg[param])
    ks_gama_pass, pval_gama_pass = sp.ks_2samp(SDSS_CARS[param], gama_pass[param])
    ks_gama_elg, pval_gama_elg = sp.ks_2samp(SDSS_CARS[param], gama_elg[param])
    
    # Compute the Anderson-Darling test for k-samples
    (ad_gama_blagn, critvals_gama_blagn,
     siglvl_gama_blagn) = sp.anderson_ksamp([SDSS_CARS[param], gama_blagn[param]])
    (ad_gama_sey, critvals_gama_sey,
     siglvl_gama_sey) = sp.anderson_ksamp([SDSS_CARS[param], gama_sey[param]])
    (ad_gama_lin, critvals_gama_lin,
     siglvl_gama_lin) = sp.anderson_ksamp([SDSS_CARS[param], gama_lin[param]])
    (ad_gama_comp, critvals_gama_comp,
     siglvl_gama_comp) = sp.anderson_ksamp([SDSS_CARS[param], gama_comp[param]])
    (ad_gama_sfg, critvals_gama_sfg,
     siglvl_gama_sfg) = sp.anderson_ksamp([SDSS_CARS[param], gama_sfg[param]])
    (ad_gama_pass, critvals_gama_pass,
     siglvl_gama_pass) = sp.anderson_ksamp([SDSS_CARS[param], gama_pass[param]])
    (ad_gama_elg, critvals_gama_elg,
     siglvl_gama_elg) = sp.anderson_ksamp([SDSS_CARS[param], gama_elg[param]])
    
    gama_blagn_tuple = (gama_blagn_len, ks_gama_blagn, pval_gama_blagn, ad_gama_blagn, 100*siglvl_gama_blagn)
    gama_sey_tuple = (gama_sey_len, ks_gama_sey, pval_gama_sey, ad_gama_sey, 100*siglvl_gama_sey)
    gama_lin_tuple = (gama_lin_len, ks_gama_lin, pval_gama_lin, ad_gama_lin, 100*siglvl_gama_lin)
    gama_comp_tuple = (gama_comp_len, ks_gama_comp, pval_gama_comp, ad_gama_comp, 100*siglvl_gama_comp)
    gama_sfg_tuple = (gama_sfg_len, ks_gama_sfg, pval_gama_sfg, ad_gama_sfg, 100*siglvl_gama_sfg)
    gama_pass_tuple = (gama_pass_len, ks_gama_pass, pval_gama_pass, ad_gama_pass, 100*siglvl_gama_pass)
    gama_elg_tuple = (gama_elg_len, ks_gama_elg, pval_gama_elg, ad_gama_elg, 100*siglvl_gama_elg)
    
    gama_blagn_label = 'GAMA BLAGN (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_blagn_tuple
    gama_sey_label = 'GAMA Seyferts (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_sey_tuple
    gama_lin_label = 'GAMA LINERs (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_lin_tuple
    gama_comp_label = 'GAMA Composites (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_comp_tuple
    gama_sfg_label = 'GAMA SFGs (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_sfg_tuple
    gama_pass_label = 'GAMA Passives (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_pass_tuple
    gama_elg_label = 'GAMA not ELGs (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % gama_elg_tuple
    
    # SDSS
    sdss_dir = 'catalogs/CARS_SDSS/'
    sdss_blagn = Table.read(sdss_dir + 'SDSS_comparison_BLAGN_logMass10_complete_0-519.fits')
    sdss_sey = Table.read(sdss_dir + 'SDSS_comparison_Seyferts_logMass10_complete_0-2196.fits')
    sdss_lin = Table.read(sdss_dir + 'SDSS_comparison_LINERs_logMass10_complete_0-12189.fits')
    sdss_comp = Table.read(sdss_dir + 'SDSS_comparison_Comps_logMass10_complete_0-11611.fits')
    sdss_sfg = Table.read(sdss_dir + 'SDSS_comparison_SFGs_logMass10_complete_0-11093.fits')
    sdss_pass = Table.read(sdss_dir + 'SDSS_comparison_Passives_logMass10_complete_0-946.fits')
    sdss_elg = Table.read(sdss_dir + 'SDSS_comparison_not-ELGs_logMass10_complete_0-13553.fits')
    
    sdss_blagn_len = len(sdss_blagn)
    sdss_sey_len = len(sdss_sey)
    sdss_lin_len = len(sdss_lin)
    sdss_comp_len = len(sdss_comp)
    sdss_sfg_len = len(sdss_sfg)
    sdss_pass_len = len(sdss_pass)
    sdss_elg_len = len(sdss_elg)
    
    sdss_blagn_weight = np.ones(sdss_blagn_len)/sdss_blagn_len
    sdss_sey_weight = np.ones(sdss_sey_len)/sdss_sey_len
    sdss_lin_weight = np.ones(sdss_lin_len)/sdss_lin_len
    sdss_comp_weight = np.ones(sdss_comp_len)/sdss_comp_len
    sdss_sfg_weight = np.ones(sdss_sfg_len)/sdss_sfg_len
    sdss_pass_weight = np.ones(sdss_pass_len)/sdss_pass_len
    sdss_elg_weight = np.ones(sdss_elg_len)/sdss_elg_len
    
    # Compute the Kolmogorov-Smirnov statistic on 2 samples
    ks_sdss_blagn, pval_sdss_blagn = sp.ks_2samp(SDSS_CARS[param], sdss_blagn[param])
    ks_sdss_sey, pval_sdss_sey = sp.ks_2samp(SDSS_CARS[param], sdss_sey[param])
    ks_sdss_lin, pval_sdss_lin = sp.ks_2samp(SDSS_CARS[param], sdss_lin[param])
    ks_sdss_comp, pval_sdss_comp = sp.ks_2samp(SDSS_CARS[param], sdss_comp[param])
    ks_sdss_sfg, pval_sdss_sfg = sp.ks_2samp(SDSS_CARS[param], sdss_sfg[param])
    ks_sdss_pass, pval_sdss_pass = sp.ks_2samp(SDSS_CARS[param], sdss_pass[param])
    ks_sdss_elg, pval_sdss_elg = sp.ks_2samp(SDSS_CARS[param], sdss_elg[param])
    
    # Compute the Anderson-Darling test for k-samples
    (ad_sdss_blagn, critvals_sdss_blagn,
     siglvl_sdss_blagn) = sp.anderson_ksamp([SDSS_CARS[param], sdss_blagn[param]])
    (ad_sdss_sey, critvals_sdss_sey,
      siglvl_sdss_sey) = sp.anderson_ksamp([SDSS_CARS[param], sdss_sey[param]])
    (ad_sdss_lin, critvals_sdss_lin,
      siglvl_sdss_lin) = sp.anderson_ksamp([SDSS_CARS[param], sdss_lin[param]])
    (ad_sdss_comp, critvals_sdss_comp,
      siglvl_sdss_comp) = sp.anderson_ksamp([SDSS_CARS[param], sdss_comp[param]])
    (ad_sdss_sfg, critvals_sdss_sfg,
      siglvl_sdss_sfg) = sp.anderson_ksamp([SDSS_CARS[param], sdss_sfg[param]])
    (ad_sdss_pass, critvals_sdss_pass,
      siglvl_sdss_pass) = sp.anderson_ksamp([SDSS_CARS[param], sdss_pass[param]])
    (ad_sdss_elg, critvals_sdss_elg,
      siglvl_sdss_elg) = sp.anderson_ksamp([SDSS_CARS[param], sdss_elg[param]])
    
    sdss_blagn_tuple = (sdss_blagn_len, ks_sdss_blagn, pval_sdss_blagn, ad_sdss_blagn, 100*siglvl_sdss_blagn)
    sdss_sey_tuple = (sdss_sey_len, ks_sdss_sey, pval_sdss_sey, ad_sdss_sey, 100*siglvl_sdss_sey)
    sdss_lin_tuple = (sdss_lin_len, ks_sdss_lin, pval_sdss_lin, ad_sdss_lin, 100*siglvl_sdss_lin)
    sdss_comp_tuple = (sdss_comp_len, ks_sdss_comp, pval_sdss_comp, ad_sdss_comp, 100*siglvl_sdss_comp)
    sdss_sfg_tuple = (sdss_sfg_len, ks_sdss_sfg, pval_sdss_sfg, ad_sdss_sfg, 100*siglvl_sdss_sfg)
    sdss_pass_tuple = (sdss_pass_len, ks_sdss_pass, pval_sdss_pass, ad_sdss_pass, 100*siglvl_sdss_pass)
    sdss_elg_tuple = (sdss_elg_len, ks_sdss_elg, pval_sdss_elg, ad_sdss_elg, 100*siglvl_sdss_elg)
    
    sdss_blagn_label = 'SDSS BLAGN (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_blagn_tuple
    sdss_sey_label = 'SDSS Seyferts (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_sey_tuple
    sdss_lin_label = 'SDSS LINERs (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_lin_tuple
    sdss_comp_label = 'SDSS Composites (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_comp_tuple
    sdss_sfg_label = 'SDSS SFGs (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_sfg_tuple
    sdss_pass_label = 'SDSS Passives (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_pass_tuple
    sdss_elg_label = 'SDSS not ELGs (%s)\nK-S=%.3g, p=%.3g\nA-D=%.3g, SL=%.3g%%' % sdss_elg_tuple
    
    multi_histo([SDSS_CARS[param], gama_blagn[param], sdss_blagn[param]],
                [CARS_label, gama_blagn_label, sdss_blagn_label],
                xlabel, ['k','cyan','darkcyan'],
                [CARS_weight, gama_blagn_weight, sdss_blagn_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_BLAGN.png', 
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    multi_histo([SDSS_CARS[param], gama_sey[param], sdss_sey[param]],
                [CARS_label, gama_sey_label, sdss_sey_label],
                xlabel, ['k','red','darkred'],
                [CARS_weight, gama_sey_weight, sdss_sey_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_Seyferts.png',
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    multi_histo([SDSS_CARS[param], gama_lin[param], sdss_lin[param]],
                [CARS_label, gama_lin_label, sdss_lin_label],
                xlabel, ['k','springgreen','darkgreen'],
                [CARS_weight, gama_lin_weight, sdss_lin_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_LINERs.png',
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    multi_histo([SDSS_CARS[param], gama_comp[param], sdss_comp[param]],
                [CARS_label, gama_comp_label, sdss_comp_label],
                xlabel, ['k','magenta','purple'],
                [CARS_weight, gama_comp_weight, sdss_comp_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_Comps.png',
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    multi_histo([SDSS_CARS[param], gama_sfg[param], sdss_sfg[param]],
                [CARS_label, gama_sfg_label, sdss_sfg_label],
                xlabel, ['k','dodgerblue','darkblue'],
                [CARS_weight, gama_sfg_weight, sdss_sfg_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_SFGs.png',
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    multi_histo([SDSS_CARS[param], gama_pass[param], sdss_pass[param]],
                [CARS_label, gama_pass_label, sdss_pass_label],
                xlabel, ['k','orange','saddlebrown'],
                [CARS_weight, gama_pass_weight, sdss_pass_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_Passives.png',
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    multi_histo([SDSS_CARS[param], gama_elg[param], sdss_elg[param]],
                [CARS_label, gama_elg_label, sdss_elg_label],
                xlabel, ['k','gold','chocolate'],
                [CARS_weight, gama_elg_weight, sdss_elg_weight], log=log, nbins=sdss_bins,
                outfile='histograms/histo_' + mass_string + param + '_not-ELGs.png',
                xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, location=loc)
    
    return

def all_scatter_plots(catalog, xparam, yparam, loc=0,
                      logx=False, logy=False, loglog=False,
                      xmin=None, xmax=None, ymin=None, ymax=None) :
    
    dictionary = {# AGN properties
                  'Z':'Redshift',
                  'logMass':r'$\log($Stellar Mass of Target$/M_\odot)$',
                  'SFR':r'Star Formation Rate ($M_\odot$ yr$^{-1}$)',
                  'OIII_FWHM':r'[O III] $\lambda5007$ FWHM (km s$^{-1}$)',
                  # best environmental parameters
                  'SurfaceDensity':r'Surface Density (Mpc$^{-2}$)',
                  'Closest_Distance':'Distance to Closest Companion (kpc)',
                  # other environmental parameters
                  'AGEPar':r'AGE Parameter (Mpc$^{-1}$)',
                  'Companions':'Number of Companions',
                  'Most_Massive_Mass':r'Most Massive Companion Mass $(\log_{10}[M_\odot])$',
                  'Overdensity':'Overdensity',
                  'd_CoM':'Projected Distance to the\nCenter of Stellar Mass (arcmin)',
                  'Most_Massive_Distance':'Distance to Most Massive Companion (kpc)',
                  'Closest_Mass':r'Closest Companion Mass $(\log_{10}[M_\odot])$'}
    
    # CARS
    main_dir = 'catalogs/CARS_SDSS/'
    SDSS_CARS = Table.read(main_dir + 'CARS_SDSS_complete_masslimited.fits')
    CARS_sfr_sigma = Table.read('catalogs/raw_cats/CARS_SFR_OIIIFWHM.fits')
    SDSS_CARS['SFR'] = CARS_sfr_sigma['SFR']
    SDSS_CARS['OIII_FWHM'] = 2*np.sqrt( 2*np.log(2) )*CARS_sfr_sigma['OIII5007_vel_dispersion']
    
    # GAMA
    gama_dir = 'catalogs/CARS_GAMA/'
    gama_blagn = Table.read(gama_dir + 'GAMA_comparison_BLAGN_logMass10_complete_0-38.fits')
    gama_sey = Table.read(gama_dir + 'GAMA_comparison_Seyferts_logMass10_complete_0-39.fits')
    gama_lin = Table.read(gama_dir + 'GAMA_comparison_LINERs_logMass10_complete_0-109.fits')
    gama_comp = Table.read(gama_dir + 'GAMA_comparison_Comps_logMass10_complete_0-125.fits')
    gama_sfg = Table.read(gama_dir + 'GAMA_comparison_SFGs_logMass10_complete_0-354.fits')
    gama_pass = Table.read(gama_dir + 'GAMA_comparison_Passives_logMass10_complete_0-52.fits')
    gama_elg = Table.read(gama_dir + 'GAMA_comparison_not-ELGs_logMass10_complete_0-641.fits')
    
    # SDSS
    sdss_dir = 'catalogs/CARS_SDSS/'
    sdss_blagn = Table.read(sdss_dir + 'SDSS_comparison_BLAGN_logMass10_complete_0-519.fits')
    sdss_sey = Table.read(sdss_dir + 'SDSS_comparison_Seyferts_logMass10_complete_0-2196.fits')
    sdss_lin = Table.read(sdss_dir + 'SDSS_comparison_LINERs_logMass10_complete_0-12189.fits')
    sdss_comp = Table.read(sdss_dir + 'SDSS_comparison_Comps_logMass10_complete_0-11611.fits')
    sdss_sfg = Table.read(sdss_dir + 'SDSS_comparison_SFGs_logMass10_complete_0-11093.fits')
    sdss_pass = Table.read(sdss_dir + 'SDSS_comparison_Passives_logMass10_complete_0-946.fits')
    sdss_elg = Table.read(sdss_dir + 'SDSS_comparison_not-ELGs_logMass10_complete_0-13553.fits')
    
    if catalog == 'SDSS' :
        # colours = ['k', 'darkcyan', 'darkred', 'darkgreen',
        #            'purple', 'darkblue', 'saddlebrown', 'chocolate']
        colours = ['k', 'cyan', 'red', 'springgreen',
                   'magenta', 'dodgerblue', 'orange', 'gold']
        xarrays = [SDSS_CARS[xparam], sdss_blagn[xparam], sdss_sey[xparam], sdss_lin[xparam],
                    sdss_comp[xparam], sdss_sfg[xparam], sdss_pass[xparam], sdss_elg[xparam]]
        yarrays = [SDSS_CARS[yparam], sdss_blagn[yparam], sdss_sey[yparam], sdss_lin[yparam],
                    sdss_comp[yparam], sdss_sfg[yparam], sdss_pass[yparam], sdss_elg[yparam]]
        labels = ['SDSS CARS', 'SDSS BLAGN', 'SDSS Seyferts', 'SDSS LINERs',
                    'SDSS Composites', 'SDSS SFGs', 'SDSS Passives', 'SDSS not ELGs']
        alpha = 0.6
        mark_every = 40 # use 2 for [O III] plots, otherwise use 30
    
    if catalog == 'GAMA' :
        colours = ['k', 'cyan', 'red', 'springgreen',
                   'magenta', 'dodgerblue', 'orange', 'gold']
        xarrays = [SDSS_CARS[xparam], gama_blagn[xparam], gama_sey[xparam], gama_lin[xparam],
                   gama_comp[xparam], gama_sfg[xparam], gama_pass[xparam], gama_elg[xparam]]
        yarrays = [SDSS_CARS[yparam], gama_blagn[yparam], gama_sey[yparam], gama_lin[yparam],
                   gama_comp[yparam], gama_sfg[yparam], gama_pass[yparam], gama_elg[yparam]]
        labels = ['SDSS CARS', 'GAMA BLAGN', 'GAMA Seyferts', 'GAMA LINERs',
                   'GAMA Composites', 'GAMA SFGs', 'GAMA Passives', 'GAMA not ELGs']
        alpha = 0.6
        mark_every = 1
    
    multi_plot(xarrays, yarrays, labels, colours,
               dictionary[xparam], dictionary[yparam], location=loc,
               outfile='correlation_plots/' + catalog + '_' + yparam + '_vs_' + xparam + '.png',
               xlog=logx, ylog=logy, log=loglog, alph=alpha, mark_nth=mark_every,
               xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    
    return

def contour3d(xx, yy, zz, xlabel, ylabel, zlabel) :
    
    from mpl_toolkits.mplot3d import Axes3D
    
    global currentFig
    
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    
#    ax.contourf(xx, yy, zz)
    ax.plot_surface(xx, yy, zz)
    
    ax.elev = 30
    ax.azim = 225
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_zlabel(zlabel, fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def density_plot(catalog) :
    
    if catalog == 'both' :
        redshifts = []
        SDSS_companions = []
        SDSS_volumes = []
        GAMA_companions = []
        GAMA_volumes = []
        for i in range(0, 12) :
            redshift, number, volume = env1.main('SDSS', i, return_vals=True)
            redshifts.append(redshift)
            SDSS_companions.append(number)
            SDSS_volumes.append(volume.value)
            if i in [5,6,8] :
                redshift, number, volume = env1.main('GAMA', i, return_vals=True)
                GAMA_companions.append(number)
                GAMA_volumes.append(volume.value)
            else :
                GAMA_companions.append(np.nan)
                GAMA_volumes.append(np.nan)
        multi(redshifts, 'Redshift',
              np.array(SDSS_companions)/np.array(SDSS_volumes),
              np.array(GAMA_companions)/np.array(GAMA_volumes),
              'Density of Companions (Mpc$^{-3}$)')
    
    redshifts = []
    number_of_companions = []
    volumes = []
    
    if catalog == 'SDSS' :
        for i in range(0, 12) :
            redshift, number, volume = env1.main('SDSS', i, return_vals=True)
            redshifts.append(redshift)
            number_of_companions.append(number)
            volumes.append(volume.value)
        plot(redshifts, 'Redshift',
             np.array(number_of_companions)/np.array(volumes),
             'Density of Companions (Mpc$^{-3}$)')
    
    if catalog == 'GAMA' :
        for i in [5,6,8] :
            redshift, number, volume = env1.main('GAMA', i, return_vals=True)
            redshifts.append(redshift)
            number_of_companions.append(number)
            volumes.append(volume.value)
        plot(redshifts, 'Redshift',
             np.array(number_of_companions)/np.array(volumes),
             'Density of Companions (Mpc$^{-3}$)')
    
    return

def diagram_BPT(x1, y1, x2, y2, x3, y3, x4, y4, standard=True,
                xmin=None, xmax=None, ymin=None, ymax=None,
                xlabel=None, ylabel=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(10.5, 6.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    # for the SDSS sample check
    ax.hist2d(x1, y1, bins=200, cmap='Reds', norm=LogNorm())
    ax.hist2d(x2, y2, bins=117, cmap='Greens', norm=LogNorm())
    ax.hist2d(x3, y3, bins=156, cmap='Purples', norm=LogNorm())
    ax.hist2d(x4, y4, bins=231, cmap='Blues', norm=LogNorm())
    
    # for the GAMA sample check
    # ax.hist2d(x1, y1, bins=5, cmap='Reds', norm=LogNorm())
    # ax.hist2d(x2, y2, bins=3, cmap='Greens', norm=LogNorm())
    # ax.hist2d(x3, y3, bins=12, cmap='Purples', norm=LogNorm())
    # ax.hist2d(x4, y4, bins=35, cmap='Blues', norm=LogNorm())
    
    # unclass_x = [0.52699786, 0.6472921, 0.51526904, 0.10082114, 0.4783239,
    #               0.64564115, 0.47279504, 0.48915058, 0.5316156, 0.52157855,
    #               0.072221704, 0.6222592, 0.061153412]
    # unclass_y = [0.20510054, 0.15230991, 0.3134719, -1.4588171, 0.8099388,
    #              0.41682592, 0.12564416, 0.5535542, 0.8346556, 0.4626364,
    #              -0.5899385, 0.4983, -0.36166635]
    # ax.plot(unclass_x, unclass_y, 'kx', label='BPT Classification Fails (13 galaxies)')
    
    # CARS_x_7p5 = [-0.2558046120106365, -0.3039520129538409, -0.2417594381357093,
    #               0.06476555533411806, 0.275461262291009, -0.3938455204465459,
    #               -0.05406315470653154, -1.0767286502681763, -0.21683621439121184,
    #               -0.2583153650781086, -0.3611262953855419, -0.041477366726339726]
    # CARS_y_7p5 = [1.0442984352563223, -0.4632845953092678, -0.17709178145221774,
    #               0.7852122081665016, 0.701362126702604, -0.3101074263261445,
    #               0.9213412646212272, 0.07531550781816444, 0.43745460649821044,
    #               0.6040263849321078, 0.02695921366978017, 0.8956928513689375]
    # ax.plot(CARS_x_7p5, CARS_y_7p5, 'ko', label='Fluxes from .eline MUSE Tables in 3" Aperture')
    
    if (standard==True) :
        xmin, xmax, ymin, ymax = -2.2, 0.7, -1.5, 1.4
        xlabel = r'$\log_{10}$([N II] $\lambda 6583 / \rm H\alpha$)'
        ylabel = r'$\log_{10}$([O III] $\lambda 5007 / \rm H\beta$)'
        
        kauff_x = np.linspace(xmin, xmax, 10000)
        kauff_y = 0.61/(kauff_x - 0.05) + 1.3 # from Kauffmann+ 2003
        #                       (http://adsabs.harvard.edu/abs/2003MNRAS.346.1055K)
        kauff_mask = ( (kauff_x > xmin) & (kauff_x < 0.05) )
        ax.plot(kauff_x[kauff_mask], kauff_y[kauff_mask], color='k', linestyle='--', lw=1.5,
                label='Kauffmann et al. (2003)', zorder=4)
        
        kewl_x = np.linspace(xmin, xmax, 10000)
        kewl_y = 0.61/(kewl_x - 0.47) + 1.19 # from Kewley+ 2001
        #                       (http://adsabs.harvard.edu/abs/2001ApJ...556..121K)
        kewl_mask = ( (kewl_x > xmin) & (kewl_x < 0.47) )
        ax.plot(kewl_x[kewl_mask], kewl_y[kewl_mask], color='k', linestyle='-', lw=1.5,
                label='Kewley et al. (2001)', zorder=4)
        
        schaw_x = np.linspace(-0.18380687748557267, xmax, 10000)
        schaw_y = 1.05*schaw_x + 0.45 # from Schawinski+ 2012
        #                       (http://adsabs.harvard.edu/abs/2007MNRAS.382.1415S)
        ax.plot(schaw_x, schaw_y, color='k', linestyle='-.', lw=1.5,
                label='Schawinski et al. (2007)', zorder=4)
        
        ax.text(-0.75, 1.1, 'Seyfert', fontsize=18)
        ax.text(0.25, -0.4, 'LINER', fontsize=18)
        ax.text(-1, -0.5, 'SFG', fontsize=18)
        ax.text(-0.14, -1, 'Comp', fontsize=18)
        
        # ax.vlines(0.47, ymin, ymax, color='b', linestyle='-.')
        
    ax.set_xlabel(xlabel, fontsize=19)
    ax.set_ylabel(ylabel, fontsize=19)
        
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()
    
    handles = [handles[1], handles[0], handles[2]]
    labels = [labels[1], labels[0], labels[2]]
    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=16,
              loc = 'lower left')
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=19)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig('correlation_plots/diagnostic_SDSS_BPT.pdf', overwrite=True)
    # plt.close()
    
    return

def diagram_WHAN(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, standard=True,
                 xmin=None, xmax=None, ymin=None, ymax=None,
                 xlabel=None, ylabel=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(10.5, 6.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)  
    
    # for the SDSS sample check
    # ax.hist2d(x1, y1, bins=40, cmap='Reds', norm=LogNorm())
    # ax.hist2d(x2, y2, bins=113, cmap='Greens', norm=LogNorm())
    # ax.hist2d(x3, y3, bins=30, cmap='Purples', norm=LogNorm())
    # ax.hist2d(x4, y4, bins=86, cmap='Blues', norm=LogNorm())
    # ax.hist2d(x5, y5, bins=80, cmap='Oranges', norm=LogNorm())
    
    # for the GAMA sample check
    ax.hist2d(x1, y1, bins=10, cmap='Reds', norm=LogNorm())
    ax.hist2d(x2, y2, bins=20, cmap='Greens', norm=LogNorm())
    ax.hist2d(x3, y3, bins=8, cmap='Purples', norm=LogNorm())
    ax.hist2d(x4, y4, bins=35, cmap='Blues', norm=LogNorm())
    ax.hist2d(x5, y5, bins=18, cmap='Oranges', norm=LogNorm())
    
    if (standard==True) :
        xmin, xmax, ymin, ymax = -2.2, 1, -0.7, 3.3
        xlabel = r'$\log_{10}$([N II] $\lambda 6583 / \rm H\alpha$)'
        ylabel = r'$\log_{10}(W_{\rm H\alpha}/ {\rm \AA})$'
        
        # regular WHAN
        ax.vlines(-0.32, ymin, ymax, color='k', linestyle='--', lw=1.5,
                  label=r'Kauffmann et al. (2003)$^{\rm T}$', zorder=2)
        ax.hlines(np.log10(6), -0.32, xmax, color='k', linestyle='-', lw=1.5,
                  label=r'Kewley et al. (2006)$^{\rm T}$', zorder=2)
        ax.vlines(-0.4, ymin, ymax, color='k', linestyle='-.', lw=1.5,
                  label=r'Stasi$\rm \'n$kska et al. (2006)$^{\rm T}$',
                  zorder=2)
        ax.hlines(np.log10(0.5), -2.2, xmax, color='k', linestyle=':', lw=1.5, zorder=2)
        
        xs = np.linspace(-3, 1.5, 1000)
        ys = -xs + np.log10(0.5)
        ax.plot(xs, ys, 'k:', lw=1.5)
        
        ax.text(-1, 1.75, 'SFG', fontsize=18)
        ax.text(-0.4, 1.75, 'Composite', rotation='vertical', fontsize=18)
        ax.text(0, 1.5, 'Seyfert', fontsize=18)
        ax.text(0.55, 0.25, 'LINER', fontsize=18)
        ax.text(-0.9, -0.1, 'Passive', fontsize=18)
    
    ax.set_xlabel(xlabel, fontsize=19)
    ax.set_ylabel(ylabel, fontsize=19)
        
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()
    
    # handles = [handles[1], handles[0], handles[2]]
    # labels = [labels[1], labels[0], labels[2]]
    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=14,
              loc = 'lower left')
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=19)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig('correlation_plots/diagnostic_SDSS_WHAN.pdf', overwrite=True)
    # plt.close()
    
    return

def emission_line_map(vals, wcs, label) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = plt.subplot(111)#, projection = wcs, slices=('x', 'y', 1, 1))
    
    # ax.coords[0].set_axislabel('Right Ascension')
    # ax.coords[1].set_axislabel('Declination')
    
    cmap = cm.plasma # or use jet
    cmap.set_bad('white', 1)
    
    # ax.set_xlim(51, 106)
    # ax.set_ylim(37, 126)
    
    frame = plt.imshow(vals, origin = 'lower', cmap = cmap,
                       norm = LogNorm(), interpolation = 'None')
    cbar = plt.colorbar()
    cbar.ax.set_title(r"$10^{-16}$ erg s$^{-1}$ cm$^{-2}$")
    cbar.set_label(label, fontsize = 15)
    
    plt.tight_layout()
    plt.show()
    
    return

def env_hists(param, log=False, xmin=None) :
    
    CARS = Table.read('catalogs/CARS_SDSS/CARS_SDSS_complete_masslimited.fits')
    SDSS = Table.read('catalogs/CARS_SDSS/SDSS-galaxies_complete_random_2500.fits')
    
    BLAGN = Table.read('catalogs/CARS_GAMA/GAMA_BLAGN_complete_0-52.fits')
    comps = Table.read('catalogs/CARS_GAMA/GAMA_Comp_complete_0-159.fits')
    LINERs = Table.read('catalogs/CARS_GAMA/GAMA_LINER_complete_0-52.fits')
    notELGs = Table.read('catalogs/CARS_GAMA/GAMA_not-ELG_complete_0-2194.fits')
    passives = Table.read('catalogs/CARS_GAMA/GAMA_Passive_complete_0-194.fits')
    seyferts = Table.read('catalogs/CARS_GAMA/GAMA_Seyfert_complete_0-39.fits')
    SFGs = Table.read('catalogs/CARS_GAMA/GAMA_SFG_complete_0-2726.fits')
    
    param_list = [
                   CARS[param],
                   SDSS[param],
                  BLAGN[param],
                  comps[param],
                  LINERs[param],
                   notELGs[param],
                  passives[param],
                  seyferts[param],
                   SFGs[param]
                  ]
    
    labels = [
               'CARS',
               'SDSS - galaxies',
              'GAMA - BLAGN',
              'GAMA - Comp.',
              'GAMA - LINERs',
               'GAMA - not_ELG',
              'GAMA - Passive',
              'GAMA - Seyferts',
               'GAMA - SFGs'
              ]
    
    multi_histo(param_list, labels, param, log, xmin)
    
    return

def filter_curves() :
    
    uu = Table.read('SDSS_response_curves/filter_curves.fits', hdu=1)
    gg = Table.read('SDSS_response_curves/filter_curves.fits', hdu=2)
    rr = Table.read('SDSS_response_curves/filter_curves.fits', hdu=3)
    ii = Table.read('SDSS_response_curves/filter_curves.fits', hdu=4)
    zz = Table.read('SDSS_response_curves/filter_curves.fits', hdu=5)
    
    # plot_multi([uu['wavelength'], gg['wavelength'], rr['wavelength'],
    #             ii['wavelength'], zz['wavelength']],
    #             [uu['respt'], gg['respt'], rr['respt'],
    #             ii['respt'], zz['respt']],
    #             [r'$u$', r'$g$', r'$r$', r'$i$', r'$z$'],
    #             ['b', 'g', 'r', 'm', 'k'],
    #             r'Wavelength $(\rm \AA)$', 'Quantum Efficiency',
    #             xmin=2800, xmax=11500, ymin=0, ymax=0.58)
    
    filters = Table.read('SDSS_response_curves/response_curves_Doi+_2010_AJ_139_1628.fits')
    
    plot_multi([filters['lambda'], filters['lambda'], filters['lambda'],
                filters['lambda'], filters['lambda']],
                [filters['u'], filters['g'], filters['r'],
                filters['i'], filters['z']],
                [r'$u$', r'$g$', r'$r$', r'$i$', r'$z$'],
                ['b', 'g', 'r', 'm', 'k'],
                r'Wavelength $(\rm \AA)$', 'Quantum Efficiency',
                outfile='correlation_plots/SDSS_filter_response_curves.pdf',
                xmin=2800, xmax=11500, ymin=0, ymax=0.58)
    
    return

def gaussian_plot(centers, sigmas, outfile) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(10,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    import functions as func
    
    xs = []
    ys = []
    colours = ['blue', 'red', 'orange', 'green']
    for i in range(len(centers)) :
        xs = np.linspace(-5, 5, 1000)
        ys = func.gaussian(xs, centers[i], sigmas[i])
        ax.plot(xs, ys, '-', color=colours[i],
                label=r'$\mu=%s$, $\sigma^2=%s$' % (centers[i], sigmas[i]))
    
    ax.set_xlabel(r'$x$', fontsize=21)
    ax.set_ylabel(r'$g(x)$', fontsize=21)
        
    ax.set_xlim(-5, 5)
    ax.set_ylim(-0.05, 1)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=17,
              loc = 'upper right')
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=16)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig(outfile, overwrite=True)
    # plt.close()
    
    return

def histo(param, label, title=None, vert_line=None, num_bins=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.hist(param) #, range=(7, 11.5)) #, bins=num_bins)
    ax.set_xlabel('%s' % label, fontsize = 15)
    
    if vert_line :
        ymin, ymax = ax.get_ylim()
        ax.vlines(vert_line, ymin, ymax,
                  color='r', linestyle='--', label=r'CARS Host')
        plt.legend()
    
    ax.set_title('%s' % title, fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def histo_with_fxn(param, label, weight, xs, ys, fit_label, loc=0,
                   num_bins=None, outfile=None,
                   xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(10,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.hist(param, weights=weight, histtype='step', linestyle='-', color='k',
            linewidth=2, bins=num_bins)
    
    ax.plot(xs, ys, 'r--', label=fit_label)
    
    ax.set_xlabel('%s' % label, fontsize=21)
    ax.set_ylabel('Fractional Frequency', fontsize=21)
        
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()
    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=17,
              loc = loc)
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=16)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig('correlation_plots/' + outfile, overwrite=True)
    # plt.close()
    
    return

def histo2d(xvals, xlab, yvals, ylab, nbins=200, xmin=None, xmax=None,
            ymin=None, ymax=None) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.hist2d(xvals, yvals, bins=nbins, norm=LogNorm())
    # cbar = plt.colorbar()
    
    ax.set_xlabel('%s' % xlab, fontsize=15)
    ax.set_ylabel('%s' % ylab, fontsize=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    # ax.legend()
    plt.tight_layout()
    plt.show()
    
    return

def multi_histo(param_list, labels, xlab, colors, weights, location=0,
                outfile=None, log=False, nbins=None, xmin=None, xmax=None, ymin=None, ymax=None) :
    
    first_param = param_list[0]
    second_param = param_list[1]
    third_param = param_list[2]
    first_weight = weights[0]
    second_weight = weights[1]
    third_weight = weights[2]
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(7,5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)    
    
    # colors = ['k', 'b', 'r', 'g', 'm', 'c', 'y', 'k', 'b']
    # styles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-']
    # styles = ['--', '-']
    
    if xlab == 'Number of Companions' :
        xlab = '1 + Number of Companions'
        first_param = 1 + param_list[0]
        second_param = 1 + param_list[1]
        third_param = 1 + param_list[2]
    
    if xlab == 'Projected Distance to the Center of Stellar Mass/arcmin' :
        xlab = r'1 + Proj. Dist. to the Center of $M_{*}$/arcmin'
        first_param = 1 + param_list[0]
        second_param = 1 + param_list[1]
        third_param = 1 + param_list[2]
    
    if xlab == 'Companions in the GAMA Cyl.' :
        xlab = '1 + Companions in the GAMA Cyl.'
        first_param = 1 + param_list[0]
        second_param = 1 + param_list[1]
        third_param = 1 + param_list[2]
    
    if log==True :
        # alternative to mask out values of 0
        # first_param = first_param[first_param > 0]
        # first_weight = np.ones(len(first_param))/len(first_param)
        # second_param = second_param[second_param > 0]
        # second_weight = np.ones(len(second_param))/len(second_param)
        
        ax.hist(np.log10(first_param), label=labels[0], weights=first_weight,
                histtype='step', linestyle='-', color=colors[0])
        ax.hist(np.log10(second_param), label=labels[1], weights=second_weight,
                histtype='step', linestyle='-', color=colors[1], linewidth=2)#, bins=int(1.5*nbins))
        ax.hist(np.log10(third_param), label=labels[2], weights=third_weight,
                histtype='step', linestyle='--', color=colors[2], linewidth=2, bins=nbins)
        xlabel = r'$\log($' + xlab + ')'
    else :
        ax.hist(first_param, label=labels[0], weights=first_weight,
                histtype='step', linestyle='-', color=colors[0], linewidth=1.5)
        ax.hist(second_param, label=labels[1], weights=second_weight,
                histtype='step', linestyle='-', color=colors[1], linewidth=2.5)#, bins=int(1.5*nbins))
        ax.hist(third_param, label=labels[2], weights=third_weight,
                histtype='step', linestyle='--', color=colors[2], linewidth=2, bins=nbins)
        xlabel = xlab
    
    ax.set_xlabel('%s' % xlabel, fontsize=19)
    ax.set_ylabel('Fractional Frequency', fontsize=19)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=12,
              loc = location)
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=14)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig(outfile, overwrite=True)
    # plt.close(fig)
    
    return

def multi(xvals, xlab, yvals1, yvals2, ylab, xmin=None, xmax=None, ymin=None,
          ymax=None, location='upper left') :
    
    global currentFig        
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.plot(xvals, yvals1, 'ko', label = "%s" % 'SDSS', zorder=2)
#    ax.axhline(np.nanmean(yvals1), color='k', linestyle='-', label='mean')
    ax.axhline(np.nanmedian(yvals1), color='k', linestyle='--', label='median',
               zorder=1)
    
    ax.plot(xvals, yvals2, 'ro', label = "%s" % 'GAMA', zorder=2)
#    ax.axhline(np.nanmean(yvals2), color='r', linestyle='-', label='mean')
    ax.axhline(np.nanmedian(yvals2), color='r', linestyle='--', label='median',
               zorder=1)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel('%s' % xlab, fontsize=15)
    ax.set_ylabel('%s' % ylab, fontsize=15)
    
    ax.legend(loc = location)
    plt.tight_layout()
    plt.show()
    
    return

def multi2(xvals1, yvals1, label1, xvals2, yvals2, label2, xlab, ylab) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig        
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.hist2d(xvals1, yvals1, bins=50, norm=LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('GAMA', fontsize = 15)
    
    ax.plot(xvals2, yvals2, 'ro', label = "%s" % label2, zorder=2)
    
    ax.set_xlim(-4, 3)
    ax.set_ylim(-1, 27)
    
    ax.set_xlabel('%s' % xlab, fontsize=15)
    ax.set_ylabel('%s' % ylab, fontsize=15)
    
    ax.legend()
    plt.tight_layout()
    plt.show()
    
    return

def multi_plot(xarrays, yarrays, labels, colours, xlab, ylab, location=0, outfile=None,
               xlog=False, ylog=False, log=False, alph=0.1, mark_nth=1,
               xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(15,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    if xlog == True :
        ax.semilogx(xarrays[0], yarrays[0], 'o', color=colours[0], label = labels[0], alpha=1.0, zorder=4, ms=9)
        i = 1
        while i < len(xarrays) : # add each dataset in a list to the plot
            ax.semilogx(xarrays[i], yarrays[i], 'o', color=colours[i], label = labels[i], alpha=alph, markevery=mark_nth)
            i += 1
    elif ylog == True :
        ax.semilogy(xarrays[0], yarrays[0], 'o', color=colours[0], label = labels[0], alpha=1.0, zorder=4, ms=9)
        i = 1
        while i < len(xarrays) : # add each dataset in a list to the plot
            ax.semilogy(xarrays[i], yarrays[i], 'o', color=colours[i], label = labels[i], alpha=alph, markevery=mark_nth)
            i += 1
    elif log == True :
        ax.loglog(xarrays[0], yarrays[0], 'o', color=colours[0], label = labels[0], alpha=1.0, zorder=4, ms=9)
        i = 1
        while i < len(xarrays) : # add each dataset in a list to the plot
            ax.loglog(xarrays[i], yarrays[i], 'o', color=colours[i], label = labels[i], alpha=alph, markevery=mark_nth)
            i += 1
    else :
        ax.plot(xarrays[0], yarrays[0], 'o', color=colours[0], label = labels[0], alpha=1.0, zorder=4, ms=9)
        i = 1
        while i < len(xarrays) : # add each dataset in a list to the plot
            ax.plot(xarrays[i], yarrays[i], 'o', color=colours[i], label = labels[i], alpha=alph, markevery=mark_nth)
            i += 1
    
    ax.set_xlabel('%s' % xlab, fontsize=21)
    ax.set_ylabel('%s' % ylab, fontsize=21)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=14,
              loc = location)
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=16)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig(outfile, overwrite=True)
    # plt.close(fig)
    
    return

def plot(xvals, xlab, yvals, ylab, cat_name='SDSS', xmin=None, xmax=None,
         ymin=None, ymax=None, hist2d=False, nbins=200, fit0=False, fit1=False, fit2=False,
         log=False, logy=False, correlation=False) :
    
    from matplotlib.colors import LogNorm
    import scipy.stats as sp
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(8,6)) # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    ax = fig.add_subplot(111) # set axes, figure location
    
    if correlation :
        spear = sp.spearmanr(xvals, yvals, nan_policy='omit')
        print(spear)
    
    if (hist2d == True) :
        plt.hist2d(xvals, yvals, bins=nbins, norm=LogNorm())
        # cbar = plt.colorbar()
        # cbar.set_label('N(%s)' % cat_name, fontsize = 15)
    elif (log == True) :
        ax.loglog(xvals, yvals, 'ko', label='data')
    elif (logy == True) :
        ax.semilogy(xvals, yvals, 'ko', label='data')
    else :
        ax.plot(xvals, yvals, 'ko', label='data')
    
    # if fit == True :
    #     xs = np.linspace(7, 11.3, 1000)
    #     xs = np.linspace(6, 13, 1000)
    #     ax.plot(xs, xs+0, 'r--', label='equality')
    #     ax.plot(xs, 0.89981547*xs + 1.5640069, 'b--', label='linear fit')
    #     ax.plot(xs, xs + 0.65219703, 'r--', label='linear offset')
    #     ax.plot(xs, 0.13952656*(xs - 5.80306785)**2 + 8.16238619, 'r--', label='parabolic fit')
    
    if fit0 == True :
        start = min(np.min(xvals), np.min(yvals))
        stop = max(np.max(xvals), np.max(yvals))
        xs = np.linspace(start, stop, 1000)
        ax.plot(xs, xs, 'k--', label='equality')
    
    if fit1 == True :
        start = min(np.min(xvals), np.min(yvals))
        stop = max(np.max(xvals), np.max(yvals))
        xs = np.linspace(start, stop, 1000)
        ax.plot(xs, xs, 'k--', label='equality')
        # ax.plot(xs, xs-1.2722659777062275, 'k-', label='linear offset fit')
        ax.plot(xs, 0.8125759463914866*xs-4.961460348702913, 'k-', label=r'fit: $0.82x-4.96$')
    
    if fit2 == True :
        start = min(np.min(xvals), np.min(yvals))
        stop = max(np.max(xvals), np.max(yvals))
        xs = np.linspace(start, stop, 1000)
        ax.plot(xs, xs, 'k--', label='equality')
        # ax.plot(xs, xs-1.2179114940900249, 'k-', label='linear offset fit')
        ax.plot(xs, 0.8492931852404171*xs-4.343099872187622, 'k-', label=r'fit: $0.85x-4.34$')
    
    ax.set_xlabel('%s' % xlab, fontsize=21)
    ax.set_ylabel('%s' % ylab, fontsize=21)
        
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()
    
    # handles = [handles[1], handles[0], handles[2]]
    # labels = [labels[1], labels[0], labels[2]]
    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=17,
               loc = 'upper left')
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=16)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig('correlation_plots/diagnostic_SDSS_WHAN.pdf', overwrite=True)
    # plt.close()
    
    return

def plot3D(xvals, yvals, zvals, xlabel, ylabel, zlabel) :
    
    from mpl_toolkits.mplot3d import Axes3D
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.plot(xvals, yvals, zvals)
    
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_zlabel(zlabel, fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_cylinder() :
    
    from matplotlib import patches
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ell1 = patches.Ellipse((0.3,0.3), width=0.25, height=0.1, angle=-45, linestyle='-',
                          edgecolor='k', facecolor='None', zorder=2)
    ax.add_patch(ell1)
    
    ell2 = patches.Ellipse((0.7,0.7), width=0.25, height=0.1, angle=-45, linestyle='-',
                          edgecolor='k', facecolor='None', zorder=2)
    ax.add_patch(ell2)
    
    xs = np.linspace(0,0.959,1000)
    ax.plot(xs, xs, 'k--')
    # ax.arrow(0, 0, 0.95, 0.95, linestyle='--', color='k')
    ax.annotate('', xy=(0.96,0.96), xytext=(0.959,0.959), arrowprops={'arrowstyle': '->'})
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    # ax.legend()
    plt.tight_layout()
    plt.show()
    
    return

def plot_decoupled_fit(file, x_array, raw_data, fit, residual, HA_NII_core, BLR1) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(8.5, 6.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.plot(x_array, raw_data, drawstyle='steps-mid', linewidth=3, color='grey')
    ax.plot(x_array, fit, 'r-', label='Model')
    resid = (x_array>=6500) & (x_array <= 6650)
    ax.plot(x_array[resid], residual[resid], 'm-.', label='Residual')
    ax.plot(x_array, HA_NII_core, 'b-', label=r'H$\rm \alpha$+[N II] core')
    # ax.plot(x_array, HA_NII_wing, 'g-', label=r'H$\rm \alpha$+[N II] wing')
    ax.plot(x_array, BLR1, 'k--', label=r'H$\rm \alpha$ BLR1')
    # ax.plot(x_array, BLR2, 'k-.', label=r'H$\rm \alpha$ BLR2')
    # ax.plot(x_array, BLR3, 'k:', label=r'H$\rm \alpha$ BLR3')
    
    ax.set_title(r'%s' % (file) ) # H$\rm \alpha$ +[N II] Spectrum
    ax.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize = 15)
    ax.set_ylabel(r'Flux Density ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\rm \AA^{-1}$)',
               fontsize = 15)
    
    ax.set_xlim(6400, 6700)
    
    ax.legend()
    plt.tight_layout()
    plt.show()
    # plt.savefig(file, overwrite=True)
    # plt.close(fig)
    
    return

def plot_multi(xarrays, yarrays, labels, colours, xlab, ylab, location=0, outfile=None,
               xlog=False, ylog=False, log=False, alph=1,
               xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(10,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    for i in range(len(xarrays)) : # add each dataset in a list to the plot
        ax.plot(xarrays[i], yarrays[i], '-', color=colours[i],
                label = labels[i], alpha=alph)
    
    ax.set_xlabel('%s' % xlab, fontsize=21)
    ax.set_ylabel('%s' % ylab, fontsize=21)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    handles, labels = ax.get_legend_handles_labels()    
    ax.legend(handles, labels, facecolor='whitesmoke', framealpha=1, fontsize=14,
              loc = location)
    
    # ax.grid(b=True, zorder=1)
    ax.tick_params(axis='both',which='both',labelsize=16)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    plt.tight_layout()
    plt.show()
    # fig.savefig(outfile, overwrite=True)
    # plt.close(fig)
    
    return

def plot_with_errors(xvals, xlab, yvals, ylab, upper_errors, lower_errors,
                     outfile=None, xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(11.54, 8))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.plot(xvals, yvals, 'k-', label='Median')
    ax.plot(xvals, upper_errors, 'k:', label=r'$1\sigma$ Limits')
    ax.plot(xvals, lower_errors, 'k:')
    
    ax.hlines(0, 6, 14, color='grey', linestyle='--')
    
    ax.set_xlabel('%s' % xlab, fontsize=22)
    ax.set_ylabel('%s' % ylab, fontsize=22)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.tick_params(axis='both',which='both',labelsize=19)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=19, loc = 'upper left')
    plt.tight_layout()
    plt.show()
    # plt.savefig(outfile, overwrite=True)
    # plt.close(fig)
    
    return

def plot_with_errorbars(xvals, x_errors, xlab, yvals, y_errors, ylab,
                        outfile=None, xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(8,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.errorbar(xvals, yvals, xerr=x_errors, yerr=y_errors,
                fmt='ko', elinewidth=0.3, capsize=1.5, errorevery=1)
    
    ax.vlines(0.35, 0, 2, color='grey', linestyle='-', label='A=0.35')
    equality = np.linspace(0, 2, 1000)
    ax.plot(equality, equality, color='grey', linestyle='--', label='Equality')
    
    ax.annotate('Merger Region', xy=(1,1.75), size=15)
    
    ax.set_xlabel('%s' % xlab, fontsize=21)
    ax.set_ylabel('%s' % ylab, fontsize=21)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.tick_params(axis='both',which='both',labelsize=17)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=15, loc = 'upper left')
    plt.tight_layout()
    plt.show()
    # plt.savefig(outfile, overwrite=True)
    # plt.close(fig)
    
    return

def population_plot(y1, y1err, y2, y2err, y3, y3err, labels,
                    outfile=None, xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(8,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    xs = [1, 2, 3, 4]
    xlab = 'Environmental Density'
    ylab = 'Fraction'
    colours = ['red', 'dodgerblue', 'gold']
    
    # for GAMA AGN
    ax.errorbar(xs, y1, yerr=y1err, fmt='-', color=colours[0], label=labels[0],
                elinewidth=1, capsize=1.5, errorevery=1)
    # for GAMA SFG
    ax.errorbar(xs, y2, yerr=y2err, fmt='-', color=colours[1], label=labels[1],
                elinewidth=1, capsize=1.5, errorevery=1)
    # for GAMA Passive
    ax.errorbar(xs, y3, yerr=y3err, fmt='-', color=colours[2], label=labels[2],
                elinewidth=1, capsize=1.5, errorevery=1)
    
    ax.set_xlabel('%s' % xlab, fontsize=21)
    ax.set_ylabel('%s' % ylab, fontsize=21)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.xticks(xs, ['Low', '', '', 'High'])
    ax.tick_params(axis='both',which='both',labelsize=17)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=15, loc = 'upper left')
    plt.tight_layout()
    plt.show()
    # plt.savefig(outfile, overwrite=True)
    # plt.close(fig)
    
    return

def sample_projection(list_of_catalogs, labels, outfile=None) :
    
    list_of_RA = []
    list_of_dec = []
    for file in list_of_catalogs :
        catalog = Table.read(file)
        sub_ra = Angle(catalog['RA'], u.deg)
        sub_dec = Angle(catalog['DEC'], u.deg)
        sub_ra = sub_ra.wrap_at(180*u.deg)
        list_of_RA.append(sub_ra)
        list_of_dec.append(sub_dec)
    
    sky_projection(list_of_RA, list_of_dec, labels, outfile)
    
    return

def sky_projection(list_of_RA, list_of_dec, labels, outfile=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(12,6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111, projection='aitoff') # aitoff or mollweide
    
    marker_style = ['ko', 'bo', 'ro', 'go', 'mo', 'co', 'yo',
                    'k^', 'b^', 'r^', 'g^', 'm^', 'c^', 'y^']
    for i in range(len(list_of_RA)) : # add each dataset to the plot
        ax.plot(list_of_RA[i].radian, list_of_dec[i].radian,
                marker_style[i], label = labels[i], ms=7)
    
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h',
                        '8h','10h'])
    ax.grid(True)
    
    ax.tick_params(axis='both',which='both',labelsize=19)
    ax.tick_params(axis='both',which='major',direction='in',width=1.5,length=10)
    ax.tick_params(axis='both',which='minor',direction='in',width=1.5,length=6)
    ax.minorticks_on()
    
    ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=19, loc = 'upper right')
    plt.tight_layout()
    # plt.show()
    plt.savefig(outfile, overwrite=True)
    plt.close(fig)
    
    return

def view_3d_plot(file) :
    
    npzfile = np.load(file)
    
    plt.contour3d(npzfile['x'], npzfile['y'], npzfile['a'],
                    'Velocity Difference (km/s)', 'Radius (Mpc)',
                    'Number of Companion Galaxies')
    
    plt.contour3d(npzfile['x'], npzfile['y'], npzfile['z'],
                    'Velocity Difference (km/s)', 'Radius (Mpc)',
                    'Projected CoM Distance to AGN (arcmin)')
    
    return

def view_all_3d_plots() :
    
    # using the GAMA catalogue
    view_3d_plot('GAMA_HE 0853+0102_10_10.npz')
    view_3d_plot('GAMA_HE 0853+0102_20_20.npz')

    # using the SDSS catalogue
    view_3d_plot('SDSS_HE 0040-1105_10_10.npz')
    view_3d_plot('SDSS_HE 0040-1105_20_20.npz')
    
    view_3d_plot('SDSS_HE 0114-0015_10_10.npz')
    view_3d_plot('SDSS_HE 0114-0015_20_20.npz')
    
    view_3d_plot('SDSS_HE 0119-0118_10_10.npz')
    view_3d_plot('SDSS_HE 0119-0118_20_20.npz')
    
    view_3d_plot('SDSS_HE 0203-0031_10_10.npz')
    view_3d_plot('SDSS_HE 0203-0031_20_20.npz')
    
    view_3d_plot('SDSS_HE 0212-0059_10_10.npz')
    view_3d_plot('SDSS_HE 0212-0059_20_20.npz')
    
    view_3d_plot('SDSS_HE 0227-0913_10_10.npz')
    view_3d_plot('SDSS_HE 0227-0913_20_20.npz')
    
    view_3d_plot('SDSS_HE 0232-0900_10_10.npz')
    view_3d_plot('SDSS_HE 0232-0900_20_20.npz')
    
    view_3d_plot('SDSS_HE 0345+0056_10_10.npz')
    view_3d_plot('SDSS_HE 0345+0056_20_20.npz')
    
    view_3d_plot('SDSS_HE 0853+0102_10_10.npz')
    view_3d_plot('SDSS_HE 0853+0102_20_20.npz')
    
    view_3d_plot('SDSS_HE 0934+0119_10_10.npz')
    view_3d_plot('SDSS_HE 0934+0119_20_20.npz')
    
    view_3d_plot('SDSS_HE 2222-0026_10_10.npz')
    view_3d_plot('SDSS_HE 2222-0026_20_20.npz')
    
    view_3d_plot('SDSS_HE 2302-0857_10_10.npz')
    view_3d_plot('SDSS_HE 2302-0857_20_20.npz')
    
    return

# env_hists('Companions', log=True)

# sample_projection(['catalogs/raw_cats/HES_sample.fits',
#                     'catalogs/raw_cats/CARS_sample.fits',
#                     'catalogs/CARS_SDSS/CARS_SDSS_base.fits'],
#                   ['HES','CARS','Thesis Sample'],
#                   'correlation_plots/HES_CARS_thesis_samples.pdf')

# sample_projection([#'catalogs/CARS_SDSS/CARS-compliments_base.fits',
                   # 'catalogs/CARS_SDSS/SDSS-galaxies_base_random_2500.fits',
                   # 'catalogs/CARS_GAMA/GAMA_BLAGN_base_0-52.fits',
                   # 'catalogs/CARS_GAMA/GAMA_Comp_base_0-159.fits',
                   # 'catalogs/CARS_GAMA/GAMA_LINER_base_0-52.fits',
                   # 'catalogs/CARS_GAMA/GAMA_not-ELG_base_0-2194.fits',
                   # 'catalogs/CARS_GAMA/GAMA_Passive_base_0-194.fits',
                   # 'catalogs/CARS_GAMA/GAMA_Seyfert_base_0-39.fits',
                   # 'catalogs/CARS_GAMA/GAMA_SFG_base_0-2726.fits'
                   # ],
                  # [#'SDSS - QSOs', 'SDSS - galaxies',
                    # 'GAMA - BLAGN',
                    # 'GAMA - Comp.',
                    # 'GAMA - LINER',
                    # 'GAMA - not_ELG',
                    # 'GAMA - Passive',
                    # 'GAMA - Seyfert',
                    # 'GAMA - SFG'
                   # ])

# to include as-is in thesis
# all_histograms('Z', 'Redshift', loc=2, xmin=0.01, xmax=0.062, ymin=0, ymax=0.44)
# all_histograms('M_g', r'Absolute $g$ Magnitude', xmin=-23, xmax=-17, ymin=0, ymax=0.45)
# all_histograms('M_i', r'Absolute $i$ Magnitude', loc=2, xmin=-24.6, xmax=-18, ymin=0, ymax=0.48)
# all_histograms('logMass', r'$\log($Stellar Mass of Target$/M_\odot)$', xmin=9.9, xmax=12.4, ymin=0, ymax=0.5)
# all_histograms('Companions', 'Number of Companions', log=True, xmin=-0.1, xmax=2.7, ymin=0, ymax=0.38)
# all_histograms('Most_Massive_Mass', r'$\log($Stellar Mass of MMG in LSE$/M_\odot)$',
#                 loc=2, xmin=9, xmax=12.4, ymin=0, ymax=0.4)
# all_histograms('Most_Massive_Distance', 'Distance to MMG in LSE/kpc',
#                 log=True, loc=2, sdss_bins=22, xmin=0, xmax=4.5, ymin=0, ymax=0.35)
# all_histograms('Closest_Mass', r'$\log($Stellar Mass of Closest Companion$/M_\odot)$',
#                 xmin=9, xmax=12, ymin=0, ymax=0.28)
# all_histograms('Closest_Distance', 'Distance to Closest Companion/kpc',
#                 log=True, loc=2, sdss_bins=22, xmin=0, xmax=4.4, ymin=0, ymax=0.5)
# all_histograms('Massive_1Mpc_Mass', r'$\log($Stellar Mass of MMC in 1 Mpc$/M_\odot)$',
#                 xmin=9, xmax=12, ymin=0, ymax=0.2)
# all_histograms('Massive_1Mpc_Distance', 'Distance to MMC in 1 Mpc/kpc',
#                 log=True, loc=2, sdss_bins=42, xmin=1.5, xmax=3.1, ymin=0, ymax=0.33)
# all_histograms('d_CoM', 'Projected Distance to the Center of Stellar Mass/arcmin',
#                 log=True, loc=2, xmin=-0.05, xmax=1.85, ymin=0, ymax=0.3)
# all_histograms('CountInCyl', 'Companions in the GAMA Cyl.', log=True, xmin=-0.1, xmax=2.25, ymin=0, ymax=0.43)
# all_histograms('SurfaceDensity', r'Surface Density/Mpc$^{-2}$',
#                 log=True, sdss_bins=18, xmin=-4.8, xmax=1, ymin=0, ymax=0.4)
# all_histograms('Overdensity', 'Overdensity', log=True, xmin=-0.9, xmax=1.45, ymin=0, ymax=0.43)
# all_histograms('AGEPar', r'AGE Parameter/Mpc$^{-1}$', # must change number of bins only for GAMA passives
#                 log=True, loc=2, sdss_bins=18, xmin=-3, xmax=2, ymin=0, ymax=0.48)

# alternatives (ie. not using a log-scale, for example)
# all_histograms('Companions', 'Number of Companions', xmin=-5, xmax=440, ymin=0, ymax=0.9)
# all_histograms('Most_Massive_Distance', 'Distance to Most Massive Companion (kpc)',
#                 xmin=-500, xmax=22500, ymin=0, ymax=0.6)
# all_histograms('Closest_Distance', 'Distance to Closest Companion (kpc)',
#                 xmin=-500, xmax=20000, ymin=0, ymax=0.95)
# all_histograms('d_CoM', 'Projected Distance to the Center of Stellar Mass (arcmin)',
#                 xmin=-2, xmax=65, ymin=0, ymax=0.35)
# all_histograms('CountInCyl', 'Companions in the GAMA Cyl.', xmin=-2, xmax=175, ymin=0, ymax=0.9)
# all_histograms('AGEPar', r'AGE Parameter (Mpc$^{-1}$)', xmin=-1, xmax=45, ymin=0, ymax=0.83)

# scatter plots showing lack of correlation between different env. parameters
# all_scatter_plots('GAMA', 'Closest_Distance', 'SurfaceDensity', loglog=True)
# all_scatter_plots('GAMA', 'Closest_Distance', 'Companions', loglog=True)
# all_scatter_plots('GAMA', 'Closest_Distance', 'Most_Massive_Mass', logx=True)
# all_scatter_plots('GAMA', 'Closest_Distance', 'Overdensity', loglog=True)
# all_scatter_plots('GAMA', 'Closest_Distance', 'AGEPar', loglog=True)

# all_scatter_plots('GAMA', 'Most_Massive_Mass', 'AGEPar', logy=True)

# scatter plots showing lack of correlation between the env. and AGN properties
# all_scatter_plots('GAMA', 'Z', 'Closest_Distance', xmin=0.009, xmax=0.061, ymin=30, ymax=2e4, logy=True)
# all_scatter_plots('GAMA', 'logMass', 'Closest_Distance', xmin=9.95, xmax=12, ymin=30, ymax=2e4, logy=True)
# all_scatter_plots('GAMA', 'SFR', 'Closest_Distance', loc=2, xmin=0.0006, xmax=100, ymin=30, ymax=2e4, loglog=True)
# all_scatter_plots('GAMA', 'OIII_FWHM', 'Closest_Distance', loc=4, xmin=100, xmax=800, ymin=30, ymax=2e4, logy=True)

# all_scatter_plots('GAMA', 'Z', 'SurfaceDensity', xmin=0.009, xmax=0.061, ymin=1e-5, ymax=10, logy=True)
# all_scatter_plots('GAMA', 'logMass', 'SurfaceDensity', xmin=9.95, xmax=12, ymin=1e-5, ymax=10, logy=True)
# all_scatter_plots('GAMA', 'SFR', 'SurfaceDensity', loc=2, xmin=0.0006, xmax=100, ymin=1e-5, ymax=10, loglog=True)
# all_scatter_plots('GAMA', 'OIII_FWHM', 'SurfaceDensity', xmin=100, xmax=800, ymin=1e-5, ymax=10, logy=True)

# all_scatter_plots('SDSS', 'Z', 'Closest_Distance', xmin=0.009, xmax=0.061, ymin=30, ymax=2e4, logy=True)
# all_scatter_plots('SDSS', 'logMass', 'Closest_Distance', xmin=9.95, xmax=12, ymin=30, ymax=2e4, logy=True)
# all_scatter_plots('SDSS', 'SFR', 'Closest_Distance', xmin=0.0006, xmax=100, ymin=30, ymax=2e4, loglog=True)
# all_scatter_plots('SDSS', 'OIII_FWHM', 'Closest_Distance', loc=4, xmin=100, xmax=800, ymin=30, ymax=2e4, logy=True)

# all_scatter_plots('SDSS', 'Z', 'SurfaceDensity', xmin=0.009, xmax=0.061, ymin=1e-5, ymax=10, logy=True)
# all_scatter_plots('SDSS', 'logMass', 'SurfaceDensity', xmin=9.95, xmax=12, ymin=1e-5, ymax=10, logy=True)
# all_scatter_plots('SDSS', 'SFR', 'SurfaceDensity', loc=2, xmin=0.0006, xmax=100, ymin=1e-5, ymax=10, loglog=True)
# all_scatter_plots('SDSS', 'OIII_FWHM', 'SurfaceDensity', xmin=100, xmax=800, ymin=1e-5, ymax=10, logy=True)

# unused
# all_scatter_plots('GAMA', 'Z', 'OIII_FWHM', xmin=0.009, xmax=0.061, ymin=100, ymax=800)
# all_scatter_plots('GAMA', 'logMass', 'OIII_FWHM', xmin=9.95, xmax=12, ymin=100, ymax=800)
# all_scatter_plots('GAMA', 'SFR', 'OIII_FWHM', xmin=1e-4, xmax=100, ymin=100, ymax=800, logx=True)
# all_scatter_plots('SDSS', 'Z', 'OIII_FWHM', xmin=0.009, xmax=0.061, ymin=100, ymax=800)
# all_scatter_plots('SDSS', 'logMass', 'OIII_FWHM', xmin=9.95, xmax=11.5, ymin=100, ymax=800)
# all_scatter_plots('SDSS', 'SFR', 'OIII_FWHM', xmin=1e-3, xmax=100, ymin=100, ymax=800, logx=True)

# other possibly interesting plots?
# all_scatter_plots('GAMA', 'Z', 'SFR', xmin=0.009, xmax=0.061, ymin=1e-4, ymax=100, logy=True)
# all_scatter_plots('SDSS', 'Z', 'SFR', xmin=0.009, xmax=0.061, ymin=1e-3, ymax=100, logy=True)

# random plots - do we need to even include these?
# all_scatter_plots('GAMA', 'Most_Massive_Mass', 'd_CoM')
# all_scatter_plots('GAMA', 'Companions', 'Most_Massive_Distance', loglog=True)
# all_scatter_plots('GAMA', 'Overdensity', 'Closest_Mass', xmin=0.1, xmax=100, logx=True)
# all_scatter_plots('GAMA', 'Closest_Mass', 'Most_Massive_Mass')

# gaussian_plot(np.array([0, 0, 0, -2]), np.array([0.2, 1, 5, 0.5]),
#               'correlation_plots/gaussian_shape.pdf')
# filter_curves()
