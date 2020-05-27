
# imports
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import Angle
from astropy.table import Table
import astropy.units as u

# constants
currentFig = 1

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
    
    plt.xlabel(xlabel, fontsize = 15 )
    plt.ylabel(ylabel, fontsize = 15 )
    ax.set_zlabel(zlabel, fontsize = 15)
    
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
    fig = plt.figure(currentFig, figsize=(8.5, 6.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    # for the SDSS sample check
    ax.hist2d(x1, y1, bins=341, cmap='Blues', norm=LogNorm())
    ax.hist2d(x2, y2, bins=156, cmap='Purples', norm=LogNorm())
    ax.hist2d(x3, y3, bins=200, cmap='Reds', norm=LogNorm())
    ax.hist2d(x4, y4, bins=117, cmap='Greens', norm=LogNorm())
    
    # for the GAMA sample check
    # ax.hist2d(x1, y1, bins=100, cmap='Blues', norm=LogNorm())
    # ax.hist2d(x2, y2, bins=30, cmap='Purples', norm=LogNorm())
    # ax.hist2d(x3, y3, bins=10, cmap='Reds', norm=LogNorm())
    # ax.hist2d(x4, y4, bins=10, cmap='Greens', norm=LogNorm())
    
    unclass_x = [0.52699786, 0.6472921, 0.51526904, 0.10082114,
               0.4783239, 0.64564115, 0.47279504, 0.48915058,
               0.5316156, 0.52157855, 0.072221704, 0.6222592, 0.061153412]
    unclass_y = [0.20510054, 0.15230991, 0.3134719, -1.4588171,
               0.8099388, 0.41682592, 0.12564416, 0.5535542,
               0.8346556, 0.4626364, -0.5899385, 0.4983, -0.36166635]
    ax.plot(unclass_x, unclass_y, 'kx', label='BPT Classification Fails (13 galaxies)')
    
    CARS_x_7p5 = [-0.2558046120106365, -0.3039520129538409, -0.2417594381357093,
                  0.06476555533411806, 0.275461262291009, -0.3938455204465459,
                  -0.05406315470653154, -1.0767286502681763, -0.21683621439121184,
                  -0.2583153650781086, -0.3611262953855419, -0.041477366726339726]
    CARS_y_7p5 = [1.0442984352563223, -0.4632845953092678, -0.17709178145221774,
                  0.7852122081665016, 0.701362126702604, -0.3101074263261445,
                  0.9213412646212272, 0.07531550781816444, 0.43745460649821044,
                  0.6040263849321078, 0.02695921366978017, 0.8956928513689375]
    ax.plot(CARS_x_7p5, CARS_y_7p5, 'ko', label='Fluxes from .eline MUSE Tables in 3" Aperture')
    
    if (standard==True) :
        xmin, xmax, ymin, ymax = -2.2, 0.7, -1.5, 1.4
        xlabel = r"$\log_{10}$([N II] $\lambda 6583 / \rm H\alpha$)"
        ylabel = r"$\log_{10}$([O III] $\lambda 5007 / \rm H\beta$)"
        
        kauff_x = np.linspace(xmin, xmax, 10000)
        kauff_y = 0.61/(kauff_x - 0.05) + 1.3 # from Kauffmann+ 2003
        #                       (http://adsabs.harvard.edu/abs/2003MNRAS.346.1055K)
        kauff_mask = ( (kauff_x > xmin) & (kauff_x < 0.05) )
        ax.plot(kauff_x[kauff_mask], kauff_y[kauff_mask], 'k--',
                label='Kauffmann et al. (2003)')
        
        kewl_x = np.linspace(xmin, xmax, 10000)
        kewl_y = 0.61/(kewl_x - 0.47) + 1.19 # from Kewley+ 2001
        #                       (http://adsabs.harvard.edu/abs/2001ApJ...556..121K)
        kewl_mask = ( (kewl_x > xmin) & (kewl_x < 0.47) )
        ax.plot(kewl_x[kewl_mask], kewl_y[kewl_mask], 'k-',
                label='Kewley et al. (2001)')
        
        schaw_x = np.linspace(-0.18380687748557267, xmax, 10000)
        schaw_y = 1.05*schaw_x + 0.45 # from Schawinski+ 2012
        #                       (http://adsabs.harvard.edu/abs/2007MNRAS.382.1415S)
        ax.plot(schaw_x, schaw_y, 'k-.', label='Schawinski et al. (2007)')
        
        ax.text(-0.8, 1.1, 'Seyfert', fontsize=15)
        ax.text(0.2, -0.3, 'LINER', fontsize=15)
        ax.text(-1.4, -0.4, 'SFG', fontsize=15)
        ax.text(-0.225, -1, 'Comp', fontsize=15)
        
        # ax.vlines(0.47, ymin, ymax, color='b', linestyle='-.')
        
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
        
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc = 'lower left')
    
    plt.tight_layout()
    plt.show()
    
    return

def diagram_WHAN(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, standard=True,
                 xmin=None, xmax=None, ymin=None, ymax=None,
                 xlabel=None, ylabel=None) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(8.5, 6.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)  
    
    # for the SDSS sample check
    ax.hist2d(x1, y1, bins=80, cmap='Oranges', norm=LogNorm())
    ax.hist2d(x2, y2, bins=86, cmap='Blues', norm=LogNorm())
    ax.hist2d(x3, y3, bins=30, cmap='Purples', norm=LogNorm())
    ax.hist2d(x4, y4, bins=40, cmap='Reds', norm=LogNorm())
    ax.hist2d(x5, y5, bins=113, cmap='Greens', norm=LogNorm())
    
    # for the GAMA sample check
    # ax.hist2d(x1, y1, bins=50, cmap='Oranges', norm=LogNorm())
    # ax.hist2d(x2, y2, bins=100, cmap='Blues', norm=LogNorm())
    # ax.hist2d(x3, y3, bins=30, cmap='Purples', norm=LogNorm())
    # ax.hist2d(x4, y4, bins=20, cmap='Reds', norm=LogNorm())
    # ax.hist2d(x5, y5, bins=20, cmap='Greens', norm=LogNorm())
    
    if (standard==True) :
        xmin, xmax, ymin, ymax = -2.2, 1, -0.7, 3.3
        xlabel = r"$\log_{10}$([N II] $\lambda 6583 / \rm H\alpha$)"
        ylabel = r'$\log_{10}(W_{\rm H\alpha}/ {\rm \AA})$'
        
        # regular WHAN
        ax.vlines(-0.4, ymin, ymax, color='k', linestyle='-.',
                  label=r'Stasi$\rm \'n$kska et al. (2006)$^{\rm T}$')
        ax.vlines(-0.32, ymin, ymax, color='k', linestyle='-',
                  label=r'Gordon et al. (2018)')
        ax.hlines(np.log10(6), -0.32, xmax, color='k', linestyle='--',
                  label=r'Kewley et al. (2006)$^{\rm T}$')
        ax.hlines(np.log10(0.5), xmin, xmax, color='k', linestyle=':')
        
        xs = np.linspace(-3, 1.5, 1000)
        ys = -xs + np.log10(0.5)
        ax.plot(xs, ys, 'k:')
        
        ax.text(-1, 2.75, 'SFG', fontsize=15)
        ax.text(0, 2.5, 'Seyfert', fontsize=15)
        ax.text(0.55, 0.5, 'LINER', fontsize=15)
        ax.text(-0.2, -0.6, 'Passive', fontsize=15)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc = 'lower left')
    
    plt.tight_layout()
    plt.show()
    
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

def histo(param, label, title=None, vert_line=None, num_bins=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    
    ax = fig.add_subplot(111)
    ax.hist(param, range=(7, 11.5)) #, bins=num_bins)
    ax.set_xlabel("%s" % label, fontsize = 15)
    
    if vert_line :
        ymin, ymax = ax.get_ylim()
        ax.vlines(vert_line, ymin, ymax,
                  color='r', linestyle='--', label=r'CARS Host')
        plt.legend()
    
    plt.title('%s' % title, fontsize = 15)
    plt.tight_layout()
    plt.show()
    
    return

def histo2d(xvals, xlab, yvals, ylab, nbins=200, xmin=None, xmax=None,
            ymin=None, ymax=None,) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    
    ax = fig.add_subplot(111)
    
    ax.hist2d(xvals, yvals, bins=nbins, norm=LogNorm())
    cbar = plt.colorbar()
    
    ax.set_xlabel("%s" % xlab, fontsize = 15 )
    ax.set_ylabel("%s" % ylab, fontsize = 15 )
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return

def multi_histo(param_list, labels, xlab, log=False, xmin=None, xmax=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    
    ax = fig.add_subplot(111)    
    
    colors = ['k', 'b', 'r', 'g', 'm', 'c', 'y', 'k', 'b']
    styles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-']
    
    if log==True :
        new_list = []
        for i in range(len(param_list)) :
            mask = param_list[i] > 0
            new_list.append( np.log10(param_list[i][mask]) )
        ax.hist(new_list, histtype='bar', density=True, label=labels,
                linestyle='--') # stacked=False
        xlabel = 'log(' + xlab + ')'
    else :
        ax.hist(param_list, histtype='step', density=True, label=labels,
                linestyle='--')
        xlabel = xlab
    
    """
    if log==True :
        for i in range(len(param_list)) :
            mask = param_list[i] > 0
            ax.hist(np.log10(param_list[i][mask]), label = labels[i],
                    density=True,
                    histtype='step', linestyle=styles[i], color=colors[i],
                    stacked=True)
        xlabel = 'log(' + xlab + ')'
    else :
        for i in range(len(param_list)) : # add each dataset to the plot
            ax.hist(param_list[i], label = labels[i], density=True,
                    histtype='step', linestyle=styles[i], color=colors[i],
                    stacked=True)
        xlabel = xlab
    """
    
    ax.set_xlabel("%s" % xlabel, fontsize = 15 )
    ax.set_xlim(xmin, xmax)
    
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return

def multi(xvals, xlab, yvals1, yvals2, ylab, xmin=None, xmax=None, ymin=None,
          ymax=None, location='upper left') :
    
    global currentFig        
    fig = plt.figure(currentFig)  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
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
    
    ax.set_xlabel("%s" % xlab, fontsize = 15 )
    ax.set_ylabel("%s" % ylab, fontsize = 15 )
    
    plt.legend(loc = location)
    
    plt.tight_layout()
    plt.show()
    
    return

def multi2(xvals1, yvals1, label1, xvals2, yvals2, label2, xlab, ylab) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig        
    fig = plt.figure(currentFig) # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111)
    
    ax.hist2d(xvals1, yvals1, bins=50, norm=LogNorm())
    cbar = plt.colorbar()
    cbar.set_label('GAMA', fontsize = 15)
    
    ax.plot(xvals2, yvals2, 'ro', label = "%s" % label2, zorder=2)
    
    ax.set_xlim(-4, 3)
    ax.set_ylim(-1, 27)
    
    ax.set_xlabel("%s" % xlab, fontsize = 15)
    ax.set_ylabel("%s" % ylab, fontsize = 15)
    
    plt.legend()
    
    plt.tight_layout()
    plt.show()
    
    return

def plot(xvals, xlab, yvals, ylab, cat_name='SDSS', xmin=None, xmax=None,
         ymin=None, ymax=None, hist2d=False, nbins=200, fit=False, log=False) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig
    fig = plt.figure(currentFig)  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111) # set axes, figure location
    
    if (hist2d == True) :
        plt.hist2d(xvals, yvals, bins=nbins)#, norm=LogNorm())
        cbar = plt.colorbar()
        cbar.set_label('N(%s)' % cat_name, fontsize = 15)
    elif (log == True) :
        ax.loglog(xvals, yvals, 'ko', label='data')
    else :
        ax.plot(xvals, yvals, 'ko', label='data')
    
    if fit == True :
        xs = np.linspace(7, 11.3, 1000)
    #    xs = np.linspace(6, 13, 1000)
        ax.plot(xs, xs+0, 'r--', label='equality')
    #    ax.plot(xs, 0.89981547*xs + 1.5640069, 'b--', label='linear fit')
    #    ax.plot(xs, xs + 0.65219703, 'r--', label='linear offset')
    #    ax.plot(xs, 0.13952656*(xs - 5.80306785)**2 + 8.16238619, 'r--', label='parabolic fit')
    
    ax.set_xlabel("%s" % xlab, fontsize = 15 )
    ax.set_ylabel("%s" % ylab, fontsize = 15 )
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.legend()
    plt.tight_layout()
    plt.show()
    
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
    
    plt.xlabel(xlabel, fontsize = 15 )
    plt.ylabel(ylabel, fontsize = 15 )
    ax.set_zlabel(zlabel, fontsize = 15)
    
    plt.show()
    
    return

def sample_projection(list_of_catalogs, labels) :
    
    list_of_RA = []
    list_of_dec = []
    for file in list_of_catalogs :
        catalog = Table.read(file)
        sub_ra = Angle(catalog['RA'], u.deg)
        sub_dec = Angle(catalog['DEC'], u.deg)
        sub_ra = sub_ra.wrap_at(180*u.deg)
        list_of_RA.append(sub_ra)
        list_of_dec.append(sub_dec)
    
    plt.sky_projection(list_of_RA, list_of_dec, labels)
    
    return

def sky_projection(list_of_RA, list_of_dec, labels) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(16,8))
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111, projection="aitoff") # aitoff or mollweide
    
    marker_style = ['ko', 'bo', 'ro', 'go', 'mo', 'co', 'yo',
                    'k^', 'b^', 'r^', 'g^', 'm^', 'c^', 'y^']
    for i in range(len(list_of_RA)) : # add each dataset to the plot
        ax.plot(list_of_RA[i].radian, list_of_dec[i].radian,
                marker_style[i], label = labels[i])
    
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h',
                        '8h','10h'])
    ax.grid(True)
    
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.show()
    
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
#                    'catalogs/raw_cats/CARS_sample.fits',
#                    'catalogs/CARS_SDSS/CARS_SDSS_base.fits'],
#                   ['HES','CARS','Thesis Sample'])

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