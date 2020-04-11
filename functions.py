
import numpy as np
import matplotlib.pyplot as plt

#import environments_part1 as env1 # will break envrionments1.py

currentFig = 1

#.....................................................................contour3d
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

#......................................................................plotting
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

#..................................................................density_plot
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

#...................................................................helper_fxns
def intercept(xx, b) :
    return xx + b

def line(xx, m, b) :
    return m*xx + b

def parabola(xx, amp, translation, offset) :
    return amp*(xx - translation)**2 + offset

#.........................................................................histo
def histo(param, label, title=None, vert_line=None, num_bins=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    
    ax = fig.add_subplot(111)
    ax.hist(param, range=(7, 11.5)) #, bins=num_bins)
    plt.xlabel("%s" % label, fontsize = 15)
    
    if vert_line :
        ymin, ymax = ax.get_ylim()
        ax.vlines(vert_line, ymin, ymax,
                  color='r', linestyle='--', label=r'CARS Host')
        plt.legend()
    
    plt.title('%s' % title, fontsize = 15)
    plt.tight_layout()
    plt.show()
    
    return

#..........................................................................plot
def histo2d(xvals, xlab, yvals, ylab, nbins=200, xmin=None, xmax=None,
            ymin=None, ymax=None,) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    
    ax = fig.add_subplot(111)
    
    plt.hist2d(xvals, yvals, bins=nbins, norm=LogNorm())
    cbar = plt.colorbar()
    
    ax.set_xlabel("%s" % xlab, fontsize = 15 )
    ax.set_ylabel("%s" % ylab, fontsize = 15 )
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return

#.........................................................................multi
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

#.........................................................................multi
def multi2(xvals1, yvals1, label1, xvals2, yvals2, label2, xlab, ylab) :
    
    from matplotlib.colors import LogNorm
    
    global currentFig        
    fig = plt.figure(currentFig)  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111)
    
    plt.hist2d(xvals1, yvals1, bins=50, norm=LogNorm())
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

#..........................................................................plot
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

#.....................................................................mollweide
def full_sky(ra1, dec1, label1, ra2, dec2, label2, ra3, dec3, label3) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(16,8))  # the current figure
    currentFig += 1
    plt.clf() # clear the figure before each run
    
    ax = fig.add_subplot(111, projection="aitoff") # aitoff or mollweide
    
    ax.plot(ra1, dec1, 'ko', label = label1)
    ax.plot(ra2, dec2, 'bo', label = label2)
    ax.plot(ra3, dec3, 'ro', label = label3)
    
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h',
                        '8h','10h'])
    ax.grid(True)
    
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.show()
    
    return

#..............................................................end of functions
