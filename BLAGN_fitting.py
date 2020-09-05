
# imports
import numpy as np
import os

import astropy.constants as const
import astropy.io.fits as fits
from astropy.table import Table
import astropy.units as u
from scipy.optimize import leastsq, least_squares
from scipy.ndimage import gaussian_filter1d

import plots as plt

# constants
currentFig = 1
# directory = 'BLAGN_candidate_spectra/'
# out_directory = 'BLAGN_spectra_plots/'
directory = 'CARS_DR7_spectra/'
# out_directory = 'CARS_DR7_spectra_plots/'
fwhm_conv = 2*np.sqrt( 2*np.log(2) )
speed_of_light = const.c.to('km/s').value

def main(method='new') :
    
    files = os.listdir(directory)
    # files = ['spSpec-51789-0398-010.fit']
    # files = ['spSpec-54484-2656-469.fit', 'spSpec-53816-2219-071.fit',
             # 'spSpec-54526-2884-160.fit']
    
    # Type II spectrum
    # files = ['Type2/spSpec-51885-0434-437.fit']
    
    BL_FWHM, plates, mjds, fiberids = [], [], [], []
    for file in files :
        with fits.open(directory + file) as hdu :
            header = hdu[0].header        
            wavestart = header['COEFF0']
            waveint = header['COEFF1']
            zz = header['Z']
            data = hdu[0].data
            # print(header)
        
        spec = data[0] # spectra
        cs_spec = data[1] # continuum subtracted spectra
        # smooth = gaussian_filter1d(cs_spec, 1)
        # cs_spec = smooth
        err = data[2]
        kk = 1 + zz        
        
        wave = np.power(10, wavestart + waveint*np.arange( len(spec) ))
        wave = wave/kk
        
        HA_window = (wave >= 6400) & (wave <= 6700)
        HA_SN_window = (wave >= 6550) & (wave <= 6575)
        HA_SN = np.sum(cs_spec[HA_SN_window] / err[HA_SN_window])
        
        x0_guesses = np.array([20, 20, 0, 100, # 0.1, 0.1, 0, 100,
                               10, 0, 1000, # 0.1, 0, 1000, # 0.1, 0, 1000,
                               -0.01, -0.01])
        
        plt.plot(wave, r'Wavelength ($\rm \AA$)', cs_spec,
                  r'Flux Density ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\rm \AA^{-1}$)')
        
        if HA_SN >= 3e8 :
            if method == 'old' :
                try :
                    popt, pcov = leastsq(complex_gauss_decoupled, x0_guesses,
                                         args = (wave[HA_window], cs_spec[HA_window], err[HA_window]),
                                         maxfev = 1000)
                except ValueError :
                    popt = np.full(9, np.nan)
            if method == 'new' :
                try :
                    result = least_squares(complex_gauss_decoupled, x0_guesses,
                                           bounds=([0, 0, -np.inf, -np.inf, 0, -np.inf,
                                                    -np.inf, -np.inf, -np.inf],
                                                   [np.inf, np.inf, np.inf, 1000, np.inf,
                                                    np.inf, 10000, np.inf, np.inf]),
                                           args = (wave[HA_window], cs_spec[HA_window], err[HA_window]),
                                           max_nfev = 1000)
                    popt = result['x']
                except ValueError :
                    popt = np.full(9, np.nan)
            
            (amp_Ha_core, amp_NII6583_core, vel_Ha_core, vel_sigma_Ha_core,
             # amp_Ha_wing, amp_NII6583_wing, vel_Ha_wing, vel_sigma_Ha_wing,
             amp_Ha_blr1, vel_blr1, vel_sigma_blr1,
             # amp_Ha_blr2, vel_blr2, vel_sigma_blr2,
             # amp_Ha_blr3, vel_blr3, vel_sigma_blr3,
             slope, intercept) = popt
            
            spSpec, MJD, PLATEID, FIBERID = file.strip('.fit').split('-')
            PLATEID = PLATEID.lstrip('0')
            FIBERID = FIBERID.lstrip('0')
            plates.append(PLATEID)
            mjds.append(MJD)
            fiberids.append(FIBERID)
            BL_FWHM.append( fwhm_conv*vel_sigma_blr1 )
            print( fwhm_conv*vel_sigma_blr1*u.km/u.s )
            k_act = 1 + (vel_Ha_core/speed_of_light)
            select = (wave >= 6400*k_act) & (wave <= 6700*k_act)
            x_array = wave[select]/k_act
            
            fit = complex_gauss_decoupled(popt, wave[select], cs_spec[select],
                                          err[select])*(err[select]) + cs_spec[select]
            residual = cs_spec[select] - fit
            
            continuum = (wave[select]/1000.0)*slope + intercept
            
            HA_NII_core = ( continuum/k_act +
                           gauss(wave[select], amp_Ha_core, vel_Ha_core,
                                 vel_sigma_Ha_core, 6562.85, 2.5) +
                           (1/3)*gauss(wave[select], amp_NII6583_core, vel_Ha_core,
                                       vel_sigma_Ha_core, 6548.05, 2.5) +
                           gauss(wave[select], amp_NII6583_core, vel_Ha_core,
                                 vel_sigma_Ha_core, 6583.45, 2.5) )
            # HA_NII_wing = ( continuum/k_act +
            #                gauss(wave[select], amp_Ha_wing, vel_Ha_wing,
            #                      vel_sigma_Ha_wing, 6562.85, 2.5) +
            #                (1/3)*gauss(wave[select], amp_NII6583_wing, vel_Ha_wing,
            #                            vel_sigma_Ha_wing, 6548.05, 2.5) +
            #                gauss(wave[select], amp_NII6583_wing, vel_Ha_wing,
            #                      vel_sigma_Ha_wing, 6583.45, 2.5) )
            
            BLR1 = continuum/k_act + gauss(wave[select], amp_Ha_blr1, vel_blr1,
                                           vel_sigma_blr1, 6562.85, 2.5)
            # BLR2 = continuum/k_act + gauss(wave[select], amp_Ha_blr2, vel_blr2,
            #                                vel_sigma_blr2, 6562.85, 2.5)
            # BLR3 = continuum/k_act + gauss(wave[select], amp_Ha_blr3, vel_blr3,
            #                                vel_sigma_blr3, 6562.85, 2.5)
            
            # plot of the Ha + [N II] region
            file_name = (out_directory + str(PLATEID) + '_' + str(MJD) + '_' +
                         str(FIBERID) + '_fit.png')
            plt.plot_decoupled_fit(file_name, x_array, cs_spec[select], fit,
                                   residual, HA_NII_core, BLR1)
            
    output = Table([plates, mjds, fiberids, BL_FWHM*u.km/u.s],
                    names=('PLATEID', 'MJD', 'FIBERID', 'BL_FWHM'))
    # output.write(out_directory + 'BL_FWHM.fits', overwrite=True)
    
    return

def complex_gauss_decoupled(input_vals, wave, data, error) :
    
    (amp_Ha_core, amp_NII6583_core, vel_Ha_core, vel_sigma_Ha_core,
     # amp_Ha_wing, amp_NII6583_wing, vel_Ha_wing, vel_sigma_Ha_wing,
     amp_Ha_blr1, vel_blr1, vel_sigma_blr1,
     # amp_Ha_blr2, vel_blr2, vel_sigma_blr2,
     # amp_Ha_blr3, vel_blr3, vel_sigma_blr3,
     slope, intercept) = input_vals
    
    HA_NII_core = (gauss(wave, amp_Ha_core, vel_Ha_core,
                         vel_sigma_Ha_core, 6562.85, 2.5) +
                   (1/3)*gauss(wave, amp_NII6583_core, vel_Ha_core,
                               vel_sigma_Ha_core, 6548.05, 2.5) +
                   gauss(wave, amp_NII6583_core, vel_Ha_core,
                         vel_sigma_Ha_core, 6583.45, 2.5) )
    # HA_NII_wing = (gauss(wave, amp_Ha_wing, vel_Ha_wing,
    #                      vel_sigma_Ha_wing, 6562.85, 2.5) +
    #                 (1/3)*gauss(wave, amp_NII6583_wing, vel_Ha_wing,
    #                             vel_sigma_Ha_wing, 6548.05, 2.5) +
    #                 gauss(wave, amp_NII6583_wing, vel_Ha_wing,
    #                       vel_sigma_Ha_wing, 6583.45, 2.5) )
    
    BLR1 = gauss(wave, amp_Ha_blr1, vel_blr1, vel_sigma_blr1, 6562.85, 2.5)
    # BLR2 = gauss(wave, amp_Ha_blr2, vel_blr2, vel_sigma_blr2, 6562.85, 2.5)
    # BLR3 = gauss(wave, amp_Ha_blr3, vel_blr3, vel_sigma_blr3, 6562.85, 2.5)
    
    continuum = (wave/1000.0)*slope + intercept
    
    return (HA_NII_core + BLR1 + continuum - data)/error

def gauss(wave, amplitude, vel, vel_sigma, rest_wave, inst_res_fwhm) :
    
    kk = 1 + vel/speed_of_light
    numerator = -(wave - rest_wave*kk)**2
    
    sigma = vel_sigma/(speed_of_light - vel_sigma)*rest_wave
    line_width_sqr = sigma**2 + (inst_res_fwhm/fwhm_conv)**2
    
    exponent = 0.5*numerator/line_width_sqr
    line = (amplitude)*np.exp( exponent )
    
    return line

main()
