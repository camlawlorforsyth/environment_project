# AGNfitter Install #

Instructions to properly install Anaconda2, Astropy, acor, and AGNfitter, as required by the SED fitting analysis.

## During CentOS 7.6 Installation ##

Choose the 'Minimal Desktop' group from the install options when using the full DVD iso.

Disable Kdump in order to use as much system memory as possible.

## Install Guest Additions ##

First, connect to the internet on the guest OS. Next, as root: 
```
yum -y update
yum -y install kernel-devel
yum -y update
```
To then ensure the correct compilers are available (again as root):
```
yum -y install gcc gcc-c++ gcc-gfortran
```
An optional restart may be required here.

Next, insert the Guest Additions disc, and follow the on-screen prompts to install the Guest Additions.

Once the Guest Additions have been installed, restart the system.

## Install Anaconda, Astropy ##

As AGNfitter uses Python 2.7, it requires an older version of Astropy. The easiest way to obtain the correct version of Python and Astropy is to use the Anaconda distribution, available at: https://www.anaconda.com/distribution/

Select the Python 2.7 version, and the 64-Bit (x86) Installer. Once the installer has completed downloading, follow the instructions for installation here: https://docs.anaconda.com/anaconda/install/linux/

## Install acor ##

Once the correct versions of Python and Astropy have been installed, acor then needs to be installed. See the installation instructions here: https://github.com/dfm/acor

## Install AGNfitter ##

Once all the required packages have been installed, follow the installation instructions here to get AGNfitter working: https://github.com/GabrielaCR/AGNfitter

Also see the paper by Calistro Rivera, describing the fitter itself: http://adsabs.harvard.edu/abs/2016ApJ...833...98C

## AGNfitter Modifications ##

### No LaTeX ###

The default AGNfitter requires LaTeX for some output files (ie. MCMC plots, SED plots) which is quite heavy if you want a lightweight system.

To disable the LaTeX requirements, we need to modify the `PLOTandWRITE_AGNfitter.py` file in the `functions/` directory.

Comment-out lines 290, 291, 549, and 550.

### Use c = 299792458 m/s ###

The defualt AGNfitter uses a value close to the actual value for the speed of light, but not quite the real value.

To correct this, we need to modify several files in the `functions/` directory, namely:

Change line 18-19 of `CONSTRUCT_modelobjects.py` to:
```
c = 2.99792458e10
c_Angst = 3.3356409519815204e-19 #(1/(c*Angstrom)
```
Change lines 423 and 625 of `DICTIONARIES_AGNfitter.py` to:
```
c=    2.99792458e8
```
Change lines 138, 169, 220, 379 of `MODEL_AGNfitter.py` to (respectively):
```
redd_x =  2.99792458 * 1e10 / (10**(bbb_x)* 1e-8)

c = 2.99792458 * 1e8

c_cm = 2.99792458e10

c = 2.99792458e8
```

Change line 580 of `PLOTandWRITE_AGNfitter.py` to:
```
x2 = (2.99792458e14/ x)[::-1] # Wavelength axis
```

### Smaller output changes ###

#### Append the galaxy ID to the traces plot filename. ####

Change line 76 of `PLOTandWRITE_AGNfitter.py` to:
```
fig.savefig(data.output_folder+str(data.name)+'/traces_mcmc_' +str(data.name) + '.' + out['plot_format'])
```

#### Comment-out the header rows for the output parameter file. ####

Change line 91 of `PLOTandWRITE_AGNfitter.py` to:
```
comments_ouput= '# Output for source ' +str(data.name) + '\n' +'# Rows are: 2.5, 16, 50, 84, 97.5 percentiles '+'\n'+ '#-----------------------------------------------------'+'\n' 
```
Change line 159 of `PLOTandWRITE_AGNfitter.py` to:
```
outputvalues_header= ' '.join([ i for i in np.hstack((P.names, 'log_Mstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', '-ln_like'))] )
```
