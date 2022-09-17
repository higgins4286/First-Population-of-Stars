# Imported Libraries
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as cn
from scipy.interpolate import interp1d
from scipy.integrate import quadrature

c = cn.c.si.value

wave, l = np.loadtxt('/Users/laurenhiggins/Desktop/CIB_Research/analysis/Jnir_spectra/santos-fig005-fesc1.txt').T #l = total spectrum emitted per wavelength
wave_units = wave * 10**-10 * u.meter #convert from angstroms to meters
nU = (c / wave_units).value #in Hz
l_units = l * u.erg / u.s / u.Hz / u.solMass
lv = (l_units.to('J s-1 Hz-1 kg-1')).value

a, Yz = np.loadtxt('/Users/laurenhiggins/Desktop/CIB_Research/analysis/Jnir_spectra/popIII_sfr_s2_skx_copy.txt', usecols=(0, 3)).T
z = (1/a) - 1
Yz_units = Yz * u.M_sun / u.yr / u.Mpc**3
Yz = Yz_units.to('kg yr-1 m-3')
'''
Now I am creating a z_array for the redshift emitted for the photon to now be at a wavelenth of 2.2 microns. 
This gives me an array of redshifts. Now I have wavelength as a function of redshift.
'''
n=1000
lamarray = np.linspace(0, 10, num=n)
lambobs = (lamarray * 10**-6 * u.meter).value 
z_arrayl2 = (np.zeros(len(lambobs))) #Change name to reflect results

for t in range(0, len(lamarray)):
    def zfcnl(lambem, lamob):
        fcn = (lamob - lambem) / lambem
        return fcn
    
    z_arrayl2 = zfcnl(wave_units.value, lambobs[t])
    
# Below is my interpolation
    SFRDint = interp1d(z, Yz, kind='nearest')
    SFRDz = np.zeros(len(z_arrayl2)) 
    for i in range(0, len(z_arrayl2)):
        if z_arrayl2[i] < np.min(z):
            SFRDz[i] = 0
        else: 
            SFRDz[i] = SFRDint(z_arrayl2[i])

# Below is my code for equation(42) from Santos et al. 2002
    Om = 0.31
    H0 = ((67.79 * u.km / u.s / u.Mpc).to('s-1')).value #convert H0 from km/s/Mpc
                                                        # by canceling out distance units --> seconds
    
    def Ihelp(z, NU, l1, Y):
        cdtdz = c / (H0 * Om**0.5 * (1+z)**2.5) #has units m/s * s
        nW = 10**9                              # This is for converting to nW
        nu = NU / (1 + z)                       #Frequency light emitted and it's new nu based on redshift 
        t = 2e6                                 # year
        jz = (l1 * t * Y) / (4*np.pi)           #sr-1 * (W * Hz-1* m-3)
        return cdtdz * nW * nu * jz             # m nW Hz Hz-1
    
    Inu = np.zeros(len(z_arrayl2))
    for k in range(0, len(z_arrayl2)-1):
        Inu[k] = -1 * quadrature(Ihelp, z_arrayl2[k], z_arrayl2[k+1], args=(nU[k], lv[k], SFRDz[k]), vec_func=False)[0]
        if z_arrayl2[k] < 7:
            Inu[k] = 0

#Below I am saving the files per observed wavelength
    tstr = str(t)
    per_lambob = np.savetxt(f'/Users/laurenhiggins/Desktop/CIB_Research/analysis/Jnir_spectra/tables/z_Inu_per_lambobs2'+tstr+'.txt', np.c_[z_arrayl2, Inu])
    #per_lambob = np.savetxt(f'/Users/laurenhiggins/Desktop/CIB_Research/Jnir_spectra/tables/z_Inu_per_lambobs2{t}.txt', np.c_[z_arrayl2, Inu]) 

from astropy.io import ascii

max_inus = []
for each_index in np.arange(0,n,1):
    file_name = '/Users/laurenhiggins/Desktop/CIB_Research/analysis/Jnir_spectra/tables/z_Inu_per_lambobs2%s.txt'%(each_index)
    loaded_file =ascii.read(file_name)
    max_Inu = np.max(loaded_file['col2'])
    max_inus.append(max_Inu)
    #print(max_inus)

np.savetxt('/Users/laurenhiggins/Desktop/CIB_Research/analysis/Jnir_spectra/tables/max_file2.txt', np.c_[lambobs, max_inus])

wave_table, Inu_table = np.loadtxt('/Users/laurenhiggins/Desktop/CIB_Research/analysis/Jnir_spectra/tables/max_file2.txt').T

fig = plt.figure(figsize=[14,8])
fig, ax1 = plt.subplots()
ax1.plot(wave_table*10**6, np.log10(Inu_table), 'g', label='Inu Max')
ax1.set_xscale('log')
ax1.set_xlim([0.7, 6])
ax1.set_xticks([0.7, 0.8, 1, 2, 4, 6])
ax1.set_ylim([-5, -2])
ax1.set_yticks([-4, -3, -2])
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.set_xlabel('$\lambda_{obs}$ [$\mu$m]', fontsize=12)
ax1.set_ylabel('$\\nu$$I$$_\\nu$ [nW m$^{-2}$ sr$^{-1}$]', fontsize=12)
ax1.legend(fontsize='8')
#fig.savefig('/Users/laurenhiggins/Desktop/CIB_Research/Jnir_spectra/maxInu_lambob.pdf', bbox_inches = "tight")
plt.show()