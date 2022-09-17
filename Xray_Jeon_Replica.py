# Imported Libraries
import numpy as np 
import matplotlib.pyplot as plt

'''
I want to compute the sum over time of the mass density for radiation form black hole accretion and from HMXB emission. This is a replica of
figure 10 from Jeon et al. 2014.
'''

'''
Importing data files of the highest lines from the top and bottom panels of figure 9.
'''

accretion_timeHMXB_dashed, accretion_rateHMXB_dashed = np.loadtxt('/Users/laurenhiggins/Desktop/Research Projects/CIB_Research/analysis/HMXB/tables/line_tables/Dashed_w_Gray_Peak.txt').T
time_Jeon_fig10, rate_Jeon_fig10 = np.loadtxt('/Users/laurenhiggins/Desktop/Research Projects/CIB_Research/analysis/HMXB/tables/Jeon_fig10_og.txt').T

fig = plt.figure(figsize=[14,8])
fig, ax1 = plt.subplots()
ax1.plot( accretion_timeHMXB_dashed, accretion_rateHMXB_dashed, 'r', label='accretion rate')
ax1.set_yscale('log')
ax1.set_xlabel('Time [Myr]')
ax1.set_ylabel('Accretion Rate $\\dot{\\rho}$ [M$\\odot$ yr$^{-1}$]')
#fig.savefig('/Users/laurenhiggins/Desktop/CIB_Research/analysis/HMXB/plots/accretion_rate_dashed.pdf', bbox_inches = "tight")
ax1.legend()
plt.show()

'''
density over cosmic time = the integration of the sum of the accretion rate densities from HMXB. This rate is from dashed lines of the top
panel of Figure 9.
'''
def rho_acc(file_sum, small_time):
    accretion_density = np.zeros(len(file_sum))
    
    def rho_help(rho, t2, t1):
        return rho * (t2 - t1) * 1e6 * (0.7)**3 / 0.3**3
    
    for i in range(0, len(file_sum)-1):
        accretion_density[i] = rho_help(file_sum[i], small_time[i], small_time[i-1])
        if i == 0:
            accretion_density[i] = rho_help(file_sum[i], small_time[i], small_time[i-1]) 
        if i > 0:
            accretion_density[i] = accretion_density[i-1] + accretion_density[i]
    return accretion_density


#dashed_density = rho_acc(interp_sum, accretion_timeHMXB_dashed)
dashed_density = rho_acc(accretion_rateHMXB_dashed, accretion_timeHMXB_dashed)

'''
Replicated plot. Note, that my input accretion rate values range from ~105 Myr to 200 Mry. This accounts for the offset in starting and ending points
in the plot.
'''

fig, ax1 = plt.subplots(figsize=[8,4])
ax1.plot( time_Jeon_fig10, np.log10(rate_Jeon_fig10), 'r-', label='Jeon fig10 w/ HMXBs')
ax1.plot( accretion_timeHMXB_dashed, np.log10(dashed_density), 'g--', label='Higgins w/ HMXB Dashed')
ax1.set_xlim([100, 207])
ax1.set_xticks([100, 120, 140, 160, 180, 200])
#ax1.set_ylim([-4, 3])
#ax1.set_yticks([-4, -2, 0, 2, 4])
#ax1.set_yscale('log')
ax1.yaxis.get_ticklocs(minor=True)
ax1.minorticks_on()
ax1.yaxis.set_ticks_position('both')
ax1.set_xlabel('t$_{H}$(z) [Myr]')
#ax1.set_ylabel('log$_{10}(\\rho_{acc})$ [M$\\odot$ yr$^{-1}$]')
ax1.set_ylabel('$\\rho_{acc}$ [M$\\odot$ yr$^{-1}$]')
ax1.legend(loc=4, fontsize=10)
#fig.savefig('/Users/laurenhiggins/Desktop/Research Projects/CIB_Research/analysis/HMXB/plots/Jeon_fig10_compairson_V5.pdf', bbox_inches = "tight")
plt.show()