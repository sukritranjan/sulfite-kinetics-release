"""
Purpose of this code is to read and plot the output of the Hu+2012 photochemical code
"""
########################
####Import useful libraries
########################
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pdb

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
deg2rad=np.pi/180.
bar2barye=1.e+6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
metricton2kg=1000. #1 metric ton in kg

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
R_earth=6371.*km2m#radius of earth in m
R_sun=69.63e9 #radius of sun in cm
AU=1.496e13#1AU in cm


########################
###Establish key
########################

#Corrected for Hu 1-indexing vs Python 0-indexing
ind_o=1-1 #O
ind_h=3-1 #H
ind_oh=4-1 #OH
ind_so2=43-1
ind_so=42-1
ind_ch4=21-1
ind_h2s=45-1
ind_h2=53-1
ind_h2o=7-1

ind_no=12-1
ind_n2o=11-1

ind_co2=52-1
ind_n2=55-1


ind_s8=79-1
ind_s8a=111-1
ind_ch4o=24-1
ind_c2h2=27-1

ind_ocs=49-1

ind_o1d=56-1
ind_co=20-1
ind_o2=54-1
ind_c2h6=31-1
ind_h2o2=6-1
ind_h2so4=73-1
ind_h2so4a=78-1
ind_ch2o=22-1
ind_o3=2-1
ind_ho2=5-1

def plot_comparison(base_file, new_file, name):
    """
    #Base file
    #New file
    #Title of plot and name of file. 
    """
    
    ###Initialize plot
    fig2, ax=plt.subplots(2, figsize=(8., 10.), sharey=True)
    markersizeval=5.

    ########################
    ###Read in base data
    ########################
    base_data=np.genfromtxt(base_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_base=base_data[:,0] # Center of altitude bins, km 
    T_z_base=base_data[:,3] # Temperature(z), in K
    P_z_base=base_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    n_z_s_base=base_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_base=np.sum(n_z_s_base,1) #sum number densities across species. This is a profile for the whole atmosphere.
    n_z_bulkatm_base=P_z_base/(k*T_z_base)
    mc_z_s_base=np.zeros(np.shape(n_z_s_base))
    mr_z_s_base=np.zeros(np.shape(n_z_s_base))

    num_s=np.shape(n_z_s_base)[1]

    for ind2 in range(0, num_s):
        mc_z_s_base[:,ind2]=n_z_s_base[:,ind2]/n_z_base#molar concentration of each species.
        mr_z_s_base[:,ind2]=n_z_s_base[:,ind2]/n_z_bulkatm_base#mixing ratio of each species.

    ########################
    ###Read in new data
    ########################    
    new_data=np.genfromtxt(new_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.

    z_center_new=new_data[:,0] # Center of altitude bins, km
    T_z_new=new_data[:,3] # Temperature(z), in K
    P_z_new=new_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    n_z_s_new=new_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)
    
    ###Get molar concentrations
    n_z_new=np.sum(n_z_s_new,1) #sum number densities across species. This is a profile for the whole atmosphere.
    n_z_bulkatm_new=P_z_new/(k*T_z_new)
    
    mc_z_s_new=np.zeros(np.shape(n_z_s_new))
    mr_z_s_new=np.zeros(np.shape(n_z_s_new))
    
    num_s=np.shape(n_z_s_new)[1]

    for ind2 in range(0, num_s):
        mc_z_s_new[:,ind2]=n_z_s_new[:,ind2]/n_z_new#molar concentration of each species.
        mr_z_s_new[:,ind2]=n_z_s_new[:,ind2]/n_z_bulkatm_new#mixing ratio

##    ########################
##    ###Print key parameters
##    ########################        
#
#    
#     print('Column-Averaged Mixing Ratios: Photochemical Products')
# ##    ##CO2-vale
#     print('SO2 (base): {0:1.1e}'.format((np.sum(mr_z_s_base[:,ind_so2]*n_z_bulkatm_base)/np.sum(n_z_bulkatm_base))))
#     print('SO2 (new): {0:1.1e}'.format((np.sum(mr_z_s_new[:,ind_so2]*n_z_bulkatm_new)/np.sum(n_z_bulkatm_new))))
#     print('H2S (base): {0:1.1e}'.format((np.sum(mr_z_s_base[:,ind_h2s]*n_z_bulkatm_base)/np.sum(n_z_bulkatm_base))))
#     print('H2S (new): {0:1.1e}'.format((np.sum(mr_z_s_new[:,ind_h2s]*n_z_bulkatm_new)/np.sum(n_z_bulkatm_new))))

    print('Surface Mixing Ratios: Photochemical Products')
    print('SO2 (base): {0:1.1e}'.format(mr_z_s_base[0,ind_so2]))
    print('SO2 (new): {0:1.1e}'.format(mr_z_s_new[0,ind_so2]))
    print('H2S (base): {0:1.1e}'.format(mr_z_s_base[0,ind_h2s]))
    print('H2S (new): {0:1.1e}'.format(mr_z_s_new[0,ind_h2s]))
    
    
    print('Surface Mixing Ratios: Photochemical Products')
    print('O2 (base): {0:1.1e}'.format(mr_z_s_base[0,ind_o2]))
    print('O2 (new): {0:1.1e}'.format(mr_z_s_new[0,ind_o2]))
   


    
    ########################
    ###Plot
    ########################
    linestyles=np.array(['-',':'])
    
    ###Top plot: Outgassed species
    ax[0].plot(mr_z_s_base[:,ind_h2], z_center_base, linewidth=2, linestyle=linestyles[0], color='purple', label='H2')
    ax[0].plot(mr_z_s_base[:,ind_n2], z_center_base, linewidth=2, linestyle=linestyles[0], color='olive', label='N2')
    ax[0].plot(mr_z_s_base[:,ind_co2], z_center_base, linewidth=2, linestyle=linestyles[0], color='magenta', label='CO2')
    ax[0].plot(mr_z_s_base[:,ind_h2o], z_center_base, linewidth=2, linestyle=linestyles[0], color='pink', label='H2O')
    ax[0].plot(mr_z_s_base[:,ind_ch4], z_center_base, linewidth=2, linestyle=linestyles[0], color='cyan', label='CH4')
    ax[0].plot(mr_z_s_base[:,ind_h2s], z_center_base, linewidth=2, linestyle=linestyles[0], color='yellow', label='H2S')
    ax[0].plot(mr_z_s_base[:,ind_so2], z_center_base, linewidth=2, linestyle=linestyles[0], color='orange', label='SO2')
    ax[0].plot(mr_z_s_base[:,ind_no], z_center_base, linewidth=2, linestyle=linestyles[0], color='black',label='NO')


    ##
    ax[0].plot(mr_z_s_new[:,ind_h2], z_center_new, linewidth=2, linestyle=linestyles[1], color='purple')
    ax[0].plot(mr_z_s_new[:,ind_n2], z_center_new, linewidth=2, linestyle=linestyles[1], color='olive')
    ax[0].plot(mr_z_s_new[:,ind_co2], z_center_new, linewidth=2, linestyle=linestyles[1], color='magenta')
    ax[0].plot(mr_z_s_new[:,ind_h2o], z_center_new, linewidth=2, linestyle=linestyles[1], color='pink')
    ax[0].plot(mr_z_s_new[:,ind_ch4], z_center_new, linewidth=2, linestyle=linestyles[1], color='cyan')
    ax[0].plot(mr_z_s_new[:,ind_h2s], z_center_new, linewidth=2, linestyle=linestyles[1], color='yellow')
    ax[0].plot(mr_z_s_new[:,ind_so2], z_center_new, linewidth=2, linestyle=linestyles[1], color='orange')
    ax[0].plot(mr_z_s_new[:,ind_no], z_center_new, linewidth=2, linestyle=linestyles[1], color='black')

    
    ax[0].set_title(name)
    ax[0].set_yscale('linear')
    ax[0].set_ylabel('Altitude (km)')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Mixing Ratio')  
    ax[0].set_xlim([1.e-14, 1.e0])
    ax[0].legend(loc=2, ncol=1, borderaxespad=0., fontsize=12)    
    
    ####Bottom Plot: Photochemical Products
    ax[1].plot(mr_z_s_base[:,ind_o], z_center_base, linewidth=2, linestyle=linestyles[0], color='green', label='O')
    ax[1].plot(mr_z_s_base[:,ind_h], z_center_base, linewidth=2, linestyle=linestyles[0], color='red', label='H')
    ax[1].plot(mr_z_s_base[:,ind_oh], z_center_base, linewidth=2, linestyle=linestyles[0], color='blue', label='OH')
    ax[1].plot(mr_z_s_base[:,ind_co], z_center_base, linewidth=2, linestyle=linestyles[0], color='grey', label='CO')
    ax[1].plot(mr_z_s_base[:,ind_o2], z_center_base, linewidth=2, linestyle=linestyles[0], color='hotpink', label='O2')
    ax[1].plot(mr_z_s_base[:,ind_o3], z_center_base, linewidth=2, linestyle=linestyles[0], color='saddlebrown', label='O3')
    ax[1].plot(mr_z_s_base[:,ind_so], z_center_base, linewidth=2, linestyle=linestyles[0], color='mediumorchid', label='SO')
    ax[1].plot(mr_z_s_base[:,ind_ho2], z_center_base, linewidth=2, linestyle=linestyles[0], color='black', label='HO2')

    
    ax[1].plot(mr_z_s_new[:,ind_o], z_center_new, linewidth=2, linestyle=linestyles[1], color='green')
    ax[1].plot(mr_z_s_new[:,ind_h], z_center_new, linewidth=2, linestyle=linestyles[1], color='red')
    ax[1].plot(mr_z_s_new[:,ind_oh], z_center_new, linewidth=2, linestyle=linestyles[1], color='blue')
    ax[1].plot(mr_z_s_new[:,ind_co], z_center_new, linewidth=2, linestyle=linestyles[1], color='grey')
    ax[1].plot(mr_z_s_new[:,ind_o2], z_center_new, linewidth=2, linestyle=linestyles[1], color='hotpink')
    ax[1].plot(mr_z_s_new[:,ind_o3], z_center_new, linewidth=2, linestyle=linestyles[1], color='saddlebrown')   
    ax[1].plot(mr_z_s_new[:,ind_so], z_center_new, linewidth=2, linestyle=linestyles[1], color='mediumorchid')   
    ax[1].plot(mr_z_s_new[:,ind_ho2], z_center_new, linewidth=2, linestyle=linestyles[1], color='black')   

    ax[1].legend(loc=2, ncol=1, borderaxespad=0., fontsize=12)    

    
    ax[1].set_yscale('linear')
    ax[1].set_ylabel('Altitude (km)')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Mixing Ratio')  
    ax[1].set_xlim([1.e-20, 1.e-0])
    
    plt.savefig('./Plots/plot'+name+'.pdf', orientation='portrait', format='pdf')
    plt.show()


###############################
###Run
###############################
plot_comparison('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','./Photochemical_Model_Outputs/Early_Earth_Depo_0.1/N2_CO2_early_earth/ConcentrationSTD.dat','Outgassing_1_0.1')
plot_comparison('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','./Photochemical_Model_Outputs/Early_Earth_Depo_0.3/N2_CO2_early_earth/ConcentrationSTD.dat','Outgassing_1_0.3')
plot_comparison('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','Outgassing_1_1')
plot_comparison('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','./Photochemical_Model_Outputs/Early_Earth_Depo_3/N2_CO2_early_earth/ConcentrationSTD.dat','Outgassing_1_3')
plot_comparison('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','./Photochemical_Model_Outputs/Early_Earth_Depo_10/N2_CO2_early_earth/ConcentrationSTD.dat','Outgassing_1_10')
plot_comparison('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat','./Photochemical_Model_Outputs/Early_Earth_Depo_30/N2_CO2_early_earth/ConcentrationSTD.dat','Outgassing_1_30')


