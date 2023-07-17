"""
Purpose of this code is to read and plot the output of the Hu+2012 photochemical code, and track our progress adding PH3. 
"""
########################
###Import useful libraries
########################
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pandas as pd
import pdb

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
m2cm=1E2 #1 m in cm
kg2g=1E3 #1 Kg in g
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
ind_cho=61-1


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
ind_n2o5=15-1
ind_hno4=70-1

def plot_comparison_rad(base_file, base_numrad, base_numz, new_file, new_numrad, new_numz, name):
    """
    #Base file
    #New file
    #Title of plot and name of file. 
    """


    ########################
    ###Build 2D grid of radiation at each wavelength and altitude
    ########################
    
    ###Base file
    base_rad=np.zeros((base_numz, base_numrad))
    base_wavs=np.zeros(base_numrad)
    for ind1 in range(0, base_numz):
        linestart=1+ind1*(base_numrad+1)
        if ind1==0:
            base_wavs, base_rad[ind1,:]=np.genfromtxt(base_file, skip_header=linestart, max_rows=base_numrad, usecols=(0,1), unpack=True) 
        else:
            base_rad[ind1,:]=np.genfromtxt(base_file, skip_header=linestart, max_rows=base_numrad, usecols=(1), unpack=True)   
    
    ###New file
    new_rad=np.zeros((new_numz, new_numrad))
    new_wavs=np.zeros(new_numrad)
    for ind1 in range(0, new_numz):
        linestart=1+ind1*(new_numrad+1)
        if ind1==0:
            new_wavs, new_rad[ind1,:]=np.genfromtxt(new_file, skip_header=linestart, max_rows=new_numrad, usecols=(0,1), unpack=True) 
        else:
            new_rad[ind1,:]=np.genfromtxt(new_file, skip_header=linestart, max_rows=new_numrad, usecols=(1), unpack=True)  
            
            
    ########################
    ###Plot
    ########################
    linestyles=np.array(['-',':'])
    
    fig, ax=plt.subplots(1, figsize=(8,6))

    
    ###Top plot: Outgassed species
    ax.plot(base_wavs, base_rad[-1,:], linewidth=2, linestyle=linestyles[0], color='blue', label='Base Rad, TOA')
    ax.plot(new_wavs, new_rad[-1,:], linewidth=2, linestyle=linestyles[0], color='red', label='New Rad, TOA')
    ax.plot(base_wavs, base_rad[0,:], linewidth=2, linestyle=linestyles[1], color='blue', label='Base Rad, BOA')
    ax.plot(new_wavs, new_rad[0,:], linewidth=2, linestyle=linestyles[1], color='red', label='New Rad, BOA')

    ax.legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    

    
    ax.set_yscale('log')
    ax.set_ylabel('Actinic Flux')
    ax.set_xscale('linear')
    ax.set_xlabel('Wavelength')  
    ax.set_xlim([100, 400])
    ax.set_ylim([1E-3, 1E+1])
   
    plt.savefig('./Plots/plot'+name+'.png', orientation='portrait', format='png')
    plt.show()

def plot_comparison_aod(h2so4_meac_file, s8_meac_file, colden_file, plotname):
    """
    Objective of this file is to plot the wavelength-dependent aerosol optical depth from MEAC files.
    h2so4_meac_file: name of H2SO4 MEAC optical parameters file 
    s8_meac_file: name of S8 MEAC optical parameters file 
    colden_file: the relevant ColumnDensity.dat file
    plotname: name of plot
    """
    ########################
    ###Aerosol parameters
    #While these can in principle vary by simulation run, all the simulations we have done here are for a single set of parameters, which tends to maximize the aerosol loading (Hu+2013). This is therefore a worst-case scenario for hazy atmospheres.
    ########################
    AERSIZE=1.0E-7*m2cm #m, converted to cm; from MEAC planet scenario file. Based on Hu+2013, this *should* be the surface area mean diameter, but in fact based on the factor of 2.0558 in RadTransfer.c, this is in fact the median particle diameter, i.e. D_0.
    AERDEN=2.0E+3*kg2g*(m2cm)**-3 #kg m^-3 converted to g cm^-3; from MEAC planet scenario file
    
    mmm_h2so4=98.0*amu2g #mass of H2SO4 molecule, in g
    mmm_s8=256*amu2g #mass of S8 molecule, in g.
    
    n_particles_h2so4=np.pi/6*AERDEN*AERSIZE**3.0*np.exp(1.5*np.log(2.0)**2)/mmm_h2so4 #number of molecules per particle
    n_particles_s8=np.pi/6*AERDEN*AERSIZE**3.0*np.exp(1.5*np.log(2.0)**2)/mmm_s8 #number of molecules per particle
        
    ########################
    ###Import optical parameters, process from per-particle parameters to per-molecule parameters (Hu+2013)
    ########################
    h2so4_wav, h2so4_xc_tot_part, h2so4_xc_scat_part, h2so4_asym_part=np.genfromtxt(h2so4_meac_file, unpack=True) #Units: nm, cm^2/particle, cm^2/particle, dimensionless

    s8_wav, s8_xc_tot_part, s8_scat_part, s8_asym_part=np.genfromtxt(s8_meac_file, unpack=True) #Units: nm, cm^2/particle, cm^2/particle, dimensionless
    
    h2so4_xc_tot=h2so4_xc_tot_part/n_particles_h2so4
    s8_xc_tot=s8_xc_tot_part/n_particles_s8
    
    
    ########################
    ###Import column densities of S8A, H2SO4A, get optical depths
    ########################
    colden_file=pd.read_csv(colden_file,names=['Species Standard Numbers', 'Column Densities', 'Mixing Ratios'], header=0,delimiter="\t")
    
    colden_s8a=colden_file.loc[colden_file['Species Standard Numbers']=='X111']['Column Densities'].values#S8A, cm^-2
    colden_h2so4a=colden_file.loc[colden_file['Species Standard Numbers']=='X78']['Column Densities'].values #H2SO4A, cm^-2
    
    tau_s8=colden_s8a*s8_xc_tot
    tau_h2so4=colden_h2so4a*h2so4_xc_tot
    
    ########################
    ###Plot
    ########################
    linestyles=np.array(['-',':'])
    fig, ax=plt.subplots(2, figsize=(8,8), sharex=True)

    
    ###Top plot: Outgassed species
    ax[0].plot(s8_wav, s8_xc_tot, linewidth=2, linestyle=linestyles[0], color='red', label='S8')
    ax[0].plot(h2so4_wav, h2so4_xc_tot, linewidth=2, linestyle=linestyles[0], color='blue', label='H2SO4')
    ax[0].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    
    ax[0].set_yscale('log')
    ax[0].set_ylim([1E-22, 1E-16])
    ax[0].set_ylabel('Cross-Section (cm2/molecule)')

    ax[1].plot(s8_wav, tau_s8, linewidth=2, linestyle=linestyles[0], color='red', label='S8')
    ax[1].plot(h2so4_wav, tau_h2so4, linewidth=2, linestyle=linestyles[0], color='blue', label='H2SO4')
    ax[1].axhline(1.0, linewidth=2, linestyle=linestyles[1], color='black')
    ax[1].axvline(550, linewidth=2, linestyle='--', color='black')
    ax[1].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)    
    ax[1].legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)  
    ax[1].set_yscale('log')
    ax[1].set_ylabel('Optical Depth')
    ax[1].set_xscale('linear')
    ax[1].set_xlabel('Wavelength')  
    ax[1].set_xlim([100, 1000])
   
    plt.savefig('./Plots/plot'+plotname+'.png', orientation='portrait', format='png')
    plt.show()
###############################
###Run
###############################
plot_comparison_rad('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/Radiation.dat', 1000, 50,'./Photochemical_Model_Outputs/Early_Earth_Depo_10/N2_CO2_early_earth/Radiation.dat', 1000, 50, 'surface_uv_base_10x')
plot_comparison_rad('./Photochemical_Model_Outputs/Early_Earth/N2_CO2_early_earth/Radiation.dat', 1000, 50,'./Photochemical_Model_Outputs/Early_Earth_Depo_30/N2_CO2_early_earth/Radiation.dat', 1000, 50, 'surface_uv_base_30x')

plot_comparison_aod('./H2SO4AER_CrossM_01.dat', './S8AER_CrossM_01.dat', './Photochemical_Model_Outputs/Early_Earth_Depo_10/N2_CO2_early_earth/ColumnDensity.dat', 'AOD_10x')
plot_comparison_aod('./H2SO4AER_CrossM_01.dat', './S8AER_CrossM_01.dat', './Photochemical_Model_Outputs/Early_Earth_Depo_30/N2_CO2_early_earth/ColumnDensity.dat', 'AOD_30x')

