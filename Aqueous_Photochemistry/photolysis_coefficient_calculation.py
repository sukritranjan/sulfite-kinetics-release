 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: sukrit, modifying DW paper 1 codes.

Purpose of these scripts is to calculate column-integrated photolysis rates for nitrate (validation case) and sulfite (science case), for use in mass balance scripts. 
"""
import numpy as np #This furnishes array operations.
import scipy as sp
from scipy import interpolate as interp
import matplotlib.pyplot as plt
import pickle
import scipy
import scipy.integrate
import pdb
import pandas as pd
import spectres
#import cookbook.py as cookbook

cm2nm=1.0E+7 #1 cm in nm
hc=1.98645e-9 #value of h*c in erg*nm, useful to convert from ergs/cm2/s/nm to photons/cm2/s/nm
day2s=86400.0 #1 day in seconds
deg2rad=np.pi/180.0
rad2deg=180.0/np.pi

N_avogadro=6.02E23
year2s=3.154e+7
m2cm=100.0

#solution specific conversion
M2cgs=(6.02e23)*1.e-3 #mol/L * particles/mol * L/cm**3 = molecules cm**-3
cgs2M=M2cgs**-1 #molecules cm**-3 to mol/L

def linear_interp(xs, x_ref, y_ref, y_ref_errs):
    """
    WARNING HUGE ASSUMPTION: ASSUMES x_ref (and xs?) INCREASE MONOTONICALLY
    MUST FORMAT INPUTS TO COMPLY.
    

    Parameters
    ----------
    xs : abscissa to project to
    x_ref : abscissa of input data
    y_ref : y-values of input data
    y_ref_errs : errors on y-values of input data
    Returns
    -------
    ys: interpolated values. Should be identical to numpy interp.
    y_errs: interpolated errors, calculated assuming uncorrelated Gaussian errors.
    
    https://stackoverflow.com/questions/24616079/how-to-interpolate-using-nearest-neighbours-for-high-dimension-numpy-python-arra
    
    POTENTIAL ALTERNATIVE: SPECTRES (https://spectres.readthedocs.io/en/latest/)
    """
    ###Assume 0 unless otherwise known.
    ys=np.zeros(np.shape(xs))
    y_errs=np.zeros(np.shape(xs))
    
    ###Get indices
    rights = np.searchsorted(x_ref, xs, 'left')
    
    for ind in range(0, len(xs)):
        right=rights[ind]
        left=right-1
        x=xs[ind]
        if right==0 and x<x_ref[0]: #If x is less than the x_ref range entirely, zero the values
            ys[ind]=0
            y_errs[ind]=0
        elif right==0 and x==x_ref[0]: #if x is exactly the beginning of the x_ref range, then it should just be that value. 
            ys[ind]=y_ref[0]
            y_errs[ind]=y_ref_errs[0]
        elif right==len(x_ref): #If x is greater than the x_ref range entirely, the it should be 0. 
            ys[ind]=0
            y_errs[ind]=0
        else:

            x1=x_ref[left]
            y1=y_ref[left]
            y1_err=y_ref_errs[left]
            
            x2=x_ref[right]
            y2=y_ref[right]
            y2_err=y_ref_errs[right]
            
            a=(x2-x)/(x2-x1)
            b=1.0-a
            ys[ind]=a*y1 + b*y2
            y_errs[ind]=np.sqrt((a*y1_err)**2.0 + (b*y2_err)**2.0)
  
    ###Compare to what np.interp returns, just to be safe.
    y_np_interp=np.interp(xs, x_ref, y_ref, left=0, right=0)
    residuals=np.abs(ys-y_np_interp)/(0.5*(ys+y_np_interp)) #residuals. Should be 0.
    
    if np.sum(residuals>10.0*np.finfo(float).eps)>0: #if the residuals exceed 10x the precision of a float ANYWHERE, kill everything dramatically so we see it. 
        ys=np.nan
        y_errs=np.nan
        print('Error in custom linear interpolation function') 
    return ys, y_errs

def get_aqueous_angle_snell(zenith_angle_air):
    """
    This function uses Snell's law to convert the air zenith angle to the aqueous zenith angle. It implements Equation 2.15 of Kirk (1994). Key assumptions are (1) unpolarized light, and that the ratio of the refractive indices of water and air are n_w/n_a=1.33. This is valid for modern air, modern sea and fresh water, and the wavelengths of light corresponding to PAR; Cockell (2000) also extend it to the UV and to early Earth. We follow Cockell (2000) in this.
    
    This formula successfully reproduces the calculations on Page 775, column 1 of Bannister 1992.
    
    Takes: 
        -Zenith angle in air, in radians
    
    Returns:
        -Zenith angle (or nadir angle if you will) in water, in radians.
    """
    if zenith_angle_air<1.0E-12:
        zenith_angle_air=1.0E-12
    
    return np.arcsin(np.sin(zenith_angle_air)/1.33)

def get_reflectance_fresnel(zenith_angle_air, zenith_angle_water):
    """
    This function uses Fresnel's Equation to calculate the reflectance of a beam traversing two mediums of different index of refraction. It implements Equation 2.12 of Kirk (1994). Key assumptions are (1) unpolarized light, and (2) flat water.
    
    Does a pretty good job reproducing Table 2.1 of Kirk (1994)
    
    Takes: 
        -Zenith angle in air, in radians
        -Zenith angle (or nadir angle if you will) in water, in radians.
       
    
    Returns:
        -Reflectance back to water.
    """
    return (0.5*(np.sin(zenith_angle_air-zenith_angle_water)**2.0/np.sin(zenith_angle_air+zenith_angle_water)**2.0) + 0.5*(np.tan(zenith_angle_air-zenith_angle_water)**2.0/np.tan(zenith_angle_air+zenith_angle_water)**2.0))

def calculate_colint_nitrate_photolysis_modernocean():
    """
    As a test of the basic approach, calculate nitrate photolysis rates for MODERN EARTH OCEANS.

    Returns
    -------
    Depth-integrated photolysis rate for modern earth ocean (cm s**-1). Multiply by concentration to obtain loss flux. 

    """
    ##############
    ###Load modern earth UV from Ranjan & Sasselov 2017, specifically when we matched Rugheimer et al. 2015 exactly. 
    ##############
    uv_file='./rugheimer_earth_modern.dat'
    # uv_file='./equatorialspec.dat'
    
    
    uv_wav_left, uv_wav_right, uv_wav, uv_toa_intensity, uv_surface_flux, uv_surface_intensity=np.genfromtxt(uv_file, skip_header=1, skip_footer=0, usecols=(0, 1, 2,3,4,6), unpack=True) #From Ranjan & Sasselov 2017, the GitHub. 
    
    ###Convert from ergs to photons
    uv_toa_intensity *= uv_wav/hc
    uv_surface_flux *= uv_wav/hc
    uv_surface_intensity *= uv_wav/hc
    
    ###Restrict data to 290 to 340 nm. Photolysis takes place up to 340 nm (Zafiriou+1979), below 290 nm blocked by ozone layer ()
    inds=np.where((uv_wav>=290.0) & (uv_wav <=340.0))
    
    uv_wav=uv_wav[inds]
    uv_wav_left=uv_wav_left[inds]
    uv_wav_right=uv_wav_right[inds]
    uv_toa_intensity=uv_toa_intensity[inds]
    uv_surface_flux=uv_surface_flux[inds]
    uv_surface_intensity=uv_surface_intensity[inds]
    
    ##############
    ###Load absorption of modern oceans, from Smith & Baker 1981. Well, in reality, the purest modern ocean waters -- highly UV-transparent. 
    ##############
    #pureoceans_smith1981_wav, pureoceans_smith1981_abs=np.genfromtxt('./Azra_Project/RSI_UV/Processed-Data/smithbaker_purest.dat', skip_header=2, skip_footer=0, unpack=True, usecols=(0,1)) #purest modern natural water; nm, cm**-1.
    
    df2=pd.read_excel('./smith_baker.xlsx', sheet_name='Sheet1', skiprows=3)
    pureoceans_smith1981_wav=np.nan_to_num(df2['lambda']) #nm 
    pureoceans_smith1981_Kd=np.nan_to_num(df2['Ksww'])*1.0E-2 #m**-1 converted to cm**-1
    pureoceans_smith1981_a=np.nan_to_num(df2['aw'])*1.0E-2 #m**-1 converted to cm**-1
    pureoceans_smith1981_b_sw=np.nan_to_num(df2['bswm'])*1.0E-2 #m**-1 converted to cm**-1, salt water
    pureoceans_smith1981_b_fw=np.nan_to_num(df2['bfwm'])*1.0E-2 #m**-1 converted to cm**-1, fresh water
    
    ##############
    ###Load absorption of nitrate (Ranjan+2022a)
    ##############
    
    data_wav={} #dict to hold wavelength scale for measured molar absorptivities (nm)
    data_abs={} #dict to hold molar absorptivities for measured molar absorptivities (M**-1 cm**-1.)
    
    ###From Gabi/Corinna's work
    df=pd.read_excel('./2020-08-23_Extinktion-coefficients-Sukrit-V7_mod.xlsx', sheet_name='Tabelle1') 
    data_wav['LK']=np.nan_to_num(df['Wavelength']) #nm LK=Lozano-Kufner
    
    data_abs['NaNO3_LK']=np.nan_to_num(df['NaNO3'])
    data_abs['NaNO3_err_LK']=np.nan_to_num(df['NaNO3-error'])
    
    ###Convert molar decadic absorption coefficients to absorption coefficients.
    nano3_xc=data_abs['NaNO3_LK']*1.0E3*np.log(10.0)/(6.02E23) #convert from M^-1 cm^-1 to cm^2 (molecule^-1)
    nano3_xc_err=data_abs['NaNO3_err_LK']*1.0E3*np.log(10.0)/(6.02E23) #convert from M^-1 cm^-1 to cm^2 (molecule^-1)
   
    ##############
    ###Get quantum yields
    ##############
    #From table 2 of Mack & Bolton. This is best for comparison to Zafiriou+1979b as it best matches the actual process they are trying to constrain. 
    nitrate_qy_wav=np.array([282.0, 305.0, 310.0, 313.0]) #nm
    nitrate_qy_qy=np.array([0.024, 0.0059, 0.0065, 0.0063]) #dimensionless
    
    qy_nitrate_func=interp.interp1d(nitrate_qy_wav,nitrate_qy_qy, kind='linear', fill_value=(nitrate_qy_qy[0], nitrate_qy_qy[-1]), bounds_error=False) #Extrapolate with last values in either direction.
    
    
    ##############
    ###SURFACE PHOTOLYSIS RATE. PRINT TEST CASES AND COMPARISONS TO SCREEN
    ##############
    print('Surface photolysis rates, different prescriptions')
    
    #####
    #FIRST, REPORT WHAT ZAFIRIOU+1979b DIRECTLY MEASURED
    #####
    print('Surface Photolysis Rate (directly inferred, Zafiriou+1979b) (s**-1; median, range): {0:1.1e} ({1:1.1e}-{2:1.1e})'.format(2.0E-3/day2s, 0.7E-3/day2s, 6.2E-3/day2s)) #Comes from conversion of rates in days**-1 on page 40 of Zafiriou+1979b.
    
    #####
    #SECOND, DO USING ZAFIRIOU+1979b VALUES, EVEN IF WRONG. THE QYs ARE DEFINITELY WRONG.
    #####
    ###For now, take solar irradiation, bins from Zafiriou+1979b.
    wav_centers_zaf=np.array([295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0])
    wav_left_zaf=wav_centers_zaf-2.5
    wav_right_zaf=wav_centers_zaf+2.5
    photflux_zaf=np.array([0.0007, 0.17, 0.96, 2.65, 5.00, 7.55, 9.95, 12.2, 14.3, 16.3])*1E12 #(photons) cm^-2 s^-1, integrated across those wavelength bins. 
    nitrat_xc_zaf=np.array([21.8, 26.8, 26.8, 22.3, 19.1, 21.6, 7.6, 3.8, 1.9, 0.9])*1.0E-21 #cm^2 (molecule^-1)
    qy_nitrate_zaf=0.04
    excit_rate_zaf=np.sum(photflux_zaf*nitrat_xc_zaf*qy_nitrate_zaf) #
    print('Surface Photolysis Rate (trusting Table V, Zafiriou+1979b) (s**-1): {0:1.1e}'.format(excit_rate_zaf))
    
    # #####
    # #THIRD, DO USING THE ZAFIRIOU+1979b SOLAR FLUX BUT INCLUDING MODERN XC AND QY. QY IN PARTICULAR HAS CHANGED.
    # #####
    # ###For now, take solar irradiation, bins from Zafiriou+1979b.
    # excit_rate_zaf_mod=np.sum(photflux_zaf*np.interp(wav_centers_zaf, data_wav['LK'], nano3_xc)*qy_nitrate_func(wav_centers_zaf)) #
    # print('Surface Photolysis Rate (Zafiriou+1979b Solar IRRAD, New XC/QY) (s**-1): {0:1.1e}'.format(excit_rate_zaf_mod))
    # ###This is really off, which makes sense given that the QY Zafiriou was using were really off. Also, we were never able to locate the source for their solar radiation, i.e. Johnston et al. 1976 -- missing from reference list also. 
    
    #####
    #FINALLY, DO IT HOW WE THINK IT SHOULD BE DONE, USING OUR MODEL RESULTS
    #####
    uv_wav_widths=uv_wav_right-uv_wav_left
    excit_rate2=np.sum(uv_wav_widths*uv_surface_flux*(np.interp(uv_wav, pureoceans_smith1981_wav, pureoceans_smith1981_Kd)/np.interp(uv_wav, pureoceans_smith1981_wav, pureoceans_smith1981_a))*np.interp(uv_wav, data_wav['LK'], nano3_xc)*qy_nitrate_func(uv_wav)) #s**-1
    print('Surface Photolysis Rate (modern solar, XCs, QY) (s**-1): {0:1.1e}'.format(excit_rate2))
    
    ##############
    ###COLUMN-INTEGRATED PHOTOLYSIS RATE. PRINT TEST CASES AND COMPARISONS TO SCREEN
    ##############
    print('Column-integrated photolysis rates, different prescriptions')

    ###Establish depth scale with depth. 
    binwidth=1.0 #cm
    depth_bin_left=np.arange(0.0, 5000.0, step=binwidth) #cm. photic depth should be 5 meters. We will go down to 50 to be safe. 
    depth_bin_right=depth_bin_left+binwidth
    depth_bin_centers=0.5*(depth_bin_left+depth_bin_right)
    
    
    #####
    #ZAFIRIOU CALCULATION UNALTERED
    #####    
    print('Column Photolysis Rate  (Zafiriou Surface Photolysis Rate, Integration Method) (cm s**-1): {0:1.1e}'.format((2.3E-8*0.5*500.0)))#Z+1979b estimate: 2.3E-8 s**-1 *0.5 (depth) * 500 cm = 5.8E-6 cm s**-1. Detailed on page 41 of Z+1979.
    
    #####
    #ZAFIRIOU SUFRACE PHOTOLYSIS, OUR INTEGRATION
    #####        
    reac_rate_depth_tot_zaf=2.3E-8*binwidth*np.sum((np.interp(305.0, pureoceans_smith1981_wav, pureoceans_smith1981_Kd)/np.interp(305.0, pureoceans_smith1981_wav, pureoceans_smith1981_a))*np.exp(-depth_bin_centers*np.interp(305.0, pureoceans_smith1981_wav, pureoceans_smith1981_Kd)))    
    print('Column Photolysis Rate (Zafiriou Surface Photolysis Rate, Our Integration Method) (cm s**-1): {0:1.1e}'.format(reac_rate_depth_tot_zaf))
    
    ####
    #OUR METHOD
    ####         
    reac_rate_depth=0.0*depth_bin_centers
    
    for ind in range(0, len(depth_bin_centers)):
        uv_flux_depth=uv_surface_flux*np.exp(-depth_bin_centers[ind]*np.interp(uv_wav, pureoceans_smith1981_wav, pureoceans_smith1981_Kd))
    
        
        spectral_reac_rate_depth=uv_flux_depth*(np.interp(uv_wav, pureoceans_smith1981_wav, pureoceans_smith1981_Kd)/np.interp(uv_wav, pureoceans_smith1981_wav, pureoceans_smith1981_a))*np.interp(uv_wav, data_wav['LK'], nano3_xc)*qy_nitrate_func(uv_wav)
        
        reac_rate_depth[ind]=np.sum(uv_wav_widths*spectral_reac_rate_depth) 
    
    reac_rate_depth_tot=np.sum(binwidth*reac_rate_depth)
    print('Column Photolysis Rate (Our Surface Photolysis Rate, Integration Method) (cm s**-1): {0:1.1e}'.format(reac_rate_depth_tot))
    
    return reac_rate_depth_tot

# calculate_colint_nitrate_photolysis_modernocean()

def sulfite_photolysis_rate(conc_s_iv, d_body, pH, I, water_abs_choice, longwaveQY0, K_d_method):
    """
    Minimum depth 1 cm. 
    If using pure water absorbance: make sure RT only goes out to 320 nm to match wavelength limits.
    """
    ##############
    ###Load early Earth UV from Ranjan & Sasselov 2017 (or should we use Ranjan et al. 2022?)
    ##############    
    uv_file='./rugheimer_earth_epoch0.dat'  
    theta_sun=60.0*deg2rad #For the source RT calculation, the solar zenith angle is 60 degrees, hard-coded.
    mu_1=1.0/np.sqrt(3.0) #For the source RT calculation, mu_1=1/sqrt(3) (quadrature closure)
    
    uv_wav_left, uv_wav_right, uv_wav, uv_toa_intensity, uv_surface_flux, uv_surface_intensity, uv_diffuse_intensity, uv_direct_intensity=np.genfromtxt(uv_file, skip_header=1, skip_footer=0, usecols=(0,1,2,3,4,6,7,8), unpack=True) #From Ranjan & Sasselov 2017, the GitHub. 
    
    ###Convert from ergs to photons
    uv_toa_intensity *= uv_wav/hc
    uv_surface_flux *= uv_wav/hc
    uv_surface_intensity *= uv_wav/hc
    uv_diffuse_intensity *= uv_wav/hc
    uv_direct_intensity *= uv_wav/hc
    
    ###Restrict data to 200 to 330 nm. 320 set by where S[IV] absorbs, 200 by what is available on early Earth. 
    inds=np.where((uv_wav>=200.0) & (uv_wav <=320.0))
    
    uv_wav=uv_wav[inds]
    uv_wav_left=uv_wav_left[inds]
    uv_wav_right=uv_wav_right[inds]
    uv_toa_intensity=uv_toa_intensity[inds]
    uv_surface_flux=uv_surface_flux[inds]
    uv_surface_intensity=uv_surface_intensity[inds]
    uv_diffuse_intensity=uv_diffuse_intensity[inds]
    uv_direct_intensity=uv_direct_intensity[inds]
    
    ###Obtain some useful secondary variables
    uv_wav_widths=uv_wav_right-uv_wav_left
    uv_diffuse_flux=mu_1*uv_diffuse_intensity
    uv_direct_flux=(np.cos(theta_sun))*uv_direct_intensity
    
    ###Calculate transmission through air/water interface, following Morel 1991
    uv_diffuse_flux_aq=uv_diffuse_flux*(1.0-0.066) #constant sky reflectance
    
    theta_sun_aq=get_aqueous_angle_snell(theta_sun) #"nadir angle" after refraction of dreict stream
    uv_direct_flux_aq=uv_direct_flux*(1.0-get_reflectance_fresnel(theta_sun, theta_sun_aq))
    uv_surface_flux_aq=uv_diffuse_flux_aq+uv_direct_flux_aq#This is the E_d(0-) term.
    
    # ##Let's make sure the things that should sum to each other do. Don't expect exactly the same due to non-perfect precision in printed values, but should be VERY close...
    # print(np.max(np.abs(uv_surface_intensity-(uv_diffuse_intensity+uv_direct_intensity))/uv_surface_intensity))
    # print(np.max(np.abs(uv_surface_flux-(uv_diffuse_flux+uv_direct_flux))/uv_surface_flux))
    # pdb.set_trace()
    # #..and so they are.
    
    # ###Plot to ensure the spectra have been properly loaded.
    # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # ax.plot(uv_wav, uv_surface_flux, linewidth=4, linestyle='-', marker='d', color='black')    
    # ax.set_yscale('log')
    # ax.set_ylabel(r'UV Surface Flux (cm$^{-2}$ s$^{-1}$ nm$^{-1}$)', fontsize=14)
    # # ax.legend(ncol=1, loc='best')    

    # ax.set_xscale('linear')
    # ax.set_xlabel('Wavelength (nm)', fontsize=14)
    # ax.set_xlim([200., 330.])
    # plt.savefig('./Plots/surface_UV.pdf', orientation='portrait', format='pdf')
    
    ##############
    ###Load absorption cross-sections of sulfite, bisulfite, from Beyad+2022 and XXXX as synthesized by Ranjan+2022
    ##############
    sulfite_xc_wav, sulfite_xc_molardecadicabsorption=np.genfromtxt('./from-ranjan2022/sulfite_spectrum.dat', skip_header=2, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, M**-1 cm**-1

    bisulfite_xc_wav, bisulfite_xc_molardecadicabsorption=np.genfromtxt('./from-ranjan2022/bisulfite_spectrum.dat', skip_header=2, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, M**-1 cm**-1
    
    # ###Plot to ensure the spectra have been properly loaded.
    # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # ax.plot(sulfite_xc_wav, sulfite_xc_molardecadicabsorption, linewidth=4, linestyle='-', marker='d', color='red', label=r'SO$_3^{2-}$')    
    # ax.plot(bisulfite_xc_wav, bisulfite_xc_molardecadicabsorption, linewidth=4, linestyle='-', marker='o', color='blue', label=r'HSO$_3^{-}$')    
    # ax.set_yscale('log')
    # ax.set_ylabel(r'Molar Decadic Absorption Coefficient (M$^{-1}$cm$^{-1}$)', fontsize=14)
    # ax.legend(ncol=1, loc='best')    

    # ax.set_xscale('linear')
    # ax.set_xlabel('Wavelength (nm)', fontsize=14)
    # ax.set_xlim([200., 330.])
    # plt.savefig('./Plots/S_IV_XCs.pdf', orientation='portrait', format='pdf')
    
    #Convert cross-sections from M**-1 cm**-1 to cm^2/molecule. Derive from two forms of Beer's law, or see p.59 of Lakowicz2006.
    sulfite_xc_xc=sulfite_xc_molardecadicabsorption*np.log(10.0)*1.0E3/(6.02E23) #cm**2 
    bisulfite_xc_xc=bisulfite_xc_molardecadicabsorption*np.log(10.0)*1.0E3/(6.02E23) #cm**2
   
    ##############
    ###Set quantum yields of sulfite, bisulfite photolysis
    ##############
    ###Sulfite
    #Data sources: Sauer, Crowell & Shkrob 2004, Lian et al. 2006, Li et al. 2012
    sulfite_qy_wav=np.array([193.0, 200.0, 248.0, 253.7])
    sulfite_qy_qy=np.array([0.39, 0.23, 0.11, 0.116])
    
    if longwaveQY0:
        qy_sulfite_func=interp.interp1d(sulfite_qy_wav,sulfite_qy_qy, kind='linear', fill_value=(sulfite_qy_qy[0], 0), bounds_error=False) #Extrapolate with last values in either direction.
    else:
        qy_sulfite_func=interp.interp1d(sulfite_qy_wav,sulfite_qy_qy, kind='linear', fill_value=(sulfite_qy_qy[0], sulfite_qy_qy[-1]), bounds_error=False) #Extrapolate with last values in either direction.

    
    
    ###Bisulfite
    #This is much harder to treat b/c much more limited data. Only one point at 214 nm (Fischer and Warneck 1996), which must be interpreted as an upper limit. 
    #Approach: (1) Convert the upper limit at 214 to an estimate under the assumption that the difference between absolute and net photolysis rates for SO3[2-] is the same as for HSO3-. (2) Assume the bisulfite photolysis curve matches the sulfite photolysis curve, just scaled to pass through the point at 214 nm. 
    bisulfite_qy_214_nm=(0.116/0.39)*0.19 #0.057
    scalefactor=bisulfite_qy_214_nm/qy_sulfite_func(213.9)
    def qy_bisulfite_func(wav):
        return scalefactor*qy_sulfite_func(wav)
    
    
    # ###Plot what these assumptions look like, together with data, as documentation 
    # ##NOTE: To make the plot for the paper, longwaveQY0=FALSE.
    # wavelengths=np.arange(180.0, 330.0, step=0.1) #wavelengths over which to plot.
    # ##The below gets the QY(>254 nm)=0 case.
    # qy_sulfite_min=np.interp(wavelengths, sulfite_qy_wav, sulfite_qy_qy, left=sulfite_qy_qy[0], right=0)
    # qy_bisulfite_min=qy_sulfite_min*scalefactor
    
    # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # ##First, plot sulfite and constraints, which are more numerous
    
    # ax.plot(wavelengths, qy_sulfite_func(wavelengths), linewidth=2, linestyle='-', color='red', label=r'SO$_3^{2-}$ (Max)')  
    # ax.plot(wavelengths, qy_bisulfite_func(wavelengths), linewidth=2, linestyle='-', color='blue', label=r'HSO$_3^{-}$ (Max)')  

    # ax.plot(wavelengths, qy_sulfite_min, linewidth=2, linestyle='--', color='red', label=r'SO$_3^{2-}$ (Min)')  
    # ax.plot(wavelengths, qy_bisulfite_min, linewidth=2, linestyle='--', color='blue', label=r'HSO$_3^{-}$ (Min)')  

    # ax.errorbar(np.array([193.0]), np.array([0.391]), yerr=0.011, color='red', marker='o', markersize=4, markeredgecolor='black', capsize=4.0, label='Sauer+2004')
    # ax.errorbar(np.array([200.0]), np.array([0.231]), yerr=0.023, color='red', marker='s', markersize=4, markeredgecolor='black', capsize=4.0, label='Lian+2006')
    # ax.errorbar(np.array([248.0]), np.array([0.108]), yerr=0.002, color='red', marker='o', markersize=4, markeredgecolor='black', capsize=4.0, label='Sauer+2004')
    # ax.errorbar(np.array([253.7]), np.array([0.39]), yerr=0.04, color='hotpink', marker='v', markersize=4, markeredgecolor='black', capsize=4.0, label='Fischer+1996')
    # ax.errorbar(np.array([253.7]), np.array([0.116]), yerr=0.002, color='red', marker='d', markersize=4, markeredgecolor='black', capsize=4.0, label='Li+2012')


    # ax.errorbar(np.array([213.9]), np.array([0.19]), yerr=0.03, color='cyan', marker='v', markersize=4, markeredgecolor='black', capsize=4.0, label='Fischer+1996')
    # ax.errorbar(np.array([213.9]), np.array([bisulfite_qy_214_nm]), yerr=bisulfite_qy_214_nm*np.sqrt((0.04/0.39)**2.0 +  (0.03/0.19)**2.0 + (0.002/0.116)**2.0), color='black', marker='v', markersize=4, markeredgecolor='black', capsize=4.0, label='Scaled estimate at 213.9 nm') #assuming uncorrelated Gaussian errors, errors add in quadrature
    # ax.set_yscale('linear')
    # ax.set_ylabel(r'$\Phi$ (Quantum Yield of Photolysis)', fontsize=14)
    # ax.set_xscale('linear')
    # ax.set_xlabel('Wavelength (nm)', fontsize=14)
    # ax.set_xlim([190., 330.])
    # ax.set_ylim([0., 1.])
    # ax.legend(ncol=1, loc='best')    
    # # plt.savefig('./Plots/S_IV_Photolysis.pdf', orientation='portrait', format='pdf')
    
    ##############
    ###Load transmission of early Earth waters, from Ranjan et al. 2022
    ##############
    prebioticocean_low_wav, prebioticocean_low_abs, prebioticocean_low_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebiotic_ocean_low.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1
    prebioticocean_high_wav, prebioticocean_high_abs, prebioticocean_high_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebiotic_ocean_high.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1

    prebioticfreshwaterlake_low_wav, prebioticfreshwaterlake_low_abs, prebioticfreshwaterlake_low_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebioticfreshwaterlake_nosIV_low.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1
    prebioticfreshwaterlake_high_wav, prebioticfreshwaterlake_high_abs, prebioticfreshwaterlake_high_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebioticfreshwaterlake_nosIV_high.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1

    prebioticcarbonatelake_low_wav, prebioticcarbonatelake_low_abs, prebioticcarbonatelake_low_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebioticcarbonatelake_nosIV_low.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1
    prebioticcarbonatelake_high_wav, prebioticcarbonatelake_high_abs, prebioticcarbonatelake_high_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebioticcarbonatelake_nosIV_high.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1

    prebioticferrouslake_low_wav, prebioticferrouslake_low_abs, prebioticferrouslake_low_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebioticferrouslake_nosIV_low.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1
    prebioticferrouslake_high_wav, prebioticferrouslake_high_abs, prebioticferrouslake_high_abs_err=np.genfromtxt('./from-ranjan2022/syntheticmolarabsorbance_prebioticferrouslake_nosIV_high.dat', skip_header=2, skip_footer=0,usecols=(0,1,2), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1

    # ###Plot to ensure the spectra have been properly loaded.   
    # ##Ocean
    # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # ax.errorbar(prebioticocean_low_wav, prebioticocean_low_abs, yerr=prebioticocean_low_abs_err, linewidth=2, linestyle='-', marker='o', color='black', label=r'Low-Absorption Endmember')    
    # ax.errorbar(prebioticocean_high_wav, prebioticocean_high_abs, yerr=prebioticocean_high_abs_err, linewidth=2, linestyle='--', marker='d', color='red', label=r'High-Absorption Endmember')
    # ax.set_yscale('log')
    # ax.set_ylabel(r'Linear Decadic Absorption Coefficient (cm$^{-1}$)', fontsize=14)
    # ax.legend(ncol=1, loc='best')  
    # ax.set_ylim([1E-3, 1.0E1])
    # ax.set_xscale('linear')
    # ax.set_xlabel('Wavelength (nm)', fontsize=14)
    # ax.set_xlim([200., 330.])
    # ax.legend(ncol=1, loc='best', fontsize=14)
    # plt.savefig('./Plots/prebioticocean_abs.pdf', orientation='portrait', format='pdf')
    # ##Looks good compared to Figure 1 of Ranjan+2022a
    
    # ##Freshwater lake
    # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # ax.errorbar(prebioticfreshwaterlake_low_wav, prebioticfreshwaterlake_low_abs, yerr=prebioticfreshwaterlake_low_abs_err, linewidth=2, linestyle='-', marker='o', color='black', label=r'Low-Absorption Endmember')    
    # ax.errorbar(prebioticfreshwaterlake_high_wav, prebioticfreshwaterlake_high_abs, yerr=prebioticfreshwaterlake_high_abs_err, linewidth=2, linestyle='--', marker='d', color='red', label=r'High-Absorption Endmember')
    # ax.set_yscale('log')
    # ax.set_ylabel(r'Linear Decadic Absorption Coefficient (cm$^{-1}$)', fontsize=14)
    # ax.legend(ncol=1, loc='best')  
    # ax.set_ylim([1E-3, 1.0E0])
    # ax.set_xscale('linear')
    # ax.set_xlabel('Wavelength (nm)', fontsize=14)
    # ax.set_xlim([200., 330.])
    # ax.legend(ncol=1, loc='best', fontsize=14)
    # plt.savefig('./Plots/prebioticfreshwaterlake_abs.pdf', orientation='portrait', format='pdf')
    
    # ##Carbonate lake
    # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # ax.errorbar(prebioticcarbonatelake_low_wav, prebioticcarbonatelake_low_abs, yerr=prebioticcarbonatelake_low_abs_err, linewidth=2, linestyle='-', marker='o', color='black', label=r'Low-Absorption Endmember')    
    # ax.errorbar(prebioticcarbonatelake_high_wav, prebioticcarbonatelake_high_abs, yerr=prebioticcarbonatelake_high_abs_err, linewidth=2, linestyle='--', marker='d', color='red', label=r'High-Absorption Endmember')
    # ax.set_yscale('log')
    # ax.set_ylabel(r'Linear Decadic Absorption Coefficient (cm$^{-1}$)', fontsize=14)
    # ax.legend(ncol=1, loc='best')  
    # ax.set_ylim([1E-3, 1.0E2])
    # ax.set_xscale('linear')
    # ax.set_xlabel('Wavelength (nm)', fontsize=14)
    # ax.set_xlim([200., 330.])
    # ax.legend(ncol=1, loc='best', fontsize=14)
    # plt.savefig('./Plots/prebioticcarbonatelake_abs.pdf', orientation='portrait', format='pdf')
    
    # # # ##Ferrous lake
    # # # fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    # # # ax.errorbar(prebioticferrouslake_low_wav, prebioticferrouslake_low_abs, yerr=prebioticferrouslake_low_abs_err, linewidth=2, linestyle='-', marker='o', color='black', label=r'Low-Absorption Endmember')    
    # # # ax.errorbar(prebioticferrouslake_high_wav, prebioticferrouslake_high_abs, yerr=prebioticferrouslake_high_abs_err, linewidth=2, linestyle='--', marker='d', color='red', label=r'High-Absorption Endmember')
    # # # ax.set_yscale('log')
    # # # ax.set_ylabel(r'Molar Decadic Absorption Coefficient (M$^{-1}$cm$^{-1}$)')
    # # # ax.legend(ncol=1, loc='best')  
    # # # ax.set_ylim([1E-3, 1.0E2])
    # # # ax.set_xscale('linear')
    # # # ax.set_xlabel('Wavelength (nm)')
    # # # ax.set_xlim([200., 300.])
    # # # ax.legend(ncol=1, loc='best')
    # # # plt.savefig('./Plots/prebioticferrouslake_abs.pdf', orientation='portrait', format='pdf')
    
    ##############
    ###Some literature spectra, unaltered, as a sensivity test.
    ##############

    ###Load purest modern ocean water optical properities, from Smith & Baker 1981. These are VERY low in organics and hence very UV transparent. Also, <300 nm the numbers are essentially made up (pure extrapolation, nothing physical).    
    ##NAPERIAN!!!
    df2=pd.read_excel('./smith_baker.xlsx', sheet_name='Sheet1', skiprows=3)
    pureoceans_smith1981_wav=np.nan_to_num(df2['lambda']) #nm 
    pureoceans_smith1981_Kd=np.nan_to_num(df2['Ksww'])*1.0E-2 #m**-1 converted to cm**-1
    pureoceans_smith1981_a=np.nan_to_num(df2['aw'])*1.0E-2 #m**-1 converted to cm**-1
    pureoceans_smith1981_b_sw=np.nan_to_num(df2['bswm'])*1.0E-2 #m**-1 converted to cm**-1, salt water
    pureoceans_smith1981_b_fw=np.nan_to_num(df2['bfwm'])*1.0E-2 #m**-1 converted to cm**-1, fresh water
    
    ###Load pure water, from Quickenden & Irving
    purewater_wav, purewater_scat, purewater_abs, purewater_abs_err=np.genfromtxt('./pure_h2o.dat', skip_header=2, skip_footer=0,usecols=(0,2,3,4), unpack=True) #units: nm, M**-1 cm**-1, M**-1 cm**-1

    ##############
    ###Select which absorbance to use
    ##############    
    if water_abs_choice=='prebiotic_ocean_low':
        solution_abs_wav=prebioticocean_low_wav
        solution_abs_abs=prebioticocean_low_abs
        solution_abs_abserr=prebioticocean_low_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_ocean_high':
        solution_abs_wav=prebioticocean_high_wav
        solution_abs_abs=prebioticocean_high_abs
        solution_abs_abserr=prebioticocean_high_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_freshwaterlake_low':
        solution_abs_wav=prebioticfreshwaterlake_low_wav
        solution_abs_abs=prebioticfreshwaterlake_low_abs
        solution_abs_abserr=prebioticfreshwaterlake_low_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_fw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_freshwaterlake_high':
        solution_abs_wav=prebioticfreshwaterlake_high_wav
        solution_abs_abs=prebioticfreshwaterlake_high_abs
        solution_abs_abserr=prebioticfreshwaterlake_high_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_fw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_carbonatelake_low':
        solution_abs_wav=prebioticcarbonatelake_low_wav
        solution_abs_abs=prebioticcarbonatelake_low_abs
        solution_abs_abserr=prebioticcarbonatelake_low_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_carbonatelake_high':
        solution_abs_wav=prebioticcarbonatelake_high_wav
        solution_abs_abs=prebioticcarbonatelake_high_abs
        solution_abs_abserr=prebioticcarbonatelake_high_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_ferrouslake_low':
        solution_abs_wav=prebioticferrouslake_low_wav
        solution_abs_abs=prebioticferrouslake_low_abs
        solution_abs_abserr=prebioticferrouslake_low_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='prebiotic_ferrouslake_high':
        solution_abs_wav=prebioticferrouslake_high_wav
        solution_abs_abs=prebioticferrouslake_high_abs
        solution_abs_abserr=prebioticferrouslake_high_abs_err
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='modern_ocean_smithbaker':
        solution_abs_wav=pureoceans_smith1981_wav
        solution_abs_abs=pureoceans_smith1981_a/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
        solution_abs_abserr=pureoceans_smith1981_a*0.25 #25% error based on page 184 of Smith & Baker 1981. Overestimated slightly for >300 nm (really +25%, -5%), underestimated for <300 nm (made up).
        solution_scat_wav=pureoceans_smith1981_wav
        solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    elif water_abs_choice=='pure_water': #Quickenden & Irving, scattering from Krocken+2014 as mediated by Ranjan+2022 ###If using this, make sure RT only goes out to 320 nm to match wavelength limits.
        solution_abs_wav=purewater_wav
        solution_abs_abs=purewater_abs
        solution_abs_abserr=purewater_abs_err
        solution_scat_wav=purewater_wav
        solution_scat_scat=purewater_scat
    else:
        print('Invalid entry for water_abs_choice')
        exit()
    
    
    ##############
    ###Prepare for depth-dependent RT calculation
    ##############
    
    ###Establish depth variables.
    #Establish depth scale with depth. 
    binwidth=1.0 #cm
    
    if d_body>1.0E+4: #cm
        d_body=1.0E+4 #No point running calculation below 100 m (1E4 cm) because water should extinguish it by then. 
    depth_bin_left=np.arange(0.0, d_body+binwidth, step=binwidth) #cm
    depth_bin_right=depth_bin_left+binwidth
    depth_bin_centers=0.5*(depth_bin_left+depth_bin_right)
    
    #Establish array to hold depth scale.
    bisulfite_reac_rate_const_depth=0.0*depth_bin_centers
    sulfite_reac_rate_const_depth=0.0*depth_bin_centers
    
    ###Finalize spectral variables
    ##Add sulfite, bisulfite into the molar absorbances.
    #Get sulfite, bisulfite concentration based on dissociation equilibria
    pKa_1=1.87 + -0.50*I**0.5 + 0.31*I #Millero+1989, p381, 25C value
    pKa_2=7.12 + -1.052*I**0.5 + 0.36*I #Millero+1989, p381. T-independent 5-25C.
    Ka_1=10**(-pKa_1)
    Ka_2=10**(-pKa_2)
    conc_H=10**(-pH)
    alpha_HSO3=(conc_H/Ka_1 + Ka_2/conc_H+1.0)**-1 #Equation 17 of Zhang & Millero
    alpha_SO3=(conc_H**2.0/(Ka_1*Ka_2)+conc_H/Ka_2+1.0)**-1  #Equation 18 of Zhang & Millero
    
    conc_HSO3=alpha_HSO3*conc_s_iv
    conc_SO3 = alpha_SO3*conc_s_iv
    
    #Add sulfite, bisulfite in. Note that we don't have any uncertainty estimate on sulfite, bisulfite absorption anyway (not reported) so might as well interpolate those ones. 
    solution_abs_abs=solution_abs_abs + conc_HSO3*np.interp(solution_abs_wav, bisulfite_xc_wav, bisulfite_xc_molardecadicabsorption, left=0., right=0.) + conc_SO3*np.interp(solution_abs_wav, sulfite_xc_wav, sulfite_xc_molardecadicabsorption, left=0., right=0.)
    #No effect on error since we don't know sulfite, bisulfite absorption errors (assume perfectly known).
    
    #Project onto common wavelength scale with spectral data. 
    solution_abs_interp, solution_abs_error_interp=linear_interp(uv_wav, solution_abs_wav, solution_abs_abs, solution_abs_abserr)
    # solution_abs_interp, solution_abs_error_interp=spectres.spectres(uv_wav, solution_abs_wav, solution_abs_abs, spec_errs=solution_abs_abserr, fill=0.0) #Try this, verify any dependency.
    solution_scat_interp=np.interp(uv_wav, solution_scat_wav, solution_scat_scat, left=0, right=0)
    

    sulfite_xc_xc_interp=np.interp(uv_wav, sulfite_xc_wav, sulfite_xc_xc, left=0, right=0)
    bisulfite_xc_xc_interp=np.interp(uv_wav, bisulfite_xc_wav, bisulfite_xc_xc, left=0, right=0)
    
    # plt.plot(uv_wav, sulfite_xc_xc_interp, color='red')
    # plt.plot(uv_wav, bisulfite_xc_xc_interp, color='black')
    # plt.yscale('log')
    
    ##############
    ###Estimate K_d
    ###3 ways, varying in complexity.
    ##############
    a=np.log(10.0)*solution_abs_interp#Need to adjust for Naperian vs Decadic absorption. 
    b=np.log(10.0)*solution_scat_interp#Need to adjust for Naperian vs Decadic absorption.
    b_b=0.5*b #Assume equal forward and backscattering
    
    mu_0_aq=(1.0/uv_surface_flux_aq)*(np.cos(theta_sun_aq)*uv_direct_flux_aq + 0.86*uv_diffuse_flux_aq)#Average cosine for downwelling radiation inside the water (follows Kirk 1984)

        
    if K_d_method=='Ranjan2022': 
        #Use method of Ranjan+2022
        #This method is very simplistic and just approximates K_d by the molar absorption, adjusted by the slant angle. Improves on Ranjan+2022 by accounting for refraction of the solar radiation. 
        
        K_d=a/np.cos(theta_sun_aq)
    elif K_d_method=='Morel2007':
        #Use approximation of Morel+2007 equation 3
        #This ignores slant angle correction...
        K_d=(1.0395/mu_0_aq)*(a+b_b)
    elif K_d_method=='Morel1991':
        ###Use method of Morel 1991 
        #Calculate K_d
        K_d=(a/mu_0_aq)*(1.0+(0.425*mu_0_aq-0.19)*(b/a))**0.5
    
    ###Do loop over depth
    for ind in range(0, len(depth_bin_centers)):
        E_d=uv_surface_flux_aq*np.exp(-depth_bin_centers[ind]*K_d) #cm**-2 s**-1 nm**-1
        
        E_dot_d=E_d*K_d/a #cm**-2 s**-1 nm**-1
        
        #Impose ceiling on E_dot_d. Scalar irradiance should not exceed 2x the deprojected planar irradiance in the two-stream approximation. But it might since our formalism assumes absorption-dominated medium and we can't rule out scattering-dominated medium for early earth at long wavelengths. So ad hoc correction.
        E_dot_d_max=(E_d/mu_0_aq)*2.5 
        inds=np.where(E_dot_d>E_dot_d_max)
        E_dot_d[inds]=E_dot_d_max[inds]
    
        bisulfite_spectral_reac_rate_const=qy_bisulfite_func(uv_wav)*bisulfite_xc_xc_interp*E_dot_d #Units: s**-1 nm**-1
        sulfite_spectral_reac_rate_const=qy_sulfite_func(uv_wav)*sulfite_xc_xc_interp*E_dot_d #Units: s**-1 nm**-1
        
        bisulfite_reac_rate_const_depth[ind]=np.sum(uv_wav_widths*bisulfite_spectral_reac_rate_const) #Units: s**-1
        sulfite_reac_rate_const_depth[ind]=np.sum(uv_wav_widths*sulfite_spectral_reac_rate_const) #Units: s**-1
    
    reac_rate_depth=conc_HSO3*bisulfite_reac_rate_const_depth+conc_SO3*sulfite_reac_rate_const_depth #Units: M s**-1
    reac_rate_depth_cgs=reac_rate_depth*M2cgs #Units: cm**-3 s**-1
    reac_rate_depth_cgs_tot=np.sum(binwidth*reac_rate_depth_cgs) #Units: cm**-2 s**-1
    
    # fig, ax=plt.subplots(1, figsize=(8, 6.), sharex=False)    
    # ax.plot(reac_rate_depth_cgs, depth_bin_centers)
    # ax.set_xscale('linear')
    # ax.set_xlim([0, np.max(reac_rate_depth_cgs)])
    # ax.invert_yaxis()

    return reac_rate_depth_cgs_tot




def nitrate_photolysis_rate_modernocean_fullcalc(K_d_method):
    """
    
    """
    ##############
    ###Load early Earth UV from Ranjan & Sasselov 2017 (or should we use Ranjan et al. 2022?)
    ##############    
    uv_file='./rugheimer_earth_modern.dat'  
    theta_sun=60.0*deg2rad #For the source RT calculation, the solar zenith angle is 60 degrees, hard-coded.
    mu_1=1.0/np.sqrt(3.0) #For the source RT calculation, mu_1=1/sqrt(3) (quadrature closure)
    
    uv_wav_left, uv_wav_right, uv_wav, uv_toa_intensity, uv_surface_flux, uv_surface_intensity, uv_diffuse_intensity, uv_direct_intensity=np.genfromtxt(uv_file, skip_header=1, skip_footer=0, usecols=(0,1,2,3,4,6,7,8), unpack=True) #From Ranjan & Sasselov 2017, the GitHub. 
    
    ###Convert from ergs to photons
    uv_toa_intensity *= uv_wav/hc
    uv_surface_flux *= uv_wav/hc
    uv_surface_intensity *= uv_wav/hc
    uv_diffuse_intensity *= uv_wav/hc
    uv_direct_intensity *= uv_wav/hc
    
    ###Restrict data to 290 to 340 nm. Photolysis takes place up to 340 nm (Zafiriou+1979), below 290 nm blocked by ozone layer ()
    inds=np.where((uv_wav>=290.0) & (uv_wav <=340.0))
    
    uv_wav=uv_wav[inds]
    uv_wav_left=uv_wav_left[inds]
    uv_wav_right=uv_wav_right[inds]
    uv_toa_intensity=uv_toa_intensity[inds]
    uv_surface_flux=uv_surface_flux[inds]
    uv_surface_intensity=uv_surface_intensity[inds]
    uv_diffuse_intensity=uv_diffuse_intensity[inds]
    uv_direct_intensity=uv_direct_intensity[inds]
    
    ###Obtain some useful secondary variables
    uv_wav_widths=uv_wav_right-uv_wav_left
    uv_diffuse_flux=mu_1*uv_diffuse_intensity
    uv_direct_flux=(np.cos(theta_sun))*uv_direct_intensity
    
    ###Calculate transmission through air/water interface, following Morel 1991
    uv_diffuse_flux_aq=uv_diffuse_flux*(1.0-0.066) #constant sky reflectance
    
    theta_sun_aq=get_aqueous_angle_snell(theta_sun) #"nadir angle" after refraction of dreict stream
    uv_direct_flux_aq=uv_direct_flux*(1.0-get_reflectance_fresnel(theta_sun, theta_sun_aq))
    uv_surface_flux_aq=uv_diffuse_flux_aq+uv_direct_flux_aq#This is the E_d(0-) term.
    
    ##############
    ###Load absorption cross-sections of nitrate (Ranjan+2022)
    ##############

    data_wav={} #dict to hold wavelength scale for measured molar absorptivities (nm)
    data_abs={} #dict to hold molar absorptivities for measured molar absorptivities (M**-1 cm**-1.)
    
    ###From Gabi/Corinna's work
    df=pd.read_excel('./2020-08-23_Extinktion-coefficients-Sukrit-V7_mod.xlsx', sheet_name='Tabelle1') 
    data_wav['LK']=np.nan_to_num(df['Wavelength']) #nm LK=Lozano-Kufner
    
    data_abs['NaNO3_LK']=np.nan_to_num(df['NaNO3']) #M**-1 cm**-1
    
    nitrate_xc_wav=data_wav['LK']
    nitrate_xc_molardecadicabsorption=data_abs['NaNO3_LK']
    nitrate_xc_xc=nitrate_xc_molardecadicabsorption*np.log(10.0)*1.0E3/(6.02E23) #cm**2 
   

    ##############
    ###Get quantum yield
    ##############
    #From table 2 of Mack & Bolton. This is best for comparison to Zafiriou+1979b as it best matches the actual process they are trying to constrain. 
    nitrate_qy_wav=np.array([282.0, 305.0, 310.0, 313.0]) #nm
    nitrate_qy_qy=np.array([0.024, 0.0059, 0.0065, 0.0063]) #dimensionless
    
    qy_nitrate_func=interp.interp1d(nitrate_qy_wav,nitrate_qy_qy, kind='linear', fill_value=(nitrate_qy_qy[0], nitrate_qy_qy[-1]), bounds_error=False) #Extrapolate with last values in either direction.

    
    ##############
    ###Load absorption of modern oceans, from Smith & Baker 1981. Well, in reality, the purest modern ocean waters -- highly UV-transparent. 
    ##############

    ###Load purest modern ocean water optical properities, from Smith & Baker 1981. These are VERY low in organics and hence very UV transparent. Also, <300 nm the numbers are essentially made up (pure extrapolation, nothing physical).    
    ##NAPERIAN!!!
    df2=pd.read_excel('./smith_baker.xlsx', sheet_name='Sheet1', skiprows=3)
    pureoceans_smith1981_wav=np.nan_to_num(df2['lambda']) #nm 
    pureoceans_smith1981_Kd=np.nan_to_num(df2['Ksww'])*1.0E-2 #m**-1 converted to cm**-1
    pureoceans_smith1981_a=np.nan_to_num(df2['aw'])*1.0E-2 #m**-1 converted to cm**-1
    pureoceans_smith1981_b_sw=np.nan_to_num(df2['bswm'])*1.0E-2 #m**-1 converted to cm**-1, salt water
    pureoceans_smith1981_b_fw=np.nan_to_num(df2['bfwm'])*1.0E-2 #m**-1 converted to cm**-1, fresh water
    
    solution_abs_wav=pureoceans_smith1981_wav
    solution_abs_abs=pureoceans_smith1981_a/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    solution_abs_abserr=pureoceans_smith1981_a*0.25 #25% error based on page 184 of Smith & Baker 1981. Overestimated slightly for >300 nm (really +25%, -5%), underestimated for <300 nm (made up).
    solution_scat_wav=pureoceans_smith1981_wav
    solution_scat_scat=pureoceans_smith1981_b_sw/np.log(10.0) #have to convert from Naperian to decadic, for consistency with rest of inputs. 
    
    
    ##############
    ###Prepare for depth-dependent RT calculation
    ##############
    
    ###Establish depth variables.
    #Establish depth scale with depth. 
    binwidth=1.0 #cm
    d_body=1.0E+4 #cm
    
    depth_bin_left=np.arange(0.0, d_body+binwidth, step=binwidth) #cm
    depth_bin_right=depth_bin_left+binwidth
    depth_bin_centers=0.5*(depth_bin_left+depth_bin_right)
    
    #Establish array to hold depth scale.
    nitrate_reac_rate_const_depth=0.0*depth_bin_centers
    
    ###Finalize spectral variables
    ##Add nitrate into the molar absorbances.
    conc_nitrate=2.0E-6
    solution_abs_abs=solution_abs_abs + conc_nitrate*np.interp(solution_abs_wav, nitrate_xc_wav, nitrate_xc_molardecadicabsorption, left=0., right=0.)
    
    #Project onto common wavelength scale with spectral data. 
    solution_abs_interp, solution_abs_error_interp=linear_interp(uv_wav, solution_abs_wav, solution_abs_abs, solution_abs_abserr)
    # solution_abs_interp, solution_abs_error_interp=spectres.spectres(uv_wav, solution_abs_wav, solution_abs_abs, spec_errs=solution_abs_abserr, fill=0.0) #Try this, verify any dependency.
    solution_scat_interp=np.interp(uv_wav, solution_scat_wav, solution_scat_scat, left=0, right=0)
    

    nitrate_xc_xc_interp=np.interp(uv_wav, nitrate_xc_wav, nitrate_xc_xc, left=0, right=0)
    
    
    ##############
    ###Estimate K_d
    ###3 ways, varying in complexity.
    ##############
    a=np.log(10.0)*solution_abs_interp#Need to adjust for Naperian vs Decadic absorption. 
    b=np.log(10.0)*solution_scat_interp#Need to adjust for Naperian vs Decadic absorption.
    b_b=0.5*b #Assume equal forward and backscattering
    
    mu_0_aq=(1.0/uv_surface_flux_aq)*(np.cos(theta_sun_aq)*uv_direct_flux_aq + 0.86*uv_diffuse_flux_aq)#Average cosine for downwelling radiation inside the water (follows Kirk 1984)

        
    if K_d_method=='Ranjan2022': 
        #Use method of Ranjan+2022
        #This method is very simplistic and just approximates K_d by the molar absorption, adjusted by the slant angle. Improves on Ranjan+2022 by accounting for refraction of the solar radiation. 
        
        K_d=a/np.cos(theta_sun_aq)
    elif K_d_method=='Morel2007':
        #Use approximation of Morel+2007 equation 3
        #This ignores slant angle correction...
        K_d=(1.0395/mu_0_aq)*(a+b_b)
    elif K_d_method=='Morel1991':
        ###Use method of Morel 1991 
        #Calculate K_d
        K_d=(a/mu_0_aq)*(1.0+(0.425*mu_0_aq-0.19)*(b/a))**0.5
    
    ###Do loop over depth
    for ind in range(0, len(depth_bin_centers)):
        E_d=uv_surface_flux_aq*np.exp(-depth_bin_centers[ind]*K_d) #cm**-2 s**-1 nm**-1
        
        E_dot_d=E_d*K_d/a #cm**-2 s**-1 nm**-1
        
        #Impose ceiling on E_dot_d. Scalar irradiance should not exceed 2x the deprojected planar irradiance in the two-stream approximation. But it might since our formalism assumes absorption-dominated medium and we can't rule out scattering-dominated medium for early earth at long wavelengths. So ad hoc correction.
        E_dot_d_max=(E_d/mu_0_aq)*2.0 
        inds=np.where(E_dot_d>E_dot_d_max)
        E_dot_d[inds]=E_dot_d_max[inds]
    
        nitrate_spectral_reac_rate_const=qy_nitrate_func(uv_wav)*nitrate_xc_xc_interp*E_dot_d #Units: s**-1 nm**-1
        
        nitrate_reac_rate_const_depth[ind]=np.sum(uv_wav_widths*nitrate_spectral_reac_rate_const) #Units: s**-1
    reac_rate_depth_zaf_tot=np.sum(binwidth*nitrate_reac_rate_const_depth)
    # reac_rate_depth=conc_nitrate*nitrate_reac_rate_const_depth#Units: M s**-1
    # reac_rate_depth_cgs=reac_rate_depth*M2cgs #Units: cm**-3 s**-1
    # reac_rate_depth_cgs_tot=np.sum(binwidth*reac_rate_depth_cgs) #Units: cm**-2 s**-1
    # reac_rate_depth_zaf_tot=reac_rate_depth_cgs_tot*1/N_avogadro*(m2cm)**2.0*year2s #convert from cm**-2 s**-1 to mol m**-2 year**-1.
    
    print(nitrate_reac_rate_const_depth[0])
    print(reac_rate_depth_zaf_tot)







# sulfite_photolysis_rate(1.0E-6, 1.0E+2, 6.0, 1.0E-3, 'prebiotic_carbonatelake_low', False, 'Morel1991')
# calculate_colint_nitrate_photolysis_modernocean()
# nitrate_photolysis_rate_modernocean_fullcalc('Morel1991')