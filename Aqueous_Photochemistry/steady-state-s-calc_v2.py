#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This version changes the way figures are plotted.
"""
###############################################################################
###Control switches
###############################################################################
singlescenariocalc=False #whether to calculate steady-stateS[IV] in particular scenario.
plotlossmechanisms=False #whether to make plot showing different loss mechanisms.DEFUNCT
plotsupplymechanisms=False #whether to make plot showing different supply mechanisms.

plotsteadystatecalc_fso2_ocean=True #whether to make plot showing steady-state calc for specific scenario as function of pSO2 for ocean
plotsteadystatecalc_fso2_lake=True #whether to make plot showing steady-state calc for specific scenario as function of pSO2 for terrestrial waters
plotsteadystatecalc_fso2_lake_maxsulfite=False #whether to make plot showing controls on sulfite concentration.
plotsteadystatecalc_fso2_lake_depths=False #whether to make plot showing steady-state calc for different depths


###############################################################################
###Constants and packages
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
import pdb
import scipy.integrate
import scipy.optimize
from scipy import interpolate as interp
from matplotlib.pyplot import cm
from scipy.optimize import fsolve
import photolysis_coefficient_calculation as photcalc

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
atm2bar=1.0/0.9869 # 1 atm in bar
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
deg2rad=np.pi/180.
bar2barye=1.e+6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
m2cm=1.0E+2 #1 m in cm
metricton2kg=1000. #1 metric ton in kg
year2s=3.154e7 #1 year in seconds
day2s=86400. #1 day in seconds

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
R_earth=6371.*km2m#radius of earth in m
R_sun=69.63e9 #radius of sun in cm
AU=1.496e13#1AU in cm
R=1.987e-3# kcal K**-1 mol**-1

#solution specific conversion
M2cgs=(6.02e23)*1.e-3 #mol/L * particles/mol * L/cm**3 = molecules cm**-3
cgs2M=M2cgs**-1 #molecules cm**-3 to mol/L


###############################################################################
###Chemical and Geological Parameters for marine, terrestrial water calculations
###############################################################################

###########
###Planet/atmosphere parametrs
###########

##General assumptions, consistent with assumptions in photochemical modeling
P_surf_0=1.0 #bar
T_surf_0=288.0 #K
mr_co2_surf_0=0.1
v_dep_so2_0=1.0 #cm/s
k_h2o_0=2.0E-6 #s**-1 (Hu+2012)

##This block comes from modeling by Sangita
phi_volc_list=([0.1, 0.3, 1.0, 3.0, 10.0, 30.0]) #volcanism flux relative to present-day
mr_so2_list=np.array([7.E-12, 2.E-11, 1.E-10, 3.E-10, 1.E-9, 3.E-9]) #surface mixing ratios
so2_wetdep_list=np.array([3.E7, 9.E7, 4.E8, 2.E9, 5.E9, 1.E10]) #cm**-2 s**-1
N_h2o_list=np.array([4.486252e+22, 4.486253e+22, 4.486257e+22, 4.486294e+22, 4.486313e+22, 4.486327e+22]) #Total column of H2O in the model runs, cm**-2
p_o2_list=np.array([6.E-11, 5.E-11, 1.E-16, 5.E-19, 4.E-20, 5.E-21]) #surface mixing ratios, bar

###########
###Chemical parameters
#These parameters not specific to any body of water, but control S[IV] destruction rate.
###########

wetdep_method_0='model' #If 'model', takes it from calculation from MEAC. If 'calc-raindrop' calculate it from pSO2, pCO2. If 'calc-std', use "standard" value for Giorgi & Chameides 1985 (different from what Hu+2012 use, which may instead be from Seinfeld & Pandis)
K_d_method_0='Morel1991'# Options are Ranjan2022, Morel2007, Morel1991

###
#Uncertain chemical parameters (range of values)
###
##Maximize sulfite
T_sulfite_disprop_exp_max=5.0*year2s #From Meyer+1979, lower bound
conc_s_iv_exp_max=1.0 #M From LBNL tech report 1980.
disprop_rxn_order_max=4.0 #This maximizes sulfite abundance because the reaction order is high, meaning reaction is suppressed at low concentrations
longwaveQY0_max=True #This maximizes sulfite abundance because it does not photolyze at longer wavelengths.

##Minimize sulfite
T_sulfite_disprop_exp_min=1.0*year2s #Meyer+1982, lower bound. 
conc_s_iv_exp_min=1.0#M From LBNL tech report 1980.
disprop_rxn_order_min=1.0 #This minimizes sulfite abundance because the reaction order is low, meaning reaction is maximized at low concentrations.
longwaveQY0_min=False #This minimizes sulfite abundance because it does photolyze at longer wavelengths.

###########
#Geological parameters: ocean
###########
d_ocean=3.8E5 #cm; from CRC, via Ranjan+2019.
effective_precip_rate_ocean_list=k_h2o_0*N_h2o_list*(18.02/6.02E23*1.0)*(year2s/m2cm) #m/year. For ocean, this should just be the model global rate. 
seepage_rate_ocean=0.0

###
#Uncertain ocean geological parameters (range of values)
###
    
##First block: maximize sulfite
water_molabs_ocean_max='prebiotic_ocean_high' #Maximizes UV shielding from oceanic constitutents

pH_ocean_max=6.25 #Minimize ocean pH to maximize sulfite because more S[IV] is in bisulfite, which absorbs less. From Kadoya+2020, Krissansen-Totton+2018.
I_ocean_max=0.3 #consistent with specified water composition

##Second block: minimize sulfite    
water_molabs_ocean_min='prebiotic_ocean_low' #Minimizes UV shielding from oceanic constitutents

pH_ocean_min=9.0 #Maximize ocean pH to minimize sulfite because more S[IV] is in sulfite, which absorbs more. From Kadoya+2020, Krissansen-Totton+2018.
I_ocean_min=0.72 #Maximum allowed by the formalism we are using. #1.2 #consistent with specified water composition
    
###########
#Geological parameters: lakes/ponds
###########
d_lake=1.0E2 #cm
d_lake_shallow=1.0E1 #cm. For the sole case of an extremely shallow system, similar to Don Juan Pond.
evaporationrate_lake=-0.12+0.06*(T_surf_0-273.15)# m/year. expression from Boyd+1985 via Pearce+2017
    
###Carbonate lake
##First block: maximize sulfite
water_molabs_carbonatelake_max='prebiotic_carbonatelake_high' #Maximizes UV shielding

pH_carbonatelake_max=6.5 #Minimize pH to maximize sulfite because more S[IV] is in bisulfite, which absorbs less. Toner & Catling 2020, via Ranjan+2022
I_carbonatelake_max=0.1 #consistent with specified water composition

seepage_rate_carbonatelake_max=0.0
effective_precip_rate_carbonatelake_max=seepage_rate_carbonatelake_max+evaporationrate_lake

##Second block: minimize sulfite
water_molabs_carbonatelake_min='prebiotic_carbonatelake_low' #Minimizes UV shielding

pH_carbonatelake_min=9.0 #Maximize pH to minimize sulfite because more S[IV] is in sulfite, which absorbs more. Toner & Catling 2020, via Ranjan+2022
I_carbonatelake_min=0.72 #Maximum allowed by the formalism we are using. #6 #consistent with specified water composition

seepage_rate_carbonatelake_min=2.0
effective_precip_rate_carbonatelake_min=seepage_rate_carbonatelake_min+evaporationrate_lake
    
###Freshwater lake
pH_freshwaterlake=6.34 #Hao+2017
I_freshwaterlake=0.001 #consistent with specified water composition

##First block: maximize sulfite
water_molabs_freshwaterlake_max='prebiotic_freshwaterlake_high' #Maximizes UV shielding

seepage_rate_freshwaterlake_max=0.0
effective_precip_rate_freshwaterlake_max=seepage_rate_freshwaterlake_max+evaporationrate_lake

##Second block: minimize sulfite
water_molabs_freshwaterlake_min='prebiotic_freshwaterlake_low' #Minimizes UV shielding

seepage_rate_freshwaterlake_min=2.0
effective_precip_rate_freshwaterlake_min=seepage_rate_freshwaterlake_min+evaporationrate_lake


###############################################################################
###Sulfite supply
###############################################################################

#####
###Sulfite dry deposition
#####
def s_iv_dry_deposition_flux(P_surf, T_surf, mr_so2_surf, v_dep_so2):
    """
    Function to calculate column supply of S[+IV] to surface from atmospheric dry deposition alone. 
    
    Could just use n_so2_0 instead of P_0, T_0, mr_so2_0, but practically the latter is what is reported in the literature so this is better from comparison perspective
    ----------
    P_0 : surface pressure (bar)
    T_0 : surface temperature (K)
    mr_so2 : mixing ratio of SO2
    v_dep_so2 : dry deposition velocity of SO2 (cm/s)

    Returns
    -------
    column flux of S[+IV] (SO2, which instantly interconverts to sulfite), in cm**-2 s**-1
    """
    n_so2_surf=mr_so2_surf*(P_surf*bar2barye)/(k*T_surf) #cm**-3
    return n_so2_surf*v_dep_so2 #cm**-2 s**-1


#####
###Sulfite wet deposition
#####

###
#Raindrop pH
###
def raindrop_pH(P_surf,  mr_co2_surf, mr_so2_surf):
    """
    To estimate pH of raindrop given (effectively) pCO2 and assuming only their equilibria.
    """
    
    pCO2=P_surf*mr_co2_surf
    pSO2=P_surf*mr_so2_surf
    # T=290.0
    
    H_CO2=3.0E-2 #M/bar. From Sander (2015) via Ranjan+2018
    H_SO2=1.34 #M/bar. From Sander (2015) via Ranjan+2018
    # H_CO2=3.11E-2*np.exp(2423.0*(1.0/T-1.0/298.0))*1.0/(atm2bar) #M/bar. From Sander (2015) via Ranjan+2018
    # H_SO2=1.23*np.exp(3120.0*(1.0/T-1.0/298.0))*1.0/(atm2bar)#M/bar. From Sander (2015) via Ranjan+2018
    
    pKa_co2_1=6.35 #From Rumble 2017 via Ranjan+2018 
    pKa_co2_2=10.33 #From Rumble 2017 via Ranjan+2018   
    pKa_so2_1=1.86 #Neta+1985 via Ranjan+2018
    pKa_so2_2=7.2  #Neta+1985 via Ranjan+2018
    
    conc_h2co3_0=pCO2*H_CO2    
    Ka_co2_1_0=10.0**-pKa_co2_1
    Ka_co2_2_0=10.0**-pKa_co2_2
    
    conc_h2so3_0=pSO2*H_SO2    
    Ka_so2_1_0=10.0**-pKa_so2_1
    Ka_so2_2_0=10.0**-pKa_so2_2
    # Ka_so2_1_0=1.7E-2*np.exp(2090.0*(1.0/T-1.0/298.0))
    # Ka_so2_2_0=6.0E-8*np.exp(1120.0*(1.0/T-1.0/298.0))
    # pdb.set_trace
    
    def system(x, conc_h2co3, Ka_co2_1, Ka_co2_2, conc_h2so3, Ka_so2_1, Ka_so2_2):
        (log_conc_hco3, log_conc_co3, log_conc_hso3, log_conc_so3, pH)=x
        conc_hco3=10.0**log_conc_hco3
        conc_co3=10.0**log_conc_co3        
        conc_hso3=10.0**log_conc_hso3        
        conc_so3=10.0**log_conc_so3    
        conc_h=10**-pH
       
        eqn1=np.log(Ka_co2_1)-np.log(conc_h)-np.log(conc_hco3)+np.log(conc_h2co3)
        eqn2=np.log(Ka_co2_2)-np.log(conc_h)-np.log(conc_co3)+np.log(conc_hco3)
        eqn3=np.log(Ka_so2_1)-np.log(conc_h)-np.log(conc_hso3)+np.log(conc_h2so3)
        eqn4=np.log(Ka_so2_2)-np.log(conc_h)-np.log(conc_so3)+np.log(conc_hso3)
        eqn5=conc_h + -1.0*conc_hco3 + -2.0*conc_co3 + -1.0*conc_hso3 + -2.0*conc_so3
        return (eqn1, eqn2, eqn3, eqn4, eqn5)
    
    #construct sensible initial guess (pH=5)
    pH_guess=5.0
    conc_h_guess=10.0**-pH_guess
    conc_hco3_guess=conc_h2co3_0*Ka_co2_1_0/conc_h_guess
    conc_co3_guess=conc_hco3_guess*Ka_co2_2_0/conc_h_guess
    conc_hso3_guess=conc_h2so3_0*Ka_so2_1_0/conc_h_guess
    conc_so3_guess=conc_hso3_guess*Ka_so2_2_0/conc_h_guess    
    initial_guess=(np.log(conc_hco3_guess), np.log(conc_co3_guess), np.log(conc_hso3_guess), np.log(conc_so3_guess), pH_guess)
    (log_conc_hco3_0, log_conc_co3_0, log_conc_hso3_0, log_conc_so3_0,pH_0)=fsolve(system, initial_guess, args=(conc_h2co3_0, Ka_co2_1_0, Ka_co2_2_0, conc_h2so3_0, Ka_so2_1_0, Ka_so2_2_0), xtol=1.0E-6, maxfev=10000)
    
    conc_hso3_0=10.0**log_conc_hso3_0
    conc_so3_0=10.0**log_conc_so3_0
    effective_H_SO2_0=(conc_h2so3_0+conc_hso3_0+conc_so3_0)/pSO2
    # print(pH_0, effective_H_SO2_0)
    return (pH_0, effective_H_SO2_0)

# raindrop_pH(1.0, 3.5E-4, 30E-12) #Modern earth case. SO2 concentration of 237E-12 taken from Hu+2012 for modern Earth overall.
# #In a pollutant-free atmosphere, raindrop has a pH=5.6, for rCO2=350 ppm (Seinfeld & Pandis 2016 page 874). 30-260 ppt SO2 measured in troposphere (Seinfeld & Pandis 2016 page 26), we take 30 to represent pre-anthropogenic SO2. This reproduces raindrop pH=5.6
# raindrop_pH(1.0, 1.0E-1, 1.0E-10)
# raindrop_pH(1.0, 0.9, 9.0E-10)

###
#Estimate wet deposition rate. 
###
def s_iv_wet_deposition_flux(method, so2_wet_dep, N_h2o, k_h2o, P_surf, mr_so2_surf, mr_co2_surf, effective_precip_rate):
    """
    Function to calculate column supply of S[+IV] to surface from atmospheric wet deposition alone. 
    
    Takes as input:
        -method. If "model", calculates based on wet deposition rate reported from model. If "calc-std", calculates using mr_so2 and standard Henry's Law coefficient from Giorgi & Chameides, which assumes raindrop pH=5. If "calc-raindrop", accounts for effect of pSO2 and pCO2 on raindrop pH.
        -model SO2 wet deposition (cm**-2 s**-1)
        -model H2O column density (cm**-2)
        -model rainout rate (s**-1)
        -Surface pressure (bar)
        -Surface mixing ratio of SO2 (dimensionless)
        -Surface mixing ratio of CO2 (dimensionless)
        -effective precipitation rate into lake, in units of m year**-1. 

    Returns:
        -column flux of S[+IV] (SO2, which instantly interconverts to sulfite), in cm**-2 s**-1
    
    The calculation methods assume equilibrium with bottommost layer of atmosphere. This is where pSO2 should be highest, so if equilibrium is not efficient may slightly overestimate supply. 
    """
    
    precip_rate_cgs=effective_precip_rate*m2cm*1.0/(year2s) #precip rate converted to cm/s from m/year

    if method=='calc-raindrop':
        pH_raindrop, H_prime_so2=raindrop_pH(P_surf,  mr_co2_surf, mr_so2_surf)
        column_flux_cgs=(P_surf*mr_so2_surf*H_prime_so2*M2cgs)*precip_rate_cgs #cm**-2 s**-1
    elif method=='calc-std':
        H_prime_so2=4E+3*1.0/(atm2bar) #effective Henry's Law constant for SO2 from Giorgi & Chameides 1985, for pH=5, T=290K, converted from M atm**-1 to M bar**-1
        column_flux_cgs=(P_surf*mr_so2_surf*H_prime_so2*M2cgs)*precip_rate_cgs #cm**-2 s**-1
    elif method=='model':
        model_precip_rate_list=k_h2o*N_h2o*(18.02/6.02E23*1.0) #s**-1 * cm**-2, use molar mass of H2O and mass density of H2O to convert to cm/s
        column_flux_cgs=so2_wet_dep/model_precip_rate_list*precip_rate_cgs
        
    return column_flux_cgs 



###############################################################################
###Define sulfite loss reactions
###############################################################################

#####
###Sulfite seepage
#Treat as semi-free parameter due to wide range of seepage rates possible. 
#Some reasonable values for seepage rates
#Pearce et al. 2017: 0.95 m/year, based on averaging a range of seepage rates, ultimately taken from Boyd et al. 1982 
#Steinmen et al. 2010 state lake seeps 0.5% volume/month, 0.1-0.2 of evap+seep. For a 100-cm lake, this works out to 100*0.005 cm/month = 0.06 m/year. 
#####
def s_iv_seepage_flux(conc_s_iv, seepage_rate):
    """
    Calculates loss rate of 

    Parameters
    ----------
    conc_s_iv : concentration of S[IV] species in M
    seepage_rate: seepage rate of water body in m year**-1 [NOTE NON-CGS UNITS]
    
    Returns
    -------
    Column-integrated sulfite seepage flux (cm**-2 s**-1)
    """
    conc_s_iv_cbgs=conc_s_iv*M2cgs #convert concentration from M=mol/L to cm**-3
    seepage_rate_cgs=seepage_rate*m2cm*1.0/(year2s)
    col_int_rate_cgs=conc_s_iv_cbgs*seepage_rate_cgs
    return col_int_rate_cgs


#####
###Sulfite disproportionation
#####
def k_s_iv_disproportionation(rxn_order, T_sulfite_disprop_exp, conc_s_iv_exp):
    """
    Calculate sulfite disproportionation reaction rate constant from sulfite disproportionation lifetime experiments, assumed reaction order.
    
    Fiducial standard: 1st-order in sulfite, following Halevy+2013 PNAS SI. 
    In truth, it is probably 2nd or 4th order given stoichiometry, but let's start with what's in the literature. 

    Parameters
    ----------
    rxn_order: assumed order of reaction
    T_sulfite_disprop_0: lifetime of sulfite from a given experiment
    conc_s_iv_0: concentration of sulfite from the same experiment

    Returns
    -------
    rate constant in units of M**(1-rxn_order) s**-1.
    """
    
    return (T_sulfite_disprop_exp**-1.0)*conc_s_iv_exp**(1.0-rxn_order)

def s_iv_disproportionation_flux(conc_s_iv, d_body, rxn_order, T_sulfite_disprop_exp, conc_s_iv_exp):
    """
    Calculate column-integrated S[IV] disproportionation flux

    Parameters
    ----------
    conc_s_iv : concentration of S[IV] species in M
    d_body : depth of water body, cm
    rxn_order: assumed order of reaction
    T_sulfite_disprop_exp: lifetime of sulfite from a given experiment
    conc_s_iv_exp: concentration of sulfite from the same experiment
    Returns
    -------
    Column-integrated sulfite disproportionation flux (cm**-2 s**-1)
    """
    k_sivdisprop=k_s_iv_disproportionation(rxn_order, T_sulfite_disprop_exp, conc_s_iv_exp) #Units: M***(1-rxn_order) s**-1
    rate=k_sivdisprop*conc_s_iv**rxn_order #M s**-1
    col_int_rate=(rate*M2cgs)*d_body #cm**-2 s**-1
    return col_int_rate


#####
###Sulfite direct oxidation. 
#####

def k_s_iv_o2_oxidation_seawater(T_surf, pH, I):
    """
    Calculate direct sulfite oxidation constant, assuming SEAWATER relative composition. 
    
    Parameters
    ----------
        T_surf is temp in K
        pH
        I is ionic strength (molal)
    Returns
    -------
        rate constant for sulfite oxidation by o2, M**-1.5 s**-1
        
    Implements Zhang & Millero 1991. 
    
    The overall rate law takes form:
        dS(IV)/dt = k [S(IV)]^2 [O2]^0.5    (((Zhang & Millero Eqn 21)))
                  = k" alpha_HSO3 alpha_SO3 [S(IV)]^2 [O2]^0.5 (((Zhang & Millero Eqn 20)))
                  = k" [HSO3-][SO3--][02]^0.5 ((Zhang & Millero Eqn 19))
        where alpha_HSO3 is the mole fraction of HSO3- and alpha_SO3 is the mole fraction of SO3(2-)
   
    We take alpha_HSO3 and alpha_SO3 via manipulating equations 17 and 18 of Zhang & Millero 1991, and via equations for pKa_1 and pKa_2 from Millero et al. 1989
 
   
    For seawater,
        log10(k in M^-1.5 min^-1) = 19.54 - 5069.47/T + 14.74I^0.5 - 2.93I -2877.0 I^0.5/T
    where T = T_surf is the temperature in K, and I is the molal ionic strength. 
    
    This seawater point can be used to calibrate
        k" = k/(alpha_HSO3 * alpha_SO3)
    We do this internally, not relying on the value quoted on pg. 677 of Zhang & Millero 1991, because it has proven impossible to recover their alpha_SO3 and alpha_HSO3 with the methods they describe. That said, the difference is practically less than 0.1 log unit, so not a bit deal (~20% effect)
        
           
    Some notes:
        -Based on measurements spanning 15-45C, pH=4-8.5, salinity=0.2-35 PSU [note that the paper is contradictory on 0.2 vs 0 as the lower limit]
        -Based on experiments conducted at middling pH and sulfite concentrations 5-10 uM, meaning that the only species of S[IV] present were HSO3-, SO3[2-]. 
        -Adddition of sulfate, calcium, magnesium cause the oxidation to decrease, possibly due to interference chemistry
        -In conditions of "excess O2" (not specified) the dependence on [O2] goes away. Thus, this prescription may overestimate the loss rate in modern-Earth like conditions.
        -There is significant pH dependence (~0.5 log units in 4.5 pH units) left in for seawater (Zhang+1991 Fig. 7). Probably due to pH-dependent complexation in the seawater. 
        
    Some important caveats:
        -The paper states k" for seawater is 6.17, and for NaCl solution a bit higher -- that appears to be a typo, it should be log10(k"), see also Figure 7. 
        -alpha_SO3, alpha_HSO3 caveat above. 
        
    """

    ###Calibrate k_dp = k" using seawater, S=35, T=25 reference case.
    ##Derive seawater parameters
    S_seawater=35.0 #salinity of seawater, as in Zhang & Millero
    I_seawater=0.0199*S_seawater/(1.0-1.0E-3*S_seawater) #See page 679 of Zhang & Millero
    pKa_1_seawater=1.87 + -0.50*I_seawater**0.5 + 0.31*I_seawater #Millero+1989, p381, 25C value
    pKa_2_seawater=7.12 + -1.052*I_seawater**0.5 + 0.36*I_seawater #Millero+1989, p381. T-independent 5-25C.
    Ka_1_seawater=10**(-pKa_1_seawater)
    Ka_2_seawater=10**(-pKa_2_seawater)
    conc_H_seawater=10**(-8.2)
    alpha_HSO3_seawater=(conc_H_seawater/Ka_1_seawater + Ka_2_seawater/conc_H_seawater+1.0)**-1 #Equation 17 of Zhang & Millero
    alpha_SO3_seawater=(conc_H_seawater**2.0/(Ka_1_seawater*Ka_2_seawater)+conc_H_seawater/Ka_2_seawater+1.0)**-1  #Equation 18
    
    
    # T_seawater=273.15+25.0 # Temperature of measurements
    # log_k_seawater_min_std=(19.54 - 5069.47/T_seawater + 14.74*I_seawater**0.5 - 2.93*I_seawater -2877.0*I_seawater**0.5/T_seawater) #units: M**-1.5 min**-1 #k for seawater at pH=8.2, S=35.
    # log_k_dp_seawater_min_std=log_k_seawater_min_std-np.log10(0.023*0.977) #This produces a log10(k") for seawater of 6.395, matching the graphclicked 6.37 from Figure 7 almost perfectly. This means that we confirm we know how the calculation is being done. However, the alpha_HSO3_seawater=0.019 and alpha_SO3_seawater=0.981 we calculate yields a log10(k") of 6.48, off by about 0.11 dex
    
    # ###Don't permit wild extrapolation
    # if np.isscalar(T_surf):
    #     if T_surf>(273.15+45):
    #         T_surf=273.15+45
    #     if T_surf<(273.15+15):
    #         T_surf=273.15+15
    # else:
    #     T_surf[T_surf>(273.15+45)]=273.15+45
    #     T_surf[T_surf<(273.15+15)]=273.15+15

    ###Calculate k for this specific case. 
    k_dp_min=10.0**(19.54 - 5069.47/T_surf + 14.74*I**0.5 - 2.93*I -2877.0*I**0.5/T_surf)/(alpha_HSO3_seawater*alpha_SO3_seawater) #units: M**-1.5 min**-1. Combination of equations 8, 22 of Zhang & Millero
   
    k_dp = k_dp_min/60. #convert to M**-1.5 s**-1
    # k_dp = 10**6.17/60. #convert to M**-1.5 s**-1 ##This is if we just take the seawater value.
 
    ##Calculate mole fractions
    pKa_1=1.87 + -0.50*I**0.5 + 0.31*I #Millero+1989, p381, 25C value
    pKa_2=7.12 + -1.052*I**0.5 + 0.36*I #Millero+1989, p381. T-independent 5-25C.
    Ka_1=10**(-pKa_1)
    Ka_2=10**(-pKa_2)
    conc_H=10**(-pH)
    alpha_HSO3=(conc_H/Ka_1 + Ka_2/conc_H+1.0)**-1 #Equation 17 of Zhang & Millero
    alpha_SO3=(conc_H**2.0/(Ka_1*Ka_2)+conc_H/Ka_2+1.0)**-1  #Equation 18 of Zhang & Millero
    
    ##Calculate and return k
    k=k_dp*alpha_HSO3*alpha_SO3
    return k

# ###We reproduce the 6.5 pH maximum. Note that to match the scale in Figure 6, we have to replace the detailed law with k_dp=10**6.17 M**-1.5 min**-1. 
# pHs=np.linspace(4.0, 8.5, num=100)
# logks=np.log10(60*k_dp_s_iv_o2_oxidation_seawater(273.15+25, pHs, 0.722))
# plt.plot(pHs, logks)

def s_iv_o2_oxidation_seawater_flux(conc_s_iv, d_body, p_O2_surf, T_surf, pH, I):
    """
    Calculate column-integrated S[IV] direct oxidation flux, assuming SEAWATER relative composition. 
    ONLY VALID WHEN HSO3-, SO32- ARE DOMINANT FORMS OF S[IV]. I.e. dilute S[IV] (<0.6 M), non-acidic pH
    
    Parameters
    ----------
        [S(IV)] in M
        d_body : depth of water body, cm
        pO2 is oxygen partial pressure in bar at the surface (assumes saturation)
        k is the rate constant (M**-1.5 s**-1)
    Returns
    -------
        Column-integrated sulfite oxidation flux (cm**-2 s**-1)
        
    
    Implements Zhang & Millero 1991. 
    
    The overall rate law takes form:
        dS(IV)/dt = k [S(IV)]^2 [O2]^0.5    (((Zhang & Millero Eqn 21)))
                  = k" alpha_HSO3 alpha_SO3 [S(IV)]^2 [O2]^0.5 (((Zhang & Millero Eqn 20)))
                  = k" [HSO3-][SO3--][02]^0.5 ((Zhang & Millero Eqn 19))
        where alpha_HSO3 is the mole fraction of HSO3- and alpha_SO3 is the mole fraction of SO3(2-)
        
    Some important caveats:
        -The use of a first-order formalism in [S[IV]] in Table S4 of Halevy+2013 is a typo (I. Halevy, personal communication, April 26 2020). 
        
    """
    ###Get rate constant, from previous function
    k=k_s_iv_o2_oxidation_seawater(T_surf, pH, I) #units: M**-1.5 s**-1

    ###Get [O2] from pO2. For time being, ignore salinity, temperature effects and just use straight Henry's Law.
    H_O2=1.3e-5*1.0e2 #Henry's Law constant for O2, converted from mol m^-3 Pa^-1 to M/bar.
    conc_O2=H_O2*p_O2_surf
    
    ###Evaluate loss rate
    rate=k*(conc_s_iv**2.0)*(conc_O2**0.5) #M s**-1
    
    ###Convert to column-integrated loss rate
    col_int_rate=(rate*M2cgs)*d_body #cm**-2 s**-1

    return col_int_rate

def get_conc_siv_tot(pSO2, pH, I):
    H_SO2=1.34 #M/bar; Ranjan+2018 Table C1.
    conc_SO2=pSO2*H_SO2
        
    pKa_1=1.87 + -0.50*I**0.5 + 0.31*I #Millero+1989, p381, 25C value
    pKa_2=7.12 + -1.052*I**0.5 + 0.36*I #Millero+1989, p381. T-independent 5-25C.
    Ka_1=10**(-pKa_1)
    Ka_2=10**(-pKa_2)
    conc_H=10**(-pH)
        
    conc_HSO3=Ka_1*conc_SO2/conc_H
    conc_SO3 = Ka_2*conc_HSO3/conc_H
        
    return (conc_SO2+conc_HSO3+conc_SO3)

###############################################################################
###Establish equation
###############################################################################
def eqn(conc_s_iv, P_surf, T_surf, mr_so2_surf, mr_co2_surf, v_dep_so2, d_body, disprop_rxn_order, T_sulfite_disprop_exp, conc_s_iv_exp, pH, I, p_O2_surf, water_molabs_choice, longwaveQY0, K_d_method, seepage_rate, effective_precip_rate, wetdep_method, so2_wet_dep, N_h2o, k_h2o):
    supply_flux=s_iv_dry_deposition_flux(P_surf, T_surf, mr_so2_surf, v_dep_so2)\
        +s_iv_wet_deposition_flux(wetdep_method, so2_wet_dep, N_h2o, k_h2o, P_surf, mr_so2_surf, mr_co2_surf, effective_precip_rate)
    
    loss_flux=s_iv_disproportionation_flux(conc_s_iv, d_body, disprop_rxn_order, T_sulfite_disprop_exp, conc_s_iv_exp) \
             +s_iv_o2_oxidation_seawater_flux(conc_s_iv, d_body, p_O2_surf, T_surf, pH, I)\
            +photcalc.sulfite_photolysis_rate(conc_s_iv, d_body, pH, I, water_molabs_choice, longwaveQY0, K_d_method)\
                +s_iv_seepage_flux(conc_s_iv, seepage_rate)
    
    return supply_flux-loss_flux

###############################################################################
###Solve equation for specific case
###############################################################################
if singlescenariocalc:
    ###Planet parameters
    P_surf_0=1.0
    T_surf_0=288.0
    p_O2_surf_0=0.2*1E-7

    mr_so2_surf_0= 1.0E-10#from Claire et al. 2014, 1E10 S flux case.
    mr_co2_surf_0= 0.01#from Claire et al. 2014, 1E10 S flux case.
    
    so2_wetdep_0=4.E8 #cm**-2 s**-1
    N_h2o_0=4.486E22 #cm**-2
    k_h2o_0=2.0E-6 #s**-1, from Hu+2012.
    wetdep_method_0='calc-raindrop' #If 'model', takes it from calculation from MEAC. If 'calc-raindrop' calculate it from pSO2, pCO2. If 'calc-std', use "standard" value for Giorgi & Chameides 1985 (different from what Hu+2012 use, which may instead be from Seinfeld & Pandis)

    ###Body-of-water parameters
    d_body_0=100.0 #cm
    water_molabs_choice_0='prebiotic_carbonatelake_low'
    pH_0=5.0
    I_0=0.001
    seepage_rate_0=0.15 #m/year
    
    evaporationrate_0=-0.12+0.06*(T_surf_0-273.15)# m/year. expression from Boyd+1985 via Pearce+2017
    effective_precip_rate_0=evaporationrate_0+seepage_rate_0 #m/year
    
    ###Chemistry parameters
    T_sulfite_disprop_exp=5.0*year2s
    conc_s_iv_exp=0.1 #M
    disprop_rxn_order=1.0
    v_dep_so2_0=1.0
    longwaveQY0_0=True #If true, QY at wavelengths past which we have constraints is set to 0. If false, QY is set to the last value. 
    K_d_method_0='Morel1991'# Options are Ranjan2022, Morel2007, Morel1991
   
   ###Numerical parameters
    initial_guess=1.0E-7 #initial guess for sulfite concentration.
    conc_sulfite=fsolve(eqn, initial_guess, args=(P_surf_0, T_surf_0, mr_so2_surf_0, mr_co2_surf_0, v_dep_so2_0, d_body_0, disprop_rxn_order, T_sulfite_disprop_exp, conc_s_iv_exp, pH_0, I_0, p_O2_surf_0, water_molabs_choice_0, longwaveQY0_0, K_d_method_0, seepage_rate_0, effective_precip_rate_0, wetdep_method_0, so2_wetdep_0, N_h2o_0, k_h2o_0), xtol=1.0E-6, maxfev=10000)
    
    print(conc_sulfite) #sulfite concentration in M?

###############################################################################
###Plot all supply mechanisms (like in Ranjan+2018)
###############################################################################

if plotsupplymechanisms:
    ###Planet/Atmosphere parameters
    #This block comes from modeling by Sangita
    phi_volc_list=np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
    phi_so2_list=np.array([3.0E+8, 9.0E+8, 3.0E+9, 9.0E+9, 3.0E+10, 9.0E+10])
    mr_so2_list=np.array([7.E-12, 2.E-11, 1.E-10, 3.E-10, 1.E-9, 3.E-9])
    so2_wetdep_list=np.array([3.E7, 9.E7, 4.E8, 2.E9, 5.E9, 1.E10])
    h2o_col_den_list=np.array([4.486252e+22, 4.486253e+22, 4.486257e+22, 4.486294e+22, 4.486313e+22, 4.486327e+22])
    model_precip_rate_list=2.0E-6*h2o_col_den_list*(18*1/(6.02E23))*(year2s/m2cm) #s**-1 * cm**-2, use molar mass of H2O and mass density of H2O to convert to cm/s, convert cm/s to m/year.
    v_dep_so2_0=1.0 #cm/s
    P_surf_0=1.0 #bar
    T_surf_0=288.0    
    k_h2o=2.0E-6 #s**-1, from Hu+2012

    ###Pond parameters
    seepage_rate_0=0#0.15 #m/year
    evaporationrate_0=-0.12+0.06*(T_surf_0-273.15)# m/year. expression from Boyd+1985 via Pearce+2017    
    
    ###LITERATURE PARAMETERS
    phi_so2_kasting=np.array([3.0E+9]) #Kasting+1989
    p_so2_kasting=np.array([2.0E-9*2.9*atm2bar]) #This is the partial pressure in bar ##Kasting+1989
    
    phi_so2_claire=np.array([3.0E+8, 9.0E+8, 3.0E+9, 9.0E+9]) #Claire+2014
    p_so2_claire=np.array([5.9E-12, 1.8E-11, 4.6E-11, 9.6E-11]) #Claire+2014
    
    ###Plot: pSO2
    fig, ax=plt.subplots(2, figsize=(8., 8.), sharex=False)
    ax[0].plot(phi_so2_list, mr_so2_list*1.0, linewidth=2, linestyle='-', markersize=10, marker='d', color='black', label='This Work')
    ax[0].plot(phi_so2_kasting, p_so2_kasting, linewidth=0, color='blue', markersize=10,  marker='s',label='Kasting+1989')
    ax[0].plot(phi_so2_claire, p_so2_claire, linewidth=0, color='red', markersize=10,  marker='o', label='Claire+2014')

    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'pSO$_2$ (bar)', fontsize=16)
    ax[0].legend(ncol=1, loc='best', fontsize=11.5)    

    ax[0].set_xscale('log')
    ax[0].set_xlabel(r'$\phi_{SO_{2}}$ (cm$^{-2}$ s$^{-1}$)',fontsize=16)
    ax[0].set_xlim([3.0E+8, 9.0E+10])
    ax[0].xaxis.tick_top()
    ax[0].xaxis.set_label_position('top')

    ###Plot: Deposition.                 
    # fig2, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    ax[1].plot(phi_volc_list, s_iv_dry_deposition_flux(P_surf_0, T_surf_0, mr_so2_list, v_dep_so2_0), linewidth=2, linestyle='-', color='black', marker='s', markersize=5, label=r'Dry Deposition')  

    ax[1].plot(phi_volc_list,s_iv_wet_deposition_flux('model', so2_wetdep_list, h2o_col_den_list, k_h2o, P_surf_0, mr_so2_list, 0.1, evaporationrate_0),linewidth=2, marker='o', markersize=5,linestyle='-', color='purple', label=r'Wet Deposition, $S=0$')  
    ax[1].plot(phi_volc_list,s_iv_wet_deposition_flux('model', so2_wetdep_list, h2o_col_den_list, k_h2o, P_surf_0, mr_so2_list, 0.1, evaporationrate_0+0.15),linewidth=2,marker='o', markersize=5, linestyle='--', color='purple', label=r'Wet Deposition, $S=0.15$ m/yr')  
    ax[1].plot(phi_volc_list,s_iv_wet_deposition_flux('model', so2_wetdep_list, h2o_col_den_list, k_h2o, P_surf_0, mr_so2_list, 0.1, evaporationrate_0+2.0),linewidth=2, marker='o', markersize=5,linestyle=':', color='purple', label=r'Wet Deposition, $S=2$ m/yr')  

    ax[1].set_yscale('log')
    ax[1].set_ylabel(r'$F_{S[IV]}$ (cm$^{-2}$ s$^{-1})$', fontsize=16)
    ax[1].legend(ncol=1, loc='best', fontsize=11.5)    

    ax[1].set_xscale('log')
    ax[1].set_xlabel(r'$\frac{\phi}{\phi_0}$', fontsize=16)
    ax[1].set_xlim([1.0E-1, 3.0E+1])

    fig.subplots_adjust(hspace=0.33)
    ax[0].yaxis.set_tick_params(labelsize=14)
    ax[0].xaxis.set_tick_params(labelsize=14)
    ax[1].yaxis.set_tick_params(labelsize=14)
    ax[1].xaxis.set_tick_params(labelsize=14)
    
    ax[0].set_title('a', loc='left', fontsize=15)
    ax[1].set_title('b', loc='left', fontsize=15)
    
    def so2_to_volc(phi_so2):
        return phi_so2/(3.0E9)
    
    def volc_to_so2(phi_volc):
        return phi_volc*3.0E9
    
    secax0=ax[0].secondary_xaxis('bottom', functions=(so2_to_volc, volc_to_so2))
    # secax0.set_xlabel(r'$\frac{\phi}{\phi_0}$', fontsize=16, labelpad=-10)
    secax0.xaxis.set_tick_params(labelsize=14)
    
    secax1=ax[1].secondary_xaxis('top', functions=(volc_to_so2, so2_to_volc))
    # secax1.set_xlabel(r'$\phi_{SO_{2}}$ (cm$^{-2}$ s$^{-1}$)', fontsize=16)    
    secax1.xaxis.set_tick_params(labelsize=14)

    plt.savefig('./Plots/S-IV_Prod.pdf')
    
    
    
###############################################################################
###Plot all loss mechanisms (like in Ranjan+2019)
###############################################################################
if plotlossmechanisms:
    
    ###Planet/Atmosphere parameters
    T_surf_0=288.0
    p_O2_surf_0=0.2*1E-7
   
    ###Body-of-water parameters
    d_body_0=100.0 #cm
    water_molabs_choice_0='prebiotic_carbonatelake_low'
    pH_0=5.0
    I_0=0.15
    
    seepage_rate_0=0.15 #m/year
    evaporationrate_0=-0.12+0.06*(T_surf_0-273.15)# m/year. expression from Boyd+1985 via Pearce+2017
    effective_precip_rate_0=evaporationrate_0+seepage_rate_0 #m/year
    
    ###Chemistry parameters
    T_sulfite_disprop_exp=5.0*year2s
    conc_s_iv_exp=0.1 #M
    disprop_rxn_order=1.0
    v_dep_so2_0=1.0
    longwaveQY0_0=True #If true, QY at wavelengths past which we have constraints is set to 0. If false, QY is set to the last value. 
    K_d_method_0='Morel1991'# Options are Ranjan2022, Morel2007, Morel1991
    
    
    ###Plot
    conc_s_iv_list=np.logspace(-9.0, -3.0, num=50, endpoint=True, base=10.0) #sulfite concentrations.
                
    fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=True)
    ax.plot(conc_s_iv_list, s_iv_o2_oxidation_seawater_flux(conc_s_iv_list, d_body_0, 2.0e-8, T_surf_0, pH_0, I_0), linewidth=2, linestyle='-', color='green', label=r'Direct Oxidation, pO$_2=2\times10^{-8}$ bar (SMIF upper limit)')  
    ax.plot(conc_s_iv_list, s_iv_o2_oxidation_seawater_flux(conc_s_iv_list, d_body_0, 3.0E-12, T_surf_0, pH_0, I_0), linewidth=2, linestyle='--', color='green', label=r'Direct Oxidation, pO$_2=3\times10^{-12}$ bar (model prediction)')  

    
    ax.plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_body_0, 1.0, 10.0*day2s, 0.09), linewidth=2, linestyle='-', color='red', label=r'Disproportionation, $n=1$, $T_{disp,0}=$10 days, [S[IV]]$_0$=0.09M (Guekezian+1997)')  
    ax.plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_body_0, 2.0, 10.0*day2s, 0.09), linewidth=2, linestyle='--', color='red', label=r'Disproportionation, $n=2$, $T_{disp,0}=$10 days, [S[IV]]$_0$=0.09M (Guekezian+1997)')  

    ax.plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_body_0, 1.0, 5.0*year2s, 1.0E-0), linewidth=2, linestyle='-', color='blue', label=r'Disproportionation, $n=1$, $T_{disp,0}=$5 years, [S[IV]]$_0$=1M (Meyer+1979)')  
    ax.plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_body_0, 2.0, 5.0*year2s, 1.0E-0), linewidth=2, linestyle='--', color='blue', label=r'Disproportionation, $n=2$, $T_{disp,0}=$5 years, [S[IV]]$_0$=1M (Meyer+1979)')  
    
    
    ax.plot(conc_s_iv_list, s_iv_seepage_flux(conc_s_iv_list, 2.0), linewidth=2, linestyle='-', color='orange', label=r'Seepage, S=2 m/yr (Alambama Fishponds)')  
    ax.plot(conc_s_iv_list, s_iv_seepage_flux(conc_s_iv_list, 0.15), linewidth=2, linestyle='--', color='orange', label=r'Seepage, S=0.15 m/yr (Lake Titicaca UL)')  
    #Need Lonar Lake
    ax.plot(conc_s_iv_list, s_iv_seepage_flux(conc_s_iv_list, 0.01), linewidth=2, linestyle=':', color='orange', label=r'Seepage, S=0.01 m/yr')
    
    photrate1=np.zeros(np.shape(conc_s_iv_list))
    photrate2=np.zeros(np.shape(conc_s_iv_list))
    for ind in range(0, len(conc_s_iv_list)):
        conc_s_iv=conc_s_iv_list[ind]
        photrate1[ind]=photcalc.sulfite_photolysis_rate(conc_s_iv, d_body_0, pH_0, I_0, water_molabs_choice_0, True, K_d_method_0)
        photrate2[ind]=photcalc.sulfite_photolysis_rate(conc_s_iv, d_body_0, pH_0, I_0, water_molabs_choice_0, False, K_d_method_0)
        

    ax.plot(conc_s_iv_list, photrate1, linewidth=2, linestyle='-', color='purple', label=r'Photolysis, $\Phi(>254nm)=0$')  
    ax.plot(conc_s_iv_list, photrate2, linewidth=2, linestyle='--', color='purple', label=r'Photolysis, $\Phi(>254nm)>0$')  

    ax.set_yscale('log')
    ax.set_ylabel(r'$F_{S[IV]}$ (cm$^{-2}$ s$^{-1}$')
    ax.legend(ncol=1, loc='best')    

    ax.set_xscale('log')
    ax.set_xlabel('[S[IV]] (M)')
    ax.set_xlim([1.0E-8, 1.0E-4])
    plt.savefig('./Plots/S-IV_Loss.pdf')
    

    
###############################################################################
###Plot steady-state concentrations as function of volcanism level for different waters
###############################################################################
if plotsteadystatecalc_fso2_ocean:
    ocean_initial_guess=1.0E-8 #low guess for ocean.
  
    ###Run calculation to solve equation
    ocean_chemmax_geomax_list=np.zeros(np.shape(phi_volc_list))
    ocean_chemmax_geomin_list=np.zeros(np.shape(phi_volc_list))
    ocean_chemmin_geomax_list=np.zeros(np.shape(phi_volc_list))
    ocean_chemmin_geomin_list=np.zeros(np.shape(phi_volc_list))    
    
    for ind in range(0, len(phi_volc_list)):
        mr_so2=mr_so2_list[ind]
        so2_wetdep=so2_wetdep_list[ind]
        N_h2o=N_h2o_list[ind]
        p_o2=p_o2_list[ind]
        
        ###Oceans
        effective_precip_rate_ocean=effective_precip_rate_ocean_list[ind]
        
        ocean_chemmax_geomax_list[ind]=fsolve(eqn, ocean_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_ocean, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_ocean_max, I_ocean_max, p_o2, water_molabs_ocean_max, longwaveQY0_max, K_d_method_0, seepage_rate_ocean, effective_precip_rate_ocean, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        ocean_chemmax_geomin_list[ind]=fsolve(eqn, ocean_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_ocean, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_ocean_min, I_ocean_min, p_o2, water_molabs_ocean_min, longwaveQY0_max, K_d_method_0, seepage_rate_ocean, effective_precip_rate_ocean, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        ocean_chemmin_geomax_list[ind]=fsolve(eqn, ocean_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_ocean, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_ocean_max, I_ocean_max, p_o2, water_molabs_ocean_max, longwaveQY0_min, K_d_method_0, seepage_rate_ocean, effective_precip_rate_ocean, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        ocean_chemmin_geomin_list[ind]=fsolve(eqn, ocean_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_ocean, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_ocean_min, I_ocean_min, p_o2, water_molabs_ocean_min, longwaveQY0_min, K_d_method_0, seepage_rate_ocean, effective_precip_rate_ocean, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        
    ###Plot                
    fig, ax=plt.subplots(2, figsize=(8., 8.), sharex=False)
    ax[0].axhline(1.0E-6, color='black', linestyle='--')
    ax[0].axhline(1.0E-9, color='black', linestyle=':')
    
    # ax[0].fill_between(phi_volc_list, get_conc_siv_tot(mr_so2_list*1.0, pH_ocean_max, I_ocean_max), y2=get_conc_siv_tot(mr_so2_list*1.0, pH_ocean_min, I_ocean_min), color='red', label='Sulfite Saturation', alpha=0.5)
    ax[0].plot(phi_volc_list, get_conc_siv_tot(mr_so2_list*1.0, pH_ocean_max, I_ocean_max), color=(213/255.0, 94/255.0, 0), linestyle='--', label='Sulfite Saturation (Lower Limit)')

    ax[0].fill_between(phi_volc_list,ocean_chemmax_geomin_list, y2=ocean_chemmax_geomax_list, color=(0, 158/255.0, 115/255.0), label='Ocean, Inefficient Chemical Loss', alpha=0.5)
    ax[0].fill_between(phi_volc_list,ocean_chemmin_geomin_list, y2=ocean_chemmin_geomax_list, color=(230/255.0, 159/255.0, 0), label='Ocean, Efficient Chemical Loss', alpha=0.5)
    
    ax[0].set_xscale('log')
    ax[0].set_xlabel(r'$\phi/\phi_0$', fontsize=16)
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'[S[IV]] (M)', fontsize=16)
    ax[0].legend(ncol=1, loc='lower right', fontsize=10)    
    ax[0].set_xlim([1.0E-1, 3.0E1])
    ax[0].set_ylim([1.0E-11, 1.0E-4])

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    
    conc_s_iv_list=np.logspace(-11.0, -5.0, num=50, endpoint=True, base=10.0) #sulfite concentrations.

    ax[1].axhline(s_iv_dry_deposition_flux(P_surf_0, T_surf_0, mr_so2_list[2], v_dep_so2_0), color='black', linestyle='-', label=r'SO$_2$ Dry Dep., $\frac{\phi}{\phi_{0}}=1$')    
    ax[1].axhline(s_iv_wet_deposition_flux(wetdep_method_0, so2_wetdep_list[2], N_h2o_list[2], k_h2o_0, P_surf_0, mr_so2_list[2], mr_co2_surf_0, effective_precip_rate_ocean_list[2]), color=(240/255.0, 228/255.0, 66/255.0), linestyle='-', label=r'SO$_2$ Wet Dep., $\frac{\phi}{\phi_{0}}=1$')    

    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_ocean, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min), linewidth=2, linestyle='-', color=(213/255.0, 94/255.0, 0), label=r'Disprop., $n=1$, $T_{disp,0}=$1 yr')  
    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_ocean, disprop_rxn_order_max, T_sulfite_disprop_exp_min, conc_s_iv_exp_min), linewidth=2, linestyle='--', color=(213/255.0, 94/255.0, 0), label=r'Disprop., $n=4$, $T_{disp,0}=$1 yr')  

    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_ocean, disprop_rxn_order_min, T_sulfite_disprop_exp_max, conc_s_iv_exp_max), linewidth=2, linestyle='-', color=(0, 114/255.0, 178/255.0), label=r'Disprop., $n=1$, $T_{disp,0}=$5 yr')  
    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_ocean, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max), linewidth=2, linestyle='--', color=(0, 114/255.0, 178/255.0), label=r'Disprop., $n=4$, $T_{disp,0}=$5 yr')  
    
    photrate1=np.zeros(np.shape(conc_s_iv_list))
    photrate2=np.zeros(np.shape(conc_s_iv_list))
    for ind in range(0, len(conc_s_iv_list)):
        conc_s_iv=conc_s_iv_list[ind]
        photrate1[ind]=photcalc.sulfite_photolysis_rate(conc_s_iv, d_ocean, pH_ocean_max, I_ocean_max, water_molabs_ocean_max, longwaveQY0_max, K_d_method_0)
        photrate2[ind]=photcalc.sulfite_photolysis_rate(conc_s_iv, d_ocean, pH_ocean_min, I_ocean_min, water_molabs_ocean_min, longwaveQY0_min, K_d_method_0)
        
    ax[1].plot(conc_s_iv_list, photrate2, linewidth=2, linestyle='-', color=(204/255.0, 121/255.0, 167/255.0), label=r'Phot., Max. Efficacy')  
    ax[1].plot(conc_s_iv_list, photrate1, linewidth=2, linestyle='--', color=(204/255.0, 121/255.0, 167/255.0), label=r'Phot., Min. Efficacy')  
    
    ax[1].plot(conc_s_iv_list, s_iv_o2_oxidation_seawater_flux(conc_s_iv_list, d_ocean, p_o2_list[2], T_surf_0, pH_ocean_min, I_ocean_min), linewidth=2, linestyle='-', color=(0, 158/255.0, 115/255.0), label=r'Dir. Ox., Max. Efficacy ($\frac{\phi}{\phi_{0}}=1$)')  


    ax[1].set_yscale('log')
    ax[1].set_ylim([1.0E-5, 1.0E+13])
    ax[1].set_ylabel(r'$F_{S[IV]}$ (cm$^{-2}$ s$^{-1})$', fontsize=16)
    # ax[1].legend(ncol=2, loc='lower left', fontsize=12)    

    ax[1].set_xscale('log')
    ax[1].set_xlabel('[S[IV]] (M)', fontsize=16)
    ax[1].set_xlim([1.0E-11, 1.0E-5])
    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)    
    
    ax[0].set_title('a', loc='left', fontsize=15)
    ax[1].set_title('b', loc='left', fontsize=15)
    
    plt.subplots_adjust(top=0.95, bottom=0.09, hspace=0.25)
    ax[1].legend(ncol=2, loc='lower left', fontsize=10, bbox_to_anchor=(0, 0, 0.5, .102))
    plt.savefig('./Plots/S-IV_conc_fvolc_ocean_loss.pdf')
    
if plotsteadystatecalc_fso2_lake:
    ###########
    ###Lakes
    ###########
    lake_initial_guess=1.0E-7 #moderate guess for lakes.
    ###Run calculation to solve equation

    carbonatelake_chemmax_geomax_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_chemmin_geomax_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_chemmax_geomin_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_chemmin_geomin_list=np.zeros(np.shape(phi_volc_list))
    
    carbonatelake_shallow_chemmax_geomax_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_shallow_chemmin_geomax_list=np.zeros(np.shape(phi_volc_list))    
    
    freshwaterlake_chemmax_geomax_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmin_geomax_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmax_geomin_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmin_geomin_list=np.zeros(np.shape(phi_volc_list))

   
    
    for ind in range(0, len(phi_volc_list)):
        mr_so2=mr_so2_list[ind]
        so2_wetdep=so2_wetdep_list[ind]
        N_h2o=N_h2o_list[ind]
        p_o2=p_o2_list[ind]
        
        ###Carbonate Lakes
        carbonatelake_chemmax_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_max, longwaveQY0_max, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        carbonatelake_chemmin_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_max, longwaveQY0_min, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)        
        carbonatelake_chemmax_geomin_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_carbonatelake_min, I_carbonatelake_min, p_o2, water_molabs_carbonatelake_min, longwaveQY0_max, K_d_method_0, seepage_rate_carbonatelake_min, effective_precip_rate_carbonatelake_min, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)        
        carbonatelake_chemmin_geomin_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_carbonatelake_min, I_carbonatelake_min, p_o2, water_molabs_carbonatelake_min, longwaveQY0_min, K_d_method_0, seepage_rate_carbonatelake_min, effective_precip_rate_carbonatelake_min, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
 
        
        ###Freshwater lakes
        freshwaterlake_chemmax_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_max, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        freshwaterlake_chemmin_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_min, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)    
        freshwaterlake_chemmax_geomin_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_min, longwaveQY0_max, K_d_method_0, seepage_rate_freshwaterlake_min, effective_precip_rate_freshwaterlake_min, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        freshwaterlake_chemmin_geomin_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_min, longwaveQY0_min, K_d_method_0, seepage_rate_freshwaterlake_min, effective_precip_rate_freshwaterlake_min, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)    
        
        
    ###Plot                

    fig, ax=plt.subplots(2, figsize=(8, 8.), sharex=False)    
    ax[0].axhline(1.0E-6, color='black', linestyle='--')

    ax[0].fill_between(phi_volc_list,carbonatelake_chemmin_geomax_list, y2=carbonatelake_chemmax_geomax_list, color=(213/255.0, 94/255.0, 0), label='Carbonate Lake, Geologically Favorable', alpha=0.5)
    ax[0].fill_between(phi_volc_list,carbonatelake_chemmin_geomin_list, y2=carbonatelake_chemmax_geomin_list, color=(204/255.0, 121/255.0, 167/255.0), label='Carbonate Lake, Geologically Unfavorable', alpha=0.5)    
    
    ax[0].fill_between(phi_volc_list,freshwaterlake_chemmin_geomax_list, y2=freshwaterlake_chemmax_geomax_list, color=(0, 114/255.0, 178/255.0), label='Freshwater Lake, Geologically Favorable', alpha=0.5)
    ax[0].fill_between(phi_volc_list,freshwaterlake_chemmin_geomin_list, y2=freshwaterlake_chemmax_geomin_list, color=(0/255.0, 158/255.0, 115/255.0), label='Freshwater Lake, Geologically Unfavorable', alpha=0.5)        
    
    
    ax[0].set_xscale('log')
    ax[0].set_xlabel(r'$\phi/\phi_0$', fontsize=16)
    ax[0].set_yscale('log')
    ax[0].set_ylabel(r'[S[IV]] (M)', fontsize=16)
    ax[0].legend(ncol=1, loc='upper left', fontsize=10)    
    ax[0].set_ylim([1.0E-8, 1.0E-4])
    ax[0].set_xlim([0.1, 30.0])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    
    conc_s_iv_list=np.logspace(-8.0, -4.0, num=50, endpoint=True, base=10.0) #sulfite concentrations.
    
    ax[1].axhline(s_iv_dry_deposition_flux(P_surf_0, T_surf_0, mr_so2_list[2], v_dep_so2_0), color='black', linestyle='-', label=r'SO$_2$ Dry Dep.,$\frac{\phi}{\phi_{0}}=1$')    
    ax[1].axhline(s_iv_wet_deposition_flux(wetdep_method_0, so2_wetdep_list[2], N_h2o_list[2], k_h2o_0, P_surf_0, mr_so2_list[2], mr_co2_surf_0, evaporationrate_lake+0.0), color=(240/255.0, 228/255.0, 66/255.0), linestyle='-', label=r'SO$_2$ Wet Dep.,$S=0$ m/y,$\frac{\phi}{\phi_{0}}=1$')    
    ax[1].axhline(s_iv_wet_deposition_flux(wetdep_method_0, so2_wetdep_list[2], N_h2o_list[2], k_h2o_0, P_surf_0, mr_so2_list[2], mr_co2_surf_0, evaporationrate_lake+2.0), color=(240/255.0, 228/255.0, 66/255.0), linestyle=':', label=r'SO$_2$ Wet Dep.,$S=2$ m/y,$\frac{\phi}{\phi_{0}}=1$')    

    
    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min), linewidth=2, linestyle='-', color=(213/255.0, 94/255.0,0), label=r'Disprop., $n=1$, $T_{disp,0}=$1 yr')  
    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_min, conc_s_iv_exp_min), linewidth=2, linestyle='--', color=(213/255.0, 94/255.0,0), label=r'Disprop., $n=4$, $T_{disp,0}=$1 yr')  

    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_max, conc_s_iv_exp_max), linewidth=2, linestyle='-', color=(0, 114/255.0, 178/255.0), label=r'Disprop., $n=1$, $T_{disp,0}=$5 yr')  
    ax[1].plot(conc_s_iv_list, s_iv_disproportionation_flux(conc_s_iv_list, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max), linewidth=2, linestyle='--', color=(0, 114/255.0, 178/255.0), label=r'Disprop., $n=4$, $T_{disp,0}=$5 yr')  
    
    photrate1=np.zeros(np.shape(conc_s_iv_list))
    photrate2=np.zeros(np.shape(conc_s_iv_list))
    for ind in range(0, len(conc_s_iv_list)):
        conc_s_iv=conc_s_iv_list[ind]
        photrate1[ind]=photcalc.sulfite_photolysis_rate(conc_s_iv, d_lake, pH_carbonatelake_max, I_carbonatelake_max, water_molabs_carbonatelake_max, longwaveQY0_max, K_d_method_0)
        photrate2[ind]=photcalc.sulfite_photolysis_rate(conc_s_iv, d_lake, pH_carbonatelake_max, I_carbonatelake_max, water_molabs_carbonatelake_max, longwaveQY0_min, K_d_method_0)
        
    ax[1].plot(conc_s_iv_list, photrate2, linewidth=2, linestyle='-', color=(204/255.0, 121/255.0, 167/255.0), label=r'Phot., Max. Efficacy')  
    ax[1].plot(conc_s_iv_list, photrate1, linewidth=2, linestyle='-.', color=(204/255.0, 121/255.0, 167/255.0), label=r'Phot., Min. Efficacy')  
    
    ax[1].plot(conc_s_iv_list, s_iv_o2_oxidation_seawater_flux(conc_s_iv_list, d_lake, p_o2_list[2], T_surf_0, pH_carbonatelake_max, I_carbonatelake_max), linewidth=2, linestyle='-', color=(0, 158/255.0, 115/255.0), label=r'Dir. Ox., Max. Efficacy ($\frac{\phi}{\phi_{0}}=1$)')  

    ax[1].plot(conc_s_iv_list, s_iv_seepage_flux(conc_s_iv_list, 2.0), linewidth=2, linestyle=':', color=(86/255.0, 180/255.0, 233/255.0), label=r'Seepage, $S=2$ m/y')  

    
    ax[1].set_yscale('log')
    ax[1].set_ylim([3.0E-3, 3.0E+12])
    ax[1].set_ylabel(r'$F_{S[IV]}$ (cm$^{-2}$ s$^{-1})$', fontsize=16)

    ax[1].set_xscale('log')
    ax[1].set_xlabel('[S[IV]] (M)', fontsize=16)
    ax[1].set_xlim([1.0E-8, 1.0E-4])
    
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)    
    
    ax[0].set_title('a', loc='left', fontsize=15)
    ax[1].set_title('b', loc='left', fontsize=15)
    # plt.subplots_adjust(top=0.95, bottom=0.09, hspace=0.9)
    # ax[1].legend(ncol=3, loc='lower left', fontsize=10, bbox_to_anchor=(-0.17, 1.1, 0.5, .102))
    plt.subplots_adjust(top=0.95, bottom=0.09, hspace=0.25)
    ax[1].legend(ncol=2, loc='lower left', fontsize=10, bbox_to_anchor=(0, 0, 0.5, .102))
    

    plt.savefig('./Plots/S-IV_conc_fvolc_lake_loss.pdf')

if plotsteadystatecalc_fso2_lake_maxsulfite:
    ###########
    ###Lakes
    ###########
    lake_initial_guess=1.0E-7 #moderate guess for lakes.
    ###Run calculation to solve equation

    carbonatelake_chemmax_geomin_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_chemmin_geomin_list=np.zeros(np.shape(phi_volc_list))
    
    carbonatelake_chemmax_1_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_chemmin_1_list=np.zeros(np.shape(phi_volc_list))
    
    carbonatelake_chemmax_geomax_list=np.zeros(np.shape(phi_volc_list))
    carbonatelake_chemmin_geomax_list=np.zeros(np.shape(phi_volc_list))
    
    for ind in range(0, len(phi_volc_list)):
        mr_so2=mr_so2_list[ind]
        so2_wetdep=so2_wetdep_list[ind]
        N_h2o=N_h2o_list[ind]
        p_o2=p_o2_list[ind]
        
        ###Carbonate Lakes
        carbonatelake_chemmax_geomin_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_carbonatelake_min, I_carbonatelake_min, p_o2, water_molabs_carbonatelake_min, longwaveQY0_max, K_d_method_0, seepage_rate_carbonatelake_min, effective_precip_rate_carbonatelake_min, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)        
        carbonatelake_chemmin_geomin_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_carbonatelake_min, I_carbonatelake_min, p_o2, water_molabs_carbonatelake_min, longwaveQY0_min, K_d_method_0, seepage_rate_carbonatelake_min, effective_precip_rate_carbonatelake_min, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        
        carbonatelake_chemmax_1_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_min, longwaveQY0_max, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)        
        carbonatelake_chemmin_1_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_min, longwaveQY0_min, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        
        carbonatelake_chemmax_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_max, longwaveQY0_max, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        carbonatelake_chemmin_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_max, longwaveQY0_min, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)   


        # ###Shallow carbonate lake
        # carbonatelake_shallow_chemmax_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake_shallow, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_max, longwaveQY0_max, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        # carbonatelake_shallow_chemmin_geomax_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, d_lake_shallow, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_carbonatelake_max, I_carbonatelake_max, p_o2, water_molabs_carbonatelake_max, longwaveQY0_min, K_d_method_0, seepage_rate_carbonatelake_max, effective_precip_rate_carbonatelake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)   
    
    ###Plot                
    fig, ax=plt.subplots(1, figsize=(8, 6.), sharex=False)    
    ax.axhline(1.0E-6, color='black', linestyle='--')


    ax.fill_between(phi_volc_list,carbonatelake_chemmin_geomax_list, y2=carbonatelake_chemmax_geomax_list, color='blue', label='Carbonate Lake, Geologically Favorable', alpha=0.5)
    ax.fill_between(phi_volc_list,carbonatelake_chemmin_1_list, y2=carbonatelake_chemmax_1_list, color='purple', label='Carbonate Lake, GF, Low $a(\lambda)$', alpha=0.5)
    ax.fill_between(phi_volc_list,carbonatelake_chemmin_geomin_list, y2=carbonatelake_chemmax_geomin_list, color='red', label='Carbonate Lake, Geologically Unfavorable', alpha=0.5)
    

    
    ax.set_xscale('log')
    ax.set_xlabel(r'$\phi/\phi_0$', fontsize=16)
    ax.set_yscale('log')
    ax.set_ylabel(r'[S[IV]] (M)', fontsize=16)
    ax.legend(ncol=1, loc='upper left', fontsize=12)    
    ax.set_ylim([1.0E-8, 1.0E-4])
    ax.set_xlim([0.1, 30.0])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig('./Plots/S-IV_conc_fvolc_lake_maxsulfite.pdf')

if plotsteadystatecalc_fso2_lake_depths:
    ###########
    ###Lakes
    ###########
    lake_initial_guess=1.0E-7 #moderate guess for lakes.
    ###Run calculation to solve equation

    
    freshwaterlake_chemmax_1_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmin_1_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmax_10_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmin_10_list=np.zeros(np.shape(phi_volc_list))   
    freshwaterlake_chemmax_100_list=np.zeros(np.shape(phi_volc_list))
    freshwaterlake_chemmin_100_list=np.zeros(np.shape(phi_volc_list))
    
    for ind in range(0, len(phi_volc_list)):
        mr_so2=mr_so2_list[ind]
        so2_wetdep=so2_wetdep_list[ind]
        N_h2o=N_h2o_list[ind]
        p_o2=p_o2_list[ind]
        
        ###Freshwater Lakes
 
        freshwaterlake_chemmax_1_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, 1.0, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_max, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        freshwaterlake_chemmin_1_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, 1.0, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_min, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)    
        
        freshwaterlake_chemmax_10_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, 10.0, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_max, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        freshwaterlake_chemmin_10_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, 10.0, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_min, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)    
        
        freshwaterlake_chemmax_100_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, 100.0, disprop_rxn_order_max, T_sulfite_disprop_exp_max, conc_s_iv_exp_max, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_max, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)
        freshwaterlake_chemmin_100_list[ind]=fsolve(eqn, lake_initial_guess, args=(P_surf_0, T_surf_0, mr_so2, mr_co2_surf_0, v_dep_so2_0, 100.0, disprop_rxn_order_min, T_sulfite_disprop_exp_min, conc_s_iv_exp_min, pH_freshwaterlake, I_freshwaterlake, p_o2, water_molabs_freshwaterlake_max, longwaveQY0_min, K_d_method_0, seepage_rate_freshwaterlake_max, effective_precip_rate_freshwaterlake_max, wetdep_method_0, so2_wetdep, N_h2o, k_h2o_0), xtol=1.0E-6, maxfev=10000)    

    
    ###Plot                
    fig, ax=plt.subplots(1, figsize=(8, 6.), sharex=False)    
    ax.axhline(1.0E-6, color='black', linestyle='--')

    ax.fill_between(phi_volc_list,freshwaterlake_chemmin_1_list, y2=freshwaterlake_chemmax_1_list, color='gold', label='1 cm', alpha=0.5)
    ax.fill_between(phi_volc_list,freshwaterlake_chemmin_10_list, y2=freshwaterlake_chemmax_10_list, color='green', label='10 cm', alpha=0.5)
    ax.fill_between(phi_volc_list,freshwaterlake_chemmin_100_list, y2=freshwaterlake_chemmax_100_list, color='blue', label='100 cm', alpha=0.5)

    
    ax.set_xscale('log')
    ax.set_xlabel(r'$\phi/\phi_0$', fontsize=16)
    ax.set_yscale('log')
    ax.set_ylabel(r'[S[IV]] (M)', fontsize=16)
    ax.legend(ncol=1, loc='upper left', fontsize=12)    
    ax.set_ylim([1.0E-8, 1.0E-4])
    ax.set_xlim([0.1, 30.0])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig('./Plots/S-IV_conc_fvolc_lake_depths.pdf')