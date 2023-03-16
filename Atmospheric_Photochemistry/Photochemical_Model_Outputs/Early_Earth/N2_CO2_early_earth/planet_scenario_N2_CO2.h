/*----------------------- planet.h --------------------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: June 5, 2011
Note: The parameters in this file can be modified to model different planets around different stars
--------------------------------------------------------------------- */

#ifndef _PLANET_H_
#define _PLANET_H_ 

/* Planet Physical Properties */
#define MASS_PLANET           5.9376E+24  /* kg */ /* Earth */
#define RADIUS_PLANET         6371000.0   /* m */ /* Earth */

/* Planet Orbital Properties */
#define ORBIT               1.00          /* AU */ /* Earth */
#define STAR_SPEC           "Data/YoungSun3d8Ga.txt" //*the change is made in coherence to the early earth scenario*//
#define FaintSun			1.0				/* Faint early Sun factor */
#define TIDELOCK			0				/* If the planet is tidelly locked */
#define STAR_TEMP			394.109	   /* Irradiance Temperature at 1 AU */
#define THETAREF			1			/* Slant Path Angle in radian */
#define PAB					0.25				/* Planet Bond Albedo */
#define FADV				0.25			/* Advection factor: 0.25=uniformly distributed, 0.6667=no Advection */
#define PSURFAB				0.0			/* Planet Surface Albedo */
#define PSURFEM				1.0			/* Planet Surface Emissivity */
#define DELADJUST           1			/* Whether use the delta adjustment in the 2-stream diffuse radiation */
#define TAUTHRESHOLD		10.0			/* Optical Depth Threshold for multi-layer diffuse radiation */
#define TAUMAX				1000.0
#define TAUMAX1				1000.0		/* Maximum optical Depth in the diffuse radiation */
#define TAUMAX2				1000000.0
#define IFDIFFUSE			1			/* Set to 1 if want to include diffuse solar radiation into the photolysis rate */

#define IFUVMULT			0			/* Whether do the UV Multiplying */
#define FUVMULT				1.0E+3		/* Multiplying factor for FUV radiation <200 nm */
#define MUVMULT				1.0E+2		/* Multiplying factor for MUV radiation 200 - 300 nm */
#define NUVMULT				1.0E+1		/* Multiplying factor for NUV radiation 300 - 400 nm */

/* Planet Temperature-Pressure Preofile*/
#define TPMODE                 1            /* 1: import data from a ZTP list; 
                                               0: calculate TP profile from the parametized formula*/
#define TPLIST                 "Data/TPProfile_N2_CO2.dat"
#define PTOP                   1.0E-10            /* Pressure at the top of atmosphere in bar(dummy variables-dont matter as TPMODE==1) */
#define TTOP				   200.0            /* Temperature at the top of atmosphere(dummy variables-dont matter as TPMODE==1)*/
#define TSTR                   200.0            /* Temperature at the top of stratosphere(dummy variables-dont matter as TPMODE==1) */
#define TINV                   0            /* set to 1 if there is a temperature inversion (dummy variables-dont matter as TPMODE==1)*/
#define PSTR                   1.0E-4            /* Pressure at the top of stratosphere(dummy variables-dont matter as TPMODE==1)*/
#define PMIDDLE				   0            /* Pressure at the bottom of stratosphere(dummy variables-dont matter as TPMODE==1) */
#define TMIDDLE				   0            /* Temperature at the bottom of stratosphere(dummy variables-dont matter as TPMODE==1)*/
#define PBOTTOM				   1.0E-1           /* Pressure at the bottom of stratosphere(dummy variables-dont matter as TPMODE==1) */
#define TBOTTOM				   200.0            /* Temperature at the bottom of stratosphere(dummy variables-dont matter as TPMODE==1)*/
#define PPOFFSET			   0.0			/* Pressure offset in log [Pa] (dummy variables-dont matter as TPMODE==1) */
#define TINTSET				   10			/* Temperature equivalent to the internal heating at 1 AU (dummy variables-dont matter as TPMODE==1) */

/* Calculation Grids*/
#define zbin 50 /*How many altitude bin?*/
#define zmax 90.0 /*Maximum altitude in km*/
#define zmin 0.0 /*Maximum altitude in km*/
#define WaveBin 9999 /*How many wavelength bin?*/
#define WaveMin 1.0 /*Minimum Wavelength in nm*/
#define WaveMax 10000.0 /*Maximum Wavelength in nm*/
#define WaveMax1 1000.0 /*Maximum Wavelength in nm for the Calculation of UV-visible radiation and photolysis rates*/
#define TDEPMAX	300.0 /* Maximum Temperature-dependence Validity for UV Cross sections */
#define TDEPMIN 200.0 /* Minimum Temperature-dependence Validity for UV Cross sections */

/* The criteria of convergence */
#define Tol1 1.0E+10
#define Tol2 1.0E-300

/* Mode of iteration */
#define	TSINI	1.0E-5	/* Initial Trial Timestep, generally 1.0E-8 */
#define FINE1 1 /* Set to one for fine iteration: Set to 2 to disregard the bottom boundary layers */
#define FINE2 1 /* Set to one for fine iteration: Set to 2 to disregard the fastest varying point */
#define TMAX 1.0E+200 /* Maximum of time step */
#define TMIN 1.0E-5 /* Minimum of time step */
#define TSPEED	1.0E+200 /* Speed up factor */
#define NMAX 1E+8 /* Maximum iteration cycles */
#define NMAXT	1.0E+202 /* Maximum iteration cumulative time in seconds */
#define MINNUM 1.0E-0 /* Minimum number density in denominator */

/* Molecular Species */
#define NSP 111 /*Number of species in the standard list*/
#define SPECIES_LIST "scenario_library/Early_Earth/N2_CO2_early_earth/species_scenario_N2r.dat"
#define AIRM 29.6 /*Average molecular mass of atmosphere, in atomic mass unit for the early earth(comprised of 90%N2 and 10%CO2)*/
#define AIRVIS	1.6E-5	/*Dynamic viscosity in SI;Dynamic viscosity in SI*(There might be a a mistake.The dynamic viscosity of air at room temparature is 1.82e-5 Pa.s but the documented dynamic viscosity at room temprature for modern earth is 1.5E-5 Pa.s.Hence rechecking is needed;ref:https://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html)*/
#define RefIdxType 3	/* Type of Refractive Index: 0=Air, 1=CO2, 2=He, 3=N2, 4=NH3, 5=CH4, 6=H2, 7=O2 */

/* Aerosol Species */
#define	AERSIZE	1.0E-7	/* diameter in m */
#define AERDEN	2.0E+3	/* density in SI */
#define	NCONDEN	1	/* Calculate the condensation every NCONDEN iterations */
#define IFGREYAER	0	/* Contribute to the grey atmosphere Temperature? 0=no, 1=yes */
#define SATURATIONREDUCTION	0.2 /* Ad hoc reduction factor for saturation pressure of water */
#define AERRADFILE1	"Data/H2SO4AER_CrossM_01.dat"	/* radiative properties of H2SO4 */
#define AERRADFILE2 "Data/S8AER_CrossM_01.dat"	/* radiative properties of S8 */

/* Initial Concentration Setting */
#define IMODE 4             /* 1: Import from SPECIES_LIST; 
                                0: Calculate initial concentrations from chemical equilibrium sub-routines (not rad);
                                3: Calculate initial concentrations from simplied chemical equilibrium formula (not rad);
                                2: Import from results of previous calculations
								4: Import from results of previous calculations in the standard form (TP import only for rad) */
#define NATOMS 23            /* Number of atoms for chemical equil     */
#define NMOLECULES 172       /* Number of molecules for chemical equil */
#define MOL_DATA_FILE "Data/molecules_all.dat" /* Data file for chemical equilibrium calculation */
#define ATOM_ABUN_FILE "Data/atom_solar.dat" /* Data file for chemical equilibrium calculation */
#define IMPORTFILEX "./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Conx.dat" /* File of concentrations X to be imported */
#define IMPORTFILEF "./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Conf.dat" /* File of concentrations F to be imported */       
#define IFIMPORTH2O 0		/* When H2O is set to constant, 1=import mixing ratios */
#define IFIMPORTCO2 0		/* When H2O is set to constant, 1=import mixing ratios */

/* Reaction Zones */
#define REACTION_LIST "Data/zone_Exoplanet_Full.dat"
#define NKin 645   /*Number of Regular Chemical Reaction in the standard list*/
#define NKinM 87  /*Number of Thermolecular Reaction in the standard list*/
#define NKinT 93  /*Number of Thermal Dissociation Reaction in the standard list*/
#define NPho 71   /*Number of Photochemical Reaction in the standard list*/
#define	THREEBODY	1.0	/* Enhancement of THREEBODY Reaction when CO2 dominant */

/* Parametization of Eddy Diffusion Coefficient */
#define EDDYPARA 2	/* =1 from Parametization, =2 from imported list */
#define KET 1.0E+5 /*unit cm2 s-1*/
#define KEH 1.0E+6
#define ZT  20.0  /*unit km*/
#define Tback 1E+4
#define KET1 1.0E+6
#define KEH1 1.0E+8
#define EDDYIMPORT	"Data/Eddyprofile_N2_CO2.dat"
#define MDIFF_H_1	4.87 /*aprroximating the values as our concerned environment is N2 dominated*/
#define MDIFF_H_2	0.698
#define MDIFF_H2_1	2.80
#define MDIFF_H2_2	0.740
#define MDIFF_H2_F	1.0

/* Parameters of rainout rates */
#define	RainF	1.0	/* Rainout factor, 0 for no rainout, 1 for earthlike normal rainout, <1 for reduced rainout */
#define	CloudDen	1.0	/* Cloud density in the unit of g m-3 */

/* Output Options */
#define OUT_FILE1			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Conx.dat"
#define OUT_FILE2			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Conf.dat"
#define OUT_HISTORY			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/History.dat"
#define OUT_PHOTORATE		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Photorate.dat"
#define OUT_CHEMICALRATE	"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ChemicalRate.dat"
#define OUT_CONVERGENCE		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Convergence.dat"
#define OUT_TIMESCALE		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Timescale.dat"
#define OUT_COLUMN			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ColumnDensity.dat"
#define OUT_STD				"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ConcentrationSTD.dat"
#define OUT_BALANCE			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/GlobalBalance.dat"
#define OUT_RADIATION		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Radiation.dat"
#define OUT_MEANOPAC		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/MeanOpacity.dat"
#define OUT_NEWTEMP			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/NewTemperature.dat"
#define NPRINT     1E+2               /* Printout results and histories every NPRINT iterations */
#define HISTORYPRINT	0			/* print out time series of chemical composition if set to 1 */

/* Input choices for the infrared opacities */
/* Must be set to the same as the opacity code */

#define CROSSHEADING		"Cross3/N2Atmos/"

#define NTEMP  13             /* Number of temperature points in grid   */
#define TLOW 100.0           /* Temperature range in K                 */
#define THIGH 400.0

#define NPRESSURE 13         /* Number of pressure points in grid      */
#define PLOW 1.0e-06         /* Pressure range in Pa                   */
#define PHIGH 1.0e+06

#define NLAMBDA 16000         /* Number of wavelength points in grid    */
#define LAMBDALOW 1.0e-07    /* Wavelength range in m                  */
#define LAMBDAHIGH 2.0e-04
#define LAMBDATYPE 1        /* LAMBDATYPE=1 -> constant resolution    */
							/* LAMBDATYPE=2 -> constant wave step     */


/* IR emission spectra output options */

#define IRLamMin	1.0		/* Minimum wavelength in the IR emission output, in microns */
#define IRLamMax	100.0	/* Maximum wavelength in the IR emission output, in microns */
#define IRLamBin	9999		/* Number of wavelength bin in the IR emission spectra */
#define Var1STD			7
#define	Var2STD			52
#define Var3STD			21
#define	Var4STD			20
#define Var1RATIO		0.0
#define	Var2RATIO		0.0
#define Var3RATIO		0.0
#define	Var4RATIO		0.0

/*  Stellar Light Reflection output options */
#define UVRFILE			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Reflection.dat"  /* Output spectrum file name */
#define UVRFILEVar1		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ReflectionVar1.dat"  /* Output spectrum file name */
#define UVRFILEVar2		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ReflectionVar2.dat"  /* Output spectrum file name */
#define UVRFILEVar3		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ReflectionVar3.dat"  /* Output spectrum file name */
#define UVRFILEVar4		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/ReflectionVar4.dat"  /* Output spectrum file name */
#define UVROPTFILE		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/UVROpt.dat" /* Output spectrum file name*/

/* Stellar Light Transmission output options */
#define UVTFILE			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Transmission.dat" /* Output spectrum file name */
#define UVTFILEVar1		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/TransmissionVar1.dat" /* Output spectrum file name */
#define UVTFILEVar2		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/TransmissionVar2.dat" /* Output spectrum file name */
#define UVTFILEVar3		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/TransmissionVar3.dat" /* Output spectrum file name */
#define UVTFILEVar4		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/TransmissionVar4.dat" /* Output spectrum file name */
#define UVTOPTFILE		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/UVTOpt.dat" /* Output spectrum file name*/

/* Thermal Emission output options */
#define IRFILE			"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/Emission.dat"	     /* Output spectrum file name */
#define IRFILEVar1		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/EmissionVar1.dat"	 /* Output spectrum file name */
#define IRFILEVar2		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/EmissionVar2.dat"	 /* Output spectrum file name */
#define IRFILEVar3		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/EmissionVar3.dat"	 /* Output spectrum file name */
#define IRFILEVar4		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/EmissionVar4.dat"	 /* Output spectrum file name */
#define IRCLOUDFILE		"./output_folder/wet_deposition_of_SO2_wrt_phi/Early_Earth/N2_CO2_early_earth/CloudTopE.dat"      /* Output emission cloud top file name */

/* Cloud Top Determination */
#define OptCloudTop	1.0	/* Optical Depth of the Cloud Top */

#endif

/* 1 Tg yr-1 = 3.7257E+9 H /cm2/s for earth */
