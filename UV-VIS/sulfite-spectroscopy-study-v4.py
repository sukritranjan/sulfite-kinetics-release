#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script plots data from the actual data study on sulfite disproportionation. Need to plot multiple samples for multiple experimental conditions. but all samples taken on same day, and common mode water blank, which also needs to be plotted.
"""

import numpy as np
import pdb
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pandas as pd
import datetime
from matplotlib.backends.backend_pdf import PdfPages

###WARNING WARNING WARNING
#TEMPORARILY suppress RunTimeWarnings for readability
#BE VERY CAREFUL ALL WARNINGS SUPPRESSED, COMMENT TO REMOVE THIS FUNCTIONALITY
import warnings
warnings.filterwarnings("ignore")
###

#####################################
###Key defined parameters
#####################################
day2s=84000. #1 day in seconds

#####################################
###Define key functions
#####################################

def plot_sulfite_spectra(date_list, condition, sample_list, central_wavelength, plotblank):
    """
    The purpose of this function is to plot the progression of the sulfite spectra with time.
    
    Inputs:
        --List of dates for which we have data
        --experimental condition we are plotting
        --list of samples for that condition
        --Central wavelength to plot zooms around (nm)
        --plotblank: if true, plots the water blanks as well.
    Returns:
        --None. Instead, generates plot and saves it as PDF
    
    ASSUMES: 15% error, following advice from Corinna based on Ranjan+2022a
    WARNING: Assumes the abscissa of all spectra files to be the same, throughout time. 
    """
    ###Initialize dates.
    date0_date=datetime.date(int(date_list[0][0:4]),int(date_list[0][5:7]), int(date_list[0][8:10]))  #0th date of data
    days=np.zeros(np.shape(date_list))
    for ind in range(0, len(date_list)):
        date_date=datetime.date(int(date_list[ind][0:4]),int(date_list[ind][5:7]), int(date_list[ind][8:10]))
        days[ind]=(date_date-date0_date).days #number of days since 0th measurement
    
    ###initialize variables
    wav={} #holds wavelength in nm, by day.
    data_ext={} #holds decadic extinction (unitless)
    water_ext={} #holds decadic extinction of water blanks (unitless)
    
    water_ext_mean={} #holds averaged water blanks per date.
    data_ext_corr={} #holds corrected decadic extinction (mean water blank subtracted off) (unitless)
    
    data_ext_corr_frac={} #holds corrected decadic extinction, fractional change across whole spectrum
    data_ext_corr_zoom=np.zeros((len(sample_list), len(date_list))) #holds corrected decadic extinction at central wavelength specifically.
    data_ext_corr_zoom_frac=np.zeros((len(sample_list), len(date_list))) #fractional change in decadic extinction at central wavelength specifically.
    
    
    ##############
    ###Loop through dates to extract data.
    ##############
    for ind in range(0, len(date_list)):
        
        date=date_list[ind]
        
        ###Water blanks
        #Construct file name
        waterfile1=directory+date+'/'+date+'_longtermsulfite_blankwater1.txt'
        waterfile2=directory+date+'/'+date+'_longtermsulfite_blankwater2.txt'

        #Load water files, assign wavelength scale.
        wav[date], water_ext[waterfile1]=np.genfromtxt(waterfile1, skip_header=2, skip_footer=0, unpack=True, delimiter=',') 
        thunk, water_ext[waterfile2]=np.genfromtxt(waterfile2, skip_header=2, skip_footer=0, unpack=True, delimiter=',') 
        
        #Form averaged water blank
        water_ext_mean[date]=0.5*(water_ext[waterfile1]+water_ext[waterfile2])
        
        ###Import data
        for ind2 in range(0, len(sample_list)):
            sample=sample_list[ind2]
            #construct file name
            datafile=directory+date+'/'+date+'_longtermsulfite_'+condition+'_'+sample+'.txt'
            
            #import data
            thunk, data_ext[datafile] = np.genfromtxt(datafile, skip_header=2, skip_footer=0, unpack=True, delimiter=',')
            
            #blank subtract
            data_ext_corr[datafile]=data_ext[datafile]-water_ext_mean[date]
            
            ###fractional change
            #Base
            date0=date_list[0]
            date0file=directory+date0+'/'+date0+'_longtermsulfite_'+condition+'_'+sample+'.txt'
            
            #Whole spectrum
            data_ext_corr_frac[datafile]=(data_ext_corr[datafile]-data_ext_corr[date0file])/data_ext_corr[date0file]
            # pdb.set_trace()
            
            ###central wavelenth
            #Raw value
            data_ext_corr_zoom[ind2, ind]=np.interp(central_wavelength, wav[date], data_ext_corr[datafile])
            
            #fractional change
            data_ext_corr_zoom_frac[ind2, ind]=data_ext_corr_zoom[ind2, ind]/data_ext_corr_zoom[ind2, 0]
            
    
    ##############
    ###Loop through dates to plot
    #Is it inefficient doing the loop through dates twice? Yes, but the time cost is trivial and the benefit of clearer code is worth it. 
    #Make multipage plot following https://www.delftstack.com/howto/matplotlib/how-to-save-plots-as-pdf-file-in-matplotlib/
    ##############
    
    ###Plot set up    
    fig, ax=plt.subplots(2,2, figsize=(8., 8.))
    markersizeval=3.
    colorseq=cm.rainbow(np.linspace(0,1,len(days)))
    #these two are for each sample
    linestylelist=np.array(['-','--',':','-.'])
    markerlist=np.array(['s','d','o','^'])
    colorlist=np.array(['black', 'maroon', 'goldenrod', 'indigo'])
    
    for ind in range(0, len(date_list)):
        date=date_list[ind] #string
        day=days[ind] #days elapsed
        
        for ind2 in range(0, len(sample_list)):
            sample=sample_list[ind2]
            
            datafile=directory+date+'/'+date+'_longtermsulfite_'+condition+'_'+sample+'.txt'
            
            #Plot data, color coded by day, broadband
            if (ind2==0 and (ind==0 or ind==(len(date_list)-1))): #only have label command if it is first iteration and 0th/last day. 0th column focuses on colors.
                ax[0,0].plot(wav[date], data_ext_corr[datafile], label='day='+str(day), color=colorseq[ind],markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])
                # ax[0,0].errorbar(wav[date], data_ext_corr[datafile], yerr=data_ext_corr[datafile]*0.15, xerr=None, label='day='+str(day), color=colorseq[ind],markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])

            else:
                ax[0,0].plot(wav[date], data_ext_corr[datafile], color=colorseq[ind], markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])
                # ax[0,0].errorbar(wav[date], data_ext_corr[datafile], yerr=data_ext_corr[datafile]*0.15, xerr=None, color=colorseq[ind],markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])

            #Plot data, color coded by day, narrowband
            ax[1,0].plot(wav[date], data_ext_corr[datafile], color=colorseq[ind],markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])
            # ax[0,1].errorbar(wav[date], data_ext_corr[datafile],yerr=data_ext_corr[datafile]*0.15, color=colorseq[ind],markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])


    
    for ind2 in range(0, len(sample_list)):
        #Plot absorbance at 260, different lines by sample number
        sample=sample_list[ind2]
        # ax[0,1].plot(days, data_ext_corr_zoom[ind2,:], color='black', markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2], label=r'$N_{samp}=$'+str(sample))   
        ax[0,1].errorbar(days, data_ext_corr_zoom[ind2,:], yerr=data_ext_corr_zoom[ind2,:]*0.15, color=colorlist[ind2], markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2], label=r'$N_{samp}=$'+str(sample))            
         
            
        #Plot fractional change in absorbance at 260, different lines by sample numbers
        # ax[1,1].plot(days, data_ext_corr_zoom_frac[ind2,:], color='black', markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])
        ax[1,1].errorbar(days, data_ext_corr_zoom_frac[ind2,:], yerr=0.15, color=colorlist[ind2], markersize=markersizeval, linestyle=linestylelist[ind2], marker=markerlist[ind2])
        print(r'Condition {0}, sample {1}: relative value {2:1.2f}$\pm${3:1.2f}'.format(condition, sample, data_ext_corr_zoom_frac[ind2,-1], data_ext_corr_zoom_frac[ind2,-1]*0.15))

    
    ax[0,0].set_xscale('linear')
    ax[0,0].set_xlim([200.0, 320.0])
    ax[0,0].set_yscale('linear')
    ax[0,0].set_ylim([0.000001, 4.0])
    ax[0,0].set_ylabel(r'Decadic Absorbance $A$', fontsize=14)

    ax[1,0].set_xscale('linear')
    ax[1,0].set_xlim([central_wavelength-5, central_wavelength+5])
    ax[1,0].set_xlabel('Wavelength (nm)', fontsize=14)
    ax[1,0].set_yscale('log')
    ax[1,0].set_ylim([0.1, 2.0])
    ax[1,0].set_ylabel(r'Decadic Absorbance $A$', fontsize=14)

    
 
    ax[0,1].set_xscale('linear')
    # ax[0,2].set_xlim([200.0, 320.0])
    ax[0,1].set_yscale('linear')
    # ax[0,2].set_ylim([0,np.max(ext_260)])
    ax[0,1].set_ylim(bottom=0)
    ax[0,1].set_ylabel(r'$A($'+str(central_wavelength)+'$)nm$', fontsize=14)
    
    ax[1,1].set_xscale('linear')
    ax[1,1].set_xlabel('Time (Days)', fontsize=14)
    ax[1,1].set_yscale('log')
    ax[1,1].set_ylabel(r'$A($'+str(central_wavelength)+'$)nm/A_{0}($'+str(central_wavelength)+'$nm)$', fontsize=14)
    ax[1,1].set_ylim(bottom=0.1)
    plt.subplots_adjust(wspace=0.35)
    
    
    # ax[0,0].legend(loc='upper left', ncol=2, borderaxespad=0., fontsize=14, bbox_to_anchor=(-0.1, 1.1, 1.1, 0.102), mode='expand') 
    ax[0,0].legend(loc='upper right', ncol=1, borderaxespad=0., fontsize=14, mode='tight') 
    
    
    # ax[0,1].set_title(condition) 

    # ax[0,1].legend(loc='upper left', ncol=3, borderaxespad=0., fontsize=10, bbox_to_anchor=(-0.1,1.2))
    ax[0,1].legend(loc='lower left', ncol=1, borderaxespad=0., fontsize=13, mode='tight')

    
    if not(plotblank): #if you DON'T want to plot the water blanks, just the data
        plt.savefig('./spectrometric/data-run/plots/long_term_sulfite_'+condition+'.pdf', orientation='portrait', format='pdf', bbox_inches = 'tight')
    else:
        fig2={}
        ax2={}
        for ind in range(0, len(date_list)):
            date=date_list[ind]
            waterfile1=directory+date+'/'+date+'_longtermsulfite_blankwater1.txt'
            waterfile2=directory+date+'/'+date+'_longtermsulfite_blankwater2.txt'
            
            fig2[date], ax2[date]=plt.subplots(1, figsize=(8., 6.))
            ax2[date].plot(wav[date], water_ext[waterfile1], marker='.', color='blue', label='file 1', linestyle='--')
            ax2[date].plot(wav[date], water_ext[waterfile2], marker='.', color='red', label='file 2', linestyle='--')
            ax2[date].plot(wav[date], water_ext_mean[date], marker='.', color='black', label='file mean', linestyle='-')
            ax2[date].set_title(date)
            ax2[date].legend(loc='best', borderaxespad=0., fontsize=10) 
    
            ax2[date].set_xscale('linear')
            ax2[date].set_xlim([200.0, 320.0])
            ax2[date].set_xlabel('Wavelength (nm)')
            ax2[date].set_yscale('linear')
            ax2[date].set_ylim([0, 0.1])
            ax2[date].set_ylabel('Decadic Extinction', fontsize=12)
            
        
        pp=PdfPages('./spectrometric/data-run/plots/long_term_sulfite_'+condition+'.pdf')
        pp.savefig(fig, orientation='portrait', bbox_inches = 'tight')
        for date in date_list:
            pp.savefig(fig2[date], orientation='portrait', bbox_inches = 'tight')
        pp.close()



#####################################
#Run code. 
#####################################
directory='./spectrometric/data-run/Long-Term-Sulfite-'

########
#Pure Water + "low concentration" carbonate lake
########
list_of_dates=np.array(['2021-10-05','2021-10-19','2021-10-27', '2021-11-02', '2021-11-09', '2021-11-16', '2021-11-23', '2021-11-29', '2021-12-07', '2021-12-17', '2022-01-11', '2022-01-18', '2022-01-26', '2022-02-07', '2022-02-15','2022-02-22','2022-03-03', '2022-03-08', '2022-03-15', '2022-03-22', '2022-03-29', '2022-04-07', '2022-04-13', '2022-04-19','2022-04-26', '2022-05-03', '2022-06-01', '2022-06-07', '2022-06-13', '2022-07-13', '2022-07-20', '2022-07-28', '2022-08-03', '2022-08-08']) #Must start with 0th date. %, ,


##Pure water
plot_sulfite_spectra(list_of_dates, '10mM_unbuffered', np.array(['1', '2', '3', '4']), 240., False)
plot_sulfite_spectra(list_of_dates, '100mM_unbuffered', np.array(['2', '3', '4']), 260., False) #missing "1" because that sample got destroyed by accident. 
plot_sulfite_spectra(list_of_dates, '100mM_pH4', np.array(['1', '2', '3', '4']), 295., False)
plot_sulfite_spectra(list_of_dates, '100mM_pH7', np.array(['1', '2', '3', '4']), 260., False)
plot_sulfite_spectra(list_of_dates, '100mM_pH13', np.array(['1', '2', '3', '4']), 260., False)


#"carbonate lake", low concentration conditions
plot_sulfite_spectra(list_of_dates[1:], 'lowconcentration_pH9', np.array(['1', '2', '3', '4']), 260., False) #This set started late. 


########
#"Searles Lake"
########
list_of_dates_searles=np.array(['2021-12-09', '2021-12-17', '2022-01-11', '2022-01-18', '2022-01-26', '2022-02-07', '2022-02-15','2022-02-22','2022-03-03', '2022-03-08', '2022-03-15', '2022-03-22', '2022-03-29', '2022-04-07', '2022-04-13', '2022-04-19', '2022-04-26', '2022-05-03', '2022-06-01', '2022-06-07', '2022-06-13', '2022-07-13', '2022-07-20', '2022-07-28', '2022-08-03', '2022-08-08']) #Must start with 0th date.
#WARNING WARNING WARNING: Dec 9 data did not have water blanks. Used 12/7 water blank instead. README in 12/9 folder as well. 

###"Searles Lake", high concentration carbonate lake approximant
plot_sulfite_spectra(list_of_dates_searles, 'searleslake_pH7', np.array(['1', '2', '3', '4']), 255., False) #This set started late. 