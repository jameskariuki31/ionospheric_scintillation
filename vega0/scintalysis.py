import yaml
from yaml import SafeLoader as SafeLoader
import csv

import numpy as np
from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy import stats
from astropy.time import Time
from datetime import datetime
from dateutil import tz

#import plotly
#import chart_studio.plotly
#from plotly.subplots import make_subplots
#import plotly.graph_objects as go
#mport plotly.express as px

import os
from os import listdir
from os.path import isfile, join

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib import gridspec
matplotlib.rcParams.update({'errorbar.capsize': 2})



yamlpath= '/home/chege/Desktop/curtin_work/vega'
from_zone = tz.gettz('UTC')
to_zone = tz.gettz('Australia/Perth')

def get_yamls():
    yaml_data = sorted(os.listdir(yamlpath))
    #print(len(yaml_data))
    #[106, 198, 288, 163, 152, 181, 158, 203, 175, 197, 268, 211, 263, 30, 256, 236, 256, 247]
    #return yaml_data[22994:23094]# + yaml_data[22797:23044]
    # works return yaml_data[22893:22942] + yaml_data[22943:22993]
    return yaml_data

def obsid_to_date(obsid): #convert obsid to date of observation
    obsid = obsid.split('.')[0]
    t = Time(int(obsid), format='gps')
    t.format = 'datetime'
    return str(t)



def make_sets(div=0, one_night=False):
    obsid_files = get_yamls()

    if one_night:
        obsidsets = [obsid_files[:div],obsid_files[div:]]
    
    else:
        dates = []
        
        for obsid in obsid_files:
            date = obsid_to_date(obsid)
            dates.append(date)
            
        d = {}
        for day in dates:
            d.setdefault(day[0:10], []).append(day)  #this chunck makes sets of the dates
        datesets = list(d.values())
        #print(len(datesets),len(datesets[1]))
        obsidsets = []
        count = 0
        for dset in datesets:
            n = len(dset)
            obsidsets.append(obsid_files[count:count+n]) #here we make sets of the yaml files
            count+=n
    print('you have', len(obsidsets), 'nights of observations. no. of obsids/night are:', [len(obsidsets[i]) for i in range(len(obsidsets))])
    return obsidsets #return the set of obsid files

def compare_set_metric(obsidsets):
    '''
    input is the list of lists containing obsids
    make a list of lists containing their metric values
    return a list of lists containing their average values.
    '''
    setmetrics = []
    for obsidlist in obsidsets:
        obsidmetrics = []
        print('................', len(obsidlist))
        for obsid in obsidlist:
            obsid = int(obsid.split('.')[0])
            print(obsid)
            with open('qadb.csv') as csvfile:
                found = False
                csv_reader = csv.reader(csvfile, delimiter=',')
                for row in csv_reader:
                    if str(row[1]) == str(obsid):
                        print(row[9])
                        obsidmetrics .append(float(row[9]))
                        found = True
                if not found:
                    obsidmetrics .append(np.nan)
                    print('obsid without metric')
        setmetrics.append(obsidmetrics)
    mean_metric = []
    for i in setmetrics:
        print(len(i))
        meanmetric = np.nanmean(i)
        mean_metric.append(meanmetric)
    print(mean_metric)
    return mean_metric

def loadfile(data_file):
	with open(data_file, 'r') as f:
                unpacked = yaml.load(f, Loader=SafeLoader)
	

	#flux = []
	#ra = []
	#dec = []
	names = []
	#ampscales = []
	for source in unpacked['sources']:
		names.append(unpacked['sources'][source]['name'])
		#dec.append(unpacked['sources'][source]['dec'])
		#ra.append(unpacked['sources'][source]['ra'])
		#flux.append(unpacked['sources'][source]['flux_density']) 
		#ampscales.append(np.nanmedian(unpacked['sources'][source]['amp_scales']))
		#print(unpacked2['sources'][source]['amp_scales'])
	#ampscales = np.nan_to_num(ampscales)
	#print(len(ampscales))
	return names #,ampscales,ra,dec

def make_sourceset(same_night_obsids): #here input a list obsids from the same night
    sourceset = []
    for obsid in same_night_obsids:
        #names,ampscales,ra,dec = loadfile(obsid)
        yaml_file = yamlpath+'/'+obsid
        names = loadfile(yaml_file)
        sourceset.append(names)
    return sourceset #A list of lists containing sources per obsid


def common_sources(sourcelists):
    sourceset = [set(i) for i in sourcelists]#turn each list in the sourceset in to a set.
    u = set.intersection(*sourceset)
    return u

def sources_to_track(obsidsets): #sorted from brightest to dimmest
    night1_obsids = obsidsets[0]
    night2_obsids = obsidsets[1]
    print('night2_obsids', night2_obsids)
    sourceset1 = make_sourceset(night1_obsids)
    print('len sourceset1', len(sourceset1))
    sourceset2 = make_sourceset(night2_obsids)
    print('len sourceset2', len(sourceset2))
    #print("sourceset2",sourceset2)
    print('getting common sources in night 1')
    common_sources_night1 = common_sources(sourceset1)
    print('common_sources_night1', common_sources_night1)
    print('getting common sources in night 2')
    common_sources_night2 = common_sources(sourceset2)
    print('common_sources_night2', common_sources_night2)
    
    print('getting common sources in both nights')
    trackable_sources = []
    for source in common_sources_night1:
        if source in common_sources_night2:
            print('Source found in both nights')
            trackable_sources.append(source)
    print('trackable_sources', trackable_sources)
    
    print('getting the brightest ones first')
    fluxes =[]
    for source in trackable_sources:
        yaml_file = yamlpath+'/'+obsidsets[0][0]
        with open(yaml_file, 'r') as f:
                    unpacked = yaml.load(f, Loader=SafeLoader)
                    for soc in unpacked['sources']:
                        if unpacked['sources'][soc]['name'] == source:
                            f_density = unpacked['sources'][soc]['flux_density']
        fluxes.append(f_density)

    fluxes, trackable_sources = zip(*sorted(zip(fluxes, trackable_sources),reverse = True))
    print(fluxes[0])
    print('We have sources we can track over both nights', trackable_sources)

    return trackable_sources

def get_values(night_obsids,source,): 
    ampscales = []
    std = []
    dates = []
    flux = []
    for obsid in night_obsids:
        print('.....looking for %s in obsid %s'%(source, obsid))
        dates.append(obsid)
        yaml_file = yamlpath+'/'+obsid
        with open(yaml_file, 'r') as f:
                    unpacked = yaml.load(f, Loader=SafeLoader)
                    for soc in unpacked['sources']:
                        if unpacked['sources'][soc]['name'] == source:
                            print('.........found it')
                            print('................getting its flux and median ampscale')
                            median_src_scal = np.nanmedian(unpacked['sources'][soc]['amp_scales'])
                            std_src_scal = np.nanstd(unpacked['sources'][soc]['amp_scales'])
                            f_density = unpacked['sources'][soc]['flux_density']
        ampscales.append(median_src_scal)
        std.append(std_src_scal)
        flux.append(f_density)

    return dates, ampscales, std, flux

def get_rms(ampscales):
    naninds = [] 
    for ind, val in enumerate(ampscales): 
        if np.isnan(val):
            naninds.append(ind)
            print('nan value encountered')
            ampscales.pop(ind)        #remove nan value
            print('nan value discarded')
    n = len(ampscales)
    model = np.ones(n)
    rms = np.sqrt(np.nanmean((ampscales - model) ** 2))
    #rms = sqrt(mean_squared_error(model,ampscales))
    return round(rms, 4), naninds

def perc_err(ampscales):
    errs = []
    for i in ampscales:
        e = (np.abs(1-i)/1)*100
        errs.append(e)

    return errs



def scalevstime_plots(night1_obsids, night2_obsids, source, plotname, night1_mdnmetric, night2_mdnmetric, plottype = 1):
    print('working on source %s in night one obsids..............................................................................'%(source))
    night1_dates, night1_ampscales, night1_stds, night1_flux = get_values(night1_obsids, source)

    print('working on source %s in night two obsids..............................................................................'%(source))
    night2_dates, night2_ampscales, night2_stds, night2_flux = get_values(night2_obsids, source)



    night1_dates = [obsid_to_date(obsid) for obsid in night1_dates]
    print(night1_dates[0])
    night2_dates = [obsid_to_date(obsid) for obsid in night2_dates]
    print(night2_dates[0])
    #print(night1_flux)
    print(len(night1_ampscales),len(night2_ampscales))
    
    #converting dates from UTC to MWA local time
    time_format = '%Y-%m-%d %H:%M:%S'
    utc_dates1 = [datetime.strptime(str(i), time_format) for i in night1_dates]
    utc_dates2 = [datetime.strptime(str(i), time_format) for i in night2_dates]
    utc_dates1 = [utc_dates1[i].replace(tzinfo=from_zone) for i in range(len(utc_dates1))]
    utc_dates2 = [utc_dates2[i].replace(tzinfo=from_zone) for i in range(len(utc_dates2))]

    mwadates1 = [utc_dates1[i].astimezone(to_zone) for i in range(len(utc_dates1))]
    mwadates2 = [utc_dates2[i].astimezone(to_zone) for i in range(len(utc_dates2))]

    #get rms values for both nights ampscales
    rms1, naninds1 = get_rms(night1_ampscales)
    rms2, naninds2 = get_rms(night2_ampscales)

    #discard dates whose ampscales were nan
    if len(naninds1)>0:
        for i in naninds1:
            print('date and std with nan ampscale discarded')
            mwadates1.pop(i)
            night1_stds.pop(i)

    if len(naninds2)>0:
        for i in naninds2:
            print('date and std with nan ampscale discarded')
            mwadates2.pop(i)
            night2_stds.pop(i)


    #make the plot 
    #make the plot
    if plottype == 1:
        night1_err = perc_err(night1_ampscales)
        night2_err = perc_err(night2_ampscales)

        fig  = plt.figure(figsize = (16,12))
        gs1 = gridspec.GridSpec(2, 2, height_ratios=[2,.5])
        gs1.update(wspace=0.0, hspace=0.0)

        ax0 = plt.subplot(gs1[0])
        ax0.plot_date(mwadates1, night1_ampscales, 'o', linestyle = ':', color='g')
        ax0.errorbar(mwadates1, night1_ampscales, yerr=night1_stds, linestyle='None',color='g')
        ax0.set_ylim([0.6, 1.4])
        ax0.set_ylabel('Median amplitude scaling factor')
        #ax0.set_xlabel('Obsid date, time in 10/2013')
        ax0.set_title('%s, %s Jy, mean qam: %s'%(source,str(round(np.nanmean(night1_flux),2)), str(round(night1_mdnmetric,2))))
        ax0.grid()

        ax1 = plt.subplot(gs1[1])#, sharey=ax0)
        ax1.plot_date(mwadates2, night2_ampscales, 'o', linestyle = ':', color='m')
        ax1.errorbar(mwadates2, night2_ampscales, yerr=night2_stds, linestyle='None',color='M')
        ax1.set_ylim([0.6, 1.4])
        #ax1.set_ylabel('Median amplitude scaling factor')
        ax1.set_xlabel('Obsid date, time in 10/2013')
        ax1.set_title('%s, %s Jy, mean qam: %s'%(source,str(round(np.nanmean(night2_flux),2)), str(round(night2_mdnmetric,2))))
        ax1.grid()

        ax2 = plt.subplot(gs1[2])
        ax2.plot(mwadates1, night1_err, color='g', linestyle='--')
        ax2.set_ylabel(r'$\epsilon(t)$')
        ax2.grid()
        ax2.grid(True, which='minor',axis='y', color='#999999', linestyle='-', alpha=0.4)
        ax2.minorticks_on()


        ax3 = plt.subplot(gs1[3])
        ax3.plot(mwadates2, night2_err, color='m', linestyle='--')
        ax3.grid()
        ax3.grid(True, which='minor',axis='y', color='#999999', linestyle='-', alpha=0.4)
        ax3.minorticks_on()
        #plt.axis('on')

        #plt.setp(ax1.get_yticklabels(), visible=False)
        #plt.setp(ax3.get_yticklabels(), visible=False)
        
        ax1.set_yticklabels([])
        ax3.set_yticklabels([])

        yticks2 = ax2.yaxis.get_major_ticks()
        yticks2[-1].label1.set_visible(False)

        xticks = ax0.xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)

        xticks1 = ax1.xaxis.get_major_ticks()
        xticks1[-1].label1.set_visible(False)

        xticks2 = ax2.xaxis.get_major_ticks()
        xticks2[-1].label1.set_visible(False)

        ax0.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)

        ax2.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)

        #ax1.yaxis.tick_left()
        #ax3.yaxis.tick_left()
        
        fig.text(0.5, 0.12, 'Date, time on October 2013', ha='center')
        fig.autofmt_xdate()                                                                                                                                                      
        plt.show()


if __name__ == '__main__':
    #obsidsets = make_sets()
    #night1_obsids = obsidsets[0]
    #night2_obsids = obsidsets[1]
    #trackable_sources = sources_to_track(obsidsets)
    #metrics = compare_set_metric(obsidsets)
    #night1_mdnmetric = metrics[0]
    #night2_mdnmetric = metrics[1]
    night1_obsids=['1065877080.yaml', '1065877200.yaml', '1065877320.yaml', '1065877440.yaml', '1065877568.yaml', '1065877688.yaml', '1065877808.yaml', '1065877928.yaml', '1065878056.yaml', '1065878176.yaml', '1065878296.yaml', '1065878664.yaml', '1065878784.yaml', '1065878904.yaml', '1065879032.yaml', '1065879152.yaml', '1065879272.yaml', '1065879392.yaml', '1065879520.yaml', '1065879640.yaml', '1065879760.yaml', '1065879880.yaml', '1065880008.yaml', '1065880128.yaml', '1065880248.yaml', '1065880368.yaml', '1065880496.yaml', '1065880616.yaml', '1065880736.yaml', '1065880856.yaml', '1065880984.yaml', '1065881104.yaml', '1065881224.yaml', '1065881344.yaml', '1065881472.yaml', '1065881592.yaml', '1065881712.yaml', '1065881832.yaml', '1065881960.yaml', '1065882080.yaml', '1065882200.yaml', '1065882448.yaml', '1065882688.yaml', '1065882936.yaml', '1065883056.yaml', '1065883176.yaml', '1065883296.yaml', '1065883424.yaml', '1065883544.yaml', '1065883664.yaml', '1065883784.yaml', '1065883912.yaml', '1065884152.yaml', '1065885248.yaml', '1065885616.yaml']

    night2_obsids=['1066566392.yaml', '1066566512.yaml', '1066566632.yaml', '1066566760.yaml', '1066566880.yaml', '1066567000.yaml', '1066567120.yaml', '1066567248.yaml', '1066567368.yaml', '1066567488.yaml', '1066567608.yaml', '1066567736.yaml', '1066567976.yaml', '1066568096.yaml', '1066568224.yaml', '1066568344.yaml', '1066568464.yaml', '1066568584.yaml', '1066568712.yaml', '1066568832.yaml', '1066568952.yaml', '1066569072.yaml', '1066569200.yaml', '1066569320.yaml', '1066569440.yaml', '1066569560.yaml', '1066569688.yaml', '1066569808.yaml', '1066569928.yaml', '1066570048.yaml', '1066570176.yaml', '1066570296.yaml', '1066570416.yaml', '1066570536.yaml', '1066570664.yaml', '1066570784.yaml', '1066570904.yaml', '1066571024.yaml', '1066571152.yaml', '1066571272.yaml', '1066571392.yaml', '1066571512.yaml', '1066571640.yaml', '1066571760.yaml', '1066571880.yaml', '1066572000.yaml', '1066572128.yaml', '1066572248.yaml', '1066572368.yaml', '1066572488.yaml', '1066572616.yaml', '1066572736.yaml', '1066572856.yaml', '1066572976.yaml', '1066573104.yaml', '1066573224.yaml', '1066573344.yaml', '1066573464.yaml', '1066573592.yaml', '1066573712.yaml', '1066573832.yaml', '1066573952.yaml', '1066574080.yaml', '1066574200.yaml', '1066574320.yaml', '1066574440.yaml', '1066574568.yaml', '1066574688.yaml', '1066574808.yaml', '1066574928.yaml', '1066575176.yaml', '1066575296.yaml', '1066575416.yaml', '1066575544.yaml', '1066575664.yaml', '1066575784.yaml', '1066575904.yaml', '1066576032.yaml', '1066576152.yaml', '1066576272.yaml', '1066576392.yaml', '1066576520.yaml', '1066576640.yaml', '1066576760.yaml', '1066576880.yaml', '1066577008.yaml', '1066577128.yaml', '1066577248.yaml', '1066577368.yaml', '1066577496.yaml', '1066577616.yaml', '1066577736.yaml', '1066577856.yaml', '1066577984.yaml', '1066578104.yaml', '1066578224.yaml', '1066578344.yaml', '1066578472.yaml', '1066578592.yaml', '1066578712.yaml', '1066578832.yaml']

    trackable_sources = ['J004616-420739', 'J233426-412520', 'J232102-162302', 'J232519-120727', 'J232634-402715', 'J231956-272735', 'J003508-200354', 'J005906-170033']
    #no_amps = ['J002549-260211', 'J235701-344532', 'J002430-292847']
    night1_mdnmetric = 25.472890624999998
    night2_mdnmetric = 3.308683050847457

    for i in range(len(trackable_sources)):
        print('.......................................................This is source %s of %s'%(str(i+1), str(len(trackable_sources))))
        #trackable_sources = sources_to_track(obsidsets)
        #print(trackable_sources[:11])
        source_to_track = trackable_sources[i]
        scalevstime_plots(night1_obsids, night2_obsids, source_to_track, '%s.png'%(source_to_track), night1_mdnmetric, night2_mdnmetric, plottype = 1)

