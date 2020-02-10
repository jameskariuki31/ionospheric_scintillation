import os
import datetime
from dateutil import tz
from pytz import timezone
from rad2hmsdms import rad2hmsdms
import urllib.request
import yaml
from yaml import SafeLoader as SafeLoader
import matplotlib.pyplot as plt
import csv
import numpy as np
from matplotlib import gridspec

from_zone = tz.gettz('UTC')
to_zone = tz.gettz('Australia/Perth')


def get_ionex(date):
    '''
    Get the name of the required ionex file derived from the date of the observation
    Download it from CODE website
    Uncompress it
    '''
    year = int(date.split('-')[0])
    month = int(date.split('-')[1])
    day = int(date.split('-')[2])

    dayofyear = datetime.datetime.strptime(''+str(year)+' '+str(month)+' '+str(day)+'', '%Y %m %d').timetuple().tm_yday

    if dayofyear < 10:
        dayofyear = '00'+str(dayofyear)
    if dayofyear < 100 and dayofyear >= 10:
        dayofyear = '0'+str(dayofyear)

    ionex = 'CODG'+str(dayofyear)+'0.'+str(list(str(year))[2])+''+str(list(str(year))[3])+'I'
    print ('file needed:', ionex)

    path = os.path.abspath('')
    found = os.listdir(path)
    if ionex not in found:
        url = 'http://ftp.aiub.unibe.ch/CODE/' + str(year)
        urllib.request.urlretrieve('%s/%s.Z' % (url, ionex), '%s.Z' % (ionex))
        print('file downloaded')

        os.system('uncompress %s.Z' % (ionex))

    return ionex


def get_radec(source, obsid):
    '''
    Gets the Ra and Dec of a source from the yaml file who's path is given.
    '''
    yaml_file = '/home/chege/Desktop/curtin_work/vega/%s.yaml' % (obsid)
    with open(yaml_file, 'r') as f:
        unpacked = yaml.load(f, Loader=SafeLoader)
        for soc in unpacked['sources']:
            if unpacked['sources'][soc]['name'] == source:
                print('.........found source')
                print('................getting its ra & dec')
                ra = unpacked['sources'][soc]['ra']
                dec = unpacked['sources'][soc]['dec']
                ampscale = np.nanmedian(unpacked['sources'][soc]['amp_scales'][1:13])
                std = np.nanstd(unpacked['sources'][soc]['amp_scales'][1:13])
                flux = unpacked['sources'][soc]['flux_density']
    return ra, dec, flux, ampscale, std


def utc_local(hr_list, day):
    '''
    Converts the utc hours obtained from IonFR output textfile into local time.
    '''
    local_hrs = []
    for hr in hr_list:
        date_str = '%s %s:00:00' % (str(day), str(int(hr)))
        datetime_obj = datetime.datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
        datetime_obj_utc = datetime_obj.replace(tzinfo=timezone('UTC'))
        # hr = datetime_obj_utc.strftime("%Y-%m-%d %H:%M:%S %Z%z")
        loc_hr = datetime_obj_utc.astimezone(to_zone)
        loc_hr = loc_hr.strftime("%Y-%m-%d %H:%M:%S %Z%z")
        # print(loc_hr[11:13])
        local_hrs.append(float(loc_hr[11:13]))
    print(local_hrs)
    return local_hrs


def farady_plot(ionrm_txtfile1, ionrm_txtfile2, plotname, source, date1, date2):
    '''
    Make a plot comparing two IonFR output textiles for a given source/LoS, and dates.
    '''
    err1 = []
    err2 = []
    farot1 = []
    farot2 = []
    hour1 = []
    hour2 = []
    try:
        with open(ionrm_txtfile1) as csvfile1:
            csv_reader1 = csv.reader(csvfile1, delimiter=' ')
            for row in csv_reader1:
                hour1.append(row[0])
                farot1.append(row[3])
                err1.append(row[4])

        with open(ionrm_txtfile2) as csvfile2:
            csv_reader2 = csv.reader(csvfile2, delimiter=' ')
            for row in csv_reader2:
                hour2.append(row[0])
                farot2.append(row[3])
                err2.append(row[4])

        err1 = [float(i) for i in err1]
        err2 = [float(i) for i in err2]
        farot1 = [float(i) for i in farot1]
        farot2 = [float(i) for i in farot2]
        hour1 = [float(i)+8 for i in hour1]
        hour2 = [float(i)+8 for i in hour2]

        # get the percentage FR differences between the two days
        perc_change = []
        common_hours = []
        for hr in hour2:
            if hr in hour1:
                common_hours.append(hr)
                ind1 = hour1.index(hr)
                ind2 = hour2.index(hr)
                print(farot2[ind2],farot1[ind1])
                change = (float(farot2[ind2]-farot1[ind1])/float(farot1[ind1]))*100
                perc_change.append(change)
        # mean_perc_change = round(np.mean(perc_change),2)
        print(perc_change)
        print(common_hours)

        # change utc time to local time
        # hour1 = utc_local(hour1, date1)
        # hour2 = utc_local(hour2, date2)

        fig = plt.figure(figsize=(10, 12))
        gs1 = gridspec.GridSpec(2, 1, height_ratios=[2, .5])
        gs1.update(wspace=0.0, hspace=0.006)

        ax0 = plt.subplot(gs1[0])
        line1, = ax0.plot(hour1, farot1, marker='o', linestyle=':', color='g')
        ax0.errorbar(hour1, farot1, yerr=err1, linestyle='None', color='g')
        line2, = ax0.plot(hour2, farot2, marker='o', linestyle=':', color='m')
        ax0.errorbar(hour2, farot2, yerr=err2, linestyle='None', color='m')
        ax0.axvline(x=21)
        ax0.axvline(x=27)
        # ax0.set_xlabel('Time (Hours)')
        ax0.set_ylabel(r'$\phi_{ion}  (rad m^{-2})$')
        ax0.set_ylim([0, 11.6])
        # ax0.set_xlim([xmin=14, xmax=32])
        ax0.set_xbound(lower=14, upper=32)
        ax0.set_title('Faraday depth values for source %s on %s and %s'%(source, date1, date2))
        ax0.grid(b=True, which='major', axis='both', color='#666666', linestyle='-')
        ax0.grid(True, which='minor', axis='both', color='#999999', linestyle='-', alpha=0.4)
        ax0.minorticks_on()
        ax0.legend((line1, line2), ('%s: qam = 11' % (date1), '%s: qam = 3' % (date2)), loc='upper right')

        '''
        ax1 = fig.add_subplot(121)
        ax1.plot(hour1, farot1, marker='o',linestyle = ':', color='r')
        ax1.errorbar(hour1, farot1, yerr=err1, linestyle='None')
        ax1.set_xlabel('Time (Hours)')
        text1 = ax1.annotate("Mean attenuation: %s %%"%(str(mean_perc_change)), xy=(1, 8.5),color='m')
        text1.set_fontsize(20)
        ax1.set_ylabel(r'$\phi_{ion}  (rad m^{-2})$')
        ax1.set_ylim([0, 11.6])
        ax1.set_xlim([14, 32])
        ax1.set_title('%s %s mean qam = 25.47'%(source, date1))
        ax1.grid(b=True, which='major',axis='both',color='#666666', linestyle='-')
        ax1.grid(True, which='minor',axis='both', color='#999999', linestyle='-', alpha=0.4)
        ax1.minorticks_on()
        '''
        ax1 = plt.subplot(gs1[1], sharex=ax0)
        line3, = ax1.plot(common_hours, perc_change, marker='o', linestyle='-', color='c')
        ax1.axvline(x=21)
        ax1.axvline(x=27)
        ax1.set_xbound(lower=14, upper=32)
        ax1.set_xlabel('Time (Hours)')
        ax1.set_ylabel('%')
        ax1.grid(b=True, which='major', axis='both', color='#666666', linestyle='-')
        ax1.grid(b=True, which='minor', axis='both', color='#999999', linestyle='-', alpha=0.4)
        ax1.minorticks_on()
        ax1.legend('%% difference', loc='upper left')

        # yticks0 = ax0.yaxis.get_major_ticks()
        # yticks0[0].label1.set_visible(False)
        xticks = ax0.xaxis.get_major_ticks()
        xticks[0].label1.set_visible(False)
        xticks[-1].label1.set_visible(False)

        # ax0.set_xticklabels([])

        '''
        ax2 = fig.add_subplot(122)
        ax2.plot(hour2, farot2, marker='o',linestyle = ':', color='r')
        ax2.errorbar(hour2, farot2, yerr=err2, linestyle='None')
        ax2.set_xlabel('Time (Hours)')
        ax2.set_ylabel(r'$\phi_{ion}  (rad m^{-2})$')
        ax2.set_ylim([0, 11.6])
        ax2.set_xlim([14, 32])
        ax2.set_title('%s %s mean qam = 3.31'%(source, date2))
        ax2.grid(b=True, which='major',axis='both',color='#666666', linestyle='-')
        ax2.grid(b=True, which='minor',axis='both',color='#999999', linestyle='-', alpha=0.4)
        ax2.minorticks_on()
        '''
        # fig.tight_layout()
        # fig.text(.5, 4.5, 'Percentage difference: %s'%(mean_perc_change), ha='center')
        # fig.savefig('met3vs11%s.png' % (plotname))
        # plt.show()

    except FileNotFoundError:
        pass


def run_ionfr(ra, dec, source, date, time='00:00:00', sign='-', convertradec=True, lat='26d42m12ss', lon='116d40m15se'):
    '''
    Runs ionfr code given the inputs
    Returns .txt files named after each source(los) and the date
    '''
    ionex = get_ionex(date)

    scrptpath = '/home/chege/Desktop/curtin_work/ionFR/ionFRM.py'

    if convertradec:
        if ra < 0:
            print('negative ra')
            ra += 360
        ra = rad2hmsdms(ra, Type="ra", deg=True)
        print(ra)
        ra = ra[:2]+'h'+ra[3:5]+'m'+ra[6:]+'s'
        print(ra)
        dec = rad2hmsdms(dec, Type="dec", deg=True)
        # if int(DEC[0:2]) > 0:
        #    sign = '+'
        # else:
        #    sign = '-'
        dec = dec[1:3]+'d'+dec[4:6]+'m'+dec[7:]+'s'
        print(dec)

    los = str(ra) + sign + str(dec)
    print(los)

    datef = str(date)+'T'+time

    runionfr = "%s %s %s %s %s %s" % (scrptpath, los, lat, lon, datef, ionex)
    os.system(runionfr)
    os.system('mv IonRM.txt %s_%s.txt' % (source, date))  # rename the textfile made


if __name__ == '__main__':
    dates = ['2013-10-10', '2013-10-18']
    trackable_sources = ['J004616-420739', 'J233426-412520', 'J232102-162302', 'J232519-120727', 'J232634-402715', 'J231956-272735', 'J003508-200354', 'J005906-170033']

    # for date in dates:
    #    for source in trackable_sources:
    #        ra, dec, flux, ampscale, std = get_radec(source, obsid)
    #        run_ionfr(ra, dec, source, date, convertradec = True)

    # make a plots
    for source in trackable_sources:
        farady_plot('%s_%s.txt' % (source, dates[0]), '%s_%s.txt' % (source, dates[1]), 'ionfr_%s_%s_%s' % (source, dates[0], dates[1]), source, dates[0], dates[1])
