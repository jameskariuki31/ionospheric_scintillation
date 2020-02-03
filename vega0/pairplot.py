import yaml
import seaborn as sns
from yaml import CSafeLoader as SafeLoader
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import csv

yamlpath = '/home/chege/Desktop/curtin_work/vega'


def get_yamls(minn, maxx):
    minv = int(minn)
    maxv = int(maxx)
    yaml_data = sorted(os.listdir(yamlpath))[minv:maxv]
    print(len(yaml_data))
    return yaml_data


def get_mean_metric(night_obsids):
    obsidmetrics = []
    for obsid in night_obsids:
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
                obsidmetrics.append(np.nan)
                print('obsid without metric')
    return np.nanmean(obsidmetrics)


def get_values(night_obsids, source, meanmetric):
    ampscales = []
    std = []
    dates = []
    # flux = []
    for obsid in night_obsids:
        print('.....looking for %s in obsid %s' % (source, obsid))
        dates.append(obsid.split('.')[0])
        yaml_file = yamlpath+'/'+obsid
        with open(yaml_file, 'r') as f:
            unpacked = yaml.load(f, Loader=SafeLoader)
            for soc in unpacked['sources']:
                if unpacked['sources'][soc]['name'] == source:
                    print('.........found it')
                    print('.........getting its std and median ampscale')
                    scales = unpacked['sources'][soc]['amp_scales']
                    median_src_scal = np.nanmedian(scales)
                    print(median_src_scal)
                    std_src_scal = np.nanstd(scales)
                    # f_density = unpacked['sources'][soc]['flux_density']
        ampscales.append(median_src_scal)
        std.append(std_src_scal)
        # flux.append(f_density)
    return ampscales, std


if __name__ == '__main__':
    trackable_sources = ['J003508-200354']
    # trackable_sources = ['J004616-420739', 'J233426-412520', 'J232102-162302', 'J232519-120727', 'J232634-402715', 'J231956-272735', 'J003508-200354', 'J005906-170033']
    nightsbymetric = [[56, 101]]  # , [56, 155]]  # [[2047, 2147], [1781, 1850], [1675, 1730], [1991, 2047]]
    # nightsbymetric = [[2047, 2067], [1781, 1801], [1675, 1695], [1991, 2011]]

    for source in trackable_sources:
        allampscales = []
        allstds = []
        allmetrics = []
        metsname = ['nyt', 'six', 'eleven', 'twenty five']
        for i in range(len(nightsbymetric)):
            obsids = get_yamls(nightsbymetric[i][0], nightsbymetric[i][1])
            meanmetric = get_mean_metric(obsids)
            print('meanmetric', meanmetric)
            ampscales, stds = get_values(obsids, source, meanmetric)
            print(len(ampscales), len(stds))
            allampscales += ampscales
            allstds += stds
            nytmet = [metsname[i]] * len(ampscales)
            allmetrics += nytmet
            print(len(allampscales), len(allstds), len(allmetrics))
        allmetrics = [str(allmetrics[i]) for i in range(len(allmetrics))]
        df = pd.DataFrame(list(zip(allampscales, allstds, allmetrics)), columns=['Median', 'std', 'QAM'])
        df2 = df[df['std'] < 0.7]
        df2.dropna()
        # g = sns.pairplot(df2, hue='QAM', dropna=True)
        # plt.show()
        # plt.savefig('%s_pplot2_4QAMs.png' % (source))

        '''
        sns.set(style="darkgrid")
        sc = df2[df2.QAM == 'df2-nyt'].plot(kind='scatter', x='Median', y='std', color='red', label='qam nyt')
        df2[df2.QAM == 'df2-six'].plot(kind='scatter', x='Median', y='std', color='green', label='qam six', ax=sc)
        df2[df2.QAM == 'df2-eleven'].plot(kind='scatter', x='Median', y='std', color='orange', label='qam eleven', ax=sc)
        df2[df2.QAM == 'df2-twenty five'].plot(kind='scatter', x='Median', y='std', color='cyan', label='qam twenty five', ax=sc)
        sc.set_xlabel('Median')
        sc.set_ylabel('Standard deviations')
        sc.set_title('%s 4 qams Median Vs Standard deviations' % (source))
        sc = plt.gcf()
        sc.set_size_inches(10, 6)
        plt.savefig('%s_metclass_4QAMs.png' % (source))


        # Create a kde plot of median vs std for all qams/nights.
        print('plotting kde')
        sub = df2[df2['QAM'] == '%s' % (metsname[3])]
        sns.kdeplot(data=sub[['Median', 'std']], dropna=True, cmap="plasma",
                    shade=True, shade_lowest=False)
        plt.title('QAM %s' % (metsname[3]))
        plt.xlabel('Median amplitude scale')
        plt.ylabel('Standard deviation')
        plt.xlim()
        plt.savefig('test-%s_metkde_QAM.png' % (source))


        colors = ['Reds', 'Blues', 'amber', 'Cyans']
        for i in range(len(nightsbymetric)):
            nyt = df2[df2['QAM'] == '%s' % (metsname[i])]
            ax = sns.kdeplot(nyt.Median, nyt.std, cmap=colors[i], shade=True,
                            shade_lowest=False)
            plt.savefig('%s_3metskde.png' % (source))
        '''

        fig = sns.lmplot(x="Median", y="std", data=df2)
        plt.show()
