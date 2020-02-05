import os
import time
from datetime import datetime, timedelta
import yaml
from yaml import SafeLoader as SafeLoader
import matplotlib.pyplot as plt
import csv
import numpy as np
import run_ionfr
import pandas as pd
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns


path = '/home/chege/Desktop/curtin_work/vega/'


def loadfile(data_file):
    with open(data_file, 'r') as f:
        unpacked = yaml.load(f, Loader=SafeLoader)

    flux = []
    ra = []
    dec = []
    names = []
    ampscales = []
    stds = []
    sourcelist = unpacked['sources']
    for source in sourcelist:
        # print(unpacked['sources'][source]['name'])
        names.append(unpacked['sources'][source]['name'])
        dec.append(unpacked['sources'][source]['dec'])
        ra.append(unpacked['sources'][source]['ra'])
        flux.append(unpacked['sources'][source]['flux_density'])
        ampscales.append(np.nanmedian(unpacked['sources'][source]['amp_scales']))
        stds.append(np.nanstd(unpacked['sources'][source]['amp_scales']))
        # print(unpacked2['sources'][source]['amp_scales'])
    # ampscales = np.nan_to_num(ampscales)
    # print(len(ampscales))
    df = pd.DataFrame(
        list(zip(ra, dec, ampscales, stds, flux, names)),
        columns=['ra', 'dec', 'ampscales', 'stds', 'flux', 'names'])
    df2 = df.dropna(axis=0)
    return df2


def get_utc_time(obsid):
    # TAI scale with 1970-01-01 00:00:10 (TAI) epoch
    os.environ['TZ'] = 'right/UTC'
    time.tzset()
    gps_epoch_as_gps = datetime(1980, 1, 6)

    gps_time_as_gps = gps_epoch_as_gps + timedelta(seconds=float(obsid))
    gps_time_as_tai = gps_time_as_gps + timedelta(seconds=19)  # const offset
    tai_epoch_as_tai = datetime(1970, 1, 1, 0, 0, 10)
    # by definition
    tai_timestamp = (gps_time_as_tai - tai_epoch_as_tai).total_seconds()
    print(datetime.utcfromtimestamp(tai_timestamp))   # "right" timezone
    utc_time = datetime.utcfromtimestamp(tai_timestamp)

    utc_date = str(utc_time)[0:10]
    print(utc_date)
    return utc_date


def get_obsid_ionfr_txts(obsid):
    data_file = '%s' % (path+'%s.yaml' % (str(obsid)))
    df = loadfile(data_file)
    date = get_utc_time(obsid)
    for index, row in df.iterrows():
        run_ionfr.run_ionfr(
            row['ra'], row['dec'], row['names'], date, convertradec=True
            )


def get_fr_value(ionfr_txt_file, hour):
    with open(ionfr_txt_file) as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=' ')
        for row in csv_reader:
            if row[0] == str(hour):
                print('found the hr, reading its fr value')
                fr_value = row[3]
                fr_value_err = row[4]
    return fr_value, fr_value_err


def make_df(hour, csvfyl=None):
    if csvfyl is not None:  # if we have a csv file with ra dec, fr, and fr_err
        df = pd.read_csv(csvfyl)
        df = df.drop(df.columns[0], axis=1)
    else:
        txtdir = '/home/chege/Desktop/curtin_work/run_ionfr/singleLOScomparisons'
        fyls = os.listdir(txtdir)
        ras = []
        decs = []
        frs = []
        frs_errs = []
        for fyl in fyls:
            # print(fyl)
            if fyl.split('.')[-1] == 'txt':
                source = fyl.split('_')[0]
                ra, dec = run_ionfr.get_radec(source)
                ras.append(float(ra))
                decs.append(float(dec))
                fylpath = txtdir + '/' + fyl
                fr_value, fr_value_err = get_fr_value(fylpath, hour)
                frs.append(float(fr_value))
                frs_errs.append(float(fr_value_err))
        df = pd.DataFrame(
            list(zip(ras, decs, frs, frs_errs)),
            columns=['ra', 'dec', 'fr', 'fr_err'])
        print('made dataframe with radec and fr values')
    df = df.dropna(axis=0)

    return df


def get_center(df):
    bulk_centre_ra = np.mean(df.ra)
    bulk_centre_dec = np.mean(df.dec)

    radius = determine_radius(df)

    # Recalculate the centre, based on the sources within the radius,
    # and specify the sources to be used for analysis.
    filtered = np.array(
                        [[a, b, c, d] for a, b, c, d
                            in zip(df.ra, df.dec, df.fr, df.fr_err)
                            if abs(a-bulk_centre_ra) < radius
                            and abs(b-bulk_centre_dec) < radius])

    fra = filtered[:, 0]
    fdec = filtered[:, 1]
    f_fr = filtered[:, 2]
    f_fr_err = filtered[:, 3]
    ra_centre = fra.mean()
    dec_centre = fdec.mean()

    return radius, fra, fdec, f_fr, f_fr_err, ra_centre, dec_centre


def cropper(matrix, crop_factor=1./np.sqrt(2)):
    length = len(matrix)
    cropped_length = int(length * crop_factor)
    i1 = int((length-cropped_length)/2)
    i2 = length-int((length-cropped_length)/2)
    return matrix[i1:i2, i1:i2]


def determine_radius(df):
    # Determine the radius of the primary beam.
    # Fit a Gaussian to the distribution of RA and Dec positions.
    # Use 2x the determined sigma as a cutoff for sources, in steps of 2.5 deg.
    ra_hist, ra_bins = np.histogram(df.ra, bins=50)
    dec_hist, dec_bins = np.histogram(df.dec, bins=50)
    ra_p0 = [max(ra_hist), np.mean(ra_bins), 8]
    dec_p0 = [max(dec_hist), np.mean(dec_bins), 8]

    def gaussian(x, a, b, c):
        return a * np.exp(-(x-b)**2 / (2*c**2))
    try:
        ra_popt, _ = curve_fit(gaussian, ra_bins[:-1], ra_hist, p0=ra_p0)
        dec_popt, _ = curve_fit(gaussian, dec_bins[:-1], dec_hist, p0=dec_p0)
        radius = np.ceil(
            (2*np.mean([abs(ra_popt[2]), abs(dec_popt[2])]))/2.5)*2.5
        # Check this radius against the extent of source available.
        if radius > max(df.ra) - min(df.ra) or radius > max(df.dec) - min(df.dec):
            radius = max([df.ra.max() - df.ra.min(), df.dec.max() - df.dec.min()])/2
            print('gaussian fit done')
    except:
        radius = max([df.ra.max() - df.ra.min(), df.dec.max() - df.dec.min()])/2
    # print(radius)
    return radius


def interpfr(
            radius, ra, dec, fr, fr_err, ra_centre, dec_centre, obsid, hr,
            interp_method='linear', resolution=200):
    grid_x, grid_y = np.meshgrid(
                                np.linspace(-radius, radius, resolution),
                                np.linspace(-radius, radius, resolution))

    grid_x += ra_centre
    grid_y += dec_centre

    grid_fr = np.flipud(
                        np.fliplr(
                                griddata(
                                        np.vstack((ra, dec)).T, fr,
                                        (grid_x, grid_y), method=interp_method,
                                        fill_value=0)))

    beam_extent = (
                    ra_centre+radius,
                    ra_centre-radius,
                    dec_centre-radius,
                    dec_centre+radius)
    print(beam_extent)
    crop_factor = 1./np.sqrt(2)
    cropped_beam_extent = (
                        ra_centre+radius*crop_factor,
                        ra_centre-radius*crop_factor,
                        dec_centre-radius*crop_factor,
                        dec_centre+radius*crop_factor)
    print(cropped_beam_extent)

    # grid_fr = cropper(grid_fr)

    return grid_fr, cropped_beam_extent


def plot_interp_fr(grid, beam_extent):
    fig, ax = plt.subplots(1, 1)
    img1 = ax.imshow(
        grid, cmap="plasma", origin="lower")
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    fig.colorbar(img1, ax=ax, format="%.2f", fraction=0.046, pad=0.04)
    fig.suptitle('Spatial Faraday depth at %shrs %s' % (hr, str(obsid)))
    plt.savefig('%s_%shrs_ionfr.png' % (obsid, hr))


def plane_fit(grid):
    print('shape', grid.shape)
    m = grid.shape[0]  # size of the matrix

    X1, X2 = np.mgrid[:m, :m]

    fig = plt.figure()
    ax = fig.add_subplot(3, 1, 1, projection='3d')
    jet = plt.get_cmap('jet')
    ax.plot_surface(X1, X2, grid, rstride=1, cstride=1, cmap=jet, linewidth=0)

    # Regression
    X = np.hstack((np.reshape(X1, (m*m, 1)), np.reshape(X2, (m*m, 1))))
    X = np.hstack((np.ones((m*m, 1)), X))
    YY = np.reshape(grid, (m*m, 1))

    theta = np.dot(np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose()), YY)

    ax.scatter(grid[:, 0], grid[:, 1], grid[:, 2], c='r', s=20)

    plane = np.reshape(np.dot(X, theta), (m, m))

    ax = fig.add_subplot(3, 1, 2, projection='3d')
    ax.plot_surface(X1, X2, plane)
    # ax.plot_surface(X1, X2, grid, rstride=1, cstride=1, cmap=jet, linewidth=0)

    # Subtraction
    grid_sub = grid - plane
    ax = fig.add_subplot(3, 1, 3, projection='3d')
    ax.plot_surface(X1, X2, grid_sub, rstride=1, cstride=1, cmap=jet, linewidth=0)

    plt.show()


def linefit(grid, beam_lim):
    y = grid[:, 90]
    print(y.shape)
    x = np.linspace(beam_lim[2], beam_lim[3], len(y))
    regplot = sns.regplot(x=x, y=y, color="g")
    # plt.show()
    fig = regplot.get_figure()
    fig.savefig('1065880128_14hrs_1ra_ionfr.png')


def line_fit(df):
    y = grid[:, 90]
    print(y.shape)
    x = np.linspace(beam_lim[2], beam_lim[3], len(y))
    regplot = sns.regplot(x=x, y=y, color="g")
    # plt.show()
    fig = regplot.get_figure()
    fig.savefig('1065880128_14hrs_1ra_ionfr.png')


if __name__ == '__main__':
    # files = sorted(os.listdir(path))
    # for i in range(len(files)):
    # obsid = files[i].split('.')[0]
    # get_obsid_ionfr_txts(1065880128)
    obsid = 1065880128
    hr = 14
    df = make_df(hr, csvfyl='1065880128_ionfr_p3.csv')
    print(df.head())
    # df.to_csv('1065880128_ionfr_p3.csv')
    radius, fra, fdec, f_fr, f_fr_err, ra_centre, dec_centre = get_center(df)
    grid, beam_extent = interpfr(
        radius, fra, fdec, f_fr, f_fr_err, ra_centre, dec_centre, obsid, hr)
    plot_interp_fr(grid, beam_extent)
    # plane_fit(grid)
    # linefit(grid, beam_extent)
