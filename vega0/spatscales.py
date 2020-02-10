import yaml
from yaml import SafeLoader as SafeLoader
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
# import plotly
# import plotly.io as pio
# pio.renderers.default = "pdf"
# import chart_studio.plotly
# from plotly.subplots import make_subplots
import plotly.graph_objects as go
# import plotly.express as px

from scipy.interpolate import griddata
from scipy.optimize import curve_fit
# from scipy import stats

import os
sns.set()
'''
Plot the median amplitude scales of each osurce in a single
yaml file as a function of ra nad dec.
Interpolate using scipy and try to obtain the background.
'''


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
        ampscales.append(
            np.nanmedian(unpacked['sources'][source]['amp_scales'][1:13]))
        stds.append(np.nanstd(unpacked['sources'][source]['amp_scales'][1:13]))
    # ampscales = np.nan_to_num(ampscales)
    # print(len(ampscales))

    df = pd.DataFrame(
        list(zip(ra, dec, ampscales, stds, flux, names)),
        columns=['ra', 'dec', 'ampscales', 'stds', 'flux', 'names'])
    df = df.dropna(axis=0)  # drop all rows with nan value in any column
    # df = df[(np.abs(stats.zscore(df.stds)) < 3).all(axis=1)
    df = df[((df.stds - df.stds.median()) / df.stds.std()).abs() < 3]  # filter
    # rows which dont have outlier stds (<3*\mu std)
    # Get the rows with N brightest sources.
    df = df.nlargest(400, 'flux', keep='all')
    print(df.flux)
    return df


def get_center(df):
    """returns filtered dimensions centered at the fitted MWA primary beam pointing center
    
    Arguments:
        df {pandas dataframe} -- dataframe with ra, dec, median amplitude scales and std
    
    Returns:
        tuple-- filtered values that are actually within the estimated MWA primary beam
    """
    bulk_centre_ra = np.median(df.ra)
    bulk_centre_dec = np.median(df.dec)

    radius = determine_radius(df)

    # Recalculate the centre, based on the sources within the radius,
    # and specify the sources to be used for analysis.
    filtered = np.array([[a, b, c, d] for a, b, c, d
                        in zip(df.ra, df.dec, df.ampscales, df.stds)
                        if abs(a-bulk_centre_ra) < radius
                        and abs(b-bulk_centre_dec) < radius])

    fra = filtered[:, 0]
    fdec = filtered[:, 1]
    f_ampscales = filtered[:, 2]
    f_stds = filtered[:, 3]
    ra_centre = fra.mean()
    dec_centre = fdec.mean()

    return radius, fra, fdec, f_ampscales, f_stds, ra_centre, dec_centre


def determine_radius(df):
    # Determine the radius of the primary beam.
    # Fit a Gaussian to the distribution of RA and Dec positions.
    # Use 2x the determined sigma as a cutoff for sources, in steps of 2.5 deg.
    ra_hist, ra_bins = np.histogram(df.ra, bins=50)
    dec_hist, dec_bins = np.histogram(df.dec, bins=50)
    ra_p0 = [max(ra_hist), np.median(ra_bins), 8]
    dec_p0 = [max(dec_hist), np.median(dec_bins), 8]

    def gaussian(x, a, b, c):
        return a * np.exp(-(x-b)**2 / (2*c**2))
    try:
        ra_popt, _ = curve_fit(gaussian, ra_bins[:-1], ra_hist, p0=ra_p0)
        dec_popt, _ = curve_fit(gaussian, dec_bins[:-1], dec_hist, p0=dec_p0)
        radius = np.ceil(
            (2*np.median([abs(ra_popt[2]), abs(dec_popt[2])]))/2.5)*2.5
        # Check this radius against the extent of source available.
        if radius > max(df.ra) - min(df.ra) or radius > max(df.dec) - min(df.dec):
            radius = max([df.ra.max() - df.ra.min(), df.dec.max() - df.dec.min()])/2
            print('gaussian fit done')
    except:
        radius = max([df.ra.max() - df.ra.min(), df.dec.max() - df.dec.min()])/2
    # print(radius)
    return radius


def spatial_scales_plot(ra, dec, ampscales, stds, obsid, beam_lim):
    print('trying plotly')
    size = stds
    fig = go.Figure(data=[go.Scatter(
        x=ra,
        y=dec,
        mode='markers',
        text=ampscales,
        hoverinfo='text',
        marker=dict(
            color=ampscales,
            colorscale='Magma_r',
            size=stds,
            # sizemode = 'diameter',
            showscale=True,
            sizeref=2. * max(size) / (6**2)
            )
    )])
    fig.update_layout(
                    title='colour=median, size=std',
                    xaxis=dict(range=[beam_lim[0], beam_lim[1]], title='Ra [deg]'),
                    yaxis=dict(range=[beam_lim[2], beam_lim[3]], title='Dec [deg]'))

    # fig.update_xaxes(autorange="reversed")
    # fig.show()
    fig.write_image('%s_median_amps1_13_400_innercrop.png' % (obsid))


def bgplot(
            radius, ra, dec, ampscales, stds, ra_centre, dec_centre, obsid,
            interp_method='linear', resolution=1000):
    grid_x, grid_y = np.meshgrid(
        np.linspace(-radius, radius, resolution),
        np.linspace(-radius, radius, resolution))

    grid_x += ra_centre
    grid_y += dec_centre

    grid_std = np.fliplr(
                        griddata(
                                np.vstack((ra, dec)).T, stds,
                                (grid_x, grid_y), method=interp_method,
                                fill_value=0))
    grid_med = np.fliplr(
                        griddata(
                                np.vstack((ra, dec)).T, ampscales,
                                (grid_x, grid_y), method=interp_method,
                                fill_value=1))

    beam_extent = (
                    ra_centre+radius,
                    ra_centre-radius,
                    dec_centre-radius,

                    dec_centre+radius
                )
    # print(beam_extent)
    crop_factor = 1./np.sqrt(2)
    cropped_beam_extent = (
                        ra_centre+radius*crop_factor,
                        ra_centre-radius*crop_factor,
                        dec_centre-radius*crop_factor,
                        dec_centre+radius*crop_factor)
    cropped_beam_extent = [10, -5, -35, -20]
    print(obsid, 'cropped beam limits', cropped_beam_extent)

    grid_med = cropper(grid_med)
    grid_std = cropper(grid_std)
    spatial_scales_plot(df.ra, df.dec, df.ampscales, df.stds, obsid, cropped_beam_extent)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))
    fig.suptitle(
        '%s amplitude scales: median (left), standard deviation \
        (right)' % (obsid))
    img1 = ax1.imshow(
        grid_med, extent=cropped_beam_extent, cmap="plasma", origin="lower")
    fig.colorbar(img1, ax=ax1, format="%.2f", fraction=0.046, pad=0.04)
    img2 = ax2.imshow(
        grid_std, vmin=0.0, vmax=0.4, extent=cropped_beam_extent, cmap="plasma", origin="lower")
    fig.colorbar(img2, ax=ax2, format="%.2f", fraction=0.046, pad=0.04)
    plt.savefig('%smedstd_lininterp400_1_13_innercrop.png' % (obsid))
    '''
    fig, ax = plt.subplots(1, 1)
    img1 = ax.imshow(grid_std, vmin=0.0, vmax=0.8, extent=cropped_beam_extent,
        cmap="plasma", origin="lower")
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    fig.colorbar(img1, ax=ax, format="%.2f",fraction=0.046, pad=0.04)
    fig.suptitle('standard deviation of amplitude scales')
    plt.savefig('%s_stdamps_400brightest.png' % (obsid))
    '''


def cropper(matrix, crop_factor=1./np.sqrt(2)):
    length = len(matrix)
    cropped_length = int(length * crop_factor)
    i1 = int((length-cropped_length)/2)
    i2 = length-int((length-cropped_length)/2)
    return matrix[i1:i2, i1:i2]


if __name__ == "__main__":
    path = '/home/chege/Desktop/curtin_work/vega/'
    files = sorted(os.listdir(path))
    for i in files[:25]:
        obsid = i.split('.')[0]
        data_file = '%s' % (path+i)
        df = loadfile(data_file)
        radius, fra, fdec, f_ampscales, f_stds, ra_centre, dec_centre \
            = get_center(df)
        # spatial_scales_plot(df.ra, df.dec, df.ampscales, df.stds, obsid)
        bgplot(
            radius, fra, fdec, f_ampscales, f_stds, ra_centre, dec_centre,
            obsid, interp_method='linear'
            )
