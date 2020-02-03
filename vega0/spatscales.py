import yaml
from yaml import SafeLoader as SafeLoader
import numpy as np
import seaborn as sns; sns.set()
import pandas as pd
import matplotlib.pyplot as plt

import plotly
#import plotly.io as pio
#pio.renderers.default = "pdf"
#import chart_studio.plotly
#from plotly.subplots import make_subplots
import plotly.graph_objects as go
#import plotly.express as px

from scipy.interpolate import griddata
from scipy import stats

import os
from os import listdir
from os.path import isfile, join

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
        #print(unpacked['sources'][source]['name'])
        names.append(unpacked['sources'][source]['name'])
        dec.append(unpacked['sources'][source]['dec'])
        ra.append(unpacked['sources'][source]['ra'])
        flux.append(unpacked['sources'][source]['flux_density']) 
        ampscales.append(np.nanmedian(unpacked['sources'][source]['amp_scales']))
        stds.append(np.nanstd(unpacked['sources'][source]['amp_scales']))
        #print(unpacked2['sources'][source]['amp_scales'])
    #ampscales = np.nan_to_num(ampscales)
    #print(len(ampscales))

    df = pd.DataFrame(list(zip(ra, dec, ampscales, stds, flux, names)), columns=['ra', 'dec', 'ampscales', 'stds', 'flux', 'names'])
    df = df.dropna(axis=0) #drop all rows with nan value in any column
    #df = df[(np.abs(stats.zscore(df.stds)) < 3).all(axis=1)] #discard rows with outlier stds (>3*\mu std)
    df = df[((df.stds - df.stds.mean()) / df.stds.std()).abs() < 3] #filter rows which dont have outlier stds (<3*\mu std)
    df = df.nlargest(500, 'flux', keep='all') #Get the rows with first n brightest sources.
    return df

def get_center(df):
    bulk_centre_ra = np.mean(df.ra)
    bulk_centre_dec = np.mean(df.dec)

    
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
    ra_p0 = [max(ra_hist), np.mean(ra_bins), 8]
    dec_p0 = [max(dec_hist), np.mean(dec_bins), 8]

    def gaussian(x, a, b, c):
        return a * np.exp(-(x-b)**2 / (2*c**2))
    try:
        ra_popt, _ = curve_fit(gaussian, ra_bins[:-1], ra_hist, p0=ra_p0)
        dec_popt, _ = curve_fit(gaussian, dec_bins[:-1], dec_hist, p0=dec_p0)
        radius = np.ceil((2*np.mean([abs(ra_popt[2]), abs(dec_popt[2])]))/2.5)*2.5
        # Check this radius against the extent of source available.
        if radius > max(ra) - min(ra) or radius > max(dec) - min(dec):
            radius = max([ra.max() - ra.min(), dec.max() - dec.min()])/2
    except:
        radius = max([df.ra.max() - df.ra.min(), df.dec.max() - df.dec.min()])/2
    print(radius)
    return radius


def spatial_scales_plot(df, obsid):
    print('trying plotly')
    size = df.stds
    fig = go.Figure(data=[go.Scatter(
        x=df.ra,
        y=df.dec,
        mode='markers',
        text = df.ampscales,
        hoverinfo='text',
        marker=dict(
            color=df.ampscales,
            colorscale='Magma_r',
            size = df.stds,
            #sizemode = 'diameter',
            showscale=True,
            sizeref = 2. * max(size) / (6**2)
            )
)])
    fig.update_layout(title = 'colour=median, size=std')
    fig.update_xaxes(autorange="reversed")
    #fig.show()
    fig.write_image('%s_amps_500.png' % (obsid))

def bgplot(radius, ra, dec, ampscales, stds, ra_centre, dec_centre, obsid, interp_method = 'linear', resolution=200):
    grid_x, grid_y = np.meshgrid(np.linspace(-radius, radius, resolution),
                                np.linspace(-radius, radius, resolution))

    grid_x += ra_centre
    grid_y += dec_centre

    grid_std = np.flipud(np.fliplr(griddata(np.vstack((ra, dec)).T, stds,
                                    (grid_x, grid_y), method=interp_method, fill_value=0)))
    grid_med = np.flipud(np.fliplr(griddata(np.vstack((ra, dec)).T, ampscales,
                                    (grid_x, grid_y), method=interp_method, fill_value=1)))

    beam_extent = (ra_centre+radius,
                    ra_centre-radius,
                    dec_centre-radius,
                    dec_centre+radius)
    print(beam_extent)
    crop_factor=1./np.sqrt(2)
    cropped_beam_extent = (ra_centre+radius*crop_factor,
                    ra_centre-radius*crop_factor,
                    dec_centre-radius*crop_factor,
                    dec_centre+radius*crop_factor)
    print(cropped_beam_extent)
    
    grid_med = cropper(grid_med)
    grid_std = cropper(grid_std)
    #grid_ddec = np.flipud(np.fliplr(griddata(np.vstack((df.ra, df.dec)).T, df.ampscales,
    #                                        (grid_x, grid_y), method=interp_method, fill_value=0)))


    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))
    fig.suptitle('%s amplitude scales: median (left), standard deviation (right)'%(obsid))
    img1 = ax1.imshow(grid_med, extent=cropped_beam_extent, cmap="plasma", origin="lower")
    fig.colorbar(img1, ax=ax1, format="%.2f",fraction=0.046, pad=0.04)
    img2 = ax2.imshow(grid_std, extent=cropped_beam_extent, cmap="plasma", origin="lower")
    fig.colorbar(img2, ax=ax2, format="%.2f",fraction=0.046, pad=0.04)
    plt.savefig('%s_cubic_interp500.png' % (obsid)) 
    '''
    fig, ax = plt.subplots(1, 1)
    img1 = ax.imshow(grid_std, vmin=0.0, vmax=0.8, extent=cropped_beam_extent, cmap="plasma", origin="lower")
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    fig.colorbar(img1, ax=ax, format="%.2f",fraction=0.046, pad=0.04)
    fig.suptitle('standard deviation of amplitude scales')
    plt.savefig('%s_stdamps_500brightest.png' % (obsid))
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
    for i in range(len(files[:3])):
        obsid = files[i].split('.')[0]
        data_file = '%s' % (path+files[i])
        df = loadfile(data_file) 
        radius, fra, fdec, f_ampscales, f_stds, ra_centre, dec_centre = get_center(df)
        #spatial_scales_plot(df, obsid)
        bgplot(radius, fra, fdec, f_ampscales, f_stds, ra_centre, dec_centre, obsid, interp_method = 'cubic')


