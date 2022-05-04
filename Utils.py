import os
import warnings

import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.stats import mad_std
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astroquery.vizier import Vizier
from numpy import arange
from photutils.detection import DAOStarFinder
from pylab import indices
from pylab import ravel
from scipy import optimize

warnings.simplefilter("ignore")


# def PSF2FWHM(PSF_model):
#     try:
#         phi = np.arctan((PSF_model[2] ** 2.0) / (PSF_model[0] ** 2.0 - PSF_model[1] ** 2.0)) / 2.0
#         alpha1 = np.sqrt(
#             2.0 / (PSF_model[0] ** 2.0 + PSF_model[1] ** 2.0 + PSF_model[2] ** 2.0 / np.sin(2.0 * phi)))
#         alpha2 = np.sqrt(np.fabs(1 / (PSF_model[0] ** 2.0 + PSF_model[1] ** 2.0 - 1.0 / (alpha1 ** 2.0))))
#         fwhm1 = 2.0 * alpha1 * np.sqrt(2.0 ** (1.0 / (PSF_model[3])) - 1.0)
#         fwhm2 = 2.0 * alpha2 * np.sqrt(2.0 ** (1.0 / (PSF_model[3])) - 1.0)
#         return fwhm1, fwhm2, phi, PSF_model[4]
#     except:
#         print('Wrong PSF model')
#         return 0., 0., 0., 0.


def get_2mass(ra, dec, r, j_lim, cat):
    my_cat = Table()
    position = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    # print('position ', position.ra.hms, position.dec.dms)
    angle = Angle(r * u.deg)
    v = Vizier(columns=['RAJ2000', 'DEJ2000',
                        'Jmag', 'e_Jmag', "+_r"],
               # column_filters={'Jmag': '>0', 'Jmag': '<' + str(j_lim)})
               column_filters={'Jmag': '<' + str(j_lim)})
    v.ROW_LIMIT = 150
    vizier_result = v.query_region(position, radius=angle, catalog=[cat])
    if len(vizier_result) != 0:
        vizier_stars = vizier_result[0]
        #        print(vizier_stars.info())
        my_cat['ID'] = arange(0, len(vizier_stars), 1, 'int16')
        my_cat['Ra'] = vizier_stars['RAJ2000']
        my_cat['Dec'] = vizier_stars['DEJ2000']
        my_cat['Dist'] = vizier_stars['_r']
        my_cat['J'] = vizier_stars['Jmag']
        my_cat['J_err'] = vizier_stars['e_Jmag']
    else:
        print('I can\'t found the stars near the target\n')
    return my_cat


def get_image_list(image_dir, image_filter):
    dir_content = os.listdir(image_dir)
    image_list = []
    # files = glob.glob(r'G:\Kislovodsk\T1b_J-20200912\*.fts')
    for file in dir_content:
        if (file.count(".fits") or file.count(".fit") or file.count(".fts")) and file.count(image_filter):
            image_list.append(file)

    return image_list


def get_header_info(header):

    ra = Angle(header['CURRA'] + ' hours')
    dec = Angle(header['CURDEC'] + ' degrees')

    x_pix = header['X'] + header['NAXIS1'] / 2.
    y_pix = header['Y'] + header['NAXIS2'] / 2.

    exp = TimeDelta(header['EXPTIME'], format='sec')
    t = Time(header['DATE-OBS'], format='fits')
    t = t + exp / 2.

    obj = coord.SkyCoord(ra.degree, dec.degree, unit=(u.deg, u.deg), frame='icrs')
    site = coord.EarthLocation.from_geodetic(lon=header['LONGITUD'],
                                             lat=header['LATITUD'], height=header['ELEVAT'])

    helio = t.light_travel_time(obj, 'heliocentric', location=site)
    hjd = t + helio

    bary = t.light_travel_time(obj, location=site)
    bjd = t + bary

    return ra.degree, dec.degree, x_pix, y_pix, t, hjd.jd, bjd.tdb.jd


def get_center(data):
    ##read file, copy data and header
    # hdulist = fits.open(path)
    # # Header = hdulist[0].header.copy()
    # data = hdulist[0].data
    # hdulist.verify('fix')

    ##gaussian convolution
    kernel = Gaussian2DKernel(x_stddev=1)
    data = convolve(data, kernel)

    ##extract background
    data -= np.median(data)
    Bkg_sigma = mad_std(data)

    ##mask bad row
    mask = np.zeros(data.shape, dtype=bool)
    mask[90:110, 0:data.shape[1]] = True

    daofind = DAOStarFinder(fwhm=4.5, threshold=5. * Bkg_sigma, sharplo=0.25)
    sources = daofind(data, mask=mask)
    # print(Sources.info)

    # plt.imshow(Data, cmap=cm.Greys_r, aspect='equal',
    #            norm=Normalize(vmin=-30, vmax=150), interpolation='nearest')
    # plt.scatter(Sources['xcentroid'], Sources['ycentroid'], s=40, facecolors='none', edgecolors='r')
    # plt.show()

    # Sort sources in ascending order
    sources.sort('flux')
    sources.reverse()
    return sources


##################################################################
# calculate center of mass
def centroid(R1, R2, R3, arr):
    # total = 0.
    Ry = arr.shape[0] / 2
    Rx = arr.shape[1] / 2

    # mask
    X_index = np.arange(0, arr.shape[1], 1)  # index array
    Y_index = np.arange(0, arr.shape[0], 1)  # index array
    distance = np.sqrt(np.power(np.ones(arr.shape) * (X_index[None, :] - Rx), 2) + np.power(
        np.ones(arr.shape) * (Y_index[:, None] - Ry), 2))  # distance array

    # mean sky
    annulus_mask = np.copy(distance)
    annulus_mask[annulus_mask < R2] = 0.
    annulus_mask[annulus_mask > R3] = 0.
    annulus_mask[annulus_mask > 0] = 1.
    masked = arr * annulus_mask
    MSky = np.median(masked[np.nonzero(masked)])
    MSky = np.nan_to_num(MSky)

    # centroid
    # aperture_mask = np.copy(distance)
    distance[distance <= R1] = 1.
    distance[distance > R1] = 0.
    masked = arr * distance
    total = np.sum(masked)

    X = np.sum(masked * X_index[None, :]) / total
    Y = np.sum(masked * Y_index[:, None]) / total
    return X - arr.shape[1] / 2, Y - arr.shape[0] / 2, MSky


##################################################################
# D2 moffat fitter
def D2_moffat_full(A, B, C, D, E, F, x0, y0):  # B=1/sigma_x^2, C=1/sigma_y^2, E=beta
    try:
        return moff(A, B, C, D, E, F, x0, y0)
    except Exception as e:
        print('except in moffat_full ', e)
        return None


# read for correction of invalid value encountered in power
# http://stackoverflow.com/questions/16990664/scipy-minimize-uses-a-nonetype


####################################################################
# def D2_gauss(A, B, C, D, x0, y0):
#     return lambda y,x: A*np.exp(-(((x0-x)/B)**2 +((y0-y)/C)**2)/2) + D
##################################################################


def D2_moffat_fitter(ROI, MSKY, x_coo, y_coo, R3):
    x0 = x_coo - np.floor(x_coo) + R3
    y0 = y_coo - np.floor(y_coo) + R3

    #     try:
    #  moffat
    params = (ROI.max(), 0.3, 0.3, 0.1, 5.0, MSKY, x0, y0)
    errorfunction = lambda p: ravel(D2_moffat_full(*p)(*indices(ROI.shape)) - ROI)
    p, success = optimize.leastsq(errorfunction, params, maxfev=1000, ftol=0.05)

    return p[1], p[2], p[3], p[4], p[5]  # , w


def moff(A, B, C, D, E, F, x0, y0):
    return lambda y, x: A * (1 + ((x - x0) * B) ** 2. +
                             ((y - y0) * C) ** 2. + ((x - x0) * (y - y0) * (D ** 2.))) ** (-E) + F


def get_PSF(data, XY_coo, fwhm):
    PSF_model = []
    R1 = np.ceil(fwhm)  # estimation of aperture radii
    R2 = np.ceil(fwhm * 3.)  # sky annulus inner radii
    R3 = np.ceil(fwhm * 6.)  # sky annulus outer radii

    for ii in range(0, len(XY_coo)):  # for every star from coo file
        x_coo = XY_coo[ii, 0]
        y_coo = XY_coo[ii, 1]
        ROI = np.copy(data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])  # copy small area
        offset = centroid(R1, R2, R3, ROI)  # search centroid, Gauss sigma and mean sky

        if np.isnan(offset[0]) is False and np.isnan(offset[1]) is False:
            x_coo = x_coo + offset[0]
            y_coo = y_coo + offset[1]
            MSKY = offset[2]
            ROI = np.copy(data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])
            param = D2_moffat_fitter(ROI, MSKY, x_coo, y_coo, R3)  # fit 2D moffat psf
            if param is not None:
                PSF_model.append(param)
        else:
            pass
    PSF_model = np.asarray(PSF_model)
    PSF_model = np.median(PSF_model, 0)
    return PSF_model
