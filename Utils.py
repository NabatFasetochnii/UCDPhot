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
from scipy.signal import medfilt

warnings.simplefilter("ignore")


def get_delta_transit(time, planet):
    from pytransit import RoadRunnerModel
    tm = RoadRunnerModel('power-2')
    tm.set_data(time=time)

    # TRAPPIST-1c
    if planet == 'c':
        k = 0.08440
        ldc = np.array([0.1, 0.0])
        p = 2.421937
        a = 28.549
        i = np.deg2rad(89.778)
        e = 0.00654
        w = np.deg2rad(282.45)
    elif planet == 'b':
        k = 0.08590
        ldc = np.array([0.1, 0.0])
        p = 1.510826
        a = 20.843
        i = np.deg2rad(89.728)
        e = 0.00622
        w = np.deg2rad(336.86)
    else:
        print('\"planet\" can only be b or c')
        return
    fake_curve = tm.evaluate_ps(t0=np.median(time), k=k, ldc=ldc, p=p, a=a, i=i, e=e, w=w)
    return 1 - fake_curve


def med_filt(D):
    return medfilt(D, [3, 1])


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
    # site = coord.EarthLocation.from_geodetic(lon=header['LONGITUD'],
    #                                          lat=header['LATITUD'], height=header['ELEVAT'])
    site = coord.EarthLocation.from_geodetic(lon=42.6675, lat=43.73611, height=2112.0)

    frame = coord.AltAz(obstime=t, location=site)
    obj_frame = obj.transform_to(frame)
    airmass = obj_frame.secz.value

    helio = t.light_travel_time(obj, 'heliocentric', location=site)
    hjd = t + helio

    bary = t.light_travel_time(obj, location=site)
    bjd = t + bary

    return ra.degree, dec.degree, x_pix, y_pix, t, hjd.jd, bjd.tdb.jd, airmass


def get_center(data, mask=None):
    #  gaussian convolution
    kernel = Gaussian2DKernel(x_stddev=1)
    data = convolve(data, kernel)

    #  extract background
    data -= np.median(data)
    Bkg_sigma = mad_std(data)

    #  mask bad row
    if not mask:
        mask = np.zeros(data.shape, dtype=bool)
        mask[90:110, 0:data.shape[1]] = True  # only for ASTRONIRCAM

    daofind = DAOStarFinder(fwhm=6.5, threshold=10. * Bkg_sigma, sharplo=0.25)
    sources = daofind(data, mask=mask)

    # Sort sources in ascending order
    sources.sort('flux')
    sources.reverse()
    return sources