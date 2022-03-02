import os

import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time, TimeDelta
from astroquery.vizier import Vizier
from numpy import arange
from photutils.detection import DAOStarFinder


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

    for file in dir_content:
        if (file.count(".fits") or file.count(".fit")) and file.count(image_filter):
            image_list.append(file)

    return image_list


def get_header_info(header):
    ra = header['CURRA'].split(' ')  # ALPHA, RA, CRVAL, OBJCTRA
    ra = (int(ra[0]) + int(ra[1]) / 60. + int(ra[2]) / 3600.) * 15.
    dec = header['CURDEC'].split(' ')  # DELTA, DEC, CRVAL, OBJCTDE
    dec = int(dec[0]) + int(dec[1]) / 60. + int(dec[2]) / 3600.

    x_pix = header['X'] + header['NAXIS1'] / 2.
    y_pix = header['Y'] + header['NAXIS2'] / 2.

    exp = TimeDelta(header['EXPTIME'], format='sec')
    t = Time(header['DATE-OBS'], format='fits')
    t = t + exp / 2.

    obj = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    site = coord.EarthLocation.from_geodetic(lon=header['LONGITUD'],
                                             lat=header['LATITUDE'], height=header['ELEVAT'])

    helio = t.light_travel_time(obj, 'heliocentric', location=site)
    hjd = t + helio

    bary = t.light_travel_time(obj, location=site)
    bjd = t + bary

    return ra, dec, header['NAXIS1'], abs(header['CD1_1']), x_pix, y_pix, t, hjd.jd, bjd.tdb.jd


def get_center(path):
    hdu = fits.open(path)[0]
    data = np.log10(hdu.data)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    # _max = median + 10. * std
    # _min = median - 1. * std
    dao_find = DAOStarFinder(fwhm=5.0, threshold=5. * std)
    sources = dao_find(data - median)
    # for col in sources.colnames:
    #     sources[col].info.format = '%.8g'  # for consistent table output
    # np.transpose((sources['xcentroid'], sources['ycentroid']))
    sources.sort('flux')
    sources.reverse()
    print(sources)
    return sources
    # apertures = CircularAperture(positions, r=6.)
    # norm = ImageNormalize(vmin=_min, vmax=_max)
    # plt.imshow(data, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
    # apertures.plot(color='blue', lw=1.5, alpha=0.3)
    # plt.show()
