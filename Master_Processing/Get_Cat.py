from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import Angle
from astropy.table import Table
from numpy import arange


# from astropy.io import ascii


def Get_UCAC(RA, DEC, R, V_lim):
    My_Cat = Table()
    # make astropy SkyCoord object
    c = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree, frame='icrs')
    # set search cone
    a = Angle(R * u.deg)
    # set columns
    V = Vizier(columns=['RAJ2000', 'DEJ2000', 'Bmag', 'e_Bmag',
                        'Vmag', 'e_Vmag', 'gmag', 'e_gmag',
                        'rmag', 'e_rmag', 'imag', 'e_imag', "+_r"],
               column_filters={'Vmag': '>0', 'Vmag': '<' + str(V_lim)})
    # set limit of rows, sort by distance default
    V.ROW_LIMIT = 150
    # get data
    Vizier_result = V.query_region(c, radius=a, catalog=['I/322A'])
    if len(Vizier_result) != 0:
        Vizier_stars = Vizier_result[0]
        My_Cat['ID'] = arange(0, len(Vizier_stars), 1, 'int16')
        My_Cat['Ra'] = Vizier_stars['RAJ2000']
        My_Cat['Dec'] = Vizier_stars['DEJ2000']
        My_Cat['B'] = Vizier_stars['Bmag']
        My_Cat['V'] = Vizier_stars['Vmag']
        mag = Vizier_stars['rmag'] - 0.272 * (Vizier_stars['rmag'] - Vizier_stars['imag']) - 0.159
        My_Cat['R'] = mag
        mag = Vizier_stars['imag'] - 0.337 * (Vizier_stars['rmag'] - Vizier_stars['imag']) - 0.37
        My_Cat['I'] = mag

    #     ascii.write(My_Cat, 'My_Cat.txt', overwrite=True, delimiter='\t', format='commented_header')
    return My_Cat


def Get_GAIA(RA, DEC, R, R_lim):
    My_Cat = Table()
    c = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree, frame='icrs')
    a = Angle(R * u.deg)
    V = Vizier(columns=['RAJ2000', 'DEJ2000',
                        'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag',
                        'RPmag', 'e_RPmag', "+_r"],
               column_filters={'RPmag': '>0', 'RPmag': '<' + str(R_lim)})
    V.ROW_LIMIT = 150
    Vizier_result = V.query_region(c, radius=a, catalog=['I/345'])
    if len(Vizier_result) != 0:
        Vizier_stars = Vizier_result[0]
        #        print(Vizier_stars.info())
        My_Cat['ID'] = arange(0, len(Vizier_stars), 1, 'int16')
        My_Cat['Ra'] = Vizier_stars['RAJ2000']
        My_Cat['Dec'] = Vizier_stars['DEJ2000']
        My_Cat['Dist'] = Vizier_stars['_r']
        My_Cat['B'] = Vizier_stars['BPmag']
        My_Cat['R'] = Vizier_stars['RPmag']
        My_Cat['G'] = Vizier_stars['Gmag']

    #     ascii.write(My_Cat, 'My_Cat.txt', overwrite=True, delimiter='\t', format='commented_header')
    return My_Cat
# Get_GAIA(290.71129, 56.61503, 256*1.85/3600., 16.)
# Get_UCAC(290.71129, 56.61503, 256*1.85/3600., 16.)
