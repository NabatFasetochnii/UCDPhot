import os

from astropy.io import fits
from astropy.time import Time, TimeDelta
from astropy import coordinates as coord, units as u


class Core(object):

    def __init__(self):
        self.path2data = 'G:\Kislovodsk\T1b_J-20200912'
        self.path2coo = 'GAIA'  # 'UCAC'
        self.saturation = 45000  # ?
        self.FWHM_e = 2.0  # ?
        self.gain = 2.24  # ?
        self.read_noise = 10.  # ?
        self.mag_lim = 15.  # ?
        self.max_mag = 9.  # ?

    def open_fits(self):

        # read directory and create list of fits-files
        file_list = []
        dir_content = os.listdir(self.path2data)
        for ii in range(0, len(dir_content)):
            if dir_content[ii].count('.fit') or dir_content[ii].count('.fits') or dir_content[ii].count('.fts'):
                file_list.append(self.path2data + '/' + dir_content[ii])
        # read first frame for object name and other information
        file_name = file_list[0]
        print('First frame: ' + file_name.split('/')[-1])
        hdulist = fits.open(file_name)
        Header = hdulist[0].header
        # Data = hdulist[0].data.copy()
        hdulist.verify('fix')
        hdulist.close()
        # Info = Get_Info(Header) TODO make get_info

    def get_info(self, Header):
        Ra = Header['ALPHA'].split(' ')      # ALPHA, RA, CRVAL, OBJCTRA
        Ra = (int(Ra[0]) + int(Ra[1])/60. + int(Ra[2])/3600.)*15.
        Dec = Header['DELTA'].split(' ')     # DELTA, DEC, CRVAL, OBJCTDE
        Dec = int(Dec[0]) + int(Dec[1])/60. + int(Dec[2])/3600.

        X_pix = Header['XSTART'] + Header['NAXIS1'] / 2.
        Y_pix = Header['YSTART'] + Header['NAXIS2'] / 2.

        Exp = TimeDelta(Header['EXPTIME'], format='sec')
        t = Time(Header['DATE-OBS'], format='fits')
        t = t + Exp/2.

        Obj = coord.SkyCoord(Ra, Dec, unit=(u.deg, u.deg), frame='icrs')
        Site = coord.EarthLocation.from_geodetic(lon=Header['LONGITUD'],
                                                 lat=Header['LATITUDE'], height=Header['ALTITUDE'])

        helio = t.light_travel_time(Obj, 'heliocentric', location=Site)
        hjd = t + helio
        print('HJD: ', hjd.jd)

        bary = t.light_travel_time(Obj, location=Site)
        bjd = t + bary
        print('BJD:', bjd.jd)
        return Ra, Dec, Header['NAXIS1'], abs(Header['CD1_1']), X_pix, Y_pix, t, hjd.jd, bjd.tdb.jd