from astropy.time import Time, TimeDelta
from astropy import coordinates as coord, units as u
from astropy.utils import iers
iers.conf.auto_download = False
import warnings
warnings.simplefilter("ignore")


def Get_Info(Header):
    Ra = Header['ALPHA'].split(' ')      # ALPHA, RA, CRVAL, OBJCTRA
    Ra = (int(Ra[0]) + int(Ra[1])/60. + int(Ra[2])/3600.)*15.
    Dec = Header['DELTA'].split(' ')     # DELTA, DEC, CRVAL, OBJCTDE
    Dec = int(Dec[0]) + int(Dec[1])/60. + int(Dec[2])/3600.
    
    X_pix = Header['XSTART'] + Header['NAXIS1'] / 2.
    Y_pix = Header['YSTART'] + Header['NAXIS2'] / 2.
    
    Exp = TimeDelta(Header['EXPTIME'], format='sec')
    t = Time(Header['DATE-OBS'], format='fits')
    t = t + Exp/2.
##    print(t.)
    
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