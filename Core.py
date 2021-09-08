import os

from astropy.io import fits


class Core(object):
    def __init__(self):
        pass


Path2Data = 'G:\Kislovodsk\T1b_J-20200912'
path2coo = 'GAIA'
# path2coo = 'UCAC'

######################################################################
# TODO set params
# set parameters for photometry
Saturation = 45000
FWHM_e = 2.0
Gain = 2.24 # ????
Rnoise = 10.
# RAper = 4.0
V_lim = 15.
MaxMag = 9.
######################################################################
# read directory and create list of fits-files
file_list = []
dir_content = os.listdir(Path2Data)
for ii in range(0, len(dir_content)):
    if dir_content[ii].count('.fit') or dir_content[ii].count('.fits') or dir_content[ii].count('.fts'):
        file_list.append(Path2Data + '/' + dir_content[ii])

# read first frame for object name and other information
file_name = file_list[0]
print('First frame: ' + file_name.split('/')[-1])
hdulist = fits.open(file_name)
Header = hdulist[0].header
# Data = hdulist[0].data.copy()
hdulist.verify('fix')
hdulist.close()
# Info = Get_Info(Header) TODO make get_info
