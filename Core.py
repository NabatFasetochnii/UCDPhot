import os
from datetime import datetime
import numpy as np
import astropy.io.fits as fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astroquery.astrometry_net import AstrometryNet
import matplotlib.pyplot as plt
import time

from photutils import CircularAperture

import Utils


class UCDPhot(object):
    pars_names = {
        "path2data": "Images directory",
        "path2save": "Report directory",
        "filter": "Images filter",
        'saturation': 'Saturation',
        "cat": "Catalog",
        "mag_lim": "Maximum magnitude limit",
        'max_mag': 'max_mag',
        # Only stars brighter than mag_lim are included
        # "image_edge": "Images edge",

        "aperture": "Aperture radius",
        # "search_box": "Centroid search box",
        "gain": "CCD gain",
        "ron": "CCD readout noise",

        "mmd": "Ensemble stars maximal magnitude difference",
        "isr": "Ensemble stars initial search radius in arcminutes",
        "msr": "Ensemble stars maximum search radius in arcminutes",
        "fop": "Field of wive in arcsec",

        "std_lim": "Ensemble stars standard deviation limit / (S/N)"
    }

    def __init__(self, path2data=r'G:\Kislovodsk\T1b_J-20200912', **kwargs):

        self.pars = {"filter": kwargs.get("filter", "J"),
                     "aperture": kwargs.get("aperture", 10),
                     'saturation': kwargs.get("saturation", 58000),
                     "path2data": path2data,
                     'mag_lim': kwargs.get("mag_lim", "20"),  # ?
                     'max_mag': kwargs.get("max_mag", "9"),  # ?
                     "gain": kwargs.get("gain", 2.18),
                     "ron": kwargs.get("ron", 0.02),
                     "mmd": kwargs.get("mmd", 2.0),
                     "isr": kwargs.get("isr", 5),
                     "msr": kwargs.get("msr", 30),
                     "std_lim": kwargs.get("std_lim", 3),
                     "fop": kwargs.get("fop", 4.6),
                     'cat': 'II/246',
                     'astrometry.net_API': ''}

        # self.path2coo = 'GAIA'  # 'UCAC'
        # self.saturation = 58000  #
        # self.FWHM_e = 2.0  # ?
        # # self.gain = 2.18  #
        # # self.read_noise = 0.02  #
        # self.mag_lim = 18.  # ?
        # self.max_mag = 9.  # ?

        # make data paths
        if path2data[-1] == '/':
            self.pars['path2save'] = path2data[:-1] + '_report'
        else:
            self.pars['path2save'] = path2data + '_report'
        if not os.path.exists(self.pars['path2save']):
            os.makedirs(self.pars['path2save'])

        # self.pars["ext_catalog"] = kwargs.get("ext_catalog", "II/336")
        # self.pars["image_edge"] = kwargs.get("image_edge", 100)

        # self.log_pars()

    def log_pars(self):
        now = datetime.now()
        with open(self.pars['path2save' + '/log_pars'], "a") as logs_file:
            logs_file.write(f"Current date and time\t{now}\n")
            logs_file.write("New set of parameters is loaded:\n")

            for key in self.pars.keys():
                logs_file.write(f"{key}\t{UCDPhot.pars_names[key]} {self.pars[key]}\n")
                print(f"{key}\t{UCDPhot.pars_names[key]} {self.pars[key]}")

            print(f"Current date and time\t{now}")
            print("New set of parameters is loaded:")

    def differential_photometry(self):
        # list of paths to fits
        images_paths = Utils.get_image_list(self.pars['path2data'], self.pars['filter'])

        # read first frame for object name and other information
        print('First frame: ' + images_paths[0].split('/')[-1])

        hdu_list = fits.open(images_paths[0])
        header = hdu_list[0].header
        hdu_list.verify('fix')
        hdu_list.close()
        target = header['TARNAME']
        info = Utils.get_header_info(header)

    def do_astrometry(self, **kwargs):
        ast = AstrometryNet()
        ast.api_key = self.pars['astrometry.net_API']
        try_again = True
        submission_id = None
        wcs_header = None

        if kwargs.get('isFITS', False):
            while try_again:
                try:
                    if not submission_id:
                        t = time.time()
                        wcs_header = ast.solve_from_image(kwargs.get('path2image', None),
                                                          submission_id=submission_id, solve_timeout=250)
                        print('time - ', time.time()-t)
                    else:
                        wcs_header = ast.monitor_submission(submission_id,
                                                            solve_timeout=250)
                except:
                    print('Timeout')
                    # submission_id = e.args[1]
                else:
                    # got a result, so terminate
                    try_again = False

            if wcs_header:
                print(wcs_header)
            # Code to execute when solve succeeds
            else:
                print('astrometry solve fails')
                # Code to execute when solve fails
        else:
            sources = Utils.get_center(kwargs.get('path2image', None))
            wcs_header = ast.solve_from_source_list(sources['xcentroid'], sources['ycentroid'],
                                                    1024, 1024, solve_timeout=350)
            print(wcs_header)








