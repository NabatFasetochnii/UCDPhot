import shutil
from datetime import datetime

import astropy.io.fits as fits
from astropy.io import ascii
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.wcs import WCS
from photutils import BasicPSFPhotometry, IntegratedGaussianPRF, MMMBackground, DAOGroup
from tqdm import tqdm

from Utils import *


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

    def __init__(self, path2data=r'G:\Kislovodsk\T1b_J-20200912\done', **kwargs):

        self.pars = {"filter": kwargs.get("filter", "J"),
                     "aperture": kwargs.get("aperture", 10),
                     'saturation': kwargs.get("saturation", 58000),
                     "path2data": path2data,
                     'mag_lim': kwargs.get("mag_lim", 20),  # ?
                     'max_mag': kwargs.get("max_mag", 9),  # ?
                     "gain": kwargs.get("gain", 2.18),
                     "ron": kwargs.get("ron", 0.02),
                     "mmd": kwargs.get("mmd", 2.0),
                     "isr": kwargs.get("isr", 5),
                     "msr": kwargs.get("msr", 30),
                     "std_lim": kwargs.get("std_lim", 3),
                     "fop": kwargs.get("fop", 4.6),  # ?
                     'cat': 'II/246',
                     "image_edge": kwargs.get("image_edge", 100),
                     "search_box": kwargs.get("search_box", 15),
                     'astrometry.net_API': 'hipfhzhlzygnlvix'}

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

    def psf_photometry(self):
        # list of paths to fits
        images_paths = get_image_list(self.pars['path2data'], self.pars['filter'])

        # read first frame for object name and other information
        print('First frame: ' + images_paths[0].split('/')[-1])
        # print(self.pars['path2data'] + '/' + images_paths[0])
        # hdul = fits.open(self.pars['path2data'] + '/' + images_paths[0])

        # header = hdul[0].header
        # hdul.verify('fix')
        # # target = header['TARNAME']
        # # print(header)
        # # info = get_header_info(header)
        # # fil = header['FILTER']
        # # image_radius = abs(header["CD1_1"]) * ((header["NAXIS1"] ** 2 + header["NAXIS2"] ** 2) ** 0.5) / 2
        #
        # hdul.close()
        # catalog = get_2mass(ra=info[0], dec=info[1], r=image_radius,
        #                     j_lim=self.pars['mag_lim'], cat=self.pars['cat'])
        # ra_fix = Angle('0:0:1.3 hours')
        # dec_fix = Angle('0:0:10 degree')
        # catalog['Ra'][0] = catalog['Ra'][0] + ra_fix.degree
        # catalog['Dec'][0] = catalog['Dec'][0] - dec_fix.degree
        #
        # ra_fix_4 = Angle('0:0:0.4 hours')
        # dec_fix_4 = Angle('0:0:5 degree')
        #
        # catalog['Ra'][4] = catalog['Ra'][4] - ra_fix_4.degree
        # catalog['Dec'][4] = catalog['Dec'][4] - dec_fix_4.degree
        #
        # catalog['Ra'][5] = catalog['Ra'][5] - ra_fix_4.degree
        # catalog['Dec'][5] = catalog['Dec'][5] - dec_fix_4.degree

        catalog = np.genfromtxt(self.pars['path2save'] + '/Cat.txt', skip_header=1,
                                names=['ID', 'Ra', 'Dec', 'Dist', 'J'])
        sigma_psf = 2.0

        # open flux files
        flux_file = open(self.pars['path2save'] + '/Flux.txt', 'w')
        mag_file = open(self.pars['path2save'] + '/Mag.txt', 'w')
        m_err_file = open(self.pars['path2save'] + '/Mag_err.txt', 'w')

        # open log file
        df = open(self.pars['path2save'] + '/Time.txt', 'w')
        df.write('DATE-OBS\tJD\tHJD\tBJD_TDB\t' +
                 'EXPTIME\tFILTER\tTEMP_SKY\t' +
                 'TEMP_AIR\tpressure\t' +
                 'X_Shift\tY_Shift\n')
        with tqdm(total=len(images_paths), desc='PSF_photometry') as bar:
            print()
            for file_path in images_paths:
                hdu_list = fits.open(self.pars['path2data'] + '/' + file_path)
                header = hdu_list[0].header.copy()
                image_data = hdu_list[0].data
                hdu_list.verify('fix')
                info = get_header_info(header)
                hdu_list.close()

                df.write(info[4].datetime.isoformat(timespec='milliseconds')  # DATE-OBS
                         + '\t' + '{:.7f}'.format(info[4].jd) + '\t')  # JD
                df.write('{:.7f}'.format(info[5]) + '\t' + '{:.7f}'.format(info[6]) + '\t')  # HJD BJD
                df.write('{:.1f}'.format(header['EXPTIME']) + '\t')  # EXPTIME
                df.write(header['FILTER'] + '\t')  # FILTER
                df.write('{:.1f}'.format(header['TEMP_SKY']) + '\t')  # TEMP_SKY
                df.write('{:.1f}'.format(header['TEMP_AIR']) + '\t')  # TEMP_AIR
                df.write('{:.1f}'.format(header['PRESS']) + '\t')  # pressure
                df.write(str(info[2]) + '\t' + str(info[3]) + '\n')  # X_Shift Y_Shift

                # make WCS object
                wcs_object = WCS(header)

                # founding xy coords of the stars on frame
                stars_xy_coords = wcs_object.all_world2pix(catalog['Ra'], catalog['Dec'], 0)

                bad_stars_x_mask = np.where((stars_xy_coords[0] < self.pars["image_edge"] / 10) |
                                            (stars_xy_coords[0] > (image_data.shape[1] -
                                                                   self.pars["image_edge"] / 10)))[0]
                bad_stars_y_mask = np.where((stars_xy_coords[1] < self.pars["image_edge"] / 10) |
                                            (stars_xy_coords[1] > (image_data.shape[0] -
                                                                   self.pars["image_edge"] / 10)))[0]
                bad_stars_mask = np.concatenate((bad_stars_x_mask, bad_stars_y_mask), axis=0)

                if len(bad_stars_mask) > 0:
                    stars_xy_coords[0][bad_stars_mask] = 0
                    stars_xy_coords[1][bad_stars_mask] = 0

                stars_xy_coords = np.vstack((stars_xy_coords[0], stars_xy_coords[1])).T

                # make subarray of PSF stars
                index = np.where((catalog['J'] > self.pars['max_mag']) & (catalog['J'] < self.pars['mag_lim']))[0]
                psf_stars = stars_xy_coords[index, :]
                index = np.where(psf_stars[:, 0] == 0)
                psf_stars = np.delete(psf_stars, index, axis=0)

                dao_group = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
                mmm_bkg = MMMBackground()
                fitter = LevMarLSQFitter()
                psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
                photometry = BasicPSFPhotometry(
                    group_maker=dao_group,
                    bkg_estimator=mmm_bkg,
                    psf_model=psf_model,
                    fitter=fitter,
                    fitshape=13)
                pos = Table(names=['x_0', 'y_0'], data=[psf_stars[:, 0],
                                                        psf_stars[:, 1]])
                result_tub = None
                result_tub = photometry.do_photometry(image=image_data, init_guesses=pos)
                f = result_tub['flux_fit']
                bkg_value = mmm_bkg.calc_background(image_data)
                for item in f:
                    flux_file.write('{:.1f}'.format(item) + '\t')
                flux_file.write('\n')
                raper = photometry.aperture_radius
                magn_image_data = - 2.5 * np.log10(f) + 2.5 * np.log10(header["EXPTIME"])
                merr_image_data = 1.0857 * np.sqrt(f * self.pars["gain"] +
                                                   raper * (bkg_value * self.pars["gain"] +
                                                            self.pars["ron"] ** 2)) / (f * self.pars["gain"])
                for item in magn_image_data:
                    mag_file.write('{:.4f}'.format(item) + '\t')
                mag_file.write('\n')
                for item in merr_image_data:
                    m_err_file.write('{:.4f}'.format(item) + '\t')
                m_err_file.write('\n')
                # print(result_tub)
                # arr = [list(result_tub['x_fit']), list(result_tub['y_fit'])]
                # ra_dec_center = wcs_object.pixel_to_world_values(*arr)
                # for i in range(len(ra_dec_center[0])):
                #     out_cat.append([list(catalog['ID'])[i], ra_dec_center[0][i], ra_dec_center[1][i],
                #                     list(catalog['Dist'])[i], list(catalog['J'])[i]])
                # out_cat = Table(names=['x_pix', 'y_pix', 'Ra', 'Dec', 'flux', 'flux_unc'],
                #                 data=[list(result_tub['x_fit']), list(result_tub['y_fit']),
                #                       *ra_dec_center, list(result_tub['flux_fit']),
                #                       list(result_tub['flux_unc'])])
                # flux.append(list(result_tub['flux_fit']))

                # # draw a picture
                # Image = np.log10(image_data)
                # X = Image.shape[1]
                # Y = Image.shape[0]
                # _mean, _median, _std = sigma_clipped_stats(Image[Y - 50:Y + 50, X - 50:X + 50])
                # _max = _median + 10. * _std
                # _min = _median - 1. * _std
                # fig = plt.figure(figsize=(7, 7), frameon=False)
                # ax = plt.subplot(projection=wcs_object, position=[0.1, 0.1, 0.8, 0.8])
                # plt.imshow(Image, vmin=_min, vmax=_max, cmap='gray_r')
                #
                # coor = np.vstack((result_tub['x_fit'], result_tub['y_fit'])).T
                # # aper_old = CircularAperture(psf_stars, r=10)
                # # aper_old.plot(color='red', lw=1.5, alpha=0.5)
                # aper = CircularAperture(coor, r=10)
                # aper.plot(color='blue', lw=1.5, alpha=0.5)
                # for i in range(0, len(psf_stars)):
                #     plt.text(psf_stars[i, 0], psf_stars[i, 1], s=catalog['ID'][i], color='blue', alpha=0.8)
                # if header['CD1_1'] > 0:
                #     plt.gca().invert_xaxis()
                # if header['CD2_2'] > 0:
                #     plt.gca().invert_yaxis()
                # plt.title(target + ', filter ' + header['FILTER'] + '\n' + header['DATE-OBS'])
                # ax.coords[1].set_ticklabel(rotation=90)
                # ax.coords[0].set_major_formatter('hh:mm:ss')
                # ax.coords[1].set_major_formatter('dd:mm:ss')
                # ax.coords[0].set_axislabel('RA')
                # ax.coords[1].set_axislabel('Dec')
                # ax.coords.grid(color='blue', ls='--', alpha=0.7)
                # # plt.show()
                # fig.savefig('field.pdf')
                print()
                bar.update(1)
            bar.close()
            df.close()
            flux_file.close()
            mag_file.close()
            m_err_file.close()

        # out_flux = np.array(flux)
        # out_mag = np.array(mag)
        # out_merr = np.array(m_err)
        # np.savetxt(self.pars['path2save'] + '/Flux.txt', out_flux, fmt='%.1f',  delimiter='\t')
        # np.savetxt(self.pars['path2save'] + '/Mag.txt', out_mag, fmt='%.4f',  delimiter='\t')
        # np.savetxt(self.pars['path2save'] + '/Mag_err.txt', out_merr, fmt='%.4f',  delimiter='\t')

    def diff_photometry(self):
        catalog = np.genfromtxt(self.pars['path2save'] + '/Cat.txt', skip_header=1,
                                names=['ID', 'Ra', 'Dec', 'Dist', 'J'])
        flux = np.genfromtxt(self.pars['path2save'] + '/Flux.txt')
        mag = np.genfromtxt(self.pars['path2save'] + '/Mag.txt')
        mag_err = np.genfromtxt(self.pars['path2save'] + '/Mag_err.txt')
        time = ascii.read(self.pars['path2save'] + '/Time.txt', delimiter='\t')
        zero = int(time['BJD_TDB'][0])

        # list of paths to fits
        images_paths = get_image_list(self.pars['path2data'], self.pars['filter'])

        images_number = len(images_paths)
        stars_number = len(catalog)
        counter_success = 0

        clr_magn = np.zeros((images_number, stars_number))
        clr_merr = np.zeros((images_number, stars_number))

        # fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=125)  # 3, 1, figsize=(6, 7), dpi=125
        # fig.suptitle('Target', fontsize=8)
        # ax.plot(time['BJD_TDB']-zero, flux[:, 0], 'b.', markersize=3, zorder=4, linewidth=0.5, label='Data')
        # ax.legend(fontsize=6)
        # locs = ax.get_xticks()
        # t = aTime(locs, format='jd')
        # x_ticks_labels = []
        # for x in t:
        #     x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
        #
        # ax.set_xticklabels(x_ticks_labels, fontsize=5)
        # ax.set_xlabel('BJD_TDB - ' + str(zero), fontsize=6)
        # ax.tick_params(axis='both', labelsize=6, direction='in')
        # ax.grid()
        # plt.savefig(self.pars['path2save'] + '/plot.pdf')

    def do_astrometry(self, **kwargs):
        from astroquery.astrometry_net import AstrometryNet
        print('working on ', os.path.split(kwargs.get('path2image', None))[1])
        hdu_list = fits.open(kwargs.get('path2image', None), 'update', memmap=False)
        header = hdu_list[0].header.copy()
        try:
            a = header['CD1_1']
            hdu_list.close()
            return
        except:
            ast = AstrometryNet()
            ast.api_key = self.pars['astrometry.net_API']
            try_again = True
            submission_id = None
            wcs_header = None
            r = Angle('0.11d')
            dec = Angle(header['CURDEC'] + ' degrees')
            ra = Angle(header['CURRA'] + ' hours')
            # print(dec.degree)
            # print(ra.degree)
            upper = 4.7
            lower = 4.5
            if kwargs.get('isFITS', False):
                while try_again:
                    try:
                        if not submission_id:
                            wcs_header = ast.solve_from_image(kwargs.get('path2image', None),
                                                              submission_id=submission_id, solve_timeout=250,
                                                              center_ra=ra.degree, center_dec=dec.degree,
                                                              radius=kwargs.get('radius', r.degree),
                                                              scale_units='arcminwidth', scale_type='ul',
                                                              scale_lower=kwargs.get('scale_lower', lower),
                                                              scale_upper=kwargs.get('scale_upper', upper),
                                                              publicly_visible='n', parity=2)
                        else:
                            wcs_header = ast.monitor_submission(submission_id, solve_timeout=250)
                    except TimeoutError as e:
                        print('Timeout')
                        submission_id = e.args[1]
                    else:
                        try_again = False
            else:
                while try_again:
                    try:
                        if not submission_id:

                            data = hdu_list[0].data.copy()
                            # print(data[0][0])
                            hdu_list.verify('fix')
                            sources = get_center(data)
                            wcs_header = ast.solve_from_source_list(sources['xcentroid'], sources['ycentroid'],
                                                                    header['NAXIS1'], header['NAXIS2'],
                                                                    center_ra=ra.degree, center_dec=dec.degree,
                                                                    radius=r.degree,
                                                                    scale_units='arcminwidth', scale_type='ul',
                                                                    scale_lower=kwargs.get('scale_lower', lower),
                                                                    scale_upper=kwargs.get('scale_upper', upper),
                                                                    publicly_visible='n', parity=2)
                        else:
                            wcs_header = ast.monitor_submission(submission_id, solve_timeout=250)
                    except TimeoutError as e:
                        print('Timeout')
                        submission_id = e.args[1]
                    else:
                        try_again = False
            if wcs_header:
                hdu_list[0].header = header + wcs_header
                hdu_list.close()
                shutil.move(kwargs.get('path2image', None), os.path.split(kwargs.get('path2image', None))[0] + '/done/')
                print('done')
            else:
                print('astrometry solve fails')
                hdu_list.close()
