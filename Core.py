import shutil
from datetime import datetime

import astropy.io.fits as fits
from astropy.io import ascii
from astropy.stats import SigmaClip
from astropy.wcs import WCS
from photutils import Background2D, MedianBackground, CircularAperture, aperture_photometry, centroids
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
        "image_edge": "Images edge",

        "aperture": "Aperture radius",
        "search_box": "Centroid search box",
        "gain": "CCD gain",
        "ron": "CCD readout noise",

        "mmd": "Ensemble stars maximal magnitude difference",
        "isr": "Ensemble stars initial search radius in arcminutes",
        "msr": "Ensemble stars maximum search radius in arcminutes",
        "fop": "Field of wive in arcsec",

        "std_lim": "Ensemble stars standard deviation limit / (S/N)"
    }

    def __init__(self, path2data, **kwargs):

        self.pars = {"filter": kwargs.get("filter", "J"),
                     "apertures": kwargs.get("apertures", [15, 20, 25]),
                     'saturation': kwargs.get("saturation", 58000),
                     "path2data": path2data,
                     'mag_lim': kwargs.get("mag_lim", 20),
                     'max_mag': kwargs.get("max_mag", 9),
                     "gain": kwargs.get("gain", 2.18),
                     "ron": kwargs.get("ron", 0.02),
                     "mmd": kwargs.get("mmd", 2.0),
                     "isr": kwargs.get("isr", 5),
                     "msr": kwargs.get("msr", 30),
                     "std_lim": kwargs.get("std_lim", 3),
                     "fop": kwargs.get("fop", 4.6),
                     'cat': 'II/246',
                     "image_edge": kwargs.get("image_edge", 100),
                     "search_box": kwargs.get("search_box", 15),
                     'astrometry.net_API': 'hipfhzhlzygnlvix'}

        # make data paths
        if path2data[-1] == '/':
            self.pars['path2save'] = path2data[:-1] + '_report'
        else:
            self.pars['path2save'] = path2data + '_report'
        if not os.path.exists(self.pars['path2save']):
            os.makedirs(self.pars['path2save'])

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

    def aperture_photometry(self):
        # list of paths to fits
        images_paths = get_image_list(self.pars['path2data'], self.pars['filter'])

        # read first frame for object name and other information
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

        flux_out_list = []
        mag_out_list = []
        mag_err_list = []

        # open log file
        df = open(self.pars['path2save'] + '/Time.txt', 'w')
        df.write('DATE-OBS\tJD\tHJD\tBJD_TDB\t' +
                 'EXPTIME\tFILTER\tTEMP_SKY\t' +
                 'TEMP_AIR\tpressure\t' +
                 'X_Shift\tY_Shift\n')
        with tqdm(total=len(images_paths), desc='aperture_photometry') as bar:
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

                index = np.where((catalog['J'] > self.pars['max_mag']) & (catalog['J'] < self.pars['mag_lim']))[0]
                stars = stars_xy_coords[index, :]
                index = np.where(stars[:, 0] == 0)
                # stars[index] = np.nan

                background_object = Background2D(image_data, (100, 100),
                                                 filter_size=(10, 10),
                                                 sigma_clip=SigmaClip(sigma=3.),
                                                 bkg_estimator=MedianBackground())

                image_data = image_data - background_object.background
                signal_sky = background_object.background_rms_median

                stars_centroids = centroids.centroid_sources(image_data,
                                                             xpos=stars[:, 0],
                                                             ypos=stars[:, 1],
                                                             box_size=self.pars["search_box"])

                apertures_objects = [CircularAperture(np.vstack((stars_centroids[0], stars_centroids[1])).T, r=r)
                                     for r in self.pars["apertures"]]
                # aper_stat = ApertureStats(image_data, apertures_objects[-1])
                # print('\nfwhm = ', aper_stat.fwhm[0])
                phot_table = aperture_photometry(image_data, apertures_objects)
                flux = []
                mag = []
                mag_err = []
                for i in range(len(self.pars["apertures"])):
                    # flux
                    f = phot_table[f'aperture_sum_{i}']
                    f = np.array(f)
                    f[index] = np.nan
                    flux.append(f)

                    # mag
                    mag_data = 24 - 2.5 * np.log10(f) + 2.5 * np.log10(header["EXPTIME"])
                    mag_data[index] = np.nan
                    mag.append(mag_data)

                    # mag error
                    mag_err_data = (1.0857 * np.sqrt(f * self.pars["gain"] +
                                                     self.pars["apertures"][i] * (signal_sky * self.pars["gain"] +
                                                                                  self.pars["ron"] ** 2)) / (
                                            f * self.pars["gain"]))
                    mag_err_data[index] = np.nan
                    mag_err.append(mag_err_data)

                flux_out_list.append(flux)
                mag_out_list.append(mag)
                mag_err_list.append(mag_err)
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
                bar.update(1)
            bar.close()
            df.close()

        out_flux = np.array(flux_out_list)
        out_mag = np.array(mag_out_list)
        out_merr = np.array(mag_err_list)
        for i in range(len(self.pars["apertures"])):
            np.savetxt(self.pars['path2save'] + f'/Flux{i}.txt', out_flux[:, i], fmt='%.1f', delimiter='\t')
            np.savetxt(self.pars['path2save'] + f'/Mag{i}.txt', out_mag[:, i], fmt='%.4f', delimiter='\t')
            np.savetxt(self.pars['path2save'] + f'/Mag_err{i}.txt', out_merr[:, i], fmt='%.4f', delimiter='\t')

    def draw_curve(self):
        from matplotlib import pyplot as plt
        from astropy.time import Time as aTime
        time = ascii.read(self.pars['path2save'] + '/Time.txt', delimiter='\t')
        zero = int(time['BJD_TDB'][0])
        for i in range(len(self.pars["apertures"])):
            raw_magn = np.genfromtxt(self.pars['path2save'] + f'/Mag{i}.txt')
            raw_merr = np.genfromtxt(self.pars['path2save'] + f'/Mag_err{i}.txt')
            std = np.nanstd(raw_magn)

            fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=125)  # 3, 1, figsize=(6, 7), dpi=125
            r = self.pars['apertures'][i]
            fig.suptitle(f'Target, radius aperture = {r}, std = {std}', fontsize=8)
            ax.errorbar(time['BJD_TDB'] - zero, raw_magn[:, 0], raw_merr[:, 0], fmt='b.', markersize=3, zorder=2,
                        linewidth=0.5, label='Raw data')
            ax.legend(fontsize=6)
            locs = ax.get_xticks()
            t = aTime(locs, format='jd')
            x_ticks_labels = []
            for x in t:
                x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))

            ax.set_xticklabels(x_ticks_labels, fontsize=5)
            ax.set_xlabel('BJD_TDB - ' + str(zero), fontsize=6)
            ax.set_ylabel('mag')
            ax.invert_yaxis()
            ax.tick_params(axis='both', labelsize=6, direction='in')
            ax.grid()
            # plt.show()
            plt.savefig(self.pars['path2save'] + f'/plot{i}.pdf')

    def do_astrometry(self, **kwargs):
        from astroquery.astrometry_net import AstrometryNet
        print('working on ', os.path.split(kwargs.get('path2image', None))[1])
        hdu_list = fits.open(kwargs.get('path2image', None), 'update', memmap=False)
        header = hdu_list[0].header.copy()
        if header['CD1_1'] is None:
            hdu_list.close()
            return
        else:
            ast = AstrometryNet()
            ast.api_key = self.pars['astrometry.net_API']
            try_again = True
            submission_id = None
            wcs_header = None
            r = Angle('0.11d')
            dec = Angle(header['CURDEC'] + ' degrees')
            ra = Angle(header['CURRA'] + ' hours')
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
