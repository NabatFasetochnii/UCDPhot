import astropy.io.fits as fits
from astropy.io import ascii
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.time import Time as aTime
from astropy.wcs import WCS
from photutils import Background2D, MedianBackground, CircularAperture, aperture_photometry, centroids
from tqdm import tqdm
from matplotlib import pyplot as plt
from astropy.stats import sigma_clip
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
                     "search_box": kwargs.get("search_box", 15)}

        # make data paths
        if path2data[-1] == '/':
            self.pars['path2save'] = path2data[:-1] + '_report'
        else:
            self.pars['path2save'] = path2data + '_report'
        # make folder
        if not os.path.exists(self.pars['path2save']):
            os.makedirs(self.pars['path2save'])

    def aperture_photometry(self, **kwargs):

        # list of paths to fits
        images_paths = get_image_list(self.pars['path2data'], self.pars['filter'])

        if kwargs['cat_path']:
            catalog = np.genfromtxt(kwargs['cat_path'], skip_header=1,
                                    names=['ID', 'Ra', 'Dec', 'Dist', 'J'])
        else:
            # read first frame for object name and other information
            print(self.pars['path2data'] + '/' + images_paths[0])
            hdul = fits.open(self.pars['path2data'] + '/' + images_paths[0])
            header = hdul[0].header
            hdul.verify('fix')
            # target = header['TARNAME']
            info = get_header_info(header)
            # fil = header['FILTER']
            image_radius = abs(header["CD1_1"]) * ((header["NAXIS1"] ** 2 + header["NAXIS2"] ** 2) ** 0.5) / 2

            hdul.close()
            catalog = get_2mass(ra=info[0], dec=info[1], r=image_radius,
                                j_lim=self.pars['mag_lim'], cat=self.pars['cat'])
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

        flux_out_list = []
        mag_out_list = []
        mag_err_list = []

        # open log file
        df = open(self.pars['path2save'] + '/Time.txt', 'w')
        df.write('DATE-OBS\tJD\tHJD\tBJD_TDB\t' +
                 'EXPTIME\tFILTER\tAIRMASS\tTEMP_SKY\t' +
                 'TEMP_AIR\tPRESSURE\t' +
                 'X_TARGET\tY_TARGET\tSKY\n')
        # from scipy import ndimage as nd
        with tqdm(total=len(images_paths), desc='aperture_photometry') as bar:
            print()
            for file_path in images_paths:
                hdu_list = fits.open(self.pars['path2data'] + '/' + file_path)
                header = hdu_list[0].header.copy()
                image_data = hdu_list[0].data
                hdu_list.verify('fix')
                info = get_header_info(header)
                hdu_list.close()

                # image_data = nd.median_filter(image_data, size=4)

                df.write(info[4].datetime.isoformat(timespec='milliseconds')  # DATE-OBS
                         + '\t' + '{:.7f}'.format(info[4].jd) + '\t')  # JD
                df.write('{:.7f}'.format(info[5]) + '\t' + '{:.7f}'.format(info[6]) + '\t')  # HJD BJD
                df.write('{:.1f}'.format(header['EXPTIME']) + '\t')  # EXPTIME
                df.write(header['FILTER'] + '\t')  # FILTER
                df.write('{:.7f}'.format(info[7]) + '\t')  # airmass
                df.write('{:.1f}'.format(header['TEMP_SKY']) + '\t')  # TEMP_SKY
                df.write('{:.1f}'.format(header['TEMP_AIR']) + '\t')  # TEMP_AIR
                df.write('{:.1f}'.format(header['PRESS']) + '\t')  # pressure
                # df.write(str(info[2]) + '\t' + str(info[3]) + '\t')  # X_Shift Y_Shift

                # make WCS object
                wcs_object = WCS(header)

                # founding xy coords of the stars on frame
                stars_xy_coords = wcs_object.all_world2pix(catalog['Ra'], catalog['Dec'], 0)

                df.write(str(stars_xy_coords[0][0]) + '\t' + str(stars_xy_coords[1][0]) + '\t')  # X_target Y_target

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
                df.write('{:.1f}'.format(background_object.background_median) + '\n')  # sky

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

    def create_humidity(self, path):
        import json
        with open(path, 'r') as read_f:
            water_data = json.load(read_f)['results'][0]['series'][0]['values']  # humidity
            water_data = np.array(water_data).T
        time = ascii.read(self.pars['path2save'] + '/Time.txt', delimiter='\t', fast_reader=False, guess=False)
        time_s = aTime(time['DATE-OBS'])
        water_time = aTime(water_data[0])
        water_data = np.float32(water_data[1])
        humidity = np.interp(time_s.tdb.jd, water_time.tdb.jd, water_data)
        np.savetxt(self.pars['path2save'] + '/humidity.txt', humidity, fmt='%.1f', delimiter='\t')

    def draw_curve(self, planet='c'):
        import numpy.ma as ma
        # import copy
        time = ascii.read(self.pars['path2save'] + '/Time.txt', delimiter='\t', fast_reader=False, guess=False)
        Info = ascii.read(f'{self.pars["path2save"]}/Info.dat', delimiter='\t')
        zero = int(time['JD'][0])
        date = time['DATE-OBS'][0].split('T')[0]
        humidity = np.genfromtxt(self.pars['path2save'] + '/humidity.txt')
        gs_kw = dict(height_ratios=[10, 10, 10, 7, 7, 7, 6, 6, 5, 10])
        delta = get_delta_transit(time['JD'] - zero, planet)
        for i in range(len(self.pars["apertures"])):
            raw_magn = np.genfromtxt(self.pars['path2save'] + f'/Mag{i}.txt')
            m = sigma_clip(raw_magn[:, 0])
            raw_magn = med_filt(raw_magn)
            raw_merr = np.genfromtxt(self.pars['path2save'] + f'/Mag_err{i}.txt')
            std = np.round(np.nanstd(m), 4)
            fig, ax = plt.subplots(10, 1, figsize=(6, 14), sharex=True, sharey=False, gridspec_kw=gs_kw, dpi=125,
                                   constrained_layout=True)
            r = self.pars['apertures'][i]

            # fig.suptitle(f'Photometry report. Radius aperture = {r} pix\nObservation date: {date}', fontsize=8)

            ax[0].errorbar(time['JD'][:] - zero, m, raw_merr[:, 0], fmt='b.',
                           markersize=3, zorder=2, linewidth=0.5, label=f'Target, std = {std}, aperture = {r} pix')
            ax[0].set_ylabel('mag', fontsize=6)
            ax[0].invert_yaxis()
            ax[0].legend(fontsize=6, loc=3)
            ax[0].tick_params(axis='both', labelsize=6, direction='in')
            # ax[0].set_xticklabels(x_ticks_labels, fontsize=5)
            ax[0].grid()

            sub_ax = ax[0].twinx()
            sub_ax.plot(time['JD'] - zero, Info['TEMP_SKY'], 'r--', label='Температура неба',
                        markersize=3, zorder=1, linewidth=0.5)
            sub_ax.set_ylabel(r'$T_{sky}$,' + ' \N{DEGREE SIGN}C', fontsize=6)
            sub_ax.legend(fontsize=6, loc=0)
            sub_ax.tick_params(axis='both', labelsize=6, direction='in')

            for item in range(1, 6):
                star_std = np.round(np.nanstd(raw_magn[:, item]), 4)
                if star_std < 0.1:
                    ax[1].plot(time['JD'] - zero, raw_magn[:, item], '.',
                               label=f'Star #{item}, std = {star_std}', markersize=3)
                    break
            ax[1].set_ylabel('mag', fontsize=6)
            ax[1].invert_yaxis()
            magic = np.genfromtxt(f"{self.pars['path2save']}/sysrem/sysrem_out{i}_0.sys_rem.txt", skip_header=True)
            magic = ma.masked_equal(magic, 1)
            magic = ma.mask_rows(magic)
            magic_std = np.round(np.nanstd(magic[:, 1]), 4)
            SN = np.round(np.mean(1.0857 / magic[:, 2]), 4)

            ax[2].errorbar(magic[:, 0] - zero, magic[:, 1], magic[:, 2], fmt='b.', markersize=3, zorder=2,
                           linewidth=0.5, label=f'Correcting linear systematic effects, std = {magic_std}, S/N = {SN}')
            ax[2].set_ylabel('mag', fontsize=6)
            ax[2].invert_yaxis()

            ax_fake_transit = ax[2].twinx()
            ax_fake_transit.errorbar(time['JD'] - zero, m +
                                     delta*np.mean(m),
                                     raw_merr[:, 0], fmt='r--', label='Модель транзита')
            ax_fake_transit.invert_yaxis()
            ax_fake_transit.set_ylabel('mag', fontsize=6)
            ax_fake_transit.legend(fontsize=6, loc=4)
            ax_fake_transit.tick_params(axis='both', labelsize=6, direction='in')

            ax[3].plot(time['JD'] - zero, Info['TEMP_AIR'], label='Температура воздуха')
            ax[3].set_ylabel(r'$T_{air}$,' + ' \N{DEGREE SIGN}C', fontsize=6)

            press_ax = ax[3].twinx()
            press_ax.plot(time['JD'] - zero, Info['PRESS'], 'r--', label='Давление воздуха')
            press_ax.set_ylabel('P, мм. рт. ст.', fontsize=6)
            press_ax.legend(fontsize=6, loc=4)
            press_ax.tick_params(axis='both', labelsize=6, direction='in')

            ax[4].plot(time['JD'] - zero, Info['RH'], label='Относительная влажность')
            ax[4].set_ylabel('Humidity, %', fontsize=6)
            ax[4].legend(fontsize=6, loc=3)
            ax[4].tick_params(axis='both', labelsize=6, direction='in')
            ax[4].grid()

            ax_pwv = ax[4].twinx()
            ax_pwv.plot(time['JD'] - zero, humidity, 'r--', label='Осаждённая влага')
            ax_pwv.set_ylabel('Precipitable Water, мм', fontsize=6)
            ax_pwv.legend(fontsize=6, loc=0)
            ax_pwv.tick_params(axis='both', labelsize=6, direction='in')

            ax[5].plot(time['JD'] - zero, time['AIRMASS'], label='Воздушная масса')
            ax[5].set_ylabel('airmass', fontsize=6)

            ax[6].plot(time['JD'] - zero, Info['WIND_DIR'], label='Направление ветра')
            ax[6].set_ylabel('WIND_DIR, \N{DEGREE SIGN}', fontsize=6)

            ax[7].plot(time['JD'] - zero, Info['WIND'], label='Скорость ветра')
            ax[7].set_ylabel('WIND, м/с', fontsize=6)

            ax[8].plot(time['JD'] - zero, time['X_TARGET'] - np.mean(time['X_TARGET']),
                       'r.', label='X Shift', markersize=3)
            ax[8].plot(time['JD'] - zero, time['Y_TARGET'] - np.mean(time['Y_TARGET']),
                       'b.', label='Y Shift', markersize=3)
            ax[8].legend(loc=0, fontsize=6)
            ax[8].set_ylabel('Shift, pix', fontsize=6)

            ax[9].plot(time['JD'] - zero, time['SKY'], '.', label='Сигнал фона неба', markersize=3)
            ax[9].set_ylabel('Sky, ADU', fontsize=6)

            locs = ax[0].get_xticks()
            t = aTime(locs, format='jd')
            x_ticks_labels = []
            for x in t:
                x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
            for a in ax[1:]:
                if a == ax[4]:
                    continue
                a.legend(fontsize=6, loc=0)
                a.tick_params(axis='both', labelsize=6, direction='in')
                a.set_xticklabels(x_ticks_labels, fontsize=5)
                a.grid()
            ax[-1].set_xlabel('JD - ' + str(zero), fontsize=6)
            # plt.show()
            plt.savefig(self.pars['path2save'] + f'/plot_{date}_{i}.pdf')

    def draw_condition(self):
        from matplotlib import pyplot as plt
        from astropy.time import Time as aTime
        # read_condition(self.pars['path2data'])
        Info = ascii.read(f'{self.pars["path2save"]}/Info.dat', delimiter='\t')
        time = ascii.read(f'{self.pars["path2save"]}/Time.txt', delimiter='\t', fast_reader=False, guess=False)
        humidity = np.genfromtxt(self.pars['path2save'] + '/humidity.txt')
        zero = int(time['JD'][0])
        date = time['DATE-OBS'][0].split('T')[0]
        gs_kw = dict(height_ratios=[10, 10, 7, 7, 7, 7, 6, 6, 5, 10])
        fig, ax = plt.subplots(10, 1, figsize=(6, 12), sharex=True, sharey=False, gridspec_kw=gs_kw, dpi=125,
                               constrained_layout=True)

        fig.suptitle(f'Condition report\nObservation date: {date}', fontsize=8)

        ax[0].plot(time['JD'] - zero, Info['TEMP_SKY'], label='Температура неба, \N{DEGREE SIGN}C')
        ax[0].set_ylabel(r'$T_{sky}$', fontsize=6)

        ax[1].plot(time['JD'] - zero, Info['TEMP_AIR'], label='Температура воздуха, \N{DEGREE SIGN}C')
        ax[1].set_ylabel(r'$T_{air}$', fontsize=6)

        ax[2].plot(time['JD'] - zero, Info['RH'], label='Относительная влажность, %')
        ax[2].set_ylabel('Humidity', fontsize=6)

        ax[3].plot(time['JD'] - zero, humidity, label='Осаждённая влага, мм')
        ax[3].set_ylabel('Precipitable Water', fontsize=6)

        ax[4].plot(time['JD'] - zero, time['AIRMASS'], label='Воздушная масса')
        ax[4].set_ylabel('airmass', fontsize=6)

        ax[5].plot(time['JD'] - zero, Info['PRESS'], label='Давление, мм. рт. ст.')
        ax[5].set_ylabel('P', fontsize=6)

        ax[6].plot(time['JD'] - zero, Info['WIND_DIR'], label='Направление ветра, \N{DEGREE SIGN}')
        ax[6].set_ylabel('WIND_DIR', fontsize=6)

        ax[7].plot(time['JD'] - zero, Info['WIND'], label='Скорость ветра, м/с')
        ax[7].set_ylabel('WIND', fontsize=6)

        ax[8].plot(time['JD'] - zero, time['X_TARGET'] - np.mean(time['X_TARGET']), 'r.', label='X Shift', markersize=3)
        ax[8].plot(time['JD'] - zero, time['Y_TARGET'] - np.mean(time['Y_TARGET']), 'b.', label='Y Shift', markersize=3)
        ax[8].legend(loc=0, fontsize=6)
        ax[8].set_ylabel('shift (pix)', fontsize=6)

        ax[9].plot(time['JD'] - zero, time['SKY'], '.', label='Sky signal (ADU)', markersize=3)
        ax[9].set_ylabel('Sky, ADU', fontsize=6)

        locs = ax[0].get_xticks()
        t = aTime(locs, format='jd')
        x_ticks_labels = []
        for x in t:
            x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
        for a in ax:
            a.legend(fontsize=6)
            a.tick_params(axis='both', labelsize=6, direction='in')
            a.set_xticklabels(x_ticks_labels, fontsize=5)
            a.grid()
        ax[-1].set_xlabel('JD - ' + str(zero), fontsize=6)
        # plt.show()
        plt.savefig(f'{self.pars["path2save"]}\condition_{date}.pdf')

    def magic(self, planet='c'):
        from PySysRem import source_lc, sysrem
        import copy
        info = ascii.read(self.pars['path2save'] + '/Time.txt', delimiter='\t', fast_reader=False, guess=False)
        time = info['JD']
        bad = np.array([1, 4, 5, 9])

        delta_transit = get_delta_transit(time, planet)

        # plt.plot(time, fake_curve)
        # plt.grid()
        # plt.gca().invert_yaxis()
        # plt.show()

        if not os.path.exists(self.pars['path2save'] + '/sysrem'):
            os.makedirs(self.pars['path2save'] + '/sysrem')
        for i in range(len(self.pars["apertures"])):
            mag = copy.deepcopy(np.genfromtxt(self.pars['path2save'] + f'/Mag{i}.txt'))
            mag[:, 0] = sigma_clip(mag[:, 0])
            mean = np.mean(mag[:, 0])
            mag[:, 0] = mag[:, 0] + delta_transit * mean
            m_err = np.genfromtxt(self.pars['path2save'] + f'/Mag_err{i}.txt')
            # plt.errorbar(time, mag[:, 0], m_err[:, 0])
            # plt.grid()
            # plt.gca().invert_yaxis()
            # plt.savefig('last_raw.pdf')
            flag = np.zeros(len(time))
            flag[bad] = 1
            source_list = []
            for k in range(len(mag[0, :])):
                source = source_lc.Source(mag[:, k], m_err[:, k], time, flag,
                                          f"{self.pars['path2save']}/sysrem/sysrem_out{i}_{k}.txt")
                source_list.append(source)
            sysrem.sys_rem(source_list, len_lin_ef=1)

    # def draw_magic_curve(self):
    #     from matplotlib import pyplot as plt
    #     info = np.genfromtxt(f"{self.pars['path2save']}/sysrem/sysrem_out0_0.sys_rem.txt", skip_header=True)
    #     fig, ax = plt.subplots(1, 1, figsize=(7, 3), dpi=125)  # 3, 1, figsize=(6, 7), dpi=125
    #     zero = int(info[0][0])
    #     date = aTime(info[0][0], format='jd').strftime("%y-%m-%d")
    #     std = np.round(np.nanstd(info[:, 1]), 4)
    #     ax.errorbar(info[:, 0] - zero, info[:, 1], info[:, 2], fmt='b.', markersize=3, zorder=2,
    #                 linewidth=0.5, label=f'Target, magic, std = {std}')
    #     ax.set_ylabel('mag', fontsize=6)
    #     ax.invert_yaxis()
    #     locs = ax.get_xticks()
    #     t = aTime(locs, format='jd')
    #     fig.suptitle(f'Correcting systematic effects\nObservation date: {date}', fontsize=8)
    #     x_ticks_labels = []
    #     for x in t:
    #         x_ticks_labels.append(str('{:.2f}'.format(x.tdb.jd)))
    #     ax.legend(fontsize=6)
    #     ax.tick_params(axis='both', labelsize=6, direction='in')
    #     ax.set_xticklabels(x_ticks_labels, fontsize=5)
    #     ax.grid()
    #     ax.set_xlabel('JD - ' + str(zero), fontsize=6)
    #     plt.savefig(self.pars['path2save'] + f'/magic_plot_{date}.pdf')

    def read_condition(self):
        # from astropy.io import ascii, fits
        dir_content = os.listdir(self.pars['path2data'])
        Info = Table()

        MJD_AVG = []
        TEMP_SKY = []
        TEMP_AIR = []
        PRESS = []
        WIND_DIR = []
        WIND = []
        RH = []
        for f in dir_content:
            if f.count('.fit') or f.count('.fits') or f.count('.fts'):
                hdulist = fits.open(self.pars['path2data'] + '/' + f)
                # Data = hdulist[0].data
                Header = hdulist[0].header
                hdulist.close()
                exp = TimeDelta(Header['EXPTIME'], format='sec')
                t = Time(Header['DATE-OBS'], format='fits')
                t = t + exp / 2.

                MJD_AVG.append(t.mjd)
                TEMP_SKY.append(Header['TEMP_SKY'])
                TEMP_AIR.append(Header['TEMP_AIR'])
                PRESS.append(Header['PRESS'])
                WIND_DIR.append(Header['WIND_DIR'])
                WIND.append(Header['WIND'])
                RH.append(Header['RH'])

        Info['MJD-AVG'] = MJD_AVG
        Info['TEMP_SKY'] = TEMP_SKY
        Info['TEMP_AIR'] = TEMP_AIR
        Info['PRESS'] = PRESS
        Info['WIND_DIR'] = WIND_DIR
        Info['WIND'] = WIND
        Info['RH'] = RH
        Info['MJD-AVG'] = Info['MJD-AVG'] + 2400000.5 - 2459105
        # print(Info.info)
        ascii.write(Info, f'{self.pars["path2save"]}\Info.dat', delimiter='\t', overwrite=True)

    def plot_field(self):
        from matplotlib import pyplot as plt
        images_paths = get_image_list(self.pars['path2data'], self.pars['filter'])
        hdul = fits.open(self.pars['path2data'] + '/' + images_paths[0])
        header = hdul[0].header
        hdul.verify('fix')
        image_data = hdul[0].data

        # make WCS object
        wcs_object = WCS(header)

        hdul.close()

        catalog = np.genfromtxt(self.pars['path2save'] + '/Cat.txt', skip_header=1,
                                names=['ID', 'Ra', 'Dec', 'Dist', 'J'])
        stars_xy_coords = wcs_object.all_world2pix(catalog['Ra'], catalog['Dec'], 0)
        # # draw a picture
        Image = np.log10(image_data)
        X = Image.shape[1]
        Y = Image.shape[0]
        _mean, _median, _std = sigma_clipped_stats(Image[Y - 50:Y + 50, X - 50:X + 50])
        _max = _median + 10. * _std
        _min = _median - 1. * _std
        fig = plt.figure(figsize=(7, 7), frameon=False)
        ax = plt.subplot(projection=wcs_object, position=[0.1, 0.1, 0.8, 0.8])
        plt.imshow(Image, vmin=_min, vmax=_max, cmap='gray_r')

        stars_xy_coords = np.vstack((stars_xy_coords[0], stars_xy_coords[1])).T

        aper = CircularAperture(stars_xy_coords, r=10)
        aper.plot(color='blue', lw=1.5, alpha=0.5)
        for i in range(0, len(catalog)):
            plt.text(x=stars_xy_coords[i, 0], y=stars_xy_coords[i, 1], s=catalog['ID'][i], color='blue', alpha=0.8)
        if header['CD1_1'] > 0:
            plt.gca().invert_xaxis()
        if header['CD2_2'] > 0:
            plt.gca().invert_yaxis()
        plt.title(header['TARNAME'] + ', filter ' + header['FILTER'] + '\n' + header['DATE-OBS'])
        ax.coords[1].set_ticklabel(rotation=90)
        ax.coords[0].set_major_formatter('hh:mm:ss')
        ax.coords[1].set_major_formatter('dd:mm:ss')
        ax.coords[0].set_axislabel('RA')
        ax.coords[1].set_axislabel('Dec')
        ax.coords.grid(color='blue', ls='--', alpha=0.7)
        # ax.plot()
        # plt.show()
        fig.savefig('field.pdf')

    @staticmethod
    def do_astrometry(astrometry_net_API: str, path2images: str, **kwargs):
        """

        :param path2images : str
            path to fits files
        :param astrometry_net_API:
        :param kwargs:
            path2image : str or Path object
                Path to the image

            radius : float
                radius of search
            isFITS : Bool
            fast_star_mask : list-like
                mask of stars, the position of which does not match the catalog
            xy_centroids : list-like
                List of xy-coordinate of source positions
            scale_lower : float
            scale_upper : float

        :return:
        """
        from astroquery.astrometry_net import AstrometryNet
        from glob import glob
        list_of_image = glob(path2images + r'\*.f*s')
        if len(list_of_image) == 0:
            list_of_image = glob(path2images + r'\*.fit')
        ast = AstrometryNet()
        ast.api_key = astrometry_net_API
        upper = 4.8
        lower = 4.3
        r = Angle('0.054d')
        # print(r.degree)
        solve_timeout = 60
        for count, image in enumerate(list_of_image):
            print('working on', os.path.split(image)[1], str(np.round(count * 100 / len(list_of_image), 2)), '%')
            hdu_list = fits.open(image, 'update', memmap=False)
            header = hdu_list[0].header.copy()

            try:
                if header['CD1_1']:
                    hdu_list.close()
                    print('Astrometry has already been done')
                    continue
            except KeyError:
                pass

            try_again = True
            submission_id = None
            wcs_header = None
            ra = Angle(header['CURRA'] + ' hours')
            dec = Angle(header['CURDEC'] + ' degrees')
            # print(ra.degree)
            # print(dec.degree)

            # ra = Angle('346.632 degrees')
            # dec = Angle('-5.054 degrees')
            if kwargs.get('isFITS', False):
                while try_again:
                    try:
                        if not submission_id:
                            wcs_header = ast.solve_from_image(image,
                                                              submission_id=submission_id, solve_timeout=solve_timeout,
                                                              center_ra=ra.degree, center_dec=dec.degree,
                                                              radius=kwargs.get('radius', r.degree),
                                                              downsample_factor=2,
                                                              scale_units='arcminwidth', scale_type='ul',
                                                              scale_lower=kwargs.get('scale_lower', lower),
                                                              scale_upper=kwargs.get('scale_upper', upper),
                                                              publicly_visible='n', parity=2)
                        else:
                            wcs_header = ast.monitor_submission(submission_id, solve_timeout=solve_timeout)
                    except:
                        print('Timeout')
                        # submission_id = e.args[1]
                        continue
                    else:
                        try_again = False
            else:
                while try_again:
                    try:
                        if not submission_id:
                            data = hdu_list[0].data.copy()
                            hdu_list.verify('fix')
                            sources = kwargs.get('xy_centroids', get_center(data))
                            # if kwargs.get('fast_star_mask'):
                            #     sources = np.delete(sources, kwargs.get('fast_star_mask'))
                            wcs_header = ast.solve_from_source_list(sources['xcentroid'],
                                                                    sources['ycentroid'],
                                                                    header['NAXIS1'], header['NAXIS2'],
                                                                    downsample_factor=2,
                                                                    center_ra=ra.degree, center_dec=dec.degree,
                                                                    radius=r.degree,
                                                                    scale_units='arcminwidth', scale_type='ul',
                                                                    scale_lower=kwargs.get('scale_lower', lower),
                                                                    scale_upper=kwargs.get('scale_upper', upper),
                                                                    publicly_visible='n', parity=2)
                        else:
                            wcs_header = ast.monitor_submission(submission_id, solve_timeout=solve_timeout)
                    except:
                        print('Timeout')
                        # submission_id = e.args[1]
                        continue
                    else:
                        try_again = False
            if wcs_header:
                hdu_list[0].header = header + wcs_header
                hdu_list.close()
                # shutil.move(image, os.path.split(image)[0] + '/done/')
                print('done')
            else:
                hdu_list.close()
                print('astrometry solve fails')
