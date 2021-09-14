import math

import numpy as np
from astropy.stats import mad_std, sigma_clip
from photutils import CircularAperture
from photutils import aperture_photometry
from scipy import optimize

##
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# from  matplotlib.colors import Normalize as Normalize

################################################################
##global variables
PSF_model = []


##################################################################
def D2_moffat(A, F, x0, y0):
    B = PSF_model[0]
    C = PSF_model[1]
    D = PSF_model[2]
    E = PSF_model[3]
    return lambda y, x: A * (1 + ((x - x0) * B) ** 2. +
                             ((y - y0) * C) ** 2. + ((x - x0) * (y - y0) * (D ** 2.))) ** (-E) + F


##################################################################
def D2_moffat_phot(ROI, x_coo, y_coo, R1, R2, R3):
    x0 = x_coo - math.floor(x_coo) + ROI.shape[1] / 2.
    y0 = y_coo - math.floor(y_coo) + ROI.shape[0] / 2.

    # make mask
    grid = np.indices(ROI.shape)
    mask = np.sqrt((grid[1] - x0) * (grid[1] - x0) + (grid[0] - y0) * (grid[0] - y0))
    MaxPix = np.max(ROI[mask < R1])
    # calc sigma clipped background
    ROI_Copy = np.copy(ROI)
    ROI_Copy[mask < R2] = np.nan
    ROI_Copy[mask > R3] = np.nan

    Sky = sigma_clip(ROI_Copy, sigma=3, maxiters=5, stdfunc=mad_std).filled(np.nan)
    NSky = np.count_nonzero(~np.isnan(Sky))
    median_Sky = np.nanmedian(Sky)
    sigma_Sky = np.nanstd(Sky)

    #     print(median_Sky)
    #     plt.imshow(Sky)
    #     plt.show()

    # recompute center of star
    params = (MaxPix, median_Sky, x0, y0)
    errorfunction = lambda p: np.ravel(D2_moffat(*p)(*np.indices(ROI.shape)) - ROI)
    p, success = optimize.leastsq(errorfunction, params, maxfev=1000, ftol=0.05)

    #    med=np.median(ROI)
    #    stdv=np.std(ROI)
    #    plt.imshow(ROI, cmap=cm.Greys_r, aspect='equal',
    #                             norm= Normalize(vmin=med-stdv*2., vmax=med+stdv*5.), interpolation='nearest')
    #
    #    plt.plot(p[2], p[3], 'ro')
    #    plt.show()

    return (p[2] + math.floor(x_coo) - (ROI.shape[1] / 2.), p[3] + math.floor(y_coo) - (ROI.shape[0] / 2.),
            median_Sky, NSky, sigma_Sky, MaxPix)


##############################
def Phot(Data, XY, PSF, RAper, FWHM, Gain, Rnoise, Saturation):
    global PSF_model
    PSF_model = PSF

    R1 = np.ceil(FWHM * 1)  # estimation of aperture radii
    R2 = np.ceil(FWHM * 4.)  # sky annulus inner radii
    R3 = np.ceil(FWHM * 7.)  # sky annulus outer radii

    # outputs empty lists
    coo = []
    sky = []
    nsky = []
    ssky = []
    maximum = []
    flux = []

    for ii in range(0, len(XY)):
        x_coo = XY[ii, 0]
        y_coo = XY[ii, 1]
        _coo = (0., 0.)
        _sky = np.nan   # sky level
        _nsky = np.nan  # number of pixel for sky
        _ssky = np.nan  # sigma of sky
        _max = np.nan   # max pixel

        # check maximum level for linearity
        if R3 < y_coo < (Data.shape[0] - R3) and R3 < x_coo < (Data.shape[1] - R3):
            # copy subarray
            ROI = np.copy(Data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])
            #             try:
            param = D2_moffat_phot(ROI, x_coo, y_coo, R1, R2, R3)
            if np.isnan(param[0]) == 0:
                _coo = (param[0], param[1])
                _sky = param[2]
                _nsky = param[3]
                _ssky = param[4]
                if param[5] < Saturation:
                    _max = param[5]
        #             except:
        #                 print ('#', ii, 'photomety failed')
        #                 pass

        # append 0 if photometry failed
        coo.append(_coo)
        maximum.append(_max)
        sky.append(_sky)
        nsky.append(_nsky)
        ssky.append(_ssky)

    coo = np.array(coo)
    coo = np.around(coo, 3)
    sky = np.array(sky)
    sky = np.around(sky, 3)
    nsky = np.array(nsky)
    ssky = np.around(ssky, 3)
    maximum = np.array(maximum)
    apertures = [CircularAperture(coo, r=r) for r in RAper]
    # aper = CircularAperture(coo, r=RAper)
    #     print(aper)
    phot_table = aperture_photometry(Data, apertures, method='exact')
    # phot_table = aperture_photometry(Data, aper, method='exact')
    bkg_sum = [sky * a.area for a in apertures]
    # bkg_sum = sky * aper.area
    #     print (phot_table)    
    f = phot_table['aperture_sum_0'] - bkg_sum[0]
    f = np.array(f)
    flux.append(f)
    f = phot_table['aperture_sum_1'] - bkg_sum[1]
    f = np.array(f)
    flux.append(f)
    f = phot_table['aperture_sum_2'] - bkg_sum[2]
    f = np.array(f)
    flux.append(f)
    # f = phot_table['aperture_sum_3'] - bkg_sum[3]
    # f = np.array(f)
    # flux.append(f)

    #     med=np.median(Data)
    #     stdv=np.std(Data)
    #     plt.imshow(Data, cmap=cm.Greys_r, aspect='equal',
    #                              norm= Normalize(vmin=med-stdv*2., vmax=med+stdv*5.), interpolation='nearest')
    #  ###    aper = CircularAperture(coo, r=R1)
    #    aper.plot(color='blue', lw=1.5, alpha=0.5)
    #    aper = CircularAperture(coo, r=R2)
    #    aper.plot(color='green', lw=1.5, alpha=0.5)
    #    aper = CircularAperture(coo, r=R3)
    #    aper.plot(color='red', lw=1.5, alpha=0.5)
    #    plt.show()

    mag = [20. - 2.5 * np.log10(_f) for _f in flux]
    err = [1.0857 * np.sqrt(_f * Gain + a.area * (ssky * ssky * Gain)) / (_f * Gain)
           for _f in flux for a in apertures]
    sn = [1.0857 / error for error in err]

    return flux, sn, mag, err, sky, maximum

################################################################################
