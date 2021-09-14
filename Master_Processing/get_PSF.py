import warnings
from math import isnan

import numpy as np
from pylab import indices
from pylab import ravel
from scipy import optimize

# import matplotlib
# import matplotlib.pyplot
# import matplotlib.cm as cm
# from  matplotlib.colors import Normalize as Normalize

warnings.simplefilter("ignore")


##################################################################
# calculate center of mass
def centroid(R1, R2, R3, arr):
    # total = 0.
    Ry = arr.shape[0] / 2
    Rx = arr.shape[1] / 2

    # mask
    X_index = np.arange(0, arr.shape[1], 1)  # index array
    Y_index = np.arange(0, arr.shape[0], 1)  # index array
    distance = np.sqrt(np.power(np.ones(arr.shape) * (X_index[None, :] - Rx), 2) + np.power(
        np.ones(arr.shape) * (Y_index[:, None] - Ry), 2))  # distance array

    # mean sky
    annulus_mask = np.copy(distance)
    annulus_mask[annulus_mask < R2] = 0.
    annulus_mask[annulus_mask > R3] = 0.
    annulus_mask[annulus_mask > 0] = 1.
    masked = arr * annulus_mask
    MSky = np.median(masked[np.nonzero(masked)])
    MSky = np.nan_to_num(MSky)

    # centroid
    # aperture_mask = np.copy(distance)
    distance[distance <= R1] = 1.
    distance[distance > R1] = 0.
    masked = arr * distance
    total = np.sum(masked)

    X = np.sum(masked * X_index[None, :]) / total
    Y = np.sum(masked * Y_index[:, None]) / total
    return X - arr.shape[1] / 2, Y - arr.shape[0] / 2, MSky


##################################################################
# D2 moffat fitter
def D2_moffat_full(A, B, C, D, E, F, x0, y0):  # B=1/sigma_x^2, C=1/sigma_y^2, E=beta
    try:
        return lambda y, x: A * (1 + ((x - x0) * B) ** 2. +
                                 ((y - y0) * C) ** 2. + ((x - x0) * (y - y0) * (D ** 2.))) ** (-E) + F
    except:
        return (None)


# read for correction of invalid value encountered in power
# http://stackoverflow.com/questions/16990664/scipy-minimize-uses-a-nonetype


####################################################################
# def D2_gauss(A, B, C, D, x0, y0):
#     return lambda y,x: A*np.exp(-(((x0-x)/B)**2 +((y0-y)/C)**2)/2) + D
##################################################################


def D2_moffat_fitter(ROI, MSKY, x_coo, y_coo, R3):
    x0 = x_coo - np.floor(x_coo) + R3
    y0 = y_coo - np.floor(y_coo) + R3

    #     try:
    #  moffat
    params = (ROI.max(), 0.3, 0.3, 0.1, 5.0, MSKY, x0, y0)
    errorfunction = lambda p: ravel(D2_moffat_full(*p)(*indices(ROI.shape)) - ROI)
    p, success = optimize.leastsq(errorfunction, params, maxfev=1000, ftol=0.05)
    #     print(p)
    #####################test#########################
    #        PSF = D2_moffat_full(*p)(*np.indices(ROI.shape))
    #        X_index = np.arange(0, ROI.shape[1], 1) ## X index array
    #        Y_index = np.arange(0, ROI.shape[0], 1) ## Y index array
    #        distance = np.sqrt(np.power(np.ones(ROI.shape)*(X_index[None, :]-p[6]), 2) +
    #                           np.power(np.ones(ROI.shape)*(Y_index[:, None]-p[7]), 2))
    #
    #        shift = distance.ravel()
    #        flux  = ROI.ravel()
    #        model = PSF.ravel()
    #        fig_resudial_graph = matplotlib.pyplot.figure(3, figsize=(5, 3))
    #        ax_resudial_graph = fig_resudial_graph.add_subplot(111)
    #        matplotlib.pyplot.cla()
    #        ax_resudial_graph.plot(shift, flux,'g.')
    #        ax_resudial_graph.plot(shift, model,'r.')
    #        matplotlib.pyplot.show()
    #
    #        matplotlib.pyplot.plot(shift, (flux-model)/np.sqrt(flux),'b.')
    #        matplotlib.pyplot.show()
    #
    #    ##        matplotlib.pyplot.cla()
    #    ##        p[5] = 0
    #    ##        print(p[5])
    #        PSF = D2_moffat_full(*p)(*np.indices(ROI.shape))
    #        ROI = ROI - PSF
    #        med=np.median(ROI)
    #        stdv=np.std(ROI)
    #        print(med)
    #        matplotlib.pyplot.imshow(ROI, cmap=cm.Greys_r, aspect='equal',
    #                             norm= Normalize(vmin=med-stdv*2., vmax=med+stdv*5.), interpolation='nearest')
    #        matplotlib.pyplot.show()
    ########################################################

    return (p[1], p[2], p[3], p[4], p[5])  # , w


#     except:
#         return (None)

##################################################################
####
def get_PSF(Data, XY_coo, FWHM):
    PSF_model = []
    R1 = np.ceil(FWHM)  # estimation of aperture radii
    R2 = np.ceil(FWHM * 3.)  # sky annulus inner radii
    R3 = np.ceil(FWHM * 6.)  # sky annulus outer radii

    for ii in range(0, len(XY_coo)):  # for every star from coo file
        x_coo = XY_coo[ii, 0]
        y_coo = XY_coo[ii, 1]
        #         print(x_coo, y_coo)
        ROI = np.copy(Data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])  # copy small area
        offset = centroid(R1, R2, R3, ROI)  # search centroid, Gauss sigma and mean sky

        if isnan(offset[0]) == False and isnan(offset[1]) == False:
            x_coo = x_coo + offset[0]
            y_coo = y_coo + offset[1]
            MSKY = offset[2]
            ROI = np.copy(Data[int(y_coo - R3):int(y_coo + R3), int(x_coo - R3):int(x_coo + R3)])
            param = D2_moffat_fitter(ROI, MSKY, x_coo, y_coo, R3)  # fit 2D moffat psf
            if param != None:
                PSF_model.append(param)
        else:
            pass

    PSF_model = np.asarray(PSF_model)

    ####focal plane tilt test here
    ##    matplotlib.pyplot.plot(XY_coo[:, 0], PSF_model[:, 0],'b.')
    ####    matplotlib.pyplot.show()
    ##    matplotlib.pyplot.plot(XY_coo[:, 0], PSF_model[:, 1],'r.')
    ##    matplotlib.pyplot.show()
    ##
    ##    matplotlib.pyplot.plot(XY_coo[:, 1], PSF_model[:, 0],'b.')
    ####    matplotlib.pyplot.show()
    ##    matplotlib.pyplot.plot(XY_coo[:, 1], PSF_model[:, 1],'r.')
    ##    matplotlib.pyplot.show()

    ##    matplotlib.pyplot.plot(XY_coo[:, 0], PSF_model[:, 2],'b.')
    ####    matplotlib.pyplot.show()
    ##    matplotlib.pyplot.plot(XY_coo[:, 1], PSF_model[:, 2],'r.')
    ##    matplotlib.pyplot.show()
    ##
    ##    matplotlib.pyplot.plot(XY_coo[:, 0], PSF_model[:, 3],'b.')
    ####    matplotlib.pyplot.show()
    ##    matplotlib.pyplot.plot(XY_coo[:, 1], PSF_model[:, 3],'r.')
    ##    matplotlib.pyplot.show()
    ##
    ##    PSF_error = np.std(PSF_model,0)
    ##    print(PSF_error)
    PSF_model = np.median(PSF_model, 0)
    ##    print(PSF_model)
    ##    print(PSF_model / PSF_error)
    return (PSF_model)
