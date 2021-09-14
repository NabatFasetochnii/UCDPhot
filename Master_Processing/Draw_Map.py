import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from photutils import CircularAperture


def Draw_Map(Image, Header, Cat, Name, RAper, TTF, Info, XY):
    wcs = WCS(Header)
    #
    # XY = wcs.all_world2pix(Cat['Ra'], Cat['Dec'], 0)
    # XY = np.vstack((XY[0], XY[1])).T

    # find fwhm
    # FWHM = find_fwhm(Image, XY=XY[0])
    # RAper.append(4.0*FWHM)

    Image = np.log10(Image)
    X = Image.shape[1]
    Y = Image.shape[0]
    _mean, _median, _std = sigma_clipped_stats(Image[Y - 50:Y + 50, X - 50:X + 50])
    _max = _median + 10. * _std
    _min = _median - 1. * _std

    fig = plt.figure(figsize=(7, 7), frameon=False)
    ax = plt.subplot(projection=wcs, position=[0.1, 0.1, 0.8, 0.8])
    plt.imshow(Image, vmin=_min, vmax=_max, cmap='gray_r')

    # for q in range(1, len(XY)):
    #     rAper = 4*find_fwhm(Image, XY[q])
    #     aper = CircularAperture(XY[q], r=rAper)  # for all stars
    #     aper.plot(color='blue', lw=1.5, alpha=0.5)
    aper = CircularAperture(XY, r=RAper)  # for all stars
    aper.plot(color='blue', lw=1.5, alpha=0.5)
    # apertures = [CircularAperture(XY[0], r=r) for r in RAper]  # only for target
    # apertures[0].plot(color='red', lw=1.5, alpha=0.5)
    # apertures[1].plot(color='green', lw=1.5, alpha=0.5)
    # apertures[2].plot(color='yellow', lw=1.5, alpha=0.5)
    # apertures[3].plot(color='orange', lw=1.5, alpha=0.5)
    aper = CircularAperture(XY[0], r=RAper)
    aper.plot(color='red', lw=1.5, alpha=0.8)
    for i in range(0, len(XY)):
        plt.text(XY[i, 0], XY[i, 1], s=Cat['ID'][i], color='blue', alpha=0.8)

    if Header['CD1_1'] > 0:
        plt.gca().invert_xaxis()
    if Header['CD2_2'] > 0:
        plt.gca().invert_yaxis()

    Title = TTF['Name'][0] + ', ' + Info[6].datetime.strftime('%Y-%m-%d\n')
    Title += 'Filter=' + Header['FILTER']
    Title += ', aperture radius =' + '{:.1f}'.format(3600 * RAper * Info[3]) + '"'
    plt.title(Title)
    ax.coords[1].set_ticklabel(rotation=90)
    ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm:ss')
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')
    ax.coords.grid(color='blue', ls='--', alpha=0.7)
    fig.savefig(Name)
    #    plt.show()
