# import glob
# import os
# import shutil
#
# import astropy.io.fits as fits
# import matplotlib.cm as cm
# import numpy as np
# from astropy.convolution import Gaussian2DKernel
# from astropy.convolution import convolve
# from astropy.stats import mad_std
# from astroquery.astrometry_net import AstrometryNet
# from matplotlib import pyplot as plt
# from matplotlib.colors import Normalize as Normalize
# from photutils.detection import DAOStarFinder
# from tqdm import tqdm
#
# ast = AstrometryNet()
# ast.api_key = 'hipfhzhlzygnlvix'
# files = glob.glob(r'G:\Kislovodsk\T1c_J-20201017\*.fts')
# with tqdm(total=len(files), desc='Astrometry.net on fits') as bar:
#     print()
#     for item in files:
#         print()
#         try:
#             ##read file, copy data and header
#             hdulist = fits.open(item, 'update', memmap=False)
#             Header = hdulist[0].header.copy()
#             Data = hdulist[0].data
#             hdulist.verify('fix')
#
#             ##gaussian convolution
#             kernel = Gaussian2DKernel(x_stddev=1)
#             Data = convolve(Data, kernel)
#
#             ##extract background
#             Data -= np.median(Data)
#             Bkg_sigma = mad_std(Data)
#             ##mask bad row
#             mask = np.zeros(Data.shape, dtype=bool)
#             mask[90:110, 0:Data.shape[1]] = True
#
#             daofind = DAOStarFinder(fwhm=4.5, threshold=5. * Bkg_sigma, sharplo=0.25)
#             Sources = daofind(Data, mask=mask)
#             # print(Sources.info)
#
#             # plt.imshow(Data, cmap=cm.Greys_r, aspect='equal',
#             #            norm=Normalize(vmin=-30, vmax=150), interpolation='nearest')
#             # plt.scatter(Sources['xcentroid'], Sources['ycentroid'], s=40, facecolors='none', edgecolors='r')
#             # plt.show()
#
#             # Sort sources in ascending order
#             Sources.sort('flux')
#             Sources.reverse()
#
#             # ast.show_allowed_settings()
#
#             image_width = Header['NAXIS2']
#             image_height = Header['NAXIS1']
#             print(Sources)
#             wcs_header = ast.solve_from_source_list(Sources['xcentroid'],
#                                                     Sources['ycentroid'],
#                                                     image_width, image_height,
#                                                     solve_timeout=120,
#                                                     center_ra=346.622,
#                                                     center_dec=-5.047,
#                                                     radius=0.1,
#                                                     downsample_factor=2,
#                                                     scale_lower=0.26,
#                                                     scale_upper=0.29,
#                                                     scale_units='arcsecperpix')
#             # print(wcs_header)
#             hdulist[0].header = Header + wcs_header
#             hdulist.close()
#             shutil.move(item, os.path.split(item)[0] + '/done/')
#             print('done')
#         except Exception as e:
#             print(e)
#         bar.update(1)
#
#     bar.close()
