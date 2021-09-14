import numpy as np


def PSF2FWHM(PSF_model):
    try:
        phi = np.arctan((PSF_model[2]**2.0)/(PSF_model[0]**2.0-PSF_model[1]**2.0))/2.0
        alpha1 = np.sqrt(2.0/(PSF_model[0]**2.0+PSF_model[1]**2.0+PSF_model[2]**2.0/np.sin(2.0*phi)))
        alpha2 = np.sqrt(np.fabs(1/(PSF_model[0]**2.0+PSF_model[1]**2.0-1.0/(alpha1**2.0))))
        FWHM1 = 2.0*alpha1*np.sqrt(2.0**(1.0/(PSF_model[3]))-1.0)
        FWHM2 = 2.0*alpha2*np.sqrt(2.0**(1.0/(PSF_model[3]))-1.0)
        return(FWHM1, FWHM2, phi, PSF_model[4])
    
    except:
        print('Wrong PSF model')
        return (0., 0., 0., 0.)
