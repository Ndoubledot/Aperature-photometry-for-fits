import astropy, photutils, glob, ccdproc,scipy
from astropy.io import fits
from astropy  import units as u
import numpy as np
from astropy.stats import sigma_clipped_stats, SigmaClip, mad_std
from astropy.wcs import WCS
from astropy.time import Time
from photutils import  CircularAperture, CircularAnnulus, aperture_photometry
from photutils.aperture import ApertureStats
from photutils.centroids import centroid_2dg
from photutils.utils import calc_total_error
from photutils.background import Background2D, MedianBackground

def AP_Photometry(science_file,wcs_bool,wcs_file,frame_filter,gain,readout_noise,ref_cood,radius):
    if type(wcs_bool) == type(True):
        pass
    else :
        raise Exception("wcs is not a booliean")
    if type(wcs_file) == type('a'):
        pass
    else :
        raise Exception("wcs_file is not a string")
    if type(science_file) == type('a'):
        pass
    else :
        raise Exception("science_file is not a string")
    if type(frame_filter) == type('a'):
        pass
    else :
        raise Exception("frame_filter is not a string")
    if type(ref_cood) == type(np.array([1])):
        pass
    else :
        raise Exception("ref_cood is not an array")
    if type(radius) == type(float(10)):
        pass
    else :
        raise Exception("radius is not a float")
    if type(readout_noise) == type(float(10)):
        pass
    else :
        raise Exception("readout_noise is not a float")
    if wcs_bool == True :
        wcs_file=fits.open(wcs_file)
        wcs_header=wcs_file[0].header
        wcs=WCS(wcs_header)
        pix_scale = np.mean(astropy.wcs.utils.proj_plane_pixel_scales(wcs))
        pix_scale = pix_scale*u.degree
        pix_scale = pix_scale.to(u.arcsec)
        pix_scale = pix_scale/pix_scale.unit
        pix_scale= float(pix_scale)
    else:
        pix_scale=1.0
    

    file = fits.open(science_file)
    data = file[0].data
    header=file[0].header
    exposure =header['EXPTIME']
    #if frame_filter in header['FILTER'] or header['ALFLTNM']:
        #pass
    #else:
        #raise Exception("mentioned filter is not same as in header")
    if len(np.array(data).shape) == 2:
        pass
    else :
        raise Exception("Shape of data is not 2 but is ",len(data.shape))
    
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data,(data.shape[0]/10, data.shape[1]/10), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    back=bkg.background
    data = data-back

    bkg = Background2D(data,(data.shape[0]/10, data.shape[1]/10), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    back=bkg.background

    if wcs_bool == True:
        pos_wcs=wcs.wcs_world2pix(ref_cood,1)
    else :
        pos_wcs = ref_cood
    r = 10 #test aperature in units of pixel around the object to calculate its fwhm

    apertures = CircularAperture(pos_wcs, r=r) 
    aperstats = ApertureStats(data, apertures)
    fwhm = np.mean(aperstats[1:].fwhm)
    fwhm = float(fwhm/(1*u.pix))
    seeing = float(fwhm*pix_scale)
    
    r = radius
    apertures = CircularAperture(pos_wcs, r=r/pix_scale)
    annulus_aperture = CircularAnnulus(pos_wcs, r_in=(r)/pix_scale, r_out=(r+2)/pix_scale)
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    aperstats = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)
    bkg_median = aperstats.median
    aperture_area = apertures.area_overlap(data)
    total_bkg = bkg_median*np.array(aperture_area)
   
    back = ccdproc.CCDData(back,unit='adu')
    back_err=ccdproc.create_deviation(back,gain=gain*u.electron/u.adu,readnoise=readout_noise*u.electron)
    error = calc_total_error(data, back_err, gain)  
    phot_table = aperture_photometry(data, apertures, error=error) #Read out total cal_error because our code is not accurate depiction of error
    mag_1=-2.5*np.log10((phot_table['aperture_sum']-total_bkg)/exposure)
    mag_err_1=1.09*(phot_table['aperture_sum_err']/(phot_table['aperture_sum']-total_bkg))
    mag=mag_1
    mag_err=mag_err_1
    flux=(phot_table['aperture_sum']-total_bkg)/exposure
    flux_err=phot_table['aperture_sum_err']
    return flux,flux_err,mag,mag_err,fwhm,seeing

                          
