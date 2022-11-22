import numpy as np
import glob
from matplotlib import pyplot as plt
from AP import AP_Photometry
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from astropy import units as u
import astropy.coordinates as coord

#following fits file is of Mrk 142 an AGN in 27 x 27 arc minute frame
file = "example.fits"
frame = fits.open(file)
header = frame[0].header
#some fits files have wcs(world co-ordinate system) information mentioned in them
#which tells us which pixel corresponds to which co-ordinate in sky
wcs = WCS(header)

Vizier.ROW_LIMIT = -1

#vizier query below gives us objects in 13.5 arc minute of radius of centre of the frame
objects=Vizier.query_region(coord.SkyCoord(ra=wcs.wcs.crval[0],dec=wcs.wcs.crval[1],unit =(u.degree,u.degree),frame='fk5'),radius=13.5*u.arcmin,catalog=["NOMAD"])
ref_cood = []
stnd_mag = []
for row in objects[0]:
    if row['Bmag']<17:
        ref_cood.append(np.array([row['RAJ2000'],row['DEJ2000']]))
        stnd_mag.append(row['Bmag'])
    else :
        pass
#we use vizier query to find out ra and dec of all the objects which lie in frame
#and to create numpy array of those co-ordinates, if we have co-ordinates of objects then we dont need to use vizier query
# we call these array of co-ordinates ref_cood
ref_cood = np.array(ref_cood)
#from the same query we create array of standard magnitude of the objects as well
stnd_mag = np.array(stnd_mag)
#AP_photometry has multiple inputs which you can explore in its main script called AP but to list the few
#science_frame : frame in which we do our photometry
#wcs_bool : True means science_frame has wcs file or wcs is mentioned in header of the file
    #in that case radius should be given in arc seconds 
    #if wcs_bool = False then it means wcs information is not available
    #in that case radius should be given in units of pixels
#wcs_file : if fits file has wcs information in its header then we can use same file as wcs_file
    #if fits file does not contain the wcs information but wcs information is on seperate file then that file should be given as wcs_fil
#Frame_filter : we can mention the filter in which the obervation was taken in the example we have metioned B filter
#gain and read out noise: these are technical information of the CCD which is nessecary to calculate errorbars accurately
    #for the example we have taken 1.0,1.0 as gain and read out noise
#radius: it is radius of the aperature which will be drawn around the object
    # if wcs information is present then give radius in arc seconds if not then give radius in pixel

#AP_Photometry(science_file,wcs_bool,wcs_file,frame_filter,gain,readout_noise,ref_cood,radius)
result = AP_Photometry(file,True,file,'B',1.0,1.0,ref_cood,4.5)
#above funtion returns flux, error in flux measurements, magnitude, error in magnitude measurements, fwhm of object in pixels
    # and seeing in arc seconds
#thus if ref_cood of 100 objects is given then it will return 6 lists of 100 objects each list with above mentioned data
mag = result[2]
error = result[3]
mag = np.array(mag)
#we plot instrumental magnitude vs standard magnitude which ideally result in straight line however some objects are variable
plt.scatter(stnd_mag,mag)
plt.errorbar(stnd_mag,mag,yerr=error,fmt=' ')
plt.xlabel('standard mag')
plt.ylabel('instrumental mag')
plt.show()
