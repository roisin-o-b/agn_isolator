import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from astropy.table import Table, QTable
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import simple_norm, SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from photutils.background import Background2D
from photutils.aperture import CircularAperture, aperture_photometry
from photutils.segmentation import SourceCatalog, detect_sources, detect_threshold

# TODO add other functions for images (show, save?)
# TODO options for background areas etc.

def GetFitsFileList():

    fitspath = input("Enter the fits files directory path: ")

    # Creates variable with list of fits files
    file_list = sorted(glob.glob(fitspath + '*.fits'))

    return file_list

def GetCounts(file_list):

    # Sets up array with observation info
    counts = Table(names = ('Day', 'Month', 'Year', 'Number of Exposure',\
                            'Decimal Year', 'Exposure Time', 'FWHM',\
                            'Airmass', 'Star 1 Counts', 'Star 2 Counts',\
                            'Star 3 Counts', 'Mrk 1018 Counts'))

    # Extracts information from files
    for i, fits_file in enumerate(file_list):
        name = str(fits_file)
        date = name[-26:-18]
        number_exp  = name[-6]

        # Separates day, month and year and converts to float
        day  = float(date[6:8])
        month= float(date[4:6])
        year = float(date[0:4])
        # Converts to decimal year
        decimal_year = year + (month-1)/12 + (day-1)/365.25

        hdu = fits.open(fits_file)

        # Specifies u' optical filter
        if hdu[0].header['FILTER'] == 'up':

            # Extracts information from file header
            exposure_time= hdu[0].header['EXPT']
            airmass = hdu[0].header['AIRMASS']
            fwhm = hdu[0].header['FWHM']

            # Focuses on area containing Mrk 1018
            data = hdu[0].data[2020:2280, 2250:2890]

            # Creates background map
            bkg = Background2D(data, box_size=(40, 35))
            # Specifies 3 sigma minimum for source detection
            threshold = bkg.background + 3.0 * bkg.background_rms

            # Specifies a FWHM of 4 - the average STELLA FWHM
            # Converts this to sigma to use in kernel
            sigma = 4.0 * gaussian_fwhm_to_sigma
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            kernel.normalize()
            # Isolates light sources
            source_segmentation = detect_sources(data, threshold, npixels=100,\
                                  connectivity = 4)

            # Fetches table with source information
            cat = SourceCatalog(data, source_segmentation)
            source_info = cat.to_table()

            # For Star 1, Star 2 and Star 3:
            # Specifies limits on location of each light source's centroid
            # Applies limits to source table and stores result in variable
            # Removes source from source table

            loc_limits1 = \
            (source_info['xcentroid'] < 200) & \
            (source_info['xcentroid'] > 100) & \
            (source_info['ycentroid'] < 100) & \
            (source_info['ycentroid'] > 50)

            star1  = source_info[loc_limits1]
            source_info = source_info[np.invert(loc_limits1)]

            loc_limits2 = \
            (source_info['xcentroid'] < 200) & \
            (source_info['xcentroid'] > 100) & \
            (source_info['ycentroid'] < 200) & \
            (source_info['ycentroid'] > 150)

            star2  = source_info[loc_limits2]
            source_info = source_info[np.invert(loc_limits2)]

            loc_limits3 = (source_info['xcentroid'] > 500)

            star3  = source_info[loc_limits3]
            source_info = source_info[np.invert(loc_limits3)]

            # The only source left in source_info must be the AGN
            # Assigns agn source and centroids variables
            agn = source_info
            agn_x   = int(np.around(agn['xcentroid'], decimals=0))
            agn_y   = int(np.around(agn['ycentroid'], decimals=0))

            # TODO specify file path for local bkgs

            # Extracts local background image for later use
            # Save in local_bkg folder
            local_bkg = bkg.background[agn_y-100:agn_y+100, agn_x-100:agn_x+100]
            fits.writeto('./local_bkg'+str(date)+'exp'+str(number_exp)+'.fits',\
                local_bkg, overwrite=True)

            #TODO do I need these backslashes?

            # Defines positions and apertures for photometric measurements
            positions = [\
                        (float(star1['xcentroid']), float(star1['ycentroid'])),\
                        (float(star2['xcentroid']), float(star2['ycentroid'])),\
                        (float(star3['xcentroid']), float(star3['ycentroid'])),\
                        (float(agn['xcentroid']), float(agn['ycentroid'])) \
                        ]
            apertures = CircularAperture(positions, r = 10/0.322)

            # Table of background-subtracted counts for each object
            phot_tbl = aperture_photometry(data - bkg.background, apertures)

    # Adds observation info and photometry results to counts array
    counts.add_row([day, month, year, number_exp, \
                    decimal_year, exposure_time, fwhm, airmass, \
                    phot_tbl['aperture_sum'][0], \
                    phot_tbl['aperture_sum'][1], \
                    phot_tbl['aperture_sum'][2], \
                    phot_tbl['aperture_sum'][3]])

    # TODO create filepath to save results
    # Saves results to .csv document
    ascii.write(counts, 'counts.csv', format = 'csv', overwrite = True)
