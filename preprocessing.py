import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from astropy.table import Table, QTable
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import simple_norm, SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from photutils import Background2D, detect_threshold, detect_sources, CircularAperture, SourceCatalog, aperture_photometry

def f(srcpath, destpath=''): #converts fitz files to fits files

	if not srcpath:
		pass
	elif srcpath[-1] != '/':
		srcpath += '/'
	if not destpath:
			destpath = srcpath
	elif destpath[-1] != '/':
		destpath += '/'
	fitzlist = glob.glob('{}*.fitz'.format(srcpath)) # returns a list of all files that match the given srcpath
	if not os.path.exists(destpath):
		os.makedirs(destpath)
	for infile in fitzlist:
		outfile = destpath+infile.split('/')[-1]
		finalfile = '{}new_{}s'.format(destpath, infile.split('/')[-1][:-1])
		infits = fits.open(outfile)
		# store extension image
		image = infits[1]
		# save extension image as primary header of new file
		newhdu = fits.PrimaryHDU(data=image.data, header=image.header)
		hdulist = fits.HDUList([newhdu])
		hdulist.writeto(finalfile, overwrite=True)
		infits.close()
		print('{} processed'.format(infile.split('/')[-1]))
	print('Finished')

f.fitz_to_fits(srcpath='/media/roisin/SP PHD U3/work/mrk1018/optical/stella_pipeline/files/2022')

#TODO automatically delete the fitz files

fls = sorted(glob.glob('*.fits'))

counts = Table(names = ('day', 'month', 'year', 'exp', 'dec_year', 'expt', 'fwhm', 'airmass', 'ref1', 'ref2', 'ref3', 'mrk1018'))

for i, fl in enumerate(fls):
    name = str(fl)
    date = name[-26:-18]
    day  = float(date[6:8]) #separate day, month and year and convert back to float
    month= float(date[4:6])
    year = float(date[0:4])
    exp  = name[-6]
    dec_year = year + (month-1)/12 + (day-1)/365.25 #convert to decimal year

    hdu = fits.open(fl)

    if hdu[0].header['FILTER'] == 'up':

        data = hdu[0].data[2020:2280, 2250:2890] #section image

        #norm = simple_norm(data, 'sqrt', percent=99.9)
        #plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
        #plt.xlim(0, data.shape[1]-1)
        #plt.ylim(0, data.shape[0]-1)
        #plt.show()

        bkg = Background2D(data, box_size=(40, 35)) #calculate background map
        threshold = bkg.background + 3.0 * bkg.background_rms #theshold for detection 3 sigma

        sigma = 4.0 * gaussian_fwhm_to_sigma   #FWHM = 4, average STELLA FWHM
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        kernel.normalize()
        segm = detect_sources(data, threshold, npixels=100, connectivity = 4) #detect sources

        #norm = ImageNormalize(stretch=SqrtStretch())
        #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
        #ax1.imshow(data, origin='lower', cmap='Greys_r', norm=norm)#
        #ax1.set_title('Background-subtracted Data')
        #ax2.imshow(segm, origin='lower', cmap=segm.cmap, interpolation='nearest')
        #ax2.set_title('Segmentation Image')
        #plt.show()

        cat = SourceCatalog(data, segm) #fetch source information
        tbl = cat.to_table()

        bool1 = (tbl['xcentroid'] < 200 ) & (tbl['xcentroid'] > 100 ) & (tbl['ycentroid'] < 100) & (tbl['ycentroid'] > 50)

        ref1  = tbl[bool1]  #isolate reference star 1
        tbl = tbl[np.invert(bool1)]

        bool2 = (tbl['xcentroid'] < 200) & (tbl['xcentroid'] > 100) & (tbl['ycentroid'] < 200) & (tbl['ycentroid'] > 150)

        ref2  = tbl[bool2]  #isolate reference star 2
        tbl = tbl[np.invert(bool2)]

        bool3 = (tbl['xcentroid'] > 500)

        ref3  = tbl[bool3]
        tbl = tbl[np.invert(bool3)]

        agn = tbl

        #ref3_x = float(ref1['xcentroid']) + 361  #find ref 3 coords from ref1
        #ref3_y = float(ref1['ycentroid']) + 62

        #ref3_x = float(ref1['xcentroid']) + 392  #find ref 3 coords from ref2
        #ref3_y = float(ref1['ycentroid']) - 20

        agn_x   = int(np.around(agn['xcentroid'], decimals=0))
        agn_y   = int(np.around(agn['ycentroid'], decimals=0))

        #extract local background for treatment of host image
        loc_bkg = bkg.background[agn_y-100:agn_y+100, agn_x-100:agn_x+100]

        fits.writeto('./local_bkg'+str(date)+'exp'+str(exp)+'.fits', loc_bkg, overwrite=True)

        # P H O T O M E T R Y

        positions = [(float(ref1['xcentroid']), float(ref1['ycentroid'])), (float(ref2['xcentroid']), float(ref2['ycentroid'])), (float(ref3['xcentroid']), float(ref3['ycentroid'])), (float(agn['xcentroid']), float(agn['ycentroid']))]
        apertures = CircularAperture(positions, r = 10/0.322)

        phot_tbl = aperture_photometry(data - bkg.background, apertures)

        expt = hdu[0].header['EXPT'] #fetch exposure time of each image
        airmass = hdu[0].header['AIRMASS']
        fwhm = hdu[0].header['FWHM']

        counts.add_row([day, month, year, exp, dec_year, expt, fwhm, airmass, phot_tbl['aperture_sum'][0], phot_tbl['aperture_sum'][1], phot_tbl['aperture_sum'][2], phot_tbl['aperture_sum'][3]])

        #norm = simple_norm(data, 'sqrt', percent=99.9)
        #plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm, interpolation='nearest')
        #apertures.plot(color='yellow', lw=1.5)
        #plt.xlim(0, data.shape[1]-1)
        #plt.ylim(0, data.shape[0]-1)
        #plt.show()

#calculate and insert photon error
counts.add_columns([np.sqrt(counts['ref1']), np.sqrt(counts['ref2']), np.sqrt(counts['ref3']), np.sqrt(counts['mrk1018'])], names=('err1', 'err2', 'err3', 'err_mrk1018'))

ascii.write(counts, 'counts_new.csv', format = 'csv', overwrite = True)
