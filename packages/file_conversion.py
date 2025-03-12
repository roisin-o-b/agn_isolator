import glob
import os
from astropy.io import fits

# .fitz files must be convert to .fits before they can be manipulated

def ask_for_paths(): # Source path of .fitz files and destination path
                   # for .fits files requested from user

    srcpath = input("Enter the source path: ")
    destpath= input("Enter the destination path: ")

    return srcpath, destpath

# Many thanks to Fabian Emmerich for providing the convert function

def convert_fitz_to_fits(srcpath, destpath):

	if not srcpath:
		pass
	elif srcpath[-1] != '/':
		srcpath += '/'
	if not destpath:
			destpath = srcpath
	elif destpath[-1] != '/':
		destpath += '/'
	fitzlist = glob.glob('{}*.fitz'.format(srcpath))
 # Returns a list of all files that match the given srcpath
 # Problem with returning empty list
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

#TODO option to delete the fitz files
