# TESS Development Notes

origins:

matplotlib scatter: 0,0 is lower left
matplotlib img: 0,0 is upper left
fits: 0,0 is upper left?
psf: 0,0 is lower left

indexing:

matplotlib gridspec: col, row; start at 0
np.array indexing: row, col; start at 0
np.random.random: col, row

nx <-> ncol
ny <-> nrow

pixels in hdulist[0].data are counting electrons (not ADU)
photons converted to electrons just prior to PSF

each pixel in hdulist[0].data should be int32 when written to disk


according to zach:
good test of PSF: shift by 1/100th of a pixel, flux should be conserved
