import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random

def axes_rotation(input_array):
	"""reversing x and y data to correct for rotation done by matshow (if not done, X and Y axes will be swapped)"""
	output_array = np.zeros_like(input_array)
	for i in range(len(input_array)):
		for j in range(len(input_array[i])):
			for k in range(len(input_array[i])):
				output_array[i][j, k] = input_array[i][k, j]
	return output_array

def psf_bilinear_interpolation(psfarray1, psfarray2, psfarray3, psfarray4):
	weighted_psf = np.zeros((13, 13))
	x1 = y1 = 0
	x2 = y2 = len(psfarray1)
	for j in np.linspace(0, x2-1, x2):
		for k in np.linspace(0, y2-1, y2):
			V_11 = psfarray1[int(j), int(k)]
			V_21 = psfarray2[int(j), int(k)]
			V_12 = psfarray3[int(j), int(k)]
			V_22 = psfarray4[int(j), int(k)]
			w_11 = (k - y2) / (y2 - y1) * (j - x2) / (x2 - x1)
			w_21 = -(k - y2) / (y2 - y1) * (j - x1) / (x2 - x1)
			w_12 = -(k - y1) / (y2 - y1) * (j - x2) / (x2 - x1)
			w_22 = (k - y1) / (y2 - y1) * (j - x1) / (x2 - x1)
			weighted_psf[int(j), int(k)] = w_11*V_11 + w_21*V_21 + w_12*V_12 + w_22*V_22
	return weighted_psf

np.random.seed(100)

datapsf1 = np.random.random((13,13))
datapsf2 = np.random.random((13,13))
#datapsf3 = np.random.random((13,13))
datapsf4 = np.random.random((13,13))
datapsf5 = np.random.random((13,13))
datapsf3 = psf_bilinear_interpolation(datapsf1, datapsf2, datapsf4, datapsf5)

dat_array = []
dat_array.extend((datapsf1, datapsf2, datapsf3, datapsf4, datapsf5))

dat_array = axes_rotation(dat_array)

plt.figure()
gs1 = gridspec.GridSpec(3, 3)
gs1.update(wspace=0.0001, hspace=0.0001)

plots = np.linspace(0, 8, 5)

for i in range(9):
	if i in plots:
		ax = plt.subplot(gs1[i])
		plt.axis('on')
		plt.tick_params(axis=u'both', which=u'both', length=0)
		ax.set_aspect('equal')
		ax.matshow(dat_array[i/2], origin='lower', cmap='Blues', alpha=1)
		#for (i, j), z in np.ndenumerate(dat_array[i/2]):
		#    ax.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
		ax.set_xticks(np.arange(0, 13, 1))
		ax.set_yticks(np.arange(0, 13, 1))
		ax.set_xticklabels(np.arange(0, 13, 1))
		ax.set_yticklabels(np.arange(0, 13, 1))
		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		title = "{} {}".format('PSF', (i/2)+1)
		if i > 4:
			title = "{} {}".format('PSF', (i/2))
		ax.set_title(title)
		if i == 4:
			ax.set_title("Weighted PSF")

mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
plt.show()





#fig = plt.figure()
#fig.add_subplot(331)
#ax_psf1 = plt.gca()
#ax_psf1.matshow(datapsf1, origin='lower', cmap='Blues', alpha=1)
##for (i, j), z in np.ndenumerate(datapsf1):
##    ax_psf1.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
#ax_psf1.set_xticks(np.arange(-.5, 12, 1))
#ax_psf1.set_yticks(np.arange(-.45, 12, 1))
#ax_psf1.set_xticklabels(np.arange(0, 13, 1))
#ax_psf1.set_yticklabels(np.arange(0, 13, 1))
##ax_psf1.grid(color='black', linestyle='-')
#ax_psf1.set_title('PSF 1')
#
#fig.add_subplot(333)
#ax_psf2 = plt.gca()
#ax_psf2.matshow(datapsf2, origin='lower', cmap='Blues', alpha=1)
##for (i, j), z in np.ndenumerate(datapsf2):
##    ax_psf2.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
#ax_psf2.set_xticks(np.arange(-.5, 12, 1))
#ax_psf2.set_yticks(np.arange(-.45, 12, 1))
#ax_psf2.set_xticklabels(np.arange(0, 13, 1))
#ax_psf2.set_yticklabels(np.arange(0, 13, 1))
##ax_psf2.grid(color='black', linestyle='-')
#ax_psf2.set_title('PSF 2')
#
#fig.add_subplot(335)
#ax_psfw = plt.gca()
#ax_psfw.set_aspect('equal', adjustable='box')
#ax_psfw.matshow(datapsfw, origin='lower', cmap='Blues', alpha=1)
##for (i, j), z in np.ndenumerate(datapsfw):
##    ax_psfw.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
#ax_psfw.set_xticks(np.arange(-.5, 12, 1))
#ax_psfw.set_yticks(np.arange(-.45, 12, 1))
#ax_psfw.set_xticklabels(np.arange(0, 13, 1))
#ax_psfw.set_yticklabels(np.arange(0, 13, 1))
##ax_psfw.grid(color='black', linestyle='-')
#ax_psfw.set_title('Weighted PSF')
#
#fig.add_subplot(337)
#ax_psf3 = plt.gca()
#ax_psf3.matshow(datapsf3, origin='lower', cmap='Blues', alpha=1)
##for (i, j), z in np.ndenumerate(datapsf3):
##    ax_psf3.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
#ax_psf3.set_xticks(np.arange(-.5, 12, 1))
#ax_psf3.set_yticks(np.arange(-.45, 12, 1))
#ax_psf3.set_xticklabels(np.arange(0, 13, 1))
#ax_psf3.set_yticklabels(np.arange(0, 13, 1))
##ax_psf3.grid(color='black', linestyle='-')
#ax_psf3.set_title('PSF 3')
#
#fig.add_subplot(339)
#ax_psf4 = plt.gca()
#ax_psf4.matshow(datapsf4, origin='lower', cmap='Blues', alpha=1)
##for (i, j), z in np.ndenumerate(datapsf4):
##    ax_psf4.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
#ax_psf4.set_xticks(np.arange(-.5, 12, 1))
#ax_psf4.set_yticks(np.arange(-.45, 12, 1))
#ax_psf4.set_xticklabels(np.arange(0, 13, 1))
#ax_psf4.set_yticklabels(np.arange(0, 13, 1))
#ax_psf4.xaxis.set_label_position('bottom')
##ax_psf4.grid(color='black', linestyle='-')
#ax_psf4.set_title('PSF 4')
#
#plt.subplots_adjust(wspace=0.001, hspace=0.001)
#
#plt.savefig('test1.png')