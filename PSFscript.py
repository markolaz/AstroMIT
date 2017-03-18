import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

datapsf1 = np.random.random((13,13))
datapsf2 = np.random.random((13,13))
datapsf3 = np.random.random((13,13))
datapsf4 = np.random.random((13,13))
datapsfw = np.random.random((13,13))

dat_array = []
dat_array.append(datapsf1)
dat_array.append(datapsf2)
dat_array.append(datapsfw)
dat_array.append(datapsf3)
dat_array.append(datapsf4)

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


plt.figure()
gs1 = gridspec.GridSpec(3, 3)
gs1.update(wspace=0.0001, hspace=0.0001) # set the spacing between axes. 

for i in range(9):
	if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
		print i
		ax1 = plt.subplot(gs1[i])
		plt.axis('on')
		ax1.set_xticklabels([])
		ax1.set_yticklabels([])
		ax1.set_aspect('equal')
		ax1.matshow(dat_array[i/2], origin='lower', cmap='Blues', alpha=1)

plt.savefig('test2.png')