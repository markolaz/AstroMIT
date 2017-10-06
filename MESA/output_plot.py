import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import matplotlib.gridspec as gridspec

#with open('LOGS/profiles.index') as profiles:
#	for line in profiles: pass
#	profile = line.split()[-1]

#with open('LOGS/profile'+profile+'.data') as source:
with open('LOGS/history.data') as source:
	data = source.readlines()

info = {str(data[1].split()[i]): data[2].split()[i] for i in range(len(data[1].split()))} 

params_labels = data[5].strip().split()

all_dat = [[] for _ in range(len(params_labels))]

for line in data[6:]:
	for i in range(len(params_labels)):
		all_dat[i].append(line.split()[i])

parameters = {str(params_labels[i]): all_dat[params_labels.index(str(params_labels[i]))] for i in range(len(params_labels))}

print sort(parameters.keys()), len(parameters.keys())

plotting = [
	'star_age',
	'star_mass',
	'log_Teff',
	'log_R',
	'log_g',
	'total_mass_h1',
	'total_mass_he4',
	'he_core_mass'
	]

gs = gridspec.GridSpec(3, 3)

fig = plt.figure()
for i in range(len(plotting)):
	if plotting[i] == 'star_age': continue
	ax = plt.subplot(gs[i])
	ax.plot(parameters.get(plotting[0]), parameters.get(plotting[i]), 'k')
	ax.set_xlabel(plotting[0])
	ax.set_ylabel(plotting[i])
	ax.set_title((plotting[i][0].capitalize() + plotting[i][1:]).replace('_', ' '))

plt.figtext(0.1, 0.9, 'M_0 = '+str(float(info.get('initial_mass')))+' M_sun')
plt.figtext(0.1, 0.88, 'initial Z = '+str(float(info.get('initial_z'))))

plt.tight_layout()
plt.show()

