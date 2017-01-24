import matplotlib.pyplot as plt
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC012557548.dat", 'r')
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC004069063.dat", 'r') # 0.5 day period
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC010905804.dat", 'r') # 0.75 day period
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC005955731.dat", 'r') # 1.0
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC009451096.dat", 'r') # 1.25
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC007700578.dat", 'r') # 1.5
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC004902030.dat", 'r') # 1.75
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC004670267.dat", 'r') # 2.0
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC011250867.dat", 'r') # 2.25
#rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC008197406.dat", 'r') # 2.5
#rawdata.next()
#rawdata.next()
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/12557548_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/4069063_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/10905804_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/5955731_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/9451096_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/7700578_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/4902030_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/4670267_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/11250867_detrended.txt", 'r')
#spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/8197406_detrended.txt", 'r')
rawdata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC011250867_detrended_1day.txt", 'r')
spldata = open("/home/marko/Desktop/MIT Research/SplineFit/data/KIC011250867_detrended_5days.txt", 'r')
rawtime = []
rawflux = []
spltime = []
splflux = []
for line in rawdata:
	parts = line.split()
	rawtime.append(float(parts[0]))
	rawflux.append(float(parts[1]))
for line in spldata:
	parts = line.split()
	spltime.append(float(parts[0]))
	splflux.append(float(parts[1]))
fig = plt.figure()
sp1 = fig.add_subplot(211)
sp1.plot(rawtime[:-1], rawflux[:-1], 'k')
sp2 = fig.add_subplot(212, sharex=sp1, sharey=sp1)
sp2.plot(spltime, splflux, 'k')
plt.show()