#!/usr/bin/env python
# under construction
from dataio import *
import numpy as np
import scipy as sp
import os
import shutil
import random
import tarfile

def mad(data):
    # median absolute deviation of 1D array
    return np.nanmedian(abs(data - np.nanmedian(data)))



def getext_base(example,ext):
        if not(example==''):
                return os.path.basename(example)

        else:
                return ext

def getext(example,ext,result):
        if(result==''):
                assert(not example=='')
                return os.path.dirname(ext)+'/'+os.path.splitext(os.path.basename(example))[0]+os.path.splitext(os.path.basename(ext))[0]
        else:
                return result

def gstat(infile,outfile,col,llim=6):
        os.system("gstater -i %s -o - -c %d -F -g -l %d | awk '$%d<999 {print}' > %s " % (infile,col,llim,col,outfile))
        return

def getblsanal():
	infile='/home/chelsea/src/list/temp16'
	outfile='/home/chelsea/Sgnhst/temp16'
	names=[];readcolumn(names,2,infile,datformat='str')
	inpath='/home/chelsea/Sgnhst/'
	os.chdir(inpath)
	for x in names:
		os.system('grep %s Sgn-16.hst | cat >> %s' % (x, outfile))


	return
def formatfile():
#	infile='/home/chelsea/src/list/EB.tab'
	infile='/home/chelsea/src/targ.ls'
#	column=1
	column=2
#	outfile='/home/chelsea/src/list/EBlist'
#	outfile='/home/chelsea/src/list/FPlist'
	outfile='/home/chelsea/src/upload/filelist'
#	inlist=[];readcolumn(inlist,column,infile)
	inlist=[];readcolumn(inlist,column,infile)
	fout=open(outfile,mode='w')
	for x in inlist:
		fout.write('%d\n' %  (int(x)))
	return

def download():
	filelist='/home/chelsea/kepler/Fullcand/candidate_listsgn'
	colid=2;colsgn=1
	inpath='/home/chelsea/KEP/'
	opath='/home/chelsea/kepler/Fullcand/LTF'
	namelist=[];readcolumn(namelist,colid,filelist,datformat='str')
	sgnlist=[];readcolumn(sgnlist,colsgn,filelist,datformat='str')
	ext='.ltf'	

	for i in range(len(namelist)):
		filei='kplr%.9d%s' % (int(namelist[i]),ext)
		dirpath=inpath+'Sgn-%s/' % (sgnlist[i])
#		os.system('cp %s%s %s' % (dirpath,filei,opath))
#		os.system('rsync phn15:%skplr%.9d%s %s' % (dirpath,int(filei),ext,opath))
		result=os.system('rsync phn15:%s%s %s' % (dirpath,filei,opath))
		if(result==256):
			os.system('echo %s >> donwload.out' % (filei))
	return
	
def getinfo():
	inlist='/home/chelsea/src/newpipe/ana.ls'
	outlist1='/home/chelsea/src/newpipe/ana.info'
	fout1=open(outlist1,mode='w')
	colid=1;colp=3;colsgn=2
	namelist=[];readcolumn(namelist,colid,inlist,datformat='str')
	plist=[];readcolumn(plist,colp,inlist)
	sgnl=[];readcolumn(sgnl,colsgn,inlist)
	print len(namelist),len(plist)
#	inpath='/home/chelsea/mntproj/newtarg/KOIrecover/'
#	os.chdir(inpath)
	for i in range(len(namelist)):	
		inpath='/home/chelsea/mntproj/Sgn-%d' % (int(sgnl[i]))
		os.chdir(inpath)
		if(abs(plist[i])>45):
			blsafile='kplr%.9d_l.blsanal' % (int(namelist[i]))
		elif(abs(plist[i])>1):	
			blsafile='kplr%.9d.blsanal' % (int(namelist[i]))
		else:
			continue
		flag=False
		if(os.path.exists(blsafile)):	
			fana=open(blsafile,mode='r')
			line=fana.readline().split()
			while(not len(line) == 0):
				if (not (line[0].startswith('#'))):
					period=float(line[2])
					if(abs(period/abs(plist[i])-1) < 0.01):
						epoch=float(line[4])
						qvar=float(line[3])
						snr=float(line[13])
						dsp=float(line[14])
						dip=float(line[7])
						epoch+=qvar*period*0.5
						newline1='kplr%.9d %s %.7f %.7f %.2e %.2e %.2f %.2f\n' % (int(namelist[i]),namelist[i],period,epoch-2454000,dip,qvar,snr,dsp)
						newline2='%d %.6f %.6f %.3f Re\n' % (int(namelist[i]),period,plist[i],abs(period/plist[i]-1))
						fout1.write(newline1)
						flag=True
						break
				line=fana.readline().split()
			fana.close()
			if (not flag):
				newline3='kplr%.9d %s %.7f -1 -1 -1 -1\n' % (int(namelist[i]),namelist[i],plist[i])
				fout1.write(newline3)
	fout1.close()
	return


def getinfo2():
	inlist='/home/chelsea/src/list/2012_Jul/KOImulti.new'
	infofile='/home/chelsea/src/list/2012_Jul/KOI.full'
	outfile='/home/chelsea/src/list/2012_Jul/KOImulti.info'
	fout=open(outfile,mode='w')
	starid=[];readcolumn(starid,1,inlist);starid=np.array(starid)
	per=[];readcolumn(per,2,inlist);per=np.array(per)
	staridf=[];readcolumn(staridf,2,infofile);staridf=np.array(staridf)
	perf=[];readcolumn(perf,3,infofile);perf=np.array(perf)
	epof=[];readcolumn(epof,4,infofile);epof=np.array(epof)
	qf=[];readcolumn(qf,6,infofile);qf=np.array(qf)
	for i in range(len(starid)):
		x=starid[i]
		px=per[i]
		inds=(abs(perf-px)<1e-5)
		ex=epof[inds]
		qx=qf[inds]
		try:
			newline='kplr%.9d %d %.6f %f %f\n' % (int(x),int(x),px,ex,qx)
		except TypeError:
			print x,px
			newline='kplr%.9d %d %.6f %f %f\n' % (int(x),int(x),px,ex[0],qx[0])
		fout.write(newline)

	fout.close()
	return

def KOItoKIC():
	KOIlist='/home/chelsea/src/list/OfirsimpleV01'
	reflist='/home/chelsea/mntproj/reffile/fulldata/Bataha.txt'
	outfile='/home/chelsea/src/list/Ofirnew'
	fout=open(outfile,mode='w')
	KOIs=[];readcolumn(KOIs,1,KOIlist);KOIs=np.array(KOIs)
	periods=[];readcolumn(periods,2,KOIlist,datformat='str')
	KOIfull=[];readcolumn(KOIfull,1,reflist);KOIfull=np.array(KOIfull)
	KICfull=[];readcolumn(KICfull,2,reflist);KICfull=np.array(KICfull)
	for i in range(len(KOIs)):
		ind=(abs(KOIfull-KOIs[i])<0.9)
		KIC=KICfull[ind]
		if(KOIs[i]-int(KOIs[i])<0.001):
			fout.write('%d %.2f %f\n' % (KIC[0],KOIs[i],float(periods[i].strip('$'))))
	fout.close()
	return
if __name__=='__main__':
#	listgen()
#	lssgn()
#	lcgen()
#	tarsgn()
#	correctext()
#	getblsanal()
#	findsgn()
#	findsgnp()
#	findpeakinfo()
#	lstemp()
#	formatfile()
#	getdepth()
	download()
#	movedir()
#	getinfo()
#	getinfo2()
#	KOItoKIC()
