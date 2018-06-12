#!/usr/bin/env python
#
#  write header and data in the ensim r2c ascii format
#  R Benoit may june 2018
#  robert.benoit.47@gmail.com
#
import numpy as np
import glob
import argparse
import cPickle as pik
from   grib2r2c import *
#
#
#
args = argparse.ArgumentParser(description='This python script converts grib2 data to the r2c format')
args.add_argument('-var',dest='var',required=True,\
                  help='Variable to be extracted. Typically either PCP or TGL')
args.add_argument('-indir',dest='indir',required=True,\
                  help='Directory path containing the input GRIB2 files')
args.add_argument('-interval',dest='interval',required=False,\
                  default=3,\
                  help='Time interval (hours) between desired samples to be saved. Default is 3')
args.add_argument('-verbose',dest='verbose',required=False,\
                  default=False,\
                  help='verbosity.  True or False. Default is False')
args.add_argument('-crop',dest='crop',required=False,\
                  default=False,\
                  help='To extract a SUBSET of the grib array.  True or False. Default is False')
args.add_argument('-lat',dest='lat',required=False,\
                  default=np.nan,\
                  help='latitude of center of SUBSET of the grib array.  degrees. Default is nan')
args.add_argument('-lon',dest='lon',required=False,\
                  default=np.nan,\
                  help='longitude of center of SUBSET of the grib array.  degrees. Default is nan')
args.add_argument('-width',dest='width',required=False,\
                  default=np.nan,\
                  help='width (km) of SUBSET of the grib array.  degrees. Default is nan')
args.add_argument('-outLALO',dest='outLALO',required=False,\
                  default='NO',\
                  help='To output latitude and longitude of selected grid array.  NO/CORNERS/ALL. Default is NO')
args.add_argument('-contour',dest='contour',required=False,\
                  default='nil',\
                  help='To plot (contour or filled-contour the gridded array.  nil/contour/contourf. Default is nil')
args.add_argument('-r2cprefix',dest='r2cprefix',required=False,\
                  default='test',\
                  help='Prefix for name of the r2c file output .  Default is test')
args.add_argument('-outdir',dest='outdir',required=False,\
						default='.',\
                  help='Directory path for the output r2c file.  Default is .')
#
args.parse_args(namespace=args)
#
var=args.var  # PCP or TGL
interval=int(args.interval)  # number of hours between desired samples to be saved
indir=args.indir
verbose=args.verbose
r2cprefix=args.r2cprefix
outdir=args.outdir
crop={}
crop['crop']=bool(args.crop)
crop['lat']=float(args.lat)
crop['lon']=float(args.lon)
crop['width']=float(args.width)
#
print "script arguments. \nvar:",var,"\ninterval:", interval, "\ndir:",indir,"\nverbose:",verbose,"\ncrop:",crop
#print type(var), type(interval), type(indir)
#
gribFiles=sorted(glob.glob(indir+"*"+var+"*grib2"))
#
#
if len(gribFiles) == 0:
	print "no GRIB2 files found. quit"
	quit()
#
r2cname=outdir + "/" + r2cprefix + "_" + var + ".r2c"   #"test_"
r2cFile=open(r2cname,'w')
r2cHeader=True
iFile=1
pickSave=False  #True
#
for gribFile in gribFiles:
	[grib,saveit]=get_grib (gribFile,interval,crop,verbose)
	#
	if saveit:
		print ("data obtained from file %s" % grib['GribName'])
		if verbose:
			print ("keys in grib dictionary")
			print ("===================")
			for key in grib.keys():
				print (key)
			print ("paramerName:%s" % grib['parameterName']+'.')
			print ("size of values in grib %d" % np.size(grib['values']))
			print ("===================")
		#
		if pickSave:
			pikfile=var+str(iFile)+".pickle"
			pikf=open(pikfile,'wb')
			pik.dump(grib['values'],pikf)
			pikf.close()
		#
		ok=put_r2c(grib,r2cFile,'R Benoit',FrameNumber=iFile,doHeader=r2cHeader,verbose=verbose)
		r2cHeader=False
		iFile+=1
		if not ok:
			print "problem in producing r2c file:",r2cFileName
#
r2cFile.close()
print "\nExtracted data saved in this r2c file:",r2cFile.name,"\n"
#
