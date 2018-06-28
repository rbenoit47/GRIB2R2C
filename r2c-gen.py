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
                  help='Variable to be extracted. Typically either PCP or TGL.')
args.add_argument('-indir',dest='indir',required=True,\
                  help='Directory path containing the input GRIB2 files')
args.add_argument('-interval',dest='interval',required=False,\
                  default=3,\
                  help='Time interval (hours) between desired samples to be saved. Default is 3')
args.add_argument('-verbose',dest='verbose',required=False,\
                  default=[],\
                  help='verbosity.  True or False. Default is False ([])')
args.add_argument('-crop',dest='crop',required=False,\
                  default=[],\
                  help='To extract a SUBSET of the grib array.  Default is False ([])')
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
args.add_argument('-LALO',dest='LALO',required=False,\
						default=False,\
                  help='Option to get lat-lon of the grib array points. Default is False ([])')
#  special case: LALO replaces record values by pairs of latitude longitude
#
args.parse_args(namespace=args)
#
var=args.var  # PCP or TGL
interval=int(args.interval)  # number of hours between desired samples to be saved
indir=args.indir
verbose=args.verbose
r2cprefix=args.r2cprefix
outdir=args.outdir
crop=structObj() #{}
crop.crop=bool(args.crop)
crop.lat=float(args.lat)
crop.lon=float(args.lon)
crop.width=float(args.width)
LALO=bool(args.LALO)
#
print "script arguments. \nvar:",var,"\ninterval:", interval, "\ndir:",indir,"\nverbose:",verbose,"\ncrop:",structPrint(crop,'crop')
"""
print type(crop.crop)
if crop.crop :
	print "crop is true"
else:
	print "crop is false"
quit()
"""
#
gribFiles=sorted(glob.glob(indir+"*"+var+"*grib2"))
#
#
if len(gribFiles) == 0:
	print "no GRIB2 files found. quit"
	quit()
#
#if not LALO:
r2cname=outdir + "/" + r2cprefix + "_" + var + ".r2c"
#else:
#	r2cname=outdir + "/" + r2cprefix + "_" + "LALO" + ".r2c"
#
r2cFile=open(r2cname,'w')
r2cHeader=True
iFile=1
pickSave=False  #True
#
for gribFile in gribFiles:
	grib,saveit=get_grib (gribFile,interval,crop,LALO=LALO,verbose=verbose)
	#
	if saveit:
		print ("data obtained from file %s" % grib.GribName)
		if verbose:
			print ("===================")
			structPrint(grib,'grib')
			structPrint(grib.geo,'geo')
			print ("===================")
		#
		if pickSave:
			pikfile=var+str(iFile)+".pickle"
			pikf=open(pikfile,'wb')
			pik.dump(grib['values'],pikf)
			pikf.close()
		#
		ok=put_r2c(grib,r2cFile,'R Benoit',FrameNumber=iFile,doHeader=r2cHeader,verbose=verbose)
		grib.var=var
		#grib_plot_it(grib,trace=1.0)
		r2cHeader=False
		iFile+=1
		if not ok:
			print "problem in producing r2c file:",r2cFileName
#
r2cFile.close()
print "\nExtracted data saved in this r2c file:",r2cFile.name,"\n"
#
if LALO :
	for LL in ['LATS','LONS']:
		#  output LATS and LONS to 2 other r2c files
		r2cname=outdir + "/" + r2cprefix + "_" + LL + ".r2c"  #"LATS"
		#
		r2cFile=open(r2cname,'w')
		if LL == 'LATS':
			grib.values=np.transpose(grib.geo.lats.copy(order='F'))
			LLtoken='Latitudes'
		elif LL == 'LONS':
			grib.values=np.transpose(grib.geo.lons.copy(order='F'))
			LLtoken='Longitudes'
		#
		grib.parameterUnits=grib.geo.latlonUnits
		ok=put_r2c(grib,r2cFile,'R Benoit',verbose=verbose)
		r2cFile.close()
		print "\n",LLtoken,"saved in this r2c file:",r2cFile.name,"\n"
#
