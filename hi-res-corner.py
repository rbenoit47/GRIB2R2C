#!/usr/bin/env python
#
import traceback
import sys
 
from eccodes import *
from llxy import *
from grib2r2c import *
import numpy as np

 
INPUT = 'input.grb' 
VERBOSE = 1  # verbose error reporting

def geoExtra(geo,cellShift=0.5):
	#
	extra=structObj()	
	#
	#  derived parameters for geo
	#
	extra.nhem=1
	extra.dgrw=- geo.orientation - 90.
	#
	[extra.x1,extra.y1]=xyfll(geo.lat1,geo.lon1,geo.Dx,extra.dgrw,extra.nhem,lat0=geo.LaD,verbose=False,radius=6371229.0)
	#
	extra.PI=1-extra.x1  #PI: I index of Pole
	extra.PJ=1-extra.y1  #PJ   idem
	# xOrigin for r2c file.  Consider grib grid as cell CENTERS
	extra.xOriginR2C=(extra.x1 - cellShift) * geo.Dx  #for cell (0.5) or grid (0)
	extra.yOriginR2C=(extra.y1 - cellShift) * geo.Dx
	#
	return extra
#
def lon180(lon):
	return  (lon + 180.0) % 360.0 -180.0
def lon360(lon):
	return  (lon + 360.0) % 360.0 
def example():
    f = open(INPUT)
 
    keys = [
        'Ni',
        'Nj',  
        'parameterName' , 
        'parameterUnits' ,
        'latitudeOfFirstGridPoint'      , 'longitudeOfFirstGridPoint'
    ]
 
    while 1:
        gid = codes_grib_new_from_file(f)
        print "gid",gid
        if gid is None:
			print ' breaking'
			break
        for key in keys:
            try:
                print '  %s: %s' % (key, codes_get(gid, key))
            except CodesInternalError as err:
                print 'Error with key="%s" : %s' % (key, err.msg)
 
	#
	#  put all the geo metadata in a structure
	geo=structObj()
	#
	geo.Nx=codes_get(gid,'Ni')
	geo.Ny=codes_get(gid,'Nj')
	geo.grtyp=codes_get(gid,'gridType')
	geo.LaD=codes_get(gid,'LaDInDegrees')
	geo.orientation=codes_get(gid,'orientationOfTheGridInDegrees')
	geo.radius=6371229.0  #from grib convention for polar stereo projection
	geo.Dx=codes_get(gid,'DxInMetres')
	geo.Dy=codes_get(gid,'DyInMetres')
	#
	geo.lat1=np.float64(codes_get(gid,'latitudeOfFirstGridPoint') )/1000000
	geo.lon1=np.float64(codes_get(gid,'longitudeOfFirstGridPoint') )/1000000
	#
	structPrint(geo,'geo')
	#	
	geo.extra=geoExtra(geo,cellShift=0.0)
	#
	structPrint(geo,'geo')
	#
	# do the 4 corners
	#
	co=0
	coNames=['LL','LR','TL','TR']
	ishift=[0,1,0,1]
	jshift=[0,0,1,1]
	for co in range(4):
		print coNames[co]
		# [extra.x1,extra.y1]=xyfll(geo.lat1,geo.lon1,geo.Dx,extra.dgrw,extra.nhem,
		#                           lat0=geo.LaD,verbose=False,radius=6371229.0)
		x=(geo.extra.xOriginR2C + ishift[co]*(geo.Nx-1)*geo.Dx)/geo.Dx
		y=(geo.extra.yOriginR2C + jshift[co]*(geo.Ny-1)*geo.Dy)/geo.Dy
		[la,lo]=llfxy(x,y,geo.Dx,geo.extra.dgrw,geo.extra.nhem,lat0=geo.LaD,radius=geo.radius)
		print 'LO,LA:',lon360(lo),la
	codes_release(gid)

    f.close()
 
 
def main():
    try:
        example()
    except CodesInternalError as err:
        if VERBOSE:
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')
 
        return 1
 
if __name__ == "__main__":
    sys.exit(main())
