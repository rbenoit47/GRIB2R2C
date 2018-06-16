# -*- coding: cp1252 -*-
# above coding allows the copyright character used below
#
from eccodes import *
import numpy as np
import datetime as dtime
import copy
from llxy import *
#
def get_grib (gribname,interval,cropSpecs,LALO=False,verbose=False):
	#
	saveit=False
	#
	GribFile= open(gribname)
	gid = codes_grib_new_from_file(GribFile)
	if verbose : print ("opening grib file %s" % gribname)
	#  the output dictionary
	#gdict={}
	#  change to a struct
	grib=structObj()
	#
	grib.GribName=gribname    #gdict['GribName'] = gribname
	grib.Ni=codes_get(gid,'Nx')
	grib.Nj=codes_get(gid,'Ny')
	grib.parameterName=codes_get(gid,'parameterName')
	if  'Total precipitation' in grib.parameterName:
		grib.lengthOfTimeRange=codes_get(gid,'lengthOfTimeRange')
		if verbose : print 'APCP.. length of accumulation:',grib.lengthOfTimeRange
	else:
		grib.lengthOfTimeRange=0
	if  'Total precipitation' in grib.parameterName:
		grib.forecastTime=grib.lengthOfTimeRange #special case of cumulative precip
	else:
		grib.forecastTime=codes_get(gid,'forecastTime')
	# check for interval condition
	t=int(grib.forecastTime)
	my_mod=t-(t/interval)*interval
	if not my_mod == 0:
		saveit=False
		return grib, saveit
	else:
		saveit=True
	#
	#  put all the geo metadata in a structure
	geo=structObj()
	#
	grib.parameterUnits=codes_get(gid,'parameterUnits')
	#
	geo.grtyp=codes_get(gid,'gridType')
	geo.LaD=codes_get(gid,'LaDInDegrees')
	geo.orientation=codes_get(gid,'orientationOfTheGridInDegrees')
	geo.Dx=codes_get(gid,'DxInMetres')
	geo.Dy=codes_get(gid,'DyInMetres')
	#
	if geo.Dx == geo.Dy and geo.Dx == 10000.0 and geo.LaD == 60.0 :
		# values are defined by the CMC forecast centre
		"""
		NOTE:
		latitudeOfFirstGridPoint longitudeOfFirstGridPoint are in the grib file header
		in millionth of degree
		
		shapeOfTheEarth = 6 
		[Earth assumed spherical with radius of 6,371,229.0 m (grib2/tables/4/3.2.table) ]
		"""
		geo.radius=6371229.0
		lat1=np.float64(codes_get(gid,'latitudeOfFirstGridPoint') )/1000000
		lon1=np.float64(codes_get(gid,'longitudeOfFirstGridPoint') )/1000000
		geo.lat1=lat1
		geo.lon1=lon1
		d60=geo.Dx
		nhem=1
		CentreLongitude = geo.orientation
		dgrw=- CentreLongitude - 90.
		[x1,y1]=xyfll(lat1,lon1,d60,dgrw,nhem,lat0=60.0,verbose=False,radius=6371229.0)
		geo.iPole= -x1 + 1
		geo.jPole= -y1 + 1
		#  get x,yOrigin for the r2c header
		geo.xOrigin=x1*d60-d60/2.  #for the cell aspect we add the d60/2
		geo.yOrigin=y1*d60-d60/2.  #(y1+grib.Nj-1)*d60+d60/2.  # top left corner, contrary to ensim doc that says bottom left
		if verbose : print "CMC.  i,j North Pole,x/yOrigin:",geo.iPole,geo.jPole,geo.xOrigin,geo.yOrigin
	else:
		print "get_grib.  Dx must be 10000. and LaDInDegrees=60.0", geo.Dx,geo.LaD," breaking"
		quit()
	#  get the data date  (formatted as YYYYMMDD) and data time as well as forecast time (both in hours)
	grib.dataDate=codes_get(gid,'dataDate')
	grib.dataTime=codes_get(gid,'dataTime')
	#
	stepUnits=codes_get(gid,'stepUnits')  # a value of 1 indicates hours as time unit.  must be
	if stepUnits <> 1:
		print "stepUnits:",stepUnits," is not 1.  breaking"
		quit()
	#
	#
	if not LALO  :
		# make it np array rightaway, with Fortran ordering. and go to 2D also
		buf=np.reshape(np.array(codes_get_values(gid),order='F'),(grib.Nj,grib.Ni))
	else:
		print "getting lat-lon coordinates for grib array geography"
		latarr,lonarr=getlalo(grib,geo)
		if verbose:
			print "LALO shapes",latarr.shape,lonarr.shape
			print "grib Ni Nj:",grib.Ni,grib.Nj
			print "LALO. min,max lat and lo:\n LA:",np.min(latarr),np.max(latarr),"\n LO:",np.min(lonarr),np.max(lonarr)
			print "LALO. LL and TR corners:\n LL ",latarr[0,0],lonarr[0,0],"\n TR ",latarr[grib.Ni-1,grib.Nj-1],lonarr[grib.Ni-1,grib.Nj-1]
			structPrint(geo,'geo')
		buf=latarr.copy(order='F')
		grib.parameterName='Latitudes'
		grib.parameterUnits='degrees'
	#
	if verbose : print buf.size,buf.shape,grib.Ni,grib.Nj,grib.Ni*grib.Nj,grib.dataDate,grib.dataTime,grib.forecastTime
	#  cropping ...
	if cropSpecs.crop :
		print "get_grib.  cropping...",structPrint(cropSpecs,'cropSpecs')
		#
		latc=cropSpecs.lat
		lonc=cropSpecs.lon
		width=cropSpecs.width
		[xc,yc]=xyfll(latc,lonc,d60,dgrw,nhem,lat0=60.0,verbose=verbose,radius=geo.radius)  #coordinates with origin at Pole
		print "xc,yc:",xc,yc
		[dum1,dum2]=llfxy(xc,yc,d60,dgrw,nhem,lat0=60.0,verbose=verbose,radius=geo.radius)
		ic=np.int(np.round(xc+geo.iPole))   #coordinates with origin at array low-left corner
		jc=np.int(np.round(yc+geo.jPole))
		iwidth=np.int(width*1000.0/d60)
		i1=ic-int(iwidth/2)
		j1=jc-int(iwidth/2)
		#
		x1crop=i1-geo.iPole   #coordinates with origin at Pole
		y1crop=j1-geo.jPole
		#
		if verbose : print "latc,lonc,dum1,dum2,iwidth,ic,jc,i1,i2,j1,j2,x1crop,y1crop:\n",latc,lonc,dum1,dum2,iwidth,ic,jc,i1,j1,x1crop,y1crop
		if verbose : print "lat lon of i1,j1 point:",llfxy(x1crop,y1crop,d60,dgrw,nhem,lat0=60.0,verbose=False,radius=6371229.0)
		#
		#cbuf,geocrop=crop_array(buf,i1,j1,iwidth,iwidth,yflip=True,verbose=True,geoIn=geo,cropSpecs=cropSpecs)
		
		#  crop(a,x1,y1,nx,ny)
		cbuf=crop(buf,i1-1,j1-1,iwidth,iwidth)
		
		#
		#  must adjust these in grib struct
		#  Ni Nj iPole jPole xOrigin yOrigin
		#
		grib.Nj,grib.Ni=cbuf.shape
		geo.iPole=geo.iPole-(i1-geo.iPole-1)
		geo.jPole=geo.jPole-(j1-geo.jPole-1)
		# new values for x1 y1
		x1=1-geo.iPole
		y1=1-geo.iPole
		geo.xOrigin=x1crop*d60  #x1*d60-d60/2.              # for the cell aspect we add the d60/2
		geo.yOrigin=y1crop*d60   #(y1+grib.Nj-1)*d60+d60/2.  # top left corner, contrary to ensim doc that says bottom left
		#
		if verbose : 
			print "cropSpecs to Ni,Nj,i/j NorthPole,x/yOrigin:",grib.Ni,grib.Nj,"\n geo for cropped array"
			structPrint(geo,'geo')
		#
		grib.values=cbuf
		del buf, cbuf
	else:
		grib.values=buf	
		del buf
	#
	if verbose : 
		for key in ('max', 'min', 'average'):
			print('%s=%.10e' % (key, codes_get(gid, key)))
	codes_release(gid)
	GribFile.close()
	#  include geo in grib structure
	grib.geo=geo
	#
	return grib, saveit
#
#
def crop(a,x1,y1,nx,ny):
    # infos:
    # https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.isfortran.html
    if a.flags.f_contiguous:
        print "crop. input is order='F'"
        return a[x1:x1+nx,y1:y1+ny]
    elif a.flags.c_contiguous:
        #  The transpose of a C-ordered array is a FORTRAN-ordered array
        print "crop. input is order='C'"
        return a.transpose()[x1:x1+nx,y1:y1+ny]
    else:
        print "input array is neither C nor Fortran order. breaking"
        quit()
#
def structObj():
	# to define a matlab-type dotted structure
	# cf.: https://stackoverflow.com/questions/34075094/python-struct-like-matlab
	return type('', (), {})()
#
def structCopy(a):
	# to copy a struct object (a) to a new one (b)
	return copy.copy(a)
#
def structPrint(a,namea):
	adict=a.__dict__
	print "\n"
	for key in adict:
		print namea+'.'+str(key)+':  ',adict[key]
	return	
def put_r2c(grib,r2cFile,writername,FrameNumber=1,doHeader=True,verbose=False):
	ok=False
	#
	if doHeader :
		ok=r2c_header(grib,r2cFile,writername,verbose)
		if not ok:
			print "problem with r2c header.  not continuing"
			return ok
	#  repeat for each Frame
	StepNumber=grib.forecastTime
	#FrameTime="2018 05 10 00:00"
	#  do the date and time for real
	dataDatef=str(grib.dataDate)
	dataTimef=str(grib.dataTime)
	dataDateTimef=dataDatef+dataTimef
	startDateTime=dtime.datetime.strptime(dataDateTimef,'%Y%m%d%H')
	forecastDateTime=startDateTime+dtime.timedelta(hours=grib.forecastTime)
	FrameTime=forecastDateTime.strftime("%Y/%m/%d %H:%M")  #add the slashes as per N Kouwen
	#
	r2cFile.write(":Frame %d %d \"%s\"\n" % (FrameNumber, StepNumber, FrameTime))
	#
	# for outputting the grib.values we'll use numpy.savetxt method
	#  rather than the .write method
	buf=np.copy(grib.values,order='F')
	#[Ni,Nj]=np.shape(buf)
	if verbose : print "dbg. Ni, Nj of grib data array, ordering: ",np.shape(buf), buf.flags.f_contiguous
	nifmt=['%-7.2f'] * grib.Nj
	#
	# data cells ought to be written to the r2c file in the Fortran order
	# i.e. for an 8x5 spatial array:
	"""
	from cell 0 to cell 39  (=8x5 values)
	.......37 38 39
	...............
	...............
	8 9 ...........
	0 1 2 3 4 5 6 7
	"""
	#
	np.savetxt(r2cFile,buf.transpose() , fmt=nifmt)
	#
	r2cFile.write(":EndFrame\n")
	#
	if verbose : print "r2c file:",r2cFile.name," Frame produced" 
	ok=True
	return ok
#
def r2c_header(grib,r2cFile,writername,verbose=False):
	ok=False
	WriterName=writername
	SourceFile=grib.GribName
	Projection='PolarStereographic'
	CentreLatitude=grib.geo.LaD
	CentreLongitude=grib.geo.orientation
	Ellipsoid='Sphere'
	AttributeUnits=grib.parameterUnits
	xCount=grib.Ni
	yCount=grib.Nj
	xDelta=grib.geo.Dx
	yDelta=grib.geo.Dy
	#
	#xOrigin=-4556441.403315   #need to be dynamically specified
	#yOrigin=  920682.141166   #need to be dynamically specified
	#  grib gives us the i,j of North Pole
	#  compute origin from that.  i.e. this is a shift
	xOrigin=grib.geo.xOrigin     #0.0-(grib.geo.iPole-0.5)*xDelta  # 0.5 due to cell geometry rather than corner
	#yOrigin=0.0-(grib.jNorthPole']-0.5)*yDelta  # 0.5 due to cell geometry rather than corner
	# top left corner, contrary to ensim doc that says bottom left
	#y1=1-grib.geo.jPole
	yOrigin=grib.geo.yOrigin    #(y1+yCount-1+0.5)*yDelta
	
	d60=1.0
	nhem=1
	CentreLongitude = grib.geo.orientation
	dgrw=- CentreLongitude - 90.
	print "lat lon of x/yOrigin:",llfxy(xOrigin,yOrigin,d60,dgrw,nhem,lat0=60.0,verbose=False,radius=6371229.0)
	
	
	
	
	#
	r2cFile.write("############################################################################\n")
	r2cFile.write(":FileType r2c  ASCII  EnSim 1.0\n")
	r2cFile.write("# File structure National Research Council Canada © 1998-2017\n")
	r2cFile.write("############################################################################\n")
	r2cFile.write("# DataType    2D Rect Cell\n")
	r2cFile.write("#\n")
	r2cFile.write(":Application  GreenKenue\n")
	r2cFile.write(":Version      3.8.1\n")
	r2cFile.write(":WrittenBy    %s\n" % (WriterName))
	now=dtime.datetime.now()
	nowfmted=now.strftime("%Y/%m/%d %H:%M:%S.000")
	r2cFile.write(":CreationTime    " + nowfmted + "\n")
	r2cFile.write("#---------------------------------------------------------------------------\n")
	r2cFile.write("#\n")
	r2cFile.write(":SourceFile   %s\n" % (SourceFile))
	r2cFile.write("#\n")
	r2cFile.write(":Projection %s\n" % (Projection))
	r2cFile.write(":CentreLatitude %d\n" % (CentreLatitude))
	r2cFile.write(":CentreLongitude %d\n" % (CentreLongitude))
	r2cFile.write(":Ellipsoid %s\n" % (Ellipsoid))
	r2cFile.write("#\n")
	r2cFile.write(":xOrigin %f\n" % (xOrigin))
	r2cFile.write(":yOrigin %f\n" % (yOrigin))
	r2cFile.write("#\n")
	r2cFile.write("#\n")
	r2cFile.write(":AttributeUnits %s\n" % (AttributeUnits))
	r2cFile.write(":xCount %d\n" % (xCount))
	r2cFile.write(":yCount %d\n" % (yCount))
	r2cFile.write(":xDelta %f\n" % (xDelta))
	r2cFile.write(":yDelta %f\n" % (yDelta))
	r2cFile.write("#\n")
	r2cFile.write(":EndHeader\n")
	if verbose : print "r2c file:",r2cFile.name," Header produced" 
	ok=True
	return ok

