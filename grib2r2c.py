# -*- coding: cp1252 -*-
# above coding allows the copyright character used below
#
from eccodes import *
import numpy as np
import datetime as dtime
from llxy import *
#
def get_grib (gribname,interval,crop,verbose=False):
	#
	saveit=False
	#
	GribFile= open(gribname)
	gid = codes_grib_new_from_file(GribFile)
	if verbose : print ("opening grib file %s" % gribname)
	#  the output dictionary
	gdict={}
	gdict['GribName'] = gribname
	gdict['Ni']=codes_get(gid,'Nx')
	gdict['Nj']=codes_get(gid,'Ny')
	gdict['parameterName']=codes_get(gid,'parameterName')
	#gdict['shortName']=codes_get(gid,'shortName')
	if  'Total precipitation' in gdict['parameterName']:
		gdict['lengthOfTimeRange']=codes_get(gid,'lengthOfTimeRange')
		if verbose : print 'APCP.. length of accumulation:',gdict['lengthOfTimeRange']
	else:
		gdict['lengthOfTimeRange']=0
	if  'Total precipitation' in gdict['parameterName']:
		gdict['forecastTime']=gdict['lengthOfTimeRange'] #special case of cumulative precip
	else:
		gdict['forecastTime']=codes_get(gid,'forecastTime')
	# check for interval condition
	t=int(gdict['forecastTime'])
	my_mod=t-(t/interval)*interval
	if not my_mod == 0:
		saveit=False
		return [gdict, saveit]
	else:
		saveit=True
	#
	gdict['parameterUnits']=codes_get(gid,'parameterUnits')
	gdict['gridType']=codes_get(gid,'gridType')
	gdict['LaDInDegrees']=codes_get(gid,'LaDInDegrees')
	gdict['orientationOfTheGridInDegrees']=codes_get(gid,'orientationOfTheGridInDegrees')
	gdict['DxInMetres']=codes_get(gid,'DxInMetres')
	gdict['DyInMetres']=codes_get(gid,'DyInMetres')
	if gdict['DxInMetres'] == gdict['DyInMetres'] and gdict['DxInMetres'] == 10000.0 and gdict['LaDInDegrees'] == 60.0 :
		# values are defined by the CMC forecast centre
		"""
		NOTE:
		latitudeOfFirstGridPoint longitudeOfFirstGridPoint are in the grib file header
		in millionth of degree
		
		shapeOfTheEarth = 6 
		[Earth assumed spherical with radius of 6,371,229.0 m (grib2/tables/4/3.2.table) ]
		"""
		lat1=np.float64(codes_get(gid,'latitudeOfFirstGridPoint') )/1000000
		lon1=np.float64(codes_get(gid,'longitudeOfFirstGridPoint') )/1000000
		d60=gdict['DxInMetres']
		nhem=1
		CentreLongitude = gdict['orientationOfTheGridInDegrees']
		dgrw=- CentreLongitude - 90.
		[x1,y1]=xyfll(lat1,lon1,d60,dgrw,nhem,lat0=60.0,verbose=verbose,radius=6371229.0)
		gdict['iNorthPole']= -x1 + 1
		gdict['jNorthPole']= -y1 + 1
		if verbose : print "CMC.  i,j North Pole:",gdict['iNorthPole'],gdict['jNorthPole']
	else:
		gdict['iNorthPole']=float('nan')
		gdict['jNorthPole']=float('nan')
		print "get_grib.  Dx must be 10000. and LaDInDegrees=60.0", gdict['DxInMetres'],gdict['LaDInDegrees']," breaking"
		quit()
	#  get the data date  (formatted as YYYYMMDD) and data time as well as forecast time (both in hours)
	gdict['dataDate']=codes_get(gid,'dataDate')
	gdict['dataTime']=codes_get(gid,'dataTime')
	#
	stepUnits=codes_get(gid,'stepUnits')  # a value of 1 indicates hours as time unit.  must be
	if stepUnits <> 1:
		print "stepUnits:",stepUnits," is not 1.  breaking"
		quit()
	#
	#buf=codes_get_values(gid)
	buf=np.reshape(np.array(codes_get_values(gid),order='F'),(gdict['Nj'],gdict['Ni']))  #make it np array rightaway
	# and go to 2D also
	#
	if verbose : print buf.size,buf.shape,gdict['Ni'],gdict['Nj'],gdict['Ni']*gdict['Nj'],gdict['dataDate'],gdict['dataTime'],gdict['forecastTime']
	#  cropping ...
	if crop['crop'] :
		print "get_grib.  cropping...",crop
		#
		latc=crop['lat']
		lonc=crop['lon']
		width=crop['width']
		[xc,yc]=xyfll(latc,lonc,d60,dgrw,nhem,lat0=60.0,verbose=verbose,radius=6371229.0)
		print "xc,yc:",xc,yc
		[dum1,dum2]=llfxy(xc,yc,d60,dgrw,nhem,lat0=60.0,verbose=verbose,radius=6371229.0)
		ic=np.int(np.round(xc+gdict['iNorthPole']))
		jc=np.int(np.round(yc+gdict['jNorthPole']))
		iwidth=np.int(width*1000.0/d60)
		i1=ic-int(iwidth/2)
		j1=jc-int(iwidth/2)
		print "i1,j1,iwidth:",i1,j1,iwidth
		if verbose : print "latc,lonc,dum1,dum2,iwidth,ic,jc,i1,i2,j1,j2:\n",latc,lonc,dum1,dum2,iwidth,ic,jc,i1,j1
		#
		cbuf=crop_array(buf,j1,i1,iwidth,iwidth,yup=True,verbose=True)
		#
		#  must adjust these in gdict
		# Ni Nj iNorthPole jNorthPole 
		gdict['Nj'],gdict['Ni']=cbuf.shape
		gdict['iNorthPole']=gdict['iNorthPole']-i1+1
		gdict['jNorthPole']=gdict['jNorthPole']-j1+1
		if verbose : print "crop to Ni,Nj,iNorthPole,jNorthPole:",gdict['Ni'],gdict['Nj'],gdict['iNorthPole'],gdict['jNorthPole']
		#
		gdict['values']=cbuf
		del buf, cbuf
		#=ici
	else:
		gdict['values']=buf	
		#gdict['values'] = np.reshape(np.array(buf),(gdict['Ni'],gdict['Nj']))  #,order='F'
		#gdict['values'] = np.reshape(buf,(gdict['Ni'],gdict['Nj']))  #,order='F'
		del buf
	#
	if verbose : 
		for key in ('max', 'min', 'average'):
			print('%s=%.10e' % (key, codes_get(gid, key)))
	codes_release(gid)
	GribFile.close()
	#
	return [gdict, saveit]
#
def crop_array(img,startx,starty_in,cropx,cropy,yup=False,verbose=False):
    ny,nx = img.shape
    if verbose:
        print "min,max,shape input:",np.min(img),np.max(img),"\n",img.shape
    #
    if not yup:
        starty=starty_in
        endy=min(starty+cropy,ny)
    else:
        starty=max(0,ny-1-(starty_in+cropy-1))
        endy=min(ny-1-(starty_in+cropy-1)+cropy,ny)
    if startx < 0 or starty < 0:
        print "starting point must not be negative. exit", startx, starty
        return
    if endy > ny:
        print "ending y point must not exceed y-size. exit", endy, ny
        return
    crop=img[starty:endy,startx:min(startx+cropx,nx)]
    if verbose:
        print "min,max,shape output:",np.min(crop),np.max(crop),"\n",crop.shape   
    return crop
#
def put_r2c(grib,r2cFile,writername,FrameNumber=1,doHeader=True,verbose=False):
	ok=False
	#
	if doHeader :
		ok=r2c_header(grib,r2cFile,writername,verbose)
		if not ok:
			print "problem with r2c header.  not continuing"
			return ok
	#  repeat for each Frame
	StepNumber=grib['forecastTime']
	#FrameTime="2018 05 10 00:00"
	#  do the date and time for real
	dataDatef=str(grib['dataDate'])
	dataTimef=str(grib['dataTime'])
	dataDateTimef=dataDatef+dataTimef
	startDateTime=dtime.datetime.strptime(dataDateTimef,'%Y%m%d%H')
	forecastDateTime=startDateTime+dtime.timedelta(hours=grib['forecastTime'])
	FrameTime=forecastDateTime.strftime("%Y/%m/%d %H:%M")  #add the slashes as per N Kouwen
	#
	r2cFile.write(":Frame %d %d \"%s\"\n" % (FrameNumber, StepNumber, FrameTime))
	#
	# for outputting the grib['values'] we'll use numpy.savetxt method
	#  rather than the .write method
	buf=grib['values']
	#[Ni,Nj]=np.shape(buf)
	if verbose : print "dbg. Ni, Nj of grib data array: ",np.shape(buf)
	nifmt=['%-7.2f'] * grib['Nj']
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
	np.savetxt(r2cFile,buf , fmt=nifmt)
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
	SourceFile=grib['GribName']
	Projection='PolarStereographic'
	CentreLatitude=grib['LaDInDegrees']
	CentreLongitude=grib['orientationOfTheGridInDegrees']
	Ellipsoid='Sphere'
	AttributeUnits=grib['parameterUnits']
	xCount=grib['Ni']
	yCount=grib['Nj']
	xDelta=grib['DxInMetres']
	yDelta=grib['DyInMetres']
	#
	#xOrigin=-4556441.403315   #need to be dynamically specified
	#yOrigin=  920682.141166   #need to be dynamically specified
	#  grib gives us the i,j of North Pole
	#  compute origin from that.  i.e. this is a shift
	xOrigin=0.0-(grib['iNorthPole']-0.5)*xDelta  # 0.5 due to cell geometry rather than corner
	#yOrigin=0.0-(grib['jNorthPole']-0.5)*yDelta  # 0.5 due to cell geometry rather than corner
	# top left corner, contrary to ensim doc that says bottom left
	y1=1-grib['jNorthPole']
	yOrigin=(y1+yCount-1+0.5)*yDelta
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

