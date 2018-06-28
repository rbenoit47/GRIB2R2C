import numpy as np
#
def llfxy(x,y,d60,dgrw,nhem,lat0=60.0,verbose=False,radius=6.371e6):
	dgtord=np.pi/np.float64(180.0)
	rdtodg=np.float64(180.0)/np.pi
	a=np.float64(radius)   #6.371e6
	re=(1+np.sin(dgtord*lat0))*a/d60
	#
	re2=re**2
	#
	#     if point is at pole set coord to (0.,90.).
	#
	dlat=90.
	dlon=0.
	if not (x == 0. and y == 0.) :
		#
		#     calculate longitude in map coordinates.
		#
		if x == 0.: dlon=np.sign(y)*90.0
		if not x == 0.: dlon=np.arctan(y/x)*rdtodg
		if x  < 0.: dlon=dlon+np.sign(y)*180.
		#
		#     adjust longitude for grid orientation.
		#
		dlon=dlon-dgrw
		if dlon > +180.: dlon=dlon-360.
		if dlon < -180.: dlon=dlon+360.
		#
		#     calculate latitude.
		#
		r2=x**2+y**2
		dlat=(re2-r2)/(re2+r2)
		dlat= np.arcsin(dlat)*rdtodg
		#
		#     change signs if in southern hemisphere.
		#
		if (x == 0. and y == 0.) :
			if nhem == 2 : dlat=-dlat
			if nhem == 2 : dlon=-dlon
	#
	if verbose : print "llfxy. ,x,y,dlat,dlon,d60,dgrw,nhem,lat0:",x,y,dlat,dlon,d60,dgrw,nhem,lat0
	#
	return [dlat,dlon]
#
def llfxya(x,y,d60,dgrw,nhem,lat0=60.0,radius=6.371e6):
	# array version of llfxy.  x and y are array as well returned lat lon
	if nhem == 2:
		print "llfxya.  nhem=2 not supported. breaking"
		quit()
	#
	dgtord=np.pi/np.float64(180.0)
	rdtodg=np.float64(180.0)/np.pi
	a=np.float64(radius)   #6.371e6
	re=(1+np.sin(dgtord*lat0))*a/d60
	#
	re2=re**2
	#
	xv,yv=np.meshgrid(x,y,indexing='ij')  #ij chosen to have x govern 1st dimension of xv yv
	r2=np.square(xv)+np.square(yv)
	ratio=np.divide(re2-r2,re2+r2)
	dlata=np.arcsin(ratio)*rdtodg
	dlona=np.arctan2(yv,xv)*rdtodg - dgrw
	#  wrap lon to +-180
	dlona180=np.fmod(dlona+180.0,360.0)-180.0
	#
	#if verbose : print "llfxy. ,x,y,dlat,dlon,d60,dgrw,nhem,lat0:",x,y,dlat,dlon,d60,dgrw,nhem,lat0
	#
	return dlata.copy(order='F'), dlona180.copy(order='F')  #np.flipud(np.transpose(dlata)),np.flipud(np.transpose(dlona))
#
def xyfll(dlat,dlon,d60,dgrw,nhem,lat0=60.0,verbose=False,radius=6.371e6):
	dgtord=np.pi/np.float64(180.0)
	a=np.float64(radius)  #6.371e6
	re=(1+np.sin(dgtord*lat0))*a/d60
	#
	glon=dlon
	if nhem == 2: glon=-dlon
	glat=dlat
	if nhem == 2: glat=-dlat
	#
	rlon=dgtord*(glon+dgrw)
	rlat=dgtord*glat
	sinlat=np.sin(rlat)
	r=re*np.sqrt((1.-sinlat)/(1.+sinlat))
	x=r*np.cos(rlon)
	y=r*np.sin(rlon)
	#
	if verbose: print "xyfll. x,y,dlat,dlon,d60,dgrw,nhem,lat0:",x,y,dlat,dlon,d60,dgrw,nhem,lat0
	#
	return [x,y]
def getlalo(grib,geo,radius=6.371e6):
	Ni=grib.Ni
	Nj=grib.Nj
	#
	d60=geo.Dx
	lat0=geo.LaD
	dgrw=-geo.orientation-90.0
	iPole=geo.iPole
	jPole=geo.jPole
	nhem=1
	#
	x=np.array(range(0,Ni))*d60+geo.xOrigin   #metric
	y=np.array(range(0,Nj))*d60+geo.yOrigin
	lats,lons=llfxya(x,y,1.0,dgrw,nhem,lat0=lat0,radius=geo.radius)
	xPole,yPole=xyfll(lats[0,0],lons[0,0],d60,dgrw,nhem,lat0=lat0,radius=geo.radius)
	print "getlalo. checking ijPole.",-xPole+1,-yPole+1
	print "getlalo. Print 4 corners\ni j lon lon360 lat\n"
	for i in [0,Ni-1]:
		for j in [0,Nj-1]:
			print i,j,lons[i,j],np.fmod(360.0+lons[i,j],360.),lats[i,j]
	#  check for fORTRAN ORDERING and enforce
	if lats.flags.f_contiguous :
		print "getlalo.  lat lon are Fortran ordered"
		return lats,lons
	else:
		print "getlalo.  lat lon are not Fortran ordered...breaking"
		quit()
