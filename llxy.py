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

