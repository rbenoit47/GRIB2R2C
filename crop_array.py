def crop_array(full,startx_in,starty_in,cropx,cropy,yflip=False,verbose=False,fullX=[],fullY=[],geoIn=[],cropSpecs=[]):
	ny,nx = full.shape
	#
	# if a geo struct is passed as geoIn, we rely on it to drive the cropping and return an adjusted geo struct
	#
	geoIsThere=False
	if not geoIn ==[] :
		if cropSpecs == []:
			print "crop_array. Must provide cropSpecs when passing geoIn.  breaking"
			quit()
		geoIsThere=True
		dX=geoIn.Dx
		lat0=geoIn.LaD
		dgrw=-geoIn.orientation-90.0
		latc=cropSpecs.lat
		lonc=cropSpecs.lon
		width=cropSpecs.width*1000/dX
		X1=geoIn.xOrigin/dX
		X=X1+np.arange(nx)
		Y1=geoIn.yOrigin/dX  # dY = dX
		Y=Y1+np.arange(ny)
		[Xc,Yc]=xyfll(latc,lonc,dX,dgrw,1,lat0,radius=6371229.0)
		#Xc+= geoIn.iPole
		#Yc+= geoIn.jPole
		Xstart=Xc-width//2
		Ystart=Yc-width//2
		print "X1,Xc,Xstart Y1,Yc,Ystart:",X1,Xc,Xstart,Y1,Yc,Ystart
		indmin=np.argmin(np.abs(X-Xstart))
		# ==================================
		print "indmin:",indmin
		startx=indmin
		indmin=np.argmin(np.abs(Y-Ystart))
		print "indmin:",indmin
		starty=indmin
		endy=starty+width
		# ==================================
		Xcrop=X[startx:startx+width]
		Ycrop=Y[starty:starty+width]
		[lat1,lon1]=llfxy(Xcrop[0],Ycrop[0],dX,dgrw,1,lat0,radius=6371229.0)
		print "crop LALO 1:",lat1,lon1
		print "crop LALO 2:",llfxy(Xcrop[width-1],Ycrop[width-1],dX,dgrw,1,lat0,radius=6371229.0)
		geoOut=structCopy(geoIn)
		#print "geoIn"
		#print geoIn.__dict__
		structPrint(geoIn,'geoIn')
		#
		geoOut.xOrigin=Xcrop[0]*dX
		geoOut.yOrigin=Ycrop[0]*dX
		geoOut.iPole=geoIn.iPole-(startx-geoIn.iPole-1)
		geoOut.jPole=geoIn.jPole-(starty-geoIn.jPole-1)
		structPrint(geoOut,'geoOut')
		quit()
	else:
		 if not yflip:
		     starty=starty_in
		     endy=min(starty+cropy,ny)
		 else:
		     starty=max(0,ny-1-(starty_in+cropy-1))
		     endy=min(ny-1-(starty_in+cropy-1)+cropy,ny)
	 
	if verbose:
		print "min,max,shape input:",np.min(full),np.max(full),"\n",full.shape
	#
	"""
	if not yflip:
		starty=starty_in
		endy=min(starty+cropy,ny)
	else:
		starty=max(0,ny-1-(starty_in+cropy-1))
		endy=min(ny-1-(starty_in+cropy-1)+cropy,ny)
	"""
	if verbose:
		  print "starty,endy:",starty,endy 	  
	if startx < 0 or starty < 0:
		print "starting point must not be negative. exit", startx, starty
		return
	if endy > ny:
		print "ending y point must not exceed y-size. exit", endy, ny
		return
	crop=full[starty:endy,startx:min(startx+cropx,nx)]
	#
	if verbose:
		print "min,max,shape output:",np.min(crop),np.max(crop),"\n",crop.shape   
	 #
	if len(fullX) * len(fullY) > 0:
		#X and Y are given for full array.  return also cropped X Y
		if fullX.size == nx:
			cropX=fullX[startx:startx+cropx]
		else:
			print "bad fullX provided.  size wrong:",fullX.size,"  quit"
			quit()

		if fullY.size == ny:
			cropY=fullY[starty:starty+cropy]
		else:
			print "bad fullY provided.  size wrong:",fullY.size,"  quit"
			quit()
		return crop,cropX,cropY
		#
	else:
		return crop		

