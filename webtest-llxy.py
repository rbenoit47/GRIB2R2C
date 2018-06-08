#!/usr/bin/env python
# -*- coding: cp1252 -*-
#
#  main
#
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from llxy import *
#
dgrw=21. #45. #249.0   21 is ok for cmc 10 km grid
nhem=1
d60=1.0  #10000.
# -111.
m=Basemap(projection='stere',lat_0=90.,lon_0=-dgrw-90.,width=40000000.,height=40000000 ,lat_ts=60.,rsphere=6371000.)   #-111.
#  coordinates of the pole
mPoleXY=np.array(m(dgrw,90.))  #
print "mPoleXY",mPoleXY
#
for k in [1,2]:
	print "k=",k
	# ligne droite lon=constant
	N=10
	plt.hold(True)
	plt.gca().set_aspect('equal')
	plt.axis('equal')
	for lon in [0.,90.,180.,270.,dgrw,-dgrw]:
		for i in range(N):
				lat=0.+(i-1)*90./N
				if k==1: [X,Y]=xyfll(lat,lon,d60,dgrw,nhem)
				if k==2: 
					[X,Y]=m(lon,lat)
					[X,Y]=np.subtract([X,Y],mPoleXY)					
				print i,X,Y,lon,lat
				if lon == 0.: 
					col='k.'
				else:
					col='r.'
				plt.plot(X,Y,col)
	#
	# cercles pour lat=constant
	N=50
	for lat in [0.,20.,40.,60.,80.]:
		for i in range(N):
				lon=0.+(i-1)*360./N
				if k==1: [X,Y]=xyfll(lat,lon,d60,dgrw,nhem)
				if k==2: 
					[X,Y]=m(lon,lat)
					if not (X==1e30 or Y==1e30): [X,Y]=np.subtract([X,Y],mPoleXY)
				print i,X,Y,lon,lat
				if not (X==1e30 or Y==1e30): plt.plot(X,Y,'b.')

	plt.hold(False)
	plt.show()
	print "mPoleXY",mPoleXY

#
print "mPoleXY",mPoleXY
#
#  cmc grid check
#
ipole=456.2
jpole=732.4
d60=10000.
#  thus  from X=Xpole+(i-ipole)*d60
X1=0.+(1-ipole)*d60
Y1=0.+(1-jpole)*d60
#
X1=X1/d60
Y1=Y1/d60
#
[lat1,lon1]=llfxy(X1,Y1,d60,dgrw,nhem)
print "point 1,1 check."
print "X,Y,lat,lon:",X1,Y1,lat1,lon1
#
print "xyfll point 1,1",xyfll(18.1429,-142.8968,d60,dgrw,nhem)
#
