#!/usr/bin/env python
# -*- coding: cp1252 -*-
#
#  main
#
from llxy import *
from grib2r2c import crop_array
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#
dgrw=21. #45. #249.0   21 is ok for cmc 10 km grid
nhem=1
d60=1.0  #10000.
#
#  cmc grid check
#
ipole=456.2
jpole=732.4
#  thus  from X=Xpole+(i-ipole)*d60
I1=0.+(1-ipole)
J1=0.+(1-jpole)
#
d60=10000.   #units for X and Y
X1=I1
Y1=J1
#
[lat1,lon1]=llfxy(X1,Y1,d60,dgrw,nhem)
print "point 1,1 check."
print "X,Y,lat,lon:",X1,Y1,lat1,lon1
#
print "xyfll point 1,1",xyfll(18.1429,-142.8968,d60,dgrw,nhem)
#
# get xOrigin and yOrigin for r2c header
#
ni=935
nj=824
#
xOrigin=X1*d60-d60/2.  #for the cell aspect we add the d60/2
yOrigin=(Y1+nj-1)*d60+d60/2.  # top left corner, contrary to ensim doc that says bottom left
#
print "x,yOrigin:",xOrigin,yOrigin," meters"
#
#   add a crop test to debug navigation of the subarray
#
fullX=np.arange(X1,X1+ni,1)
fullY=np.arange(Y1,Y1+nj,1)   #Y is upward
fulla=np.empty((nj,ni))
print "full X,Y,a shapes:",fullX.shape,fullY.shape,fulla.shape
# fill a
for j in range(nj):
	for i in range(ni):
		fulla[j,i]=fullX[i]-fullY[j]
#
cs1=plt.contour(fullX,fullY,fulla,colors='k')
plt.clabel(cs1)
del cs1
#  draw grid axes
plt.hold(True)
plt.plot(fullX,np.zeros(fullX.size),'r')
plt.plot(np.zeros(fullY.size),fullY,'b')
#
i1=250
j1=200
iwidth=350
cropa=crop_array(fulla,j1,i1,iwidth,iwidth,yup=False,verbose=True)
cropnj,cropni=cropa.shape
cropX1=X1+(j1-1)
cropY1=Y1+(i1-1)
cropX=np.arange(cropX1,cropX1+cropni,1)
cropY=np.arange(cropY1,cropY1+cropnj,1)
# overlay the subarray
cs2=plt.contour(cropX,cropY,cropa,colors='g')
plt.clabel(cs2)
plt.plot(cropX1,cropY1,'k+')
plt.plot(cropX1+cropni,cropY1+cropnj,'k+')
plt.show()
# do it also with geography
