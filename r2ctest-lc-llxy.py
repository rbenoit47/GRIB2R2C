#!/usr/bin/env python
# -*- coding: cp1252 -*-
#
#  main
#
from llxy import *
from grib2r2c import crop_array
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cPickle as pik
#
dgrw=21. #21 is ok for cmc 10 km grid
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
nc=ni  # nb columns and lines
nl=nj
#
fullX=np.arange(X1,X1+nc,1)
fullY=np.arange(Y1,Y1+nl,1)   #Y is upward
# fill a
#
pikf=open("TMP1.pickle",'rb')  #PCP1
fulla=np.array(pik.load(pikf))
pikf.close()
print "from pickle"
# convert 9999.0 to NaN
imiss=np.argwhere(fulla == 9999.0)
print "len miss:",imiss.shape
fulla[imiss[:,0],imiss[:,1]]=np.nan
#
"""
fulla=np.empty((nl,nc))
for j in range(nl):
	for i in range(nc):
		fulla[j,i]=fullX[i]-fullY[j]
"""
print "full X,Y,a shapes:",fullX.shape,fullY.shape,fulla.shape
#
cvalues=np.arange(250.,295.,5.)   #(0.,170.,10.)
cs1=plt.contour(fullX,fullY,fulla,cvalues,colors='k',linewidths=3)
#plt.clabel(cs1)
del cs1
#  draw grid axes
plt.hold(True)
plt.plot(fullX,np.zeros(fullX.size),'r')
plt.plot(np.zeros(fullY.size),fullY,'b')
#
c1=550   #i1
l1=200   #j1
cwidth=350
lwidth=200
#
cropa,cropX,cropY=crop_array(fulla,c1,l1,cwidth,lwidth,yflip=False,verbose=True,fullX=fullX,fullY=fullY)
cropnl,cropnc=cropa.shape
#cropX1=X1+(c1-0)
#cropY1=Y1+(l1-0)
#cropX=np.arange(cropX1,cropX1+cropnc,1)
#cropY=np.arange(cropY1,cropY1+cropnl,1)
cropX1=cropX[0]
cropY1=cropY[0]
# overlay the subarray
cs2=plt.contour(cropX,cropY,cropa,cvalues,colors='r')   #
#plt.clabel(cs2)
plt.plot(cropX1,cropY1,'r+')
plt.plot(cropX1+cropnc,cropY1+cropnl,'r+')
plt.show()
"""
del cs2
cs3=plt.contour(fulla)
plt.clabel(cs3)
plt.show()
del cs3
cs4=plt.contour(cropa)
plt.clabel(cs4)
plt.show()
"""
# do it also with geography
