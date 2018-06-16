#!/usr/bin/env python
from   grib2r2c import *
import numpy as np
#
radius=6371229.0
a=np.float64(radius)
dgtord=np.pi/np.float64(180.0)
#
grib=structObj()
grib.Ni=5
grib.Nj=5
#
geo=structObj()
geo.Dx=1.0
geo.LaD=60.0
geo.orientation=-90.
geo.iPole=3
geo.jPole=3
#
re=(1+np.sin(dgtord*geo.LaD))*a/geo.Dx
#
geo.xOrigin=-re
geo.Dx=re/(geo.iPole-1)  #2
geo.yOrigin=-geo.Dx*(geo.jPole-1)
#
structPrint(grib,'grib')
structPrint(geo,'geo')
#
lats,lons=getlalo(grib,geo,radius=a)
#
print lats.shape
print lons.shape
#
print "lats\n",np.flipud(np.transpose(lats))
print "lons\n",np.flipud(np.transpose(lons))
