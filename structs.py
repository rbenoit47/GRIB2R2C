#!/usr/bin/env python
import copy
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
	for key in adict:
		print namea+'.'+str(key)+':  ',adict[key]
	return	
#
a=structObj()
#
a.x=1
a.y=1.2
a.z='text'
#
structPrint(a,'a')
#
b=structCopy(a)
b.z='change'
#
structPrint(b,'b')
