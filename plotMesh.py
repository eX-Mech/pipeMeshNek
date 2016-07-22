#!/usr/bin/env python
#
# load and plot 2D mesh created by pipeMeshNek
# jcanton@mech.kth.se

#------------------------------------------------------------------------------
# import pymech
#
import os
try:
	assert os.path.exists('./pymech')
except AssertionError:
	import subprocess
	subprocess.call('git clone https://github.com/jcanton/pymech.git', shell=True)

import sys
sys.path.append('./pymech/src/')
import neksuite as ns
import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------
# plotting function
#
def arc_2points(x, y, r):
	#
	# plot arc given 2 points and radius
	#
	import numpy as np
	#
	x1 = x[0]; x2 = x[1]
	y1 = y[0]; y2 = y[1]
	#
	alpha = (x1 - x2)/(y2 - y1)
	beta  = (x2**2 - x1**2 + y2**2 - y1**2)/(2*(y2 - y1))
	a = 1 + alpha**2
	b = -2*x1 - 2*y1*alpha + 2*alpha*beta
	c = x1**2 + y1**2 - 2*y1*beta + beta**2 - r**2
	#
	delta = b**2 - 4*a*c
	#
	xc1 = (-b + np.sqrt(delta))/(2*a)
	xc2 = (-b - np.sqrt(delta))/(2*a)
	#
	yc1 = alpha*xc1 + beta
	yc2 = alpha*xc2 + beta
	#
	d1 = np.sqrt(xc1**2 + yc1**2)
	d2 = np.sqrt(xc2**2 + yc2**2)
	#
	if ( d1<d2 ):
		xc = xc1; yc = yc1
	else:
		xc = xc2; yc = yc2
	#
	tt1 = np.arctan2( (y1-yc), (x1-xc) )
	tt2 = np.arctan2( (y2-yc), (x2-xc) )
	start_t = min(tt1,tt2)
	end_t   = max(tt1,tt2)
	#
	tt = np.linspace(start_t, end_t, 100)
	xx = xc + np.abs(r)*np.cos(tt)
	yy = yc + np.abs(r)*np.sin(tt)
	#
	plt.plot(xx, yy, '-k')

#------------------------------------------------------------------------------
# load mesh
#
fname = 'base2d.rea'
field = ns.readrea(fname)

#------------------------------------------------------------------------------
# plot
#
plt.figure(1)
plt.clf()

nedges = 4
for iel in range(field.nel):
	xv = np.reshape(field.elem[iel].pos[0, 0,:,:], (4,1) )
	yv = np.reshape(field.elem[iel].pos[1, 0,:,:], (4,1) )
	xl = (np.min(xv)+np.max(xv))/2
	yl = (np.min(yv)+np.max(yv))/2
	plt.text(xl, yl, '%d' % (iel+1))
	for iedge in range(nedges):
		xe = np.roll(xv, -iedge)[0:2]
		ye = np.roll(yv, -iedge)[0:2]
		if (field.elem[iel].curv[iedge] == 0):
			plt.plot(xe, ye, '-k')
		else:
			arc_2points(xe, ye, field.elem[iel].curv[iedge])

plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title(r'2D section of the mesh')
plt.grid(True)
plt.axis('equal')
plt.draw()
plt.show(block=True)
