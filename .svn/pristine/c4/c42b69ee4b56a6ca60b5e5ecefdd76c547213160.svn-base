#!/usr/bin/env python

import numpy as nm

class BoundElem:
	'''Create boundary elements on various domains. Unit tested on 05/06/2010 using circle with 4, 8, and 256 nodes with 10 points per chord. Visualized in Matlab.'''
	def __init__(self):
		pass
		
	def LinearElements2D(self):
		'''The *closed* boundary curve in a 2D domain will be divided into linear chords discretized with M points apiece.'''
		zx =[]
		zy=[]
		M = self.M
		y=self.nodes
		cspacing = 1./M
		t = nm.arange(0,1,cspacing)
		for k in range(y.shape[0]):
			if k == y.shape[0]-1:
				pos = y[0,:]
			else:
				pos = y[k+1,:]
			xcoord=y[k,0] + t*(pos[0]-y[k,0])
			ycoord=y[k,1] + t*(pos[1]-y[k,1])
			zx.append(xcoord)
			zy.append(ycoord)
		zx=nm.asarray(zx).flatten()
		zy=nm.asarray(zy).flatten()
		self.quadpts = nm.column_stack([zx,zy])
		self.quadspacing = self.nodespacing*cspacing
		
	
		
class Circle(BoundElem):
	'''Create boundary elements on a circle.'''
	def __init__(self,numNodes,numPtsPerElem,circleRadius=1):
		BoundElem.__init__(self)
		self.N=numNodes
		self.M = numPtsPerElem
		self.R=circleRadius
		angles = nm.arange(0,2*nm.pi,2*nm.pi/self.N)
		self.nodes=nm.column_stack([self.R*nm.cos(angles),self.R*nm.sin(angles)])
		self.nodespacing=nm.sqrt(nm.sum( (self.nodes[0,:]-self.nodes[1,:])**2 ))
		
