#!/usr/bin/env python

import numpy as nm
import mien.parsers.fileIO as io
import mien.nmpml.data as mdat
import mien.parsers.nmpml as nmpml


def writeFile(f,sps,fname='CubicFunctions',start=0):
	ndoc = nmpml.blankDocument()
	tf = mdat.newData(f, {'Name':fname, 'SampleType':'timeseries', 'SamplesPerSecond':sps, "StartTime":start}) 
	ndoc.newElement(tf)
	io.write(ndoc, fname+".mdat")		
	
	
def CubicBasisFunctions(y,x,j,N):
	#find the values of basis function j for all points in the array y (line geometry); x are interval end points 
	y = nm.sort(y)
	phi_j =nm.empty((0,))
	if j ==1:
		print('j='); print(j)
		for a in y[ nm.logical_and(y>=x[0],y<x[1]) ]:
			val=1
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[1], y<x[2]) ]:
			val = (a-x[1])*(a-x[2])*(a-x[3]) / ( (x[0] - x[1])*(x[0] - x[2])*(x[0] - x[3]) ) + (a-x[0])*(a-x[2])*(a-x[3]) / ( (x[1] - x[0])*(x[1] - x[2])*(x[1] - x[3]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[2],y<x[3]) ]:
			val = (a-x[2])*(a-x[3])*(a-x[4]) / ( (x[1] - x[2])*(x[1] - x[3])*(x[1] - x[4]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[3], y<=x[-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
	elif j ==2:
		print('j='); print(j)
		for a in y[ nm.logical_and(y>=x[0], y<x[1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)			
		for a in y[ nm.logical_and(y>=x[1], y<x[2]) ]:
			val = (a-x[0])*(a-x[1])*(a-x[3]) / ( (x[2] - x[0])*(x[2] - x[1])*(x[2] - x[3]) )
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[2],y<x[3]) ]:
			val = (a-x[1])*(a-x[3])*(a-x[4]) / ( (x[2] - x[1])*(x[2] - x[3])*(x[2] - x[4]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[3], y<x[4]) ]:
			val = (a-x[3])*(a-x[4])*(a-x[5]) / ( (x[2] - x[3])*(x[2] - x[4])*(x[2] - x[5]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[4],y<=x[-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
	elif j < N-1:
		print('j='); print(j)
		for a in y[ nm.logical_and(y>=x[0], y<x[j-2]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j-2], y<x[j-1]) ]:
			val = (a-x[j-3])*(a-x[j-2])*(a-x[j-1]) / ( (x[j] - x[j-3])*(x[j] - x[j-2])*(x[j] - x[j-1]) )
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j-1], y<x[j]) ]:
			val = (a-x[j-2])*(a-x[j-1])*(a-x[j+1]) / ( (x[j] - x[j-2])*(x[j] - x[j-1])*(x[j] - x[j+1]) )
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j],y<x[j+1]) ]:
			val = (a-x[j-1])*(a-x[j+1])*(a-x[j+2]) / ( (x[j] - x[j-1])*(x[j] - x[j+1])*(x[j] - x[j+2]) )
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j+1], y<x[j+2]) ]:
			val = (a-x[j+1])*(a-x[j+2])*(a-x[j+3]) / ( (x[j] - x[j+1])*(x[j] - x[j+2])*(x[j] - x[j+3]) )
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j+2], y<=x[-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
	elif j ==N-1:
		print('j='); print(j)
		for a in y[ nm.logical_and(y>=x[0], y<x[-5]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-5], y<x[-4]) ]:
			val = (a-x[-4])*(a-x[-5])*(a-x[-6]) / ( (x[-3] - x[-4])*(x[-3] - x[-5])*(x[-3] - x[-6]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-4], y<x[-3]) ]:
			val = (a-x[-2])*(a-x[-4])*(a-x[-5]) / ( (x[-3] - x[-2])*(x[-3] - x[-4])*(x[-3] - x[-5]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-3], y<x[-2]) ]:
			val = (a-x[-1])*(a-x[-2])*(a-x[-4]) / ( (x[-3] - x[-1])*(x[-3] - x[-2])*(x[-3] - x[-4]) )
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-2], y<=x[-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)			
	elif j == N:
		print('j='); print(j)
		for a in y[ nm.logical_and(y>=x[0], y<x[-4]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-4], y<x[-3])]:
			val = (a-x[-3])*(a-x[-4])*(a-x[-5]) / ( (x[-2] - x[-3])*(x[-2] - x[-4])*(x[-2] - x[-5]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-3], y<x[-2]) ]:
			val = (a-x[-2])*(a-x[-3])*(a-x[-4]) / ( (x[-1] - x[-2])*(x[-1] - x[-3])*(x[-1] - x[-4]) ) + (a-x[-1])*(a-x[-3])*(a-x[-4]) / ( (x[-2] - x[-1])*(x[-2] - x[-3])*(x[-2] - x[-4]) ) 
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-2], y<=x[-1]) ]:
			val=1
			phi_j=nm.append(phi_j,val)
	else:
		print('Skipping '); print(j)
	return phi_j
				
def HatBasisFunctions(y,x,j,N):
	#find the values of basis function j for all points in the array y (line geometry)
	y = nm.sort(y)
	phi_j =nm.empty((0,))
	if j ==1:
		for a in y[ nm.logical_and(y>=x[0],y<x[1]) ]:
			val=1
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[1], y<x[2]) ]:
			val = 1 - N * (a-x[1])
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[2], y<=x[-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
	elif j < N:
		for a in y[ nm.logical_and(y>=x[0], y<x[j-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)			
		for a in y[ nm.logical_and(y>=x[j-1], y<x[j]) ]:
			val = N * (a-x[j-1])
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j],y<x[j+1]) ]:
			val = 1 - N * (a-x[j])
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[j+1],y<=x[-1]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
	elif j == N:
		for a in y[ nm.logical_and(y>=x[0], y<=x[-3]) ]:
			val=0
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-3], y<x[-2]) ]:
			val = N * (a-x[-3])
			phi_j=nm.append(phi_j,val)
		for a in y[ nm.logical_and(y>=x[-2],y<=x[-1]) ]:
			val=1
			phi_j=nm.append(phi_j,val)
	else:
		print('Skipping '); print(j)
	return phi_j
		
	
def GaussianKernel(x,y):
	K = nm.exp(-(x-y)**2)
	return K
	
def LameKernel(x,y):
	K = nm.ones(y.shape)
	return K
	
def makeQuadVals(QuadPts,ptOfInterest,Poles,N,basisFunc,kernFunc):
	KernVals = kernFunc(ptOfInterest,QuadPts)
	QuadVals = []
	for s in range(N):
		p = basisFunc(QuadPts,Poles,s+1,N)
		#FIXME -- I think quadrature rule has to go in here
		v = nm.dot(p,KernVals)
		QuadVals.append(v)
	return nm.asarray(QuadVals)
	
def TrapezoidalRule(QuadVals,qs):
	quad = qs*( (1./2) * QuadVals[0] + (1./2) * QuadVals[-1] + nm.sum(QuadVals[1:-1]) )
	return quad



if __name__ == '__main__':
	N = 20 #number of force pts and the number of basis functions; N+1 is number of intervals; N+2 is final number of pts (the end pts)
	h = 1./N
	x = nm.arange(h/2,1,h)
	x = nm.hstack((0,x,1))
	sps = N*20
	y = nm.arange(0,1,1./sps)
	# f = nm.ones(N)
	QuadVals = makeQuadVals(y,0,x,N,HatBasisFunctions,LameKernel)
	Quad = TrapezoidalRule(QuadVals,1./sps)
	print(Quad)
	# HF = nm.empty((len(y),N))
	# for k in range(N):
	# 	HF[:,k] = HatBasisFunctions(y,x,k+1,N)
	writeFile(QuadVals,1,'QuadVals',1)
	