#!/usr/bin/env python

import numpy as nm
import mien.parsers.fileIO as io
import mien.nmpml.data as mdat
import mien.parsers.nmpml as nmpml


def writeFile(f,fname='TestMe',sps=None,start=0):
	ndoc = nmpml.blankDocument()
	if sps:
		tf = mdat.newData(f, {'Name':fname, 'SampleType':'timeseries', 'SamplesPerSecond':sps, "StartTime":start}) 
	else:			
		tf = mdat.newData(f, {'Name':fname, 'SampleType':'function'}) 
	ndoc.newElement(tf)
	io.write(ndoc, fname+".mdat")		
		
def DiscretizeCircle(N,R=1):
	# N pt discretization of a circle; x1 in first row, x2 in second row
	angles = nm.arange(0,2*nm.pi,2*nm.pi/N)
	y=[R*nm.cos(angles),R*nm.sin(angles)]
	return nm.asarray(y)
	
def DiscretizeChords(y,M,R=1):
	# N2 pt discretization of each of N chords on a circle given by the points y; x1 in first row, x2 in second row
	zx =[]
	zy=[]
	cspacing = 1./M
	t = nm.arange(0,1,cspacing)
	for k in range(y.shape[1]):
		if k == y.shape[1]-1:
			pos = nm.array([R,0])
		else:
			pos = y[:,k+1]
		xcoord=y[0,k] + t*(pos[0]-y[0,k])
		ycoord=y[1,k] + t*(pos[1]-y[1,k])
		zx.append(xcoord)
		zy.append(ycoord)
	zx=nm.asarray(zx).flatten()
	zy=nm.asarray(zy).flatten()
	z = nm.row_stack([zx,zy])
	N = y.shape[1]
	chord = 2*R*nm.sin(nm.pi/N)
	spacing = chord*cspacing
	return z, spacing

	
def CubicBasisFunctions(z,y,j,N):
	#find the values of basis function j for all points in the array z (quad pts, circle geometry); y are interval end points
	pass
				
def HatBasisFunctions(z,y,j,N,R=1):
	#find the values of basis function j for all points in the array z (quad pts, circle geometry); y are interval end points
	#start counting basis functions at 1, not 0
	# each phi_j is associated with pole y[j-1]
	chord = 2*R*nm.sin(nm.pi/N) #chord length between adjacent y points (needed for linear interpolation)
	angz=nm.arctan2(z[1,:],z[0,:]) #find the angular positions of both y and z to figure out where the z points fall
	angz[angz<0]+=2*nm.pi			# change interval from (-pi,pi] to [0,2*pi)
	angy=nm.arctan2(y[1,:],y[0,:]) 
	angy[angy<0]+=2*nm.pi
	phi_j =nm.empty((0,)) #get ready to store function values
	
	#pad with zeros below interval [y_{j-2}, y_{j}), except for phi_1, which is a special case
	if j==1:
		a = z[ :, nm.logical_and(angz >= angy[0], angz < angy[1]) ]
		val = 1 - (1./chord) * nm.sqrt(nm.sum((a-nm.array([y[:,0]]).transpose())**2,0))
		phi_j=nm.append(phi_j,val)
		a = z[ :, nm.logical_and(angz >= angy[1], angz < angy[N-1]) ]
		val = nm.zeros(a.shape[1])
		phi_j=nm.append(phi_j,val)		
		a = z[ :, nm.logical_and(angz >= angy[N-1], angz < 2*nm.pi) ]
		val=(1./chord) * nm.sqrt(nm.sum((a-nm.array([y[:,N-1]]).transpose())**2,0))
		phi_j=nm.append(phi_j,val)
		return phi_j
	elif j >2 and j <= N:
		a = z[ :, nm.logical_and(angz >= angy[0], angz < angy[j-2]) ]
		val = nm.zeros(a.shape[1])
		phi_j=nm.append(phi_j,val)
	elif j ==2:
		pass
	else:
		print('Skipping ' + str(j))
		return None
	#upward going linear interpolation on [y_{j-2}, y_{j-1})
	a = z[ :, nm.logical_and(angz >= angy[j-2], angz < angy[j-1]) ]
	val=(1./chord) * nm.sqrt(nm.sum((a-nm.array([y[:,j-2]]).transpose())**2,0))
	phi_j=nm.append(phi_j,val)
	#downward going linear interpolation on [y_{j-1}, y_{j}); fill in 2*pi for the "N+1"st point
	if j < N:
		ang = angy[j]
	elif j == N:
		ang = 2*nm.pi
	a = z[ :, nm.logical_and(angz >= angy[j-1], angz < ang) ]
	val = 1 - (1./chord) * nm.sqrt(nm.sum((a-nm.array([y[:,j-1]]).transpose())**2,0))
	phi_j=nm.append(phi_j,val)
	#pad with zeros above interval [y_{j-2}, y_{j}) if j < N
	if j < N:
		a = z[ :, nm.logical_and(angz >= ang, angz < 2*nm.pi) ]
		val=nm.zeros(a.shape[1])
		phi_j=nm.append(phi_j,val)
	return phi_j
	
def ConstantBasisFunctions(z,y,j,N):
	#find the values of basis function j for all points in the array z (quad pts, circle geometry); y are interval end points 
	#start counting basis functions at 1, not 0
	# each phi_j is associated with pole y[j-1]
	angz=nm.arctan2(z[1,:],z[0,:]) #find the angular positions of both y and z to figure out where the z points fall
	angz[angz<0]+=2*nm.pi			# change interval from (-pi,pi] to [0,2*pi)
	angy=nm.arctan2(y[1,:],y[0,:]) 
	angy[angy<0]+=2*nm.pi
	phi_j =nm.empty((0,)) #get ready to store function values
	
	#pad with zeros below interval [y_{j-1}, y_{j})
	a = z[ :, nm.logical_and(angz >= angy[0], angz < angy[j-1]) ]
	val = nm.zeros(a.shape[1])
	phi_j=nm.append(phi_j,val)
	#piecewise constant 1 on [y_{j-1}, y_{j}); fill in 2*pi for the "N+1"st point
	if j < N and j > 0:
		ang = angy[j]
	elif j == N:
		ang = 2*nm.pi
	else:
		print('Skipping ' + str(j))
		return None
	a = z[ :, nm.logical_and(angz >= angy[j-1], angz < ang) ]
	val=nm.ones(a.shape[1])
	phi_j=nm.append(phi_j,val)
	#pad with zeros above interval [y_{j-1}, y_{j})
	a = z[ :, nm.logical_and(angz >= ang, angz < 2*nm.pi) ]
	val=nm.zeros(a.shape[1])
	phi_j=nm.append(phi_j,val)
	return phi_j
		
def HatBasisFunctionSingle(M):
	#get ready to store function values
	phi_j =nm.empty((0,)) 
	#upward line
	val=nm.arange(0,1,1./M)
	phi_j=nm.append(phi_j,val)
	#downward line
	val=nm.arange(1,-1./M,-1./M)
	phi_j=nm.append(phi_j,val)
	return phi_j

def CstBasisFunctionSingle(M):
	return nm.ones(M+1)
			
def GaussianKernel(x,z):
	# z are integration pts, x is point of interest (evaluation point)
	K = nm.exp(-(z-x)**2)
	return K
	
def LameKernel(x,z):
	# z are integration pts, x is point of interest (evaluation point)
	K = nm.ones(z.shape[1])
	return K

def LaplacianKernel(x,z):
	# z are integration pts ON A CIRCLE (will need different functions for the chords of a circle or other domain), x is point of interest (evaluation point)
	n = CircleNormal(z) # find the outward pointing normal 
	dif  = nm.asarray([z[0,:] - x[0], z[1,:] - x[1]])
	dprod = nm.asarray([nm.dot(dif[:,k],n[:,k]) for k in range(dif.shape[1])])
	K = dprod/(2*nm.pi*(nm.sum(dif**2,0)))
	return K
	
def LameDensity(y):
	# y are poles
	D = nm.ones(y.shape[1])
	return D
	
def CosDensity(y,R=1):
	# y are poles
	angy=nm.arctan2(y[1,:],y[0,:]) 
	D = 2*R*nm.cos(3*angy)
	return D

def SinDensity(y,R=1):
	# y are poles
	angy=nm.arctan2(y[1,:],y[0,:]) 
	D = 2*R*nm.sin(3*angy)
	return D

def CircleNormal(z):
	#calculates the normal to a circle at all of the points z ON the circle (x1 in first row, x2 in second row)
	a = nm.sqrt(nm.sum(z**2,0)) #circle radius
	n = z/a
	return n
	
def ChordNormal(z):
	pass
		
def findIndexHat(j,N,M):
	if nm.mod(M,2):
		print('M must be even or there will be trouble!')	
	if j ==0:
		ind = range(N*M-M,N*M)
		ind.extend(range(M+1))
	elif j == N-1:
		ind = range(N*M-2*M,N*M)
		ind.append(0)
	elif j >0 and j < N-1:
		ind = range((j-1)*M,(j+1)*M+1)
	else:
		print('Skipping ' + str(j))
		return None
	return ind		

def findIndexCst(j,N,M):
	if nm.mod(M,2):
		print('M must be even or there will be trouble!')	
	if j == N-1:
		ind = range(N*M-M,N*M)
		ind.append(0)
	elif j >-1 and j < N-1:
		ind = range(j*M,(j+1)*M+1)
	else:
		print('Skipping ' + str(j))
		return None
	return ind		

def findIntVal(QuadPts,ptsOfInterest,Poles,N,spacing,basisFunc,kernFunc,intFunc,densityFunc):
	'''Integrate against a kernel around a domain. QuadPts are a set of points where the actual quadrature is to take place, ptsOfInterest is where the evaluation is taking place, Poles are where the boundary condition (or density function) values are given, N is the number of poles, spacing is the spacing between QuadPts, basisFunc is a function handle for calculating basis elements, kernFunc is a function handle for calculating kernel values, intFunc is a function handle for which integration rule to use, and densityFunc is a function handle that returns the values of the B.C.s (or density function) at Poles. Note all coordinate matrices have x-values in the first ROW, and y-values in the second ROW.'''
	if len(ptsOfInterest.shape) < 2:
		num = 1
	else:
		num = ptsOfInterest.shape[1]
	DensVals = densityFunc(Poles)
	IntVal = nm.zeros((num,))
	for k in range(num):
		if num == 1:
			x = ptsOfInterest[:]
		else:
			x = ptsOfInterest[:,k]
		KernVals = kernFunc(x,QuadPts)
		for s in range(N):
			p = basisFunc(QuadPts,Poles,s+1,N) #make basis function s+1
			v = p*KernVals #broadcast multiply the basis function against the kernel
			s_int_val = intFunc(v,spacing) #apply quadrature rule to the product
			IntVal[k]+=s_int_val*DensVals[s] #multiply the integral of the s+1-th basis function and kernel by the density function value at the s-th pole
	return IntVal

def findIntValSingle(QuadPts,spacing,ptsOfInterest,Poles,N,M,kernFunc,intFunc,densityFunc,elementType='hat'):
	'''Integrate against a kernel around a domain. QuadPts are a set of points where the actual quadrature is to take place, ptsOfInterest is where the evaluation is taking place, Poles are where the boundary condition (or density function) values are given, N is the number of poles, spacing is the spacing between QuadPts, basisFunc is a function handle for calculating basis elements, kernFunc is a function handle for calculating kernel values, intFunc is a function handle for which integration rule to use, and densityFunc is a function handle that returns the values of the B.C.s (or density function) at Poles. Note all coordinate matrices have x-values in the first ROW, and y-values in the second ROW.'''
	num = ptsOfInterest.shape[1]
	DensVals = densityFunc(Poles)
	if elementType == 'hat':
		Element = HatBasisFunctionSingle(M)
		findIndex = findIndexHat
	elif elementType == 'cst':
		Element = CstBasisFunctionSingle(M)
		findIndex = findIndexCst
	else:
		print('Element type not recognized. Ending program.')
		return None
	IntVal = nm.zeros((num,))
	for k in range(num):
		if num == 1:
			x = ptsOfInterest[:]
		else:
			x = ptsOfInterest[:,k]
		KernVals = kernFunc(x,QuadPts)
		for s in range(N):
			ind = findIndex(s,N,M) #make the indices for basis function s
			v = Element*KernVals[ind] #broadcast multiply the basis function against the kernel
			s_int_val = intFunc(v,spacing) #apply quadrature rule to the product
			IntVal[k]+=s_int_val*DensVals[s] #multiply the integral of the s+1-th basis function and kernel by the density function value at the s-th pole
	return IntVal

def TrapezoidalRule(ptvals,quadspacing):
	quad = quadspacing*( 0.5*ptvals[0] + 0.5*ptvals[-1] + nm.sum(ptvals[1:-1]) )
	return quad

def SimpsonsRule(ptvals,quadspacing):
	# Count the first point twice because we are on a circle. 
	quad = (quadspacing/3.)*( 2*ptvals[0] + 2*nm.sum(ptvals[2::2]) + 4*nm.sum(ptvals[1::2]) )
	return quad



if __name__ == '__main__':
	#test the single basis functions and the findIndex routine
	# r = nm.arange(0,1,0.1)
	# a = nm.arange(0,2*nm.pi, nm.pi/6)
	# xc = []
	# yc = []
	# for k in range(r.shape[0]):
	# 	xc.append(r[k]*nm.cos(a))
	# 	yc.append(r[k]*nm.sin(a))
	# xc=nm.asarray(xc).flatten()
	# yc=nm.asarray(yc).flatten()
	# ptsOfInterest = nm.row_stack([xc,yc]) 
	ptsOfInterest = 1.01*nm.array([[nm.cos(0.7),nm.sin(0.7)]]).transpose()
	ActualVals = -( (nm.sum(ptsOfInterest**2,0))**(-3./2) ) * nm.sin(3*nm.arctan2(ptsOfInterest[1,:],ptsOfInterest[0,:]))
	ErrValsHat = []
	ErrValsCst = []
	numchords=[2**k for k in range(2,11)]	
	M = 1000
	for N in numchords:
		Poles = DiscretizeCircle(N)
		QuadPts, spacing = DiscretizeChords(Poles,M)
		IntVals=findIntValSingle(QuadPts,spacing,ptsOfInterest,Poles,N,M,LaplacianKernel,TrapezoidalRule,SinDensity,'hat')
		errhat = nm.sqrt((IntVals-ActualVals)**2).max()
		IntVals=findIntValSingle(QuadPts,spacing,ptsOfInterest,Poles,N,M,LaplacianKernel,TrapezoidalRule,SinDensity,'cst')
		errcst = nm.sqrt((IntVals-ActualVals)**2).max()
		ErrValsHat.append(errhat)
		ErrValsCst.append(errcst)
		print(errhat)
		print(errcst)
	ErrValsCst=nm.asarray(ErrValsCst).transpose()
	ErrValsHat=nm.asarray(ErrValsHat).transpose()
	numchords=nm.asarray(numchords).transpose()
	writeFile(nm.column_stack([numchords,ErrValsCst,ErrValsHat]),'ErrorValsM1000')
	
	# N = 10 #number of pts on the circle and the number of basis functions and intervals
	# y = DiscretizeCircle(N) #2D points
	# # make quadpts
	# M = 20
	# z, sp = DiscretizeChords(y,M) 
	# phic = []
	# phih = []
	# for k in range(N):
	# 	ck=ConstantBasisFunctions(z,y,k+1,N)
	# 	phic.append(ck)
	# 	hk=HatBasisFunctions(z,y,k+1,N)
	# 	phih.append(hk)
	# phic=nm.asarray(phic).transpose()	
	# phih=nm.asarray(phih).transpose()	
	# writeFile(nm.column_stack([phic,phih]),'TestBases',sp)	
	
	# #whole test case
	# # r = nm.arange(0,1,0.1)
	# # a = nm.arange(0, 2*nm.pi, nm.pi/6)
	# r = nm.arange(0,1,0.05)
	# a = nm.arange(0,2*nm.pi, nm.pi/6)
	# xc = []
	# yc = []
	# for k in range(r.shape[0]):
	# 	xc.append(r[k]*nm.cos(a))
	# 	yc.append(r[k]*nm.sin(a))
	# xc=nm.asarray(xc).flatten()
	# yc=nm.asarray(yc).flatten()
	# xgrid = nm.row_stack([xc,yc]) 
	# actualIntVals = ( (nm.sum(xgrid**2,0))**(3./2) ) * nm.cos(3*nm.arctan2(xgrid[1,:],xgrid[0,:]))
	# IntValsHat = []
	# IntValsCst = []
	# ErrValsHat = []
	# ErrValsCst = []
	# numchords=[]
	# numptsperchord = []
	# for N in [2**k for k in range(2,11)]: #number of pts on the circle and the number of basis functions and intervals
	# 	y = DiscretizeCircle(N) #2D points
	# 	for M in 1000: 
	# 		print((N,M))
	# 		z, sp = DiscretizeChords(y,M) 
	# 		iv = findIntVal(z,xgrid,y,N,sp,ConstantBasisFunctions,LaplacianKernel,TrapezoidalRule,CosDensity)
	# 		IntValsCst.append(iv) 
	# 		err = nm.sqrt((iv - actualIntVals)**2)
	# 		ErrValsCst.append(err)
	# 		iv = findIntVal(z,xgrid,y,N,sp,HatBasisFunctions,LaplacianKernel,TrapezoidalRule,CosDensity)
	# 		IntValsHat.append(iv) 
	# 		err = nm.sqrt((iv - actualIntVals)**2)
	# 		ErrValsHat.append(err)
	# 		numchords.append(N)
	# 		numptsperchord.append(M)
	# IntValsCst=nm.asarray(IntValsCst).transpose()
	# IntValsHat=nm.asarray(IntValsHat).transpose()
	# ErrValsCst=nm.asarray(ErrValsCst).transpose()
	# ErrValsHat=nm.asarray(ErrValsHat).transpose()
	# numchords=nm.asarray(numchords)
	# numptsperchord=nm.asarray(numptsperchord)
	# writeFile(nm.column_stack([xgrid.transpose(),actualIntVals.transpose()]),'ActualValsOnGrid_N16_run1')
	# writeFile(IntValsCst,'IntValsCst_N16_run1')
	# writeFile(IntValsHat,'IntValsHat_N16_run1')
	# writeFile(ErrValsCst,'ErrValsCst_N16_run1')
	# writeFile(ErrValsHat,'ErrValsHat_N16_run1')
	# numpts=nm.column_stack([numchords,numptsperchord,nm.amax(ErrValsCst,0),nm.amax(ErrValsHat,0)])
	# writeFile(numpts,'ResultSummary_N16_run1')
	# print('Summary')
	# print(numpts)
	# # print('Cst Err')
	# # print(ErrValsCst)
	# # print('Hat Err')
	# # print(ErrValsHat)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	