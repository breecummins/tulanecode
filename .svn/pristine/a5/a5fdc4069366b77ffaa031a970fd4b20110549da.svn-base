#!/usr/bin/env python

import numpy as nm
import sys

class BasisFunc:
	'''Create basis functions over the boundary elements.'''
	def __init__(self,M):
		self.M = M
		
		
class HatBasisFunc(BasisFunc):
	'''Bilinear basis functions on a *closed* curve.'''
	def __init__(self,M,quadtype):
		BasisFunc.__init__(self,M)
		#make generic hat. Quadrature type should not affect hat basis function, but it (quadtype) is included as as argument to support a uniform API.
		upline=nm.arange(0,1,1./M)
		downline=nm.arange(1,-1./(2*M),-1./M)
		self.phi=nm.append(upline,downline)
		self.quadtype = quadtype
		
	def findIndex(self,j,N):
		'''Find the location of hat j when there are N chords with M points apiece on a closed curve.'''
		M = self.M
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

class CstBasisFunc(BasisFunc):
	'''Constant basis functions on a *closed* curve.'''
	def __init__(self,M,quadtype):
		BasisFunc.__init__(self,M)
		if quadtype == 'r':  #Riemann sum. Want to avoid overlapping points.
			self.phi=nm.ones(M)
		elif quadtype == 't' or quadtype == 's': #trapezoid or simpson's rule. Need both endpts of the subinterval.
			self.phi=nm.ones(M+1) 
		else:
			print('Quadrature type not recognized.')
			sys.exit()	
		self.quadtype = quadtype
			
	def findIndex(self,j,N):
		'''Find the location of constant basis j when there are N chords with M points apiece on a closed curve.'''
		M = self.M
		if self.quadtype == 'r':
			if j ==0: 
				ind = range(int(N*M-nm.floor(M/2.)),N*M)
				ind.extend(range(int(nm.ceil(M/2.))))
			elif j > 0 and j < N:
				ind = range(int(j*M-nm.floor(M/2.)), int(j*M + nm.ceil(M/2.)))
			else:
				print('Skipping ' + str(j))
				return None
		else:	
			if j ==0: 
				ind = range(int(N*M-nm.floor(M/2.)),N*M)
				ind.extend(range(int(nm.ceil(M/2.)+1)))
			elif j > 0 and j < N:
				ind = range(int(j*M-nm.floor(M/2.)), int(j*M + nm.ceil(M/2.)+1))
			else:
				print('Skipping ' + str(j))
				return None
		return ind		
				

if __name__ == '__main__':
	#unit test for basis functions completed 05/07/2010. The findIndex method is not completely tested in conjunction with IntegrationKernels.py, although I believe it is doing the right thing in isolation. The results here showed h^2 convergence to the real solution for all combinations of quadrature method and basis function. The loglog plots and convergence rates were computed in Matlab. M was tried at 1, 4, or 20, although not exhaustively on all combinations.
	data = []
	numchords = [2**k for k in range(2,12)] 
	quadtype = 'r' #'r' or 't' or 's'
	basistype = 'c'	#'c' or 'h'
	for N in numchords:
		ang = nm.arange(0,2*nm.pi,2*nm.pi/N)
		f = -ang**2 + 2*nm.pi*ang #function to integrate at nodes
		import BoundaryElements as be
		import QuadratureMethods as qm
		M = 1
		circ = be.Circle(N,M)
		circ.LinearElements2D()
		if quadtype == 't':
			quad = qm.CompTrapRule(circ.quadspacing)
			print('Trapezoid Rule')
		elif quadtype == 's':
			quad = qm.CompSimpRule(circ.quadspacing)
			print("Simpson's Rule")
		elif quadtype == 'r':
			quad = qm.RiemannSum(circ.quadspacing)
			print('Riemann Sum')
		if basistype == 'c':
			basfun = CstBasisFunc(M,quadtype)
			print(basfun.phi)
		elif basistype == 'h':
			basfun = HatBasisFunc(M,quadtype)
			print(basfun.phi)
		P = basfun.phi.shape[0]
		coeffs = quad.makeCoeffs(P)
		subintegral = nm.sum(coeffs*basfun.phi) 
		wholeintegral = nm.sum(f*subintegral)
		print(wholeintegral)
		data.append(wholeintegral)
	data = nm.asarray(data)
	exactans = (4*nm.pi**3)/3.
	import mat2py as mf
	mf.write('/Users/bcummins/basistest.mat',{'N':numchords,'intvals':data,'exactans':exactans})
	print('Exact answer is ' + str(exactans))
	
	