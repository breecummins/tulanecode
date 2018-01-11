#!/usr/bin/env python

import numpy as nm

class Stokeslet2DReg2():
	'''Regularized 2D Stokeslet, cube power.'''
	def __init__(self,blobParam,fluidViscosity):
		self.eps = blobParam		
		self.mu = fluidViscosity	

	def singleKernVal(self,obsPt=nm.array([1.,0.]),poleStrength=nm.array([1.,0.]),poleLocation=nm.array([0.,0.])):
		'''For testing purposes.'''
		dif = obsPt - poleLocation
		r2 = nm.sum(dif**2) + self.eps**2
		H1 = 2*self.eps**2/r2 - nm.log(r2)
		H2 = 2/r2
		term1 = poleStrength*H1
		term2 = nm.dot(poleStrength,dif)*dif*H2
		val = (term1+term2)/(8*nm.pi*self.mu)
		return val

	def kernVals2D(self,obsPt,nodes):
		'''Builds the values for the original matrix.'''
		dif = obsPt - nodes
		r2 = nm.sum(dif**2,1) + self.eps**2
		H1 = 2*self.eps**2/r2 - nm.log(r2)
		H2 = 2/r2
		xterm =  H1 + (dif[:,0]**2)*H2
		mixedterm = (dif[:,0]*dif[:,1])*H2
		yterm = H1 + (dif[:,1]**2)*H2
		return xterm/(8*nm.pi*self.mu), mixedterm/(8*nm.pi*self.mu), yterm/(8*nm.pi*self.mu)

	def makeMatrix2DOriginal(self,obspts,nodes):
		if len(nodes.shape) == 1:
			if len(nodes) != 2:
				raise ValueError('Requires two dimensional inputs.')
			N = 1
		elif nodes.shape[1] != 2:
			raise ValueError('Requires two dimensional inputs.')
		else:
			N = nodes.shape[0]
		if len(obspts.shape) == 2:
			npts = obspts.shape[0]
		elif len(obspts.shape) == 1:
			npts = 1
		else:
			raise ValueError('Observation points must be an n by 2 matrix.')
		mat = nm.zeros((2*npts,2*N))
		for k in range(npts):
			if npts == 1:
				xterm, mixedterm, yterm = self.kernVals2D(obspts,nodes)
			else:	
				xterm, mixedterm, yterm = self.kernVals2D(obspts[k,:],nodes)
			mat[2*k,2*nm.arange(N)] = xterm
			mat[2*k+1,2*nm.arange(N)+1]= yterm
			mat[2*k,2*nm.arange(N)+1] = mixedterm
			mat[2*k+1,2*nm.arange(N)] = mixedterm
		return mat
		
		
if __name__ == "__main__":
	#make a blob
	blob = Stokeslet2DReg2(2.e-2,3.e-1)
	#discretize boundary to make nodes
	angles = nm.arange(0,2*nm.pi,2*nm.pi/32)
	nodes=nm.column_stack([0.25*nm.cos(angles),0.25*nm.sin(angles)])
	#choose observation points at which to calculate fluid velocity
	spacing=0.01
	r = nm.arange(-0.5,0.5+spacing/2.,spacing)
	n=r.shape[0]
	obspts=nm.zeros((n**2,2))
	for j in range(n):
		obspts[j*n:(j+1)*n,0] = r[j]*nm.ones(n)
		obspts[j*n:(j+1)*n,1] = r
	#build the matrix to find the forces
	M1 = blob.makeMatrix2DOriginal(nodes,nodes)
	#Now, solve U0 = M1*F for F, where U0 are the boundary conditions and F = [f11, f12, f21, f22, f31, f32,...] where (fn1, fn2) is the force at the n-th node. 
	#Then make the matrix to find the velocity at obspts
	M2 = blob.makeMatrix2DOriginal(obspts,nodes)
	#At last find u = M2*F, where u is the sought after velocity and F are the forces that you just found
	
	
	
	