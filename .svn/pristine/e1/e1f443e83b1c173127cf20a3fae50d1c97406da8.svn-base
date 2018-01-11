#!/usr/bin/env python

import numpy as nm

class Kernel:
	'''This is the base class for some integration kernels.'''
	def __init__(self):
		pass
		
	def makeMatrix2DBEM(self,obspts,BEobj,BFobj,QMobj):
		quadpts = BEobj.quadpts
		if quadpts.shape[1] != 2:
			raise ValueError('Requires two dimensional inputs.')
		N = BEobj.N
		P = BFobj.phi.shape[0]
		coeffs = QMobj.makeCoeffs(P)
		intcoeffs = coeffs*BFobj.phi #what actually multiplies the kernel in the quad rule
		if len(obspts.shape) == 2:
			npts = obspts.shape[0]
		elif len(obspts.shape) == 1:
			npts = 1
		else:
			raise ValueError('Observation points must be an n by 2 matrix.')
		mat = nm.zeros((2*npts,2*N))
		for k in range(npts):
			if npts == 1:
				xterm, mixedterm, yterm = self.kernVals2D(obspts,quadpts)
			else:	
				xterm, mixedterm, yterm = self.kernVals2D(obspts[k,:],quadpts)
			for j in range(N):
				ind=BFobj.findIndex(j,N)
				mat[2*k,2*j] = nm.sum(intcoeffs*xterm[ind])
				mat[2*k+1,2*j+1]= nm.sum(intcoeffs*yterm[ind])
				vm = nm.sum(intcoeffs*mixedterm[ind])
				mat[2*k,2*j+1] = vm
				mat[2*k+1,2*j] = vm
		return mat

	def makeMatrix2D(self,obspts,poles):
		if len(poles.shape) == 1:
			if len(poles) != 2:
				raise ValueError('Requires two dimensional inputs.')
			N = 1
		elif poles.shape[1] != 2:
			raise ValueError('Requires two dimensional inputs.')
		else:
			N = poles.shape[0]
		if len(obspts.shape) == 2:
			npts = obspts.shape[0]
		elif len(obspts.shape) == 1:
			npts = 1
		else:
			raise ValueError('Observation points must be an n by 2 matrix.')
		mat = nm.zeros((2*npts,2*N))
		for k in range(npts):
			if npts == 1:
				xterm, mixedterm, yterm = self.kernVals2D(obspts,poles)
			else:	
				xterm, mixedterm, yterm = self.kernVals2D(obspts[k,:],poles)
			for j in range(N):
				mat[2*k,2*j] = xterm[j]
				mat[2*k+1,2*j+1]= yterm[j]
				vm = mixedterm[j]
				mat[2*k,2*j+1] = vm
				mat[2*k+1,2*j] = vm
		return mat

	
class Stokeslet2DReg(Kernel):
	'''Regularized 2D Stokeslet from 2001 paper.'''
	def __init__(self,blobParam,fluidViscosity):
		Kernel.__init__(self)
		self.eps = blobParam		
		self.mu = fluidViscosity	

	def singleKernVal(self,obsPt=nm.array([1.,0.]),poleStrength=nm.array([1.,0.]),poleLocation=nm.array([0.,0.])):
		'''For testing purposes.'''
		dif = obsPt - poleLocation
		rr = nm.sqrt(nm.sum(dif**2) + self.eps**2)
		rre = rr + self.eps
		rr2e = rr + 2*self.eps
		term1 = -poleStrength*( nm.log(rre) - self.eps*rr2e/(rre*rr) ) 
		term2 = nm.dot(poleStrength,dif)*dif*rr2e/((rre**2)*rr)
		val = (term1+term2)/(4*nm.pi*self.mu)
		return val
		
	def kernVals2D(self,obsPt,quadpts):
		'''Builds the values for the matrix.'''
		dif = obsPt - quadpts
		rr = nm.sqrt(nm.sum(dif**2,1) + self.eps**2)
		rre = rr + self.eps
		rr2e = rr + 2*self.eps
		lr = -nm.log(rre) + self.eps*rr2e/(rre*rr)
		mf = rr2e/((rre**2)*rr)
		xterm =  lr + (dif[:,0]**2)*mf
		mixedterm = (dif[:,0]*dif[:,1])*mf
		yterm = lr + (dif[:,1]**2)*mf
		return xterm/(4*nm.pi*self.mu), mixedterm/(4*nm.pi*self.mu), yterm/(4*nm.pi*self.mu)


class Stokeslet2DReg2(Kernel):
	'''Regularized 2D Stokeslet, cube power.'''
	def __init__(self,blobParam,fluidViscosity):
		Kernel.__init__(self)
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

	def kernVals2D(self,obsPt,quadpts):
		'''Builds the values for the matrix.'''
		dif = obsPt - quadpts
		r2 = nm.sum(dif**2,1) + self.eps**2
		H1 = 2*self.eps**2/r2 - nm.log(r2)
		H2 = 2/r2
		xterm =  H1 + (dif[:,0]**2)*H2
		mixedterm = (dif[:,0]*dif[:,1])*H2
		yterm = H1 + (dif[:,1]**2)*H2
		return xterm/(8*nm.pi*self.mu), mixedterm/(8*nm.pi*self.mu), yterm/(8*nm.pi*self.mu)


class Stokeslet2D(Kernel):
	'''Singular 2D Stokeslet.'''
	def __init__(self,fluidViscosity):
		Kernel.__init__(self)
		self.mu = fluidViscosity	

	def singleKernVal(self,obsPt=nm.array([1.,0.]),poleStrength=nm.array([1.,0.]),poleLocation=nm.array([0.,0.])):
		dif = obsPt - poleLocation
		r = nm.sqrt(nm.sum(dif**2))
		term1 = -poleStrength*nm.log(r) 
		term2 = (poleStrength*dif)*dif/(r**2)
		val = (term1+term2)/(4*nm.pi*self.mu)
		return val

	def kernVals2D(self,obsPt,quadpts):
		dif = obsPt - quadpts
		r = nm.sqrt(nm.sum(dif**2,1))
		lr = -nm.log(r)
		xterm =  lr + (dif[:,0]/r)**2
		mixedterm = dif[:,0]*dif[:,1]/r**2
		yterm = lr + (dif[:,1]/r)**2
		return xterm/(4*nm.pi*self.mu), mixedterm/(4*nm.pi*self.mu), yterm/(4*nm.pi*self.mu)


class Identity(Kernel):
	'''Identity.'''
	def __init__(self):
		Kernel.__init__(self)
		
	def singleKernVal(self,obsPt=nm.array([1.,0.]),poleStrength=nm.array([1.,0.]),poleLocation=nm.array([0.,0.])):
		return poleStrength
		
	def kernVals2D(self,obspts,quadpts):
		xterm = nm.ones(quadpts.shape[0])
		yterm = nm.ones(quadpts.shape[0])
		mixedterm = nm.zeros(quadpts.shape[0])
		return xterm, mixedterm, yterm


if __name__ == '__main__':
	#test scenario for matrix building for the original method
	mu = 1
	eps = 0.25
	blob = Stokeslet2DReg(eps,mu)

	f1 = nm.array([1,0])
	x1 = nm.array([0,1])
	f2 = nm.array([-1,0])
	x2 = nm.array([0,-1])
	x0 = nm.array([0,0])
	poles = nm.row_stack([x1,x2])
	forces = nm.row_stack([f1,f2])
	
	#one by one
	velh = nm.array([0,0])
	for k in range(2):
		v=blob.singleKernVal(x0,forces[k,:],poles[k,:])
		velh = velh + v
	
	#by matrix
	mat = blob.makeMatrix2D(x0,poles)
	fvec = nm.append(f1,f2)
	velm = nm.dot(mat,fvec)
	
	print(velh)
	print(velm)
	print(nm.sqrt(nm.sum(velh-velm)**2))








	