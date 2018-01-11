#!/usr/bin/env python

import numpy as nm
import sys
import PeriodicBCDerivatives as P

def makeGrid(N,h):
	gp = nm.arange(0,1-h/2,h)
	ux = nm.tile(gp,(1,N))
	uy = nm.tile(gp,(N,1))
	return ux,uy
	
def makeH(u,v,h):
	H1 = u*P.xDerivWide(u,h) + P.vInterp(v)*P.yDerivWide(u,h)
	H2 = P.uInterp(u)*P.xDerivWide(v,h) + v*P.yDerivWide(v,h)
	return H1, H2
	
def timeExtrap(Hn,Hnm1):
	return 1.5*Hn - 0.5*Hnm1
	
def pressureGrad(q,h):
	'''Getting pressure at cell edges. See MAC grid for derivative choice.'''
	dqx = P.xDerivSlimMinus(q,h)
	dqy = P.yDerivSlimMinus(q,h)
	return dqx, dqy
	
def makeFFTSine(N):
	'''Make sine factor for FFt method of solving equations. numpy's FFT starts indexing at 0.'''
	k = nm.asarray([[i]*N for i in range(N)])
	l = k.transpose()
	sindenom = (nm.sin(nm.pi*k/N))^2 + (nm.sin(nm.pi*l/N))^2
	return sindenom

def makeFFTRHS(RHS):
	Rhat = nm.fft.fft2(RHS)
	sindenom = makeFFTSine(RHS.shape[0])
	return Rhat, sindenom
	
		
def phiSolver(m1,m2,h):
	'''Use FFTs to solve the Laplacian for phi.'''
	RHS = P.divEdge2Cen(m1,m2,h)
	Rhat, sindenom = makeFFTRHS(RHS)
	Rhat[0,0] = 0 # set dc element to 0
	sindenom[0,0] = 1 #We're going to divide by this, so it needs to be nonzero.
	phihat = h^2*Rhat/(-4*sindenom)
	return numpy.fft.ifft2(phihat)
	
def mSolver(RHS,h,nu,dt):
	'''Use FFTs to solve for m.'''
	Rhat, sindenom = makeFFTRHS(RHS) #No bad points, so no correction needed
	mhat = Rhat/(1 + (2*nu*dt/h^2)*sindenom)
	return numpy.fft.ifft2(mhat)
	
def timeLoop(dt, Tfinal, nu, N, h, H0, Hm1, u0, v0, qmethod, q0):
	pass
	
def exactSoln():
	pass
	
if __name__=='__main__':
	#First do q = 0.
	nu = 0.1; N = 32; h = 1./N; dt = h; Tfinal = 0.5	
	