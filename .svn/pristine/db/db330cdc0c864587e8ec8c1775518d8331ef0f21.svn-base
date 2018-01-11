#!/usr/bin/env python

import numpy as nm
import os
import sys

from scipy.sparse.linalg import gmres
import mat2py as mf

class Stokeslet2DReg2():
	'''Regularized 2D Stokeslet, cube power.'''
	def __init__(self,blobParam,fluidViscosity):
		self.eps = blobParam		
		self.mu = fluidViscosity	

	def singleKernVal(self,obsPt=(1.,0.),poleStrength=(1.,0.),poleLocation=(0.,0.)):
		'''For testing purposes.'''
		obsPt = nm.array(obsPt)
		poleStrength = nm.array(poleStrength)
		poleLocation = nm.array(poleLocation)
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

	def makeMatrix2DExactIntegrals(self,obspts,nodes):
		if nodes.shape[1] != 2:
			raise ValueError('Requires two dimensional inputs.')
		h2, nodesminus1, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1 = self.makeNodeCoeffs(nodes)
		if len(obspts.shape) == 2:
			npts = obspts.shape[0]
		elif len(obspts.shape) == 1:
			npts = 1
		else:
			raise ValueError('Observation points must be an n by 2 matrix.')
		N = nodes.shape[0]
		mat = nm.zeros((2*npts,2*N))
		for k in range(0,npts):
			if npts == 1:
				op = obspts
			else:	
				op=obspts[k,:]
			#calculate the 'k' integrals
			a0, a1, b0, b1, c0, c1, d1, D, r2 = self.makeQuadCoeffs(op, nodesminus1,difnod,h2)
			int2 = self.Int2(d1,h2,1,D,r2) - self.Int2(d1,h2,0,D,r2)
			int4 = self.Int4(d1,h2,1,D,r2) - self.Int4(d1,h2,0,D,r2)
			int5 = self.Int5(d1,h2,1,D,r2) - self.Int5(d1,h2,0,D,r2)
			int6 = self.Int6(d1,h2,1,D,r2) - self.Int6(d1,h2,0,D,r2)
			# calculate the 'k+1' integrals
			a0p1, a1p1, b0p1, b1p1, c0p1, c1p1, d1p1, Dp1, r2p1 = self.makeQuadCoeffs(op, nodes,difnodp1,h2)
			int1p1 = self.Int1(d1p1,h2,1,Dp1,r2p1) - self.Int1(d1p1,h2,0,Dp1,r2p1)
			int2p1 = self.Int2(d1p1,h2,1,Dp1,r2p1) - self.Int2(d1p1,h2,0,Dp1,r2p1)
			int3p1 = self.Int3(d1p1,h2,1,Dp1) - self.Int3(d1p1,h2,0,Dp1)
			int4p1 = self.Int4(d1p1,h2,1,Dp1,r2p1) - self.Int4(d1p1,h2,0,Dp1,r2p1)
			int5p1 = self.Int5(d1p1,h2,1,Dp1,r2p1) - self.Int5(d1p1,h2,0,Dp1,r2p1)
			int6p1 = self.Int6(d1p1,h2,1,Dp1,r2p1) - self.Int6(d1p1,h2,0,Dp1,r2p1)
			#put it altogether to make the matrix
			mat[2*k,range(0,2*N,2)] = -int2 + int4*(self.eps**2 + a0) + a1*int5 + a2*int6 -int1p1 + int2p1 + int3p1*(self.eps**2 + a0p1) + int4p1*(a1p1 -self.eps**2 - a0p1) + int5p1*(a2p1-a1p1) - a2p1*int6p1
			mat[2*k+1,range(1,2*N,2)] = -int2 + int4*(self.eps**2 + b0) + b1*int5 + b2*int6 -int1p1 + int2p1 + int3p1*(self.eps**2 + b0p1) + int4p1*(b1p1 -self.eps**2 - b0p1) + int5p1*(b2p1-b1p1) - b2p1*int6p1
			mat[2*k,range(1,2*N,2)] = c0*int4 + c1*int5 + c2*int6 + c0p1*int3p1 + int4p1*(c1p1 - c0p1) + int5p1*(c2p1 - c1p1) - c2p1*int6p1
			mat[2*k+1,range(0,2*N,2)] = mat[2*k,range(1,2*N,2)]
		return mat/(4*nm.pi*self.mu)

	def makeNodeCoeffs(self,nodes):
		h2 = nm.sum((nodes[0,:] - nodes[1,:])**2)
		nodesminus1 = nm.row_stack([nodes[-1,:],nodes[:-1,:]])
		difnod = nodesminus1 - nodes
		difnodp1 = nm.row_stack([difnod[1:,:],difnod[0,:]])
		a2, b2, c2 = self.makeNChelper(difnod)
		a2p1, b2p1, c2p1 = self.makeNChelper(difnodp1)
		return h2, nodesminus1, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1

	def makeNChelper(self,difnod):
		a2 = difnod[:,0]**2
		b2 = difnod[:,1]**2
		c2 = difnod[:,0]*difnod[:,1]
		return a2, b2, c2

	def makeQuadCoeffs(self,obspts,nodes,difnod,h2):
		difobs = obspts - nodes
		r2 = nm.sum(difobs**2,1)
		a0 = difobs[:,0]**2
		a1 = 2*difobs[:,0]*difnod[:,0]
		b0 = difobs[:,1]**2
		b1 = 2*difobs[:,1]*difnod[:,1]
		c0 = difobs[:,0]*difobs[:,1]
		c1 = difobs[:,0]*difnod[:,1] + difobs[:,1]*difnod[:,0] 
		d1 = a1 + b1 
		D = d1**2 - 4*h2*(r2+self.eps**2)
		return a0, a1, b0, b1, c0, c1, d1, D, r2

	def Afunc(self,d1,h2,s,D):
		return nm.arctan2(d1+2*h2*s,nm.sqrt(-D))

	def Lfunc(self,d1,h2,s,r2):
		return nm.log(self.eps**2 + r2 + d1*s + h2*s**2)

	def Int1(self,d1,h2,s,D,r2):
		return -s + nm.sqrt(-D)*self.Afunc(d1,h2,s,D)/(2*h2) + (d1 + 2*h2*s)*self.Lfunc(d1,h2,s,r2)/(4*h2)

	def Int2(self,d1,h2,s,D,r2):
		return 1./(8*h2**2)*( 2*h2*s*(d1-h2*s) - 2*d1*nm.sqrt(-D)*self.Afunc(d1,h2,s,D) + (-d1**2 + 2*h2*(self.eps**2+r2+h2*s**2))*self.Lfunc(d1,h2,s,r2) )

	def Int3(self,d1,h2,s,D):
		return 2*self.Afunc(d1,h2,s,D)/nm.sqrt(-D)

	def Int4(self,d1,h2,s,D,r2):
		return self.Lfunc(d1,h2,s,r2)/(2*h2) - d1*self.Afunc(d1,h2,s,D)/(h2*nm.sqrt(-D))

	def Int5(self,d1,h2,s,D,r2):
		return s/h2 - d1*self.Lfunc(d1,h2,s,r2)/(2*h2**2) + (d1**2 - 2*h2*(r2+self.eps**2))*self.Afunc(d1,h2,s,D)/(h2**2*nm.sqrt(-D))

	def Int6(self,d1,h2,s,D,r2):
		return -d1*s/h2**2 + s**2/(2*h2) + (d1**2-h2*(r2+self.eps**2))*self.Lfunc(d1,h2,s,r2)/(2*h2**3) - d1*(d1**2 - 3*h2*(r2 + self.eps**2))*self.Afunc(d1,h2,s,D) / (h2**3*nm.sqrt(-D))


def discretizeCircle(numNodes,circleRadius=1):
	'''Create boundary elements on a circle.'''
	angles = nm.arange(0,2*nm.pi,2*nm.pi/numNodes)
	nodes=nm.column_stack([circleRadius*nm.cos(angles),circleRadius*nm.sin(angles)])
	return nodes
	

def useGMRES(A,RHS):
	forces, info = gmres(A,RHS)
	if info == 0:
		pass
	elif info > 0:
		print('Maximum iterations reached: '+str(info))
	else:
		print('gmres encountered an error.')
	return forces
	
def constructpatches(spacing=0.005):
	farr = nm.arange(0.48,0.5+spacing/2.,spacing)
	nearr = nm.arange(0.255/nm.sqrt(2),0.255/nm.sqrt(2) + 0.0200001,spacing)
	n=farr.shape[0]
	ofar=nm.zeros((n**2,2))
	onear=nm.zeros((n**2,2))
	for j in range(n):
		ofar[j*n:(j+1)*n,0] = farr[j]*nm.ones(n)
		ofar[j*n:(j+1)*n,1] = farr
		onear[j*n:(j+1)*n,0] = nearr[j]*nm.ones(n)
		onear[j*n:(j+1)*n,1] = nearr
	return nm.row_stack([ofar,onear])

def constructfullgrid(spacing=0.01):
	r = nm.arange(-0.5,0.5+spacing/2.,spacing)
	n=r.shape[0]
	pts=nm.zeros((n**2,2))
	for j in range(n):
		pts[j*n:(j+1)*n,0] = r[j]*nm.ones(n)
		pts[j*n:(j+1)*n,1] = r
	return pts
	
	
def exactSoln(obspts,rad=0.25,mu=1.0,BC=(1,0)):
	'''x-values in col 0, y-vals in col 1'''
	BC = nm.array(BC)
	flag = len(obspts.shape)
	if flag ==2:
		r2 = nm.sum(obspts**2,1)
	elif flag ==1:
		r2 = nm.sum(obspts**2)
	else:
		raise ValueError('2D only bro.')
	r = nm.sqrt(r2)
	f0 = (8*nm.pi/(1-2*nm.log(rad)))*BC
	t1 = (2*nm.log(r) - rad**2/r2)/(8*nm.pi*mu)
	dp = nm.dot(obspts,f0)
	t2 = dp*(1 - rad**2/r2)/(4*nm.pi*mu*r2)
	if flag == 2:
		u = -f0[0]*t1 + obspts[:,0]*t2
		v = -f0[1]*t1 + obspts[:,1]*t2
	else:
		u = -f0[0]*t1 + obspts[0]*t2
		v = -f0[1]*t1 + obspts[1]*t2
	return u, v
		

if __name__ == '__main__':
	
	# single blob for poster, varying eps
	mu = 1.0
	rad = 0.25
	farp = 0.5*nm.array([nm.cos(0.7),nm.sin(0.7)])
	nearp = 0.251*nm.array([nm.cos(0.7),nm.sin(0.7)])
	obspts = nm.asarray([farp,nearp])
	exactu, exactv = exactSoln(obspts,rad,mu)
	numchords = [2**k for k in range(4,8)]
	for N in numchords:
		nodes = discretizeCircle(N,rad)
		BCs = nm.asarray(N*[1,0])
		origblobparams = [rad*2*nm.sin(nm.pi/N)/4]
		eintblobparams = [(rad*2*nm.sin(nm.pi/N)/(3./2))**2]
		u = nm.empty((0,))
		v = nm.empty((0,))
		origu = nm.empty((0,))
		origv = nm.empty((0,))
		for k in range(len(origblobparams)):
			print('N = ' + str(N)+ ', Blob number '+str(k))
			blob = Stokeslet2DReg2(eintblobparams[k],mu)
			forcemat = blob.makeMatrix2DWithIntegrals(nodes,nodes)
			forces = useGMRES(forcemat,BCs)
			velmat = blob.makeMatrix2DWithIntegrals(obspts,nodes)
			vel = nm.dot(velmat,forces)
			origblob = Stokeslet2DReg2(origblobparams[k],mu)
			origforcemat = origblob.makeMatrix2DOriginal(nodes,nodes)
			origforces = useGMRES(origforcemat,BCs)
			origvelmat = origblob.makeMatrix2DOriginal(obspts,nodes)
			origvel = nm.dot(origvelmat,origforces)
			if len(u) == 0:
				u = vel[0::2]                   	
				v = vel[1::2]                    	
				origu = origvel[0::2]                   	
				origv = origvel[1::2]                    	
			else:                                                     	
				u = nm.append(u,vel[0::2])	
				v = nm.append(v,vel[1::2])	
				origu = nm.append(origu,origvel[0::2])	
				origv = nm.append(origv,origvel[1::2])	
		fname = 'StokesCylTest_exactintegralsANDorig_2points_singleblob4poster_varyingeps_N'+ (4-len(str(N)))*'0' + str(N) + '.mat'                	
		mf.write(fname,{'u':u,'v':v,'origu':origu,'origv':origv,'exactu':exactu,'exactv':exactv,'obspts':obspts,'origbps':origblobparams,'eintbps':eintblobparams})
		
	
	# # many blobs
	# mu = 1.0
	# rad = 0.25
	# farp = 0.5*nm.array([nm.cos(0.7),nm.sin(0.7)])
	# nearp = 0.251*nm.array([nm.cos(0.7),nm.sin(0.7)])
	# obspts = nm.asarray([farp,nearp])
	# exactu, exactv = exactSoln(obspts,rad,mu)
	# numchords = [2**k for k in range(3,9)]
	# for N in numchords:
	# 	nodes = discretizeCircle(N,rad)
	# 	BCs = nm.asarray(N*[1,0])
	# 	origblobparams = nm.linspace(2*nm.pi*rad/(50*N),2*(2*nm.pi*rad/(N)),100)
	# 	origu = nm.empty((0,))
	# 	origv = nm.empty((0,))
	# 	for k in range(len(origblobparams)):
	# 		print('N = ' + str(N)+ ', Blob number '+str(k))
	# 		origblob = Stokeslet2DReg2(origblobparams[k],mu)
	# 		origforcemat = origblob.makeMatrix2DOriginal(nodes,nodes)
	# 		origforces = useGMRES(origforcemat,BCs)
	# 		origvelmat = origblob.makeMatrix2DOriginal(obspts,nodes)
	# 		origvel = nm.dot(origvelmat,origforces)
	# 		if len(origu) == 0:
	# 			origu = origvel[0::2]                   	
	# 			origv = origvel[1::2]                    	
	# 		else:                                                     	
	# 			origu = nm.append(origu,origvel[0::2])	
	# 			origv = nm.append(origv,origvel[1::2])	
	# 	fname = 'StokesCylTest_origmeth_2points_newblobs_N'+ (4-len(str(N)))*'0' + str(N) + '.mat'                	
	# 	mf.write(fname,{'origu':origu,'origv':origv,'exactu':exactu,'exactv':exactv,'obspts':obspts,'origbps':origblobparams})
	# 	# eintblobparams = nm.linspace((2*nm.pi*rad/(75*N))**2,(2*nm.pi*rad/(N))**1.5,500)
	# 	# u = nm.empty((0,))
	# 	# v = nm.empty((0,))
	# 	# for k in range(len(eintblobparams)):
	# 	# 	print('N = ' + str(N)+ ', Blob number '+str(k))
	# 	# 	blob = Stokeslet2DReg2(eintblobparams[k],mu)
	# 	# 	forcemat = blob.makeMatrix2DWithIntegrals(nodes,nodes)
	# 	# 	forces = useGMRES(forcemat,BCs)
	# 	# 	velmat = blob.makeMatrix2DWithIntegrals(obspts,nodes)
	# 	# 	vel = nm.dot(velmat,forces)
	# 	# 	if len(u) == 0:
	# 	# 		u = vel[0::2]                   	
	# 	# 		v = vel[1::2]                    	
	# 	# 	else:                                                     	
	# 	# 		u = nm.append(u,vel[0::2])	
	# 	# 		v = nm.append(v,vel[1::2])				
	# 	# fname = 'StokesCylTest_exactintegrals_2points_newblobs_N'+ (4-len(str(N)))*'0' + str(N) + '.mat'                	
	# 	# mf.write(fname,{'u':u,'v':v,'exactu':exactu,'exactv':exactv,'obspts':obspts,'eintbps':eintblobparams})
	
