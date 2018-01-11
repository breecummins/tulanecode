#!/usr/bin/env python

import numpy as nm
import os
import sys

from scipy.sparse.linalg import gmres, LinearOperator
import mat2py as mf

class CubicRegStokeslet():
	'''Regularized 2D Stokeslet, cube power. eps is the blob parameter, mu is dynamic fluid viscosity, nodes is the nx2 matrix of points representing the discretized boundary, and obspts is the mx2 matrix of points where we want to know the velocity. bdrytype = "closed", "open", or "surface" and specifies the boundary type: a closed curve, an open curve, or a 2D surface (may need to be differences between open and closed surfaces, not sure yet). This is needed for the linear BEM numerical scheme.'''
	def __init__(self,eps,mu,obspts,nodes,bdrytype):
		self.eps = eps		
		self.mu = mu	
		self.obspts = obspts
		self.nodes = nodes
		self.M = self.checkSize(obspts)
		self.N = self.checkSize(nodes)
		self.bdrytype = bdrytype
		# if self.M == self.N and nm.all(obspts == nodes):
		# 	self.matrixflag = 'integraleqn'
		# else:
		# 	self.matrixflag = 'evaluation'		
		
	def checkSize(self,mat):
		if len(mat.shape) == 1:
			if len(mat) != 2:
				raise ValueError('Coordinates must be 2D.')
			num = 1
		elif mat.shape[1] != 2:
			raise ValueError('Coordinates must be 2D.')
		else:
			num = mat.shape[0]
		return num
		
	def getPoint(self,k,pts):
		if pts.shape[0] == 1:
			return pts
		else:
			return pts[k,:]
		
	def kernValsOriginal(self,pt,flag=''):
		'''Build the matrix rows for the original Riemann sum method.'''
		dif = pt - self.nodes
		r2 = nm.sum(dif**2,1) + self.eps**2
		H1 = 2*self.eps**2/r2 - nm.log(r2)
		H2 = 2/r2
		xterm =  H1 + (dif[:,0]**2)*H2
		mixedterm = (dif[:,0]*dif[:,1])*H2
		yterm = H1 + (dif[:,1]**2)*H2
		row1 = nm.zeros((2*self.N,))
		row2 = nm.zeros((2*self.N,))
		row1[2*nm.arange(self.N)] = xterm/(8*nm.pi*self.mu)
		row1[2*nm.arange(self.N)+1] = mixedterm/(8*nm.pi*self.mu)
		row2[2*nm.arange(self.N)+1]= yterm/(8*nm.pi*self.mu)
		row2[2*nm.arange(self.N)] = mixedterm/(8*nm.pi*self.mu)
		if flag == 'trap' and self.bdrytype == 'open':
			row1[0:2] = 0.5*row1[0:2]; row1[-2:] = 0.5*row1[-2:]
			row2[0:2] = 0.5*row2[0:2]; row2[-2:] = 0.5*row2[-2:]
		return row1, row2
						
	def makeMatrixOriginal(self,flag=''):
		'''Build the regularized Stokeslet matrix explicitly. flag = "trap" means use the trapezoidal rule.'''
		mat = nm.zeros((2*self.M,2*self.N))
		for k in range(self.M):
			pt = self.getPoint(k,self.obspts)
			mat[2*k,:], mat[2*k+1,:] = self.kernValsOriginal(pt,flag)
		return mat
		
	def linOpOriginal(self,dt):
		'''Do not explicitly build the matrix; construct a linear operater for use with GMRES. dt is data type.'''
		return LinearOperator( (2*self.M, 2*self.N), matvec=self.linOpHelperOriginal,dtype=dt)
		
	def linOpHelperOriginal(self,vec):
		if len(vec.shape) > 1:
			vec = vec.flatten()
		output = nm.zeros((2*self.M,))
		for k in range(self.M):
			pt = self.getPoint(k,self.obspts)
			row1, row2 = self.kernValsOriginal(pt)
			output[2*k] = nm.sum(row1*vec)
			output[2*k+1] = nm.sum(row2*vec)
		return output	
		
##################################################################################################################################################################	
	# #The following algorithm is broken for open curves.  Not sure about closed curves.
	def kernValsExactIntegrals(self, k, h2, nodesminus1, nodesminus1sshift, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1):
		'''Build the matrix rows for the exact integral BEM. The first input is the point of interest, the others are the outputs from makeNodeCoeffs.'''
		pt = self.getPoint(k,self.obspts)
		N=self.N
		#calculate the 'k' integrals
		a0, a1, b0, b1, c0, c1, d1, D, r2 = self.makeQuadCoeffs(pt,nodesminus1,difnod,h2)
		int2 = self.Int2(d1,h2,1,D,r2) - self.Int2(d1,h2,0,D,r2)
		int4 = self.Int4(d1,h2,1,D,r2) - self.Int4(d1,h2,0,D,r2)
		int5 = self.Int5(d1,h2,1,D,r2) - self.Int5(d1,h2,0,D,r2)
		int6 = self.Int6(d1,h2,1,D,r2) - self.Int6(d1,h2,0,D,r2)
		# calculate the 'k+1' integrals
		a0p1, a1p1, b0p1, b1p1, c0p1, c1p1, d1p1, Dp1, r2p1 = self.makeQuadCoeffs(pt,nodesminus1sshift,difnodp1,h2)
		int1p1 = self.Int1(d1p1,h2,1,Dp1,r2p1) - self.Int1(d1p1,h2,0,Dp1,r2p1)
		int2p1 = self.Int2(d1p1,h2,1,Dp1,r2p1) - self.Int2(d1p1,h2,0,Dp1,r2p1)
		int3p1 = self.Int3(d1p1,h2,1,Dp1) - self.Int3(d1p1,h2,0,Dp1)
		int4p1 = self.Int4(d1p1,h2,1,Dp1,r2p1) - self.Int4(d1p1,h2,0,Dp1,r2p1)
		int5p1 = self.Int5(d1p1,h2,1,Dp1,r2p1) - self.Int5(d1p1,h2,0,Dp1,r2p1)
		int6p1 = self.Int6(d1p1,h2,1,Dp1,r2p1) - self.Int6(d1p1,h2,0,Dp1,r2p1)
		row1 = nm.zeros((2*N,))
		row2 = nm.zeros((2*N,))
		if self.bdrytype == 'closed':
			row1[2*nm.arange(N)] = -int2 + int4*(self.eps**2 + a0) + a1*int5 + a2*int6 -int1p1 + int2p1 + int3p1*(self.eps**2 + a0p1) + int4p1*(a1p1 -self.eps**2 - a0p1) + int5p1*(a2p1-a1p1) - a2p1*int6p1
			row1[2*nm.arange(N)+1] = c0*int4 + c1*int5 + c2*int6 + c0p1*int3p1 + int4p1*(c1p1 - c0p1) + int5p1*(c2p1 - c1p1) - c2p1*int6p1
			row2[2*nm.arange(N)+1]= -int2 + int4*(self.eps**2 + b0) + b1*int5 + b2*int6 -int1p1 + int2p1 + int3p1*(self.eps**2 + b0p1) + int4p1*(b1p1 -self.eps**2 - b0p1) + int5p1*(b2p1-b1p1) - b2p1*int6p1
			row2[2*nm.arange(N)] = row1[2*nm.arange(N)+1]
		elif self.bdrytype == 'open':
			row1[2*nm.arange(1,N)] = row1[2*nm.arange(1,N)] + -int2 + int4*(self.eps**2 + a0) + a1*int5 + a2*int6 
			row1[2*nm.arange(N-1)] = row1[2*nm.arange(N-1)] + -int1p1 + int2p1 + int3p1*(self.eps**2 + a0p1) + int4p1*(a1p1 -self.eps**2 - a0p1) + int5p1*(a2p1-a1p1) - a2p1*int6p1
			row1[2*nm.arange(1,N)+1] = row1[2*nm.arange(1,N)+1] + c0*int4 + c1*int5 + c2*int6 
			row1[2*nm.arange(N-1)+1] = row1[2*nm.arange(N-1)+1] + c0p1*int3p1 + int4p1*(c1p1 - c0p1) + int5p1*(c2p1 - c1p1) - c2p1*int6p1
			row2[2*nm.arange(1,N)+1] = row2[2*nm.arange(1,N)+1] -int2 + int4*(self.eps**2 + b0) + b1*int5 + b2*int6 
			row2[2*nm.arange(N-1)+1] = row2[2*nm.arange(N-1)+1] -int1p1 + int2p1 + int3p1*(self.eps**2 + b0p1) + int4p1*(b1p1 -self.eps**2 - b0p1) + int5p1*(b2p1-b1p1) - b2p1*int6p1
			row2[2*nm.arange(N)] = row1[2*nm.arange(N)+1]
		return row1/(4*nm.pi*self.mu), row2/(4*nm.pi*self.mu)
				
	def makeNodeCoeffs(self):
		h2 = nm.sum((self.nodes[0,:] - self.nodes[1,:])**2) #this assumes equally spaced nodes
		if self.bdrytype == 'closed':
			nodesminus1 = nm.row_stack([self.nodes[-1,:],self.nodes[:-1,:]])
			nodesminus1shift = self.nodes
			difnod = nodesminus1 - nodesminus1shift
			difnodp1 = nm.row_stack([difnod[1:,:],difnod[0,:]])
		elif self.bdrytype == 'open':
			nodesminus1 = self.nodes[:-1,:]
			nodesminus1shift = self.nodes[1:,:]
			difnod = nodesminus1 - nodesminus1shift #CHECK ME
			difnodp1 = difnod #CHECK ME
		a2, b2, c2 = self.makeNChelper(difnod)
		a2p1, b2p1, c2p1 = self.makeNChelper(difnodp1)			
		return h2, nodesminus1, nodesminus1shift, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1
	
	def makeNChelper(self,difnod):
		a2 = difnod[:,0]**2
		b2 = difnod[:,1]**2
		c2 = difnod[:,0]*difnod[:,1]
		return a2, b2, c2
	
	def makeQuadCoeffs(self,pt,ns,difnod,h2):
		difobs = pt - ns
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

#####################################################################################################################################################################	
	#Unbroken code resumes here
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

	def MultipleInts(self,d1,h2,D,r2,indexlist=range(1,7)):
		ints = nm.empty((0,))
		if 1 in indexlist:
			ints = nm.append(ints, self.Int1(d1,h2,1,D,r2) - self.Int1(d1,h2,0,D,r2) )
		if 2 in indexlist:
			ints = nm.append(ints, self.Int2(d1,h2,1,D,r2) - self.Int2(d1,h2,0,D,r2) )
		if 3 in indexlist:
			ints = nm.append(ints, self.Int3(d1,h2,1,D) - self.Int3(d1,h2,0,D)  )
		if 4 in indexlist:
			ints = nm.append(ints, self.Int4(d1,h2,1,D,r2) - self.Int4(d1,h2,0,D,r2) )
		if 5 in indexlist:
			ints = nm.append(ints, self.Int5(d1,h2,1,D,r2) - self.Int5(d1,h2,0,D,r2) )
		if 6 in indexlist:
			ints = nm.append(ints, self.Int6(d1,h2,1,D,r2) - self.Int6(d1,h2,0,D,r2) )
		return ints
	
	def makeMatrixExactIntegrals(self):
		'''Build the exact integral matrix explicitly.'''
		h2, nodesminus1, nodesminus1sshift, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1 = self.makeNodeCoeffs()
		mat = nm.zeros((2*self.M,2*self.N))
		for k in range(self.M):
			mat[2*k,:], mat[2*k+1,:] = self.kernValsExactIntegrals(k,h2, nodesminus1, nodesminus1sshift, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1)
		return mat
		
	def makeMatrixExactIntegralsSlow(self):
		'''Build one entry at a time.'''
		mat = nm.zeros((2*self.M,2*self.N))
		if self.bdrytype == 'open':
			for j in range(self.M):
				pt = self.getPoint(j,self.obspts)
				for k in range(self.N):
					if k > 0 and k < N-1:
						nodekm1 = self.getPoint(k-1,self.nodes)
						nodek = self.getPoint(k,self.nodes)
						nodekp1 = self.getPoint(k+1,self.nodes)
						#calculate the 'k-1' integrals
						a0m1, a1m1, a2m1, b0m1, b1m1, b2m1, c0m1, c1m1, c2m1, d1m1, Dm1, r2m1, h2m1 = self.makeQuadCoeffsSlow(pt,nodekm1,nodek)				
						intsm1 = self.MultipleInts(d1m1,h2m1,Dm1,r2m1,indexlist=[2,4,5,6])
						# calculate the 'k' integrals
						a0, a1, a2, b0, b1, b2, c0, c1, c2, d1, D, r2, h2 = self.makeQuadCoeffsSlow(pt,nodek,nodekp1)
						ints = self.MultipleInts(d1,h2,D,r2)
						#build block matrix
						M11 = -intsm1[0] + intsm1[1]*(self.eps**2 + a0m1) + a1m1*intsm1[2] + a2m1*intsm1[3] -ints[0] + ints[1] + ints[2]*(self.eps**2 + a0) + ints[3]*(a1 -self.eps**2 - a0) + ints[4]*(a2-a1) - a2*ints[5]
						M12 = c0m1*intsm1[1] + c1m1*intsm1[2] + c2m1*intsm1[3] + c0*ints[2] + ints[3]*(c1 - c0) + ints[4]*(c2 - c1) - c2*ints[5]
						M22 = -intsm1[0] + intsm1[1]*(self.eps**2 + b0m1) + b1m1*intsm1[2] + b2m1*intsm1[3] -ints[0] + ints[1] + ints[2]*(self.eps**2 + b0) + ints[3]*(b1 -self.eps**2 - b0) + ints[4]*(b2-b1) - b2*ints[5]
						M21 = M12
					elif k ==0: #when the boundary curve is open, we do not have the upward going integral before the first point
						nodek = self.getPoint(k,self.nodes)
						nodekp1 = self.getPoint(k+1,self.nodes)
						# calculate the 'k' integrals
						a0, a1, a2, b0, b1, b2, c0, c1, c2, d1, D, r2, h2 = self.makeQuadCoeffsSlow(pt,nodek,nodekp1)
						ints = self.MultipleInts(d1,h2,D,r2)
						#build block matrix
						M11 =  -ints[0] + ints[1] + ints[2]*(self.eps**2 + a0) + ints[3]*(a1 -self.eps**2 - a0) + ints[4]*(a2-a1) - a2*ints[5]
						M12 =  c0*ints[2] + ints[3]*(c1 - c0) + ints[4]*(c2 - c1) - c2*ints[5]
						M22 =  -ints[0] + ints[1] + ints[2]*(self.eps**2 + b0) + ints[3]*(b1 -self.eps**2 - b0) + ints[4]*(b2-b1) - b2*ints[5]
						M21 = M12
					elif k == N-1: #when the boundary curve is open, we do not have the downward going integral after the last point
						nodekm1 = self.getPoint(k-1,self.nodes)
						nodek = self.getPoint(k,self.nodes)
						#calculate the 'k-1' integrals
						a0m1, a1m1, a2m1, b0m1, b1m1, b2m1, c0m1, c1m1, c2m1, d1m1, Dm1, r2m1, h2m1 = self.makeQuadCoeffsSlow(pt,nodekm1,nodek)				
						intsm1 = self.MultipleInts(d1m1,h2m1,Dm1,r2m1,indexlist=[2,4,5,6])
						#build block matrix
						M11 = -intsm1[0] + intsm1[1]*(self.eps**2 + a0m1) + a1m1*intsm1[2] + a2m1*intsm1[3]
						M12 = c0m1*intsm1[1] + c1m1*intsm1[2] + c2m1*intsm1[3]
						M22 = -intsm1[0] + intsm1[1]*(self.eps**2 + b0m1) + b1m1*intsm1[2] + b2m1*intsm1[3]
						M21 = M12
					mat[nm.ix_([2*j,(2*j+1)],[2*k,(2*k+1)])] = nm.array([[M11, M12], [M21, M22]])
		elif self.bdrytype == 'closed':
			for j in range(self.M):
				pt = self.getPoint(j,self.obspts)
				for k in range(self.N):
					if k > 0 and k < N-1:
						nodekm1 = self.getPoint(k-1,self.nodes)
						nodek = self.getPoint(k,self.nodes)
						nodekp1 = self.getPoint(k+1,self.nodes)
					elif k == 0:  #when the boundary curve is closed, the -1st point = the N-1st point
						nodekm1 = self.getPoint(self.N-1,self.nodes)
						nodek = self.getPoint(k,self.nodes)
						nodekp1 = self.getPoint(k+1,self.nodes)
					elif k == N-1: #when the boundary curve is closed, the N-th point = the 0-th point
						nodekm1 = self.getPoint(k-1,self.nodes)
						nodek = self.getPoint(k,self.nodes)
						nodekp1 = self.getPoint(0,self.nodes)
					#calculate the 'k-1' integrals
					a0m1, a1m1, a2m1, b0m1, b1m1, b2m1, c0m1, c1m1, c2m1, d1m1, Dm1, r2m1, h2m1 = self.makeQuadCoeffsSlow(pt,nodekm1,nodek)				
					intsm1 = self.MultipleInts(d1m1,h2m1,Dm1,r2m1,indexlist=[2,4,5,6])
					# calculate the 'k' integrals
					a0, a1, a2, b0, b1, b2, c0, c1, c2, d1, D, r2, h2 = self.makeQuadCoeffsSlow(pt,nodek,nodekp1)
					ints = self.MultipleInts(d1,h2,D,r2)
					#build block matrix
					M11 = -intsm1[0] + intsm1[1]*(self.eps**2 + a0m1) + a1m1*intsm1[2] + a2m1*intsm1[3] -ints[0] + ints[1] + ints[2]*(self.eps**2 + a0) + ints[3]*(a1 -self.eps**2 - a0) + ints[4]*(a2-a1) - a2*ints[5]
					M12 = c0m1*intsm1[1] + c1m1*intsm1[2] + c2m1*intsm1[3] + c0*ints[2] + ints[3]*(c1 - c0) + ints[4]*(c2 - c1) - c2*ints[5]
					M22 = -intsm1[0] + intsm1[1]*(self.eps**2 + b0m1) + b1m1*intsm1[2] + b2m1*intsm1[3] -ints[0] + ints[1] + ints[2]*(self.eps**2 + b0) + ints[3]*(b1 -self.eps**2 - b0) + ints[4]*(b2-b1) - b2*ints[5]
					M21 = M12
					mat[nm.ix_([2*j,(2*j+1)],[2*k,(2*k+1)])] = nm.array([[M11, M12], [M21, M22]])
		return mat
				
	def makeQuadCoeffsSlow(self,pt,nodek,nodekp1):
		h2 = nm.sum((nodek-nodekp1)**2)
		r2 = nm.sum((pt-nodek)**2)
		a0 = (pt[0] - nodek[0])**2
		a1 = 2*(pt[0] - nodek[0])*(nodek[0]-nodekp1[0])
		a2 = (nodek[0]-nodekp1[0])**2
		b0 = (pt[1] - nodek[1])**2
		b1 = 2*(pt[1] - nodek[1])*(nodek[1]-nodekp1[1])
		b2 = (nodek[1]-nodekp1[1])**2
		c0 = (pt[0] - nodek[0])*(pt[1] - nodek[1])
		c1 = (pt[0] - nodek[0])*(nodek[1]-nodekp1[1]) + (pt[1] - nodek[1])*(nodek[0]-nodekp1[0])
		c2 = (nodek[0]-nodekp1[0])*(nodek[1]-nodekp1[1])
		d1 = a1 + b1 
		D = d1**2 - 4*h2*(r2+self.eps**2)
		return a0, a1, a2, b0, b1, b2, c0, c1, c2, d1, D, r2, h2
					

	def linOpExactIntegrals(self,dt):
		'''Do not explicitly build the matrix; construct a linear operater for use with GMRES.'''
		return LinearOperator( (2*self.M, 2*self.N), matvec=self.linOpHelperExactIntegrals,dtype=dt)

	def linOpHelperExactIntegrals(self,vec):
		if len(vec.shape) > 1:
			vec = vec.flatten()
		h2, nodesminus1, nodesminus1sshift, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1 = self.makeNodeCoeffs()
		output = nm.zeros((2*self.M,))
		for k in range(self.M):
			row1, row2 = self.kernValsExactIntegrals(k,h2, nodesminus1, nodesminus1sshift, difnod, a2, b2, c2, difnodp1, a2p1, b2p1, c2p1)
			output[2*k] = nm.sum(row1*vec)
			output[2*k+1] = nm.sum(row2*vec)
		return output			

def useGMRES(A,RHS,rs=None):
	forces, info = gmres(A,RHS,restart=rs)
	if info == 0:
		pass
	elif info > 0:
		print('Maximum iterations reached: '+str(info))
	else:
		print('gmres encountered an error.')
	return forces
	
	
	

if __name__ == '__main__':
	
	nm.set_printoptions(linewidth=200)
	
	# ##closed curves##
	# from StokesletCylinderExample_ExactIntegrals import discretizeCircle
	# N = 5
	# nodes = discretizeCircle(N,0.25)
	# dt=nodes.dtype
	# theta = (2*nm.pi/N)*nm.arange(N)
	# df = nm.column_stack([nm.cos(2*theta),nm.sin(2*theta)])
	# RHS = df.flatten()
	# blob = CubicRegStokeslet(theta[1]/4.,1,nodes,nodes,'closed') #to solve matrix equation, use nodes as both sets of points 
	# 
	# ##Original Method##
	# #see if both matrix and linear operator versions of gmres work
	# A = blob.makeMatrixOriginal()
	# print(A)
	# test1 = useGMRES(A,RHS)
	# A = blob.linOpOriginal(dt)
	# test2 = useGMRES(A,RHS)
	# print(nm.max(nm.abs(test1-test2)))
	# #now test against old code
	# import StokesletCylinderExample_ExactIntegrals as SCE
	# blob2 = SCE.Stokeslet2DReg2(theta[1]/4.,1)
	# A = blob2.makeMatrix2DOriginal(nodes,nodes)
	# test3 = useGMRES(A,RHS)
	# print(nm.max(nm.abs(test1-test3)))
	# print(nm.max(nm.abs(test2-test3)))
	# print(test1)
	# # print(abs(test1 - test2))
	# # print(abs(test2 - test3))
	# 
	# ##Success within roundoff##
	# 
	# ##Exact Integrals##
	# #see if both matrix and linear operator versions of gmres work in the original method
	# A = blob.makeMatrixExactIntegrals()
	# print(A)
	# test4 = useGMRES(A,RHS)
	# A = blob.linOpExactIntegrals(dt)
	# test5 = useGMRES(A,RHS)
	# print(nm.max(nm.abs(test4-test5)))
	# #now test against old code
	# import StokesletCylinderExample_ExactIntegrals as SCE
	# blob2 = SCE.Stokeslet2DReg2(theta[1]/4.,1)
	# A = blob2.makeMatrix2DExactIntegrals(nodes,nodes)
	# print(A)
	# test6 = useGMRES(A,RHS)
	# print(nm.max(nm.abs(test4-test6)))
	# print(nm.max(nm.abs(test5-test6)))
	# print(test4)
	# # print(abs(test4 - test5))
	# # print(abs(test5 - test6))
	# 
	# ##Success within roundoff##

	##open curves##
	from StokesletCylinderExample_ExactIntegrals import discretizeCircle
	dist=1.;
	num = 10;
	step = 2*dist/num
	for d in nm.arange(-dist,dist+step,step):
		if d == -dist:
			obspts = nm.column_stack([nm.linspace(-dist,dist,num), d*nm.ones(num)])
		else:
			obspts=nm.row_stack([obspts,nm.column_stack([nm.linspace(-dist,dist,num), d*nm.ones(num)])])
	# print(obspts)
	# import sys
	# sys.exit()
	
	N=10
	print('N = 10')
	nodes = nm.column_stack([nm.linspace(-N*0.05,N*0.05,N), nm.zeros(N)])
	dt=nodes.dtype
	df = nm.column_stack([nm.zeros(N),nm.ones(N)])
	RHS = df.flatten()
	eps = 0.005
	blob = CubicRegStokeslet(eps,1,nodes,nodes,'open') #to solve matrix equation, use nodes as both sets of points 
	Ao = blob.makeMatrixOriginal('trap')
	fo = useGMRES(Ao,RHS)
	blobeval = CubicRegStokeslet(eps,1,obspts,nodes,'open')
	Ao2 = blobeval.makeMatrixOriginal('trap')
	uo = nm.dot(Ao2,fo)
	uo = nm.reshape(uo,(uo.shape[0]/2,2))
	Aes = blob.makeMatrixExactIntegralsSlow()
	fes = useGMRES(Aes,RHS)
	Aes2 = blobeval.makeMatrixExactIntegralsSlow()
	ues = nm.dot(Aes2,fes)
	ue = nm.reshape(ues,(ues.shape[0]/2,2))
	# print(Ao)
	# print(fo)
	# print(Aes)
	# print(fes)
	print( nm.max( nm.sum((uo - ue)**2,1) / nm.sum(ue**2,1)) )
	print( nm.sum( nm.sum((uo - ue)**2,1) / nm.sum(ue**2,1)) )
	import mat2py
	mat2py.write('testopencurve.mat',{'obspts':obspts,'uo':uo,'ue':ue})
	
