#!/usr/bin/env python

import numpy as nm
import os
import sys
from scipy.sparse.linalg import gmres

####################################################################
# Parent class for regularized Stokeslets						   #
####################################################################
class genericStokeslet:
	def __init__(self,eps,mu,dim,bdrytype):
		'''Takes a spread parameter, a fluid viscosity, a domain dimension (2 or 3), and a boundary type ("closed curve" or "open curve").'''
		self.eps = eps
		self.mu = mu
		self.dim = dim
		if bdrytype == 'closed curve' or bdrytype == 'open curve':
			self.bdrytype = bdrytype
		else:
			raise ValueError('Boundary type not recognized. Use "closed curve" or "open curve".')
		
	def _checkSize(self,mat):
		if len(mat.shape) == 1:
			if len(mat) != self.dim:
				raise ValueError('Array is the wrong size for dimension %d.' % dim)
			N = 1
		elif mat.shape[1] != self.dim:
			raise ValueError('Array is the wrong size for dimension %d.' % dim)
		else:
			N = mat.shape[0]
		return N
		
	def _acceptArgs(self,velpts,forcepts,vels,forces):
		M = self._checkSize(velpts)
		N = self._checkSize(forcepts)
		if vels != None:
			Q = self._checkSize(vels)
			if Q != M:
				raise ValueError('Size mismatch between velocities and their locations.')
		if forces != None:
			Q = self._checkSize(forces)
			if Q != N:
				raise ValueError('Size mismatch between forces and their locations.')
		return M, N
		
	def _getPoint(self,k,pts):
		if len(pts.shape) == 1 and k==0:
			return pts
		elif k == -1:
			if self.bdrytype == 'closed curve':
				return pts[-1,:]
			else:
				return None
		elif k == len(pts):
			if self.bdrytype == 'closed curve':
				return pts[0,:]
			else:
				return None
		else:
			return pts[k,:]
			
	def _useGMRES(self,A,RHS):
		forces, info = gmres(A,RHS)
		if info == 0:
			pass
		elif info > 0:
			print('Maximum iterations reached: '+str(info))
		else:
			print('gmres encountered an error.')
		return forces
			
	def Forces2Vels(self,velpts,forcepts,forces):
		M, N = self._acceptArgs(velpts,forcepts,None,forces)
		vels = nm.zeros((self.dim*M,))
		for k in range(M):
			pt = self._getPoint(k,velpts)
			rows = self._buildRows(pt,forcepts,N)
			vels[self.dim*k:self.dim*(k+1)] = (1./self.mu)*nm.sum(rows*forces.flatten(),1)
		return vels.reshape((M,self.dim))
							
	def Vels2Forces(self,forcepts,velpts,vels):
		M, N = self._acceptArgs(velpts,forcepts,vels,None)
		kernel = nm.zeros((self.dim*M,self.dim*N))
		for k in range(M):
			pt = self._getPoint(k,velpts)
			kernel[self.dim*k:self.dim*(k+1),:] = self._buildRows(pt,forcepts,N)
		forces = self.mu*self._useGMRES(kernel,vels.flatten())
		return forces.reshape((N,self.dim))


		
####################################################################
# 2D cubic Stokeslet, linear BEM								   #
####################################################################
class cubic2DStokeslet_linearBEM(genericStokeslet):	
	def __init__(self,eps,mu,bdrytype):
		genericStokeslet.__init__(self,eps,mu,2,bdrytype)
		
	def _buildRows(self,pt,forcepts,N):
		rows = nm.zeros((self.dim,self.dim*N))
		for k in range(0,N):
			ptk = self._getPoint(k,forcepts)
			ptkm = self._getPoint(k-1,forcepts)
			ptkp = self._getPoint(k+1,forcepts)
			twobytwo = nm.zeros((2,2))
			if ptkm != None:
				h = nm.sqrt(nm.dot(ptk-ptkm,ptk-ptkm))
				twobytwo = twobytwo + h*self._rowHelper(pt,ptk,ptkm)
			if ptkp != None:
				h = nm.sqrt(nm.dot(ptk-ptkp,ptk-ptkp))
				twobytwo = twobytwo + h*self._rowHelper(pt,ptk,ptkp)
			rows[:,self.dim*k:self.dim*(k+1)] = twobytwo
		rows = (1./(4*nm.pi))*rows
		return rows
		
	def _rowHelper(self,pt,pt1,pt2):
		A, B, C, a, b, c = self._coeffs(pt,pt1,pt2)
		# print('A = %s' % str(A)); print('B = %s' % str(B)); print('C = %s' % str(C))
		I1, I2, I3, I4 = self._integrals(a,b,c)
		# if nm.isnan(I1):
		# 	print('a = %f' % a); print('b = %f' % b); print('c = %f' % c) 
		# 	print('I1 = %f' % I1); print('I2 = %f' % I2); print('I3 = %f' % I3); print('I4 = %f' % I4)
		# 	sys.exit()
		return (-I1 + (self.eps**2)*I2)*nm.eye(self.dim) + C*I2 + B*I3 + A*I4
		
	def _coeffs(self,xj,xk,xk1):
		adiff = xk1 - xk
		cdiff = xj - xk1
		A = self._coeffmat(adiff)
		C = self._coeffmat(cdiff)
		B = nm.zeros((2,2))
		B[0,0] = 2*adiff[0]*cdiff[0]
		B[0,1] = cdiff[0]*adiff[1] + cdiff[1]*adiff[0]
		B[1,0] = B[0,1]
		B[1,1] = 2*adiff[1]*cdiff[1]
		a = nm.sum(adiff**2)
		b = 2*nm.dot(adiff,cdiff)
		c = nm.sum(cdiff**2) + self.eps**2
		return A, B, C, a, b, c
		
	def _coeffmat(self,dif):
		mat = nm.zeros((2,2))
		mat[0,0] = dif[0]**2
		mat[0,1] = dif[0]*dif[1]
		mat[1,0] = mat[0,1]
		mat[1,1] = dif[1]**2
		return mat
		
	def _integrals(self,a,b,c):
		sD = nm.sqrt(-b**2 + 4*a*c)
		L0 = nm.log(c)
		L1 = nm.log(a+b+c)
		A0 = nm.arctan2(b,sD)
		A1 = nm.arctan2(b+2*a,sD)
		I1 = self._I1(1,A1,L1,a,b,c,sD) - self._I1(0,A0,L0,a,b,c,sD)
		I2 = self._I2(1,A1,L1,a,b,c,sD) - self._I2(0,A0,L0,a,b,c,sD)
		I3 = self._I3(1,A1,L1,a,b,c,sD) - self._I3(0,A0,L0,a,b,c,sD)
		I4 = self._I4(1,A1,L1,a,b,c,sD) - self._I4(0,A0,L0,a,b,c,sD)
		return I1, I2, I3, I4
		
	def _I1(self,s,A,L,a,b,c,sD):
		return -b*sD*A/(4*(a**2)) + (-b**2 + 2*a*c + 2*(a**2)*(s**2))*L/(8*(a**2)) + s*(b-a*s)/(4*a)
		
	def _I2(self,s,A,L,a,b,c,sD):
		return -b*A/(a*sD) + L/(2*a)
		
	def _I3(self,s,A,L,a,b,c,sD):
		return (b**2 - 2*a*c)*A/((a**2)*sD) - b*L/(2*(a**2)) + s/a
		
	def _I4(self,s,A,L,a,b,c,sD):
		return -b*(b**2 - 3*a*c)*A/(a**3*sD) + (b**2 - a*c)*L/(2*(a**3)) + s*(-2*b + a*s)/(2*(a**2))
	


		
####################################################################
# 2D cubic Stokeslet, original method							   #
####################################################################
class cubic2DStokeslet_orig(genericStokeslet):
	def __init__(self,eps,mu,bdrytype):
		genericStokeslet.__init__(self,eps,mu,2,bdrytype)
		
	def _buildRows(self,pt,forcepts,N):
		rows = nm.zeros((self.dim,self.dim*N))
		dif = pt - forcepts
		if dif.ndim > 1:
			r2 = nm.sum(dif**2,1) + self.eps**2
			xdiff = dif[:,0]
			ydiff = dif[:,1]
		else:
			r2 = nm.sum(dif**2) + self.eps**2
			xdiff = dif[0]
			ydiff = dif[1]
		H1 = 1./(8*nm.pi*self.mu)*(2*self.eps**2/r2 - nm.log(r2))
		H2 = 1./(8*nm.pi*self.mu)*(2/r2)
		rows[0,2*nm.arange(N)] = H1 + (xdiff**2)*H2
		rows[0,2*nm.arange(N)+1] = (xdiff*ydiff)*H2
		rows[1,2*nm.arange(N)+1]= H1 + (ydiff**2)*H2
		rows[1,2*nm.arange(N)] = rows[0,2*nm.arange(N)+1]
		return rows
	
	
	
	

	
####################################################################
# Unit tests							  						   #
####################################################################	
			
def unittest_2Dcircle_exact(eps,mu,a,testpoints):
	#calculate exact solution (for BCs of (1,0))
	f = nm.array([1./(1 - 2*nm.log(a)),0])
	x = testpoints
	if len(x.shape) > 1:
		r = nm.sqrt(nm.sum(x**2,1))
		r = r[:,nm.newaxis]
		fdxx = nm.dot(f,x)[:,nm.newaxis]*x
	else:
		r = nm.sqrt(nm.sum(x**2))
		fdxx = nm.dot(f,x)*x
	return -f*(2*nm.log(r) - a**2/r**2) + 2*(fdxx/r**2)*(1-a**2/r**2)

def unittest_2Dcircle(eps,mu,a,Nc,testpoints,blobtype):
	#make circle (a = radius, Nc = number of discretization points)
	angs = nm.arange(0,2*nm.pi,2*nm.pi/Nc)
	forcepts=nm.column_stack([a*nm.cos(angs),a*nm.sin(angs)])
	h = nm.sqrt(2*a**2*(1-nm.cos(2*nm.pi/Nc)))
	BCs = nm.zeros(forcepts.shape)
	BCs[:,0] = 1
	if eps > 0:
		print('The ratio of discretization to regularization parameters is %f' % (h/eps))
	else:
		print('Epsilon = 0 -- Stokeslet method.')
	exactsoln = unittest_2Dcircle_exact(eps,mu,a,testpoints)
	print('Exact solution is ')
	print(exactsoln)
	#approximate with regularized Stokeslets 
	blob = blobtype(eps,mu,'closed curve')
	forces = blob.Vels2Forces(forcepts,forcepts,BCs)
	uapprox = blob.Forces2Vels(testpoints,forcepts,forces)
	print('Approximate solution is ')
	print(uapprox)	
	return exactsoln, uapprox

	
def unittest_2Dcircle_forwardonly(eps,mu,a,Nc,testpoints):
	#make circle (a = radius, Nc = number of discretization points)
	angs = nm.arange(0,2*nm.pi,2*nm.pi/Nc)
	forcepts=nm.column_stack([a*nm.cos(angs),a*nm.sin(angs)])
	forces=nm.zeros(forcepts.shape)
	forces[:,0] = 4
	exactsoln = unittest_2Dcircle_exact(eps,mu,a,testpoints)
	print('Exact solution is ')
	print(exactsoln)
	#approximate with regularized BEM Stokeslets 
	blob = cubic2DStokeslet_linearBEM(eps,mu,'closed curve')
	uapprox = blob.Forces2Vels(testpoints,forcepts,forces)
	print('Approximate solution is ')
	print(uapprox)	
	return exactsoln, uapprox

	
	
if __name__=='__main__':
	# #testing....
	# mu = 1; N = 4; a = 1; h = a*2*nm.pi/N; eps = 0.1
	# unittest_2Dcircle(eps,mu,a,N,nm.array([1,0]))	
	
	#compare with Mike
	mu = 1; a = 1; eps = 0
	pts = 1.01*nm.array([nm.cos(0.7),nm.sin(0.7)])
	N = [32,64,128,256];
	xerr = [];
	yerr = [];
	for n in N:
		es, ua = unittest_2Dcircle_forwardonly(eps,mu,a,n,pts)
		xerr.append( nm.abs(es[range(0,len(es),2)]-ua.flatten()[range(0,len(es),2)]) )
		yerr.append( nm.abs(es[range(1,len(es),2)]-ua.flatten()[range(1,len(es),2)]) )
	print(xerr)
	print(yerr)	
	
	# #convergence....
	# mu = 1; a = 1; eps = 1.e-5
	# pt0 = 1.000*nm.array([nm.cos(nm.pi/4),nm.sin(nm.pi/4)])
	# pt1 = 1.001*nm.array([nm.cos(nm.pi/4),nm.sin(nm.pi/4)])
	# pt2 = 1.010*nm.array([nm.cos(nm.pi/4),nm.sin(nm.pi/4)])
	# pt3 = 1.100*nm.array([nm.cos(nm.pi/4),nm.sin(nm.pi/4)])
	# pt4 = 1.250*nm.array([nm.cos(nm.pi/4),nm.sin(nm.pi/4)])
	# pts = [pt0,pt1,pt2,pt3,pt4]
	# N = [16,32,64,128,256,512];
	# err = [];
	# for pt in pts:
	# 	print(pt)
	# 	for n in N:
	# 		es, ua = unittest_2Dcircle(eps,mu,a,n,pt)
	# 		err.append(nm.sqrt(nm.sum((es-ua.flatten())**2)))
	# print(err)
	# err = nm.asarray(err)
	# err = nm.reshape(err,(len(pts),len(N)))
	# print(err)
	# 
	# import mat2py
	# mat2py.write('/Users/bcummins/scratch/testcircle.mat',{'pts':pts,'N':N,'err':err,'mu':mu,'a':a,'eps':eps})
		
	

		