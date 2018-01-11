#!/usr/bin/env python

import numpy as np
import sys
from scipy.sparse.linalg import gmres
import mat2py
import matplotlib.pyplot as plt



def calcVelSinglePoints(x,xf,f,eps,mu):
	r = np.sqrt(np.sum((x-xf)**2))
	dx = x-xf
	H2 = 1./(8*np.pi*(r**2+eps**2)**(3/2.))
	H1 = (r**2 + 2*eps**2)*H2
	M = np.array([dx[0]*dx,dx[1]*dx,dx[2]*dx]) 
	M = H2*M
	M = M + H1*np.eye(3)
	u = (1./mu) * np.dot(M,f)
	return u
	
def calcVelMultiplePointsOneEps(xu,xf,f,eps,mu):
	
	# deal with 1D vs 2D arrays
	xu, P = arrayDims(xu)
	xf, N = arrayDims(xf)
			
	# build each row
	R = np.zeros((3,3*N))
	u = np.zeros((P,3))
	for k in range(P):
		d = xu[k,:] - xf
		d2 = d**2
		r2 = np.sum(d2,1)
		H2 = 1./(8*np.pi*(r2+eps**2)**(3/2.))
		H1 = (r2 + 2*eps**2)*H2
		R[0,0::3] = H1 + d2[:,0]*H2
		R[0,1::3] = d[:,0]*d[:,1]*H2
		R[0,2::3] = d[:,0]*d[:,2]*H2
		R[1,0::3] = d[:,1]*d[:,0]*H2
		R[1,1::3] = H1 + d2[:,1]*H2
		R[1,2::3] = d[:,1]*d[:,2]*H2
		R[2,0::3] = d[:,2]*d[:,0]*H2
		R[2,1::3] = d[:,2]*d[:,1]*H2
		R[2,2::3] = H1 + d2[:,2]*H2
		u[k,:] = np.dot(R,f.flatten())/mu
	return u

def calcVelMultiplePointsTwoEps(xu,xf,f,eps,mu,epsind):
	# deal with 1D vs 2D arrays
	xu, P = arrayDims(xu)
	xf, N = arrayDims(xf)
			
	# build each row
	R = np.zeros((3,3*N))
	u = np.zeros((P,3))
	for k in range(P):
		d = xu[k,:] - xf
		d2 = d**2
		r2 = np.sum(d2,1)
		H2eps0 = 1./(8*np.pi*(r2[:epsind]+eps[0]**2)**(3/2.))
		H2eps1 = 1./(8*np.pi*(r2[epsind:]+eps[1]**2)**(3/2.))
		H1 = np.append((r2[:epsind] + 2*eps[0]**2)*H2eps0,(r2[epsind:] + 2*eps[1]**2)*H2eps1)
		H2 = np.append(H2eps0,H2eps1)
		R[0,0::3] = H1 + d2[:,0]*H2
		R[0,1::3] = d[:,0]*d[:,1]*H2
		R[0,2::3] = d[:,0]*d[:,2]*H2
		R[1,0::3] = d[:,1]*d[:,0]*H2
		R[1,1::3] = H1 + d2[:,1]*H2
		R[1,2::3] = d[:,1]*d[:,2]*H2
		R[2,0::3] = d[:,2]*d[:,0]*H2
		R[2,1::3] = d[:,2]*d[:,1]*H2
		R[2,2::3] = H1 + d2[:,2]*H2
		u[k,:] = np.dot(R,f.flatten())/mu
	return u
	
def calcForceMultiplePoints(xu,u,xf,eps,mu,epsind=-1):
	'''Calculate the velocity at the 3D points in the matrix x, caused by the forces f at the points xf. x is px3, xf and f are nx3, eps and mu are scalars. If eps contains two values, then epsind is the index at which to switch from eps[0] to eps[1]. Default value is no switch (epsind = -1).'''
	if epsind == -1:
		M = makeMatrixOneEps(xu,xf,eps)
	else:
		M = makeMatrixTwoEps(xu,xf,eps,epsind)
	f = mu*gmres(M,u.flatten())
	return f.reshape(xf.shape)

def makeMatrixOneEps(x,xf,eps):
	'''Make the Stokeslet matrix.'''
	
	# deal with 1D vs 2D arrays
	x, P = arrayDims(x)
	xf, N = arrayDims(xf)
			
	M = np.zeros((3*P,3*N))
			
	# build the matrix	
	for k in range(P):
		ind = 3*k
		d = xu[k,:] - xf
		d2 = d**2
		r2 = np.sum(d2,1)
		H2 = 1./(8*np.pi*(r2+eps**2)**(3/2.))
		H1 = (r2 + 2*eps**2)*H2
		M[ind,0::3] =   H1 + d2[:,0]*H2
		M[ind,1::3] =   d[:,0]*d[:,1]*H2
		M[ind,2::3] =   d[:,0]*d[:,2]*H2
		M[ind+1,0::3] = d[:,1]*d[:,0]*H2
		M[ind+1,1::3] = H1 + d2[:,1]*H2
		M[ind+1,2::3] = d[:,1]*d[:,2]*H2
		M[ind+2,0::3] = d[:,2]*d[:,0]*H2
		M[ind+2,1::3] = d[:,2]*d[:,1]*H2
		M[ind+2,2::3] = H1 + d2[:,2]*H2
	return M

def makeMatrixTwoEps(x,xf,eps,epsind):
	'''Make the Stokeslet matrix, allowing for two values of the spread parameter epsilon. eps is a list of two elements and epsind is the index where we switch eps.'''
	
	# deal with 1D vs 2D arrays
	x, P = arrayDims(x)
	xf, N = arrayDims(xf)
	M = np.zeros((3*P,3*N))
	
	# build the matrix	
	for k in range(P):
		ind = 3*k
		d = xu[k,:] - xf
		d2 = d**2
		r2 = np.sum(d2,1)
		H2eps0 = 1./(8*np.pi*(r2[:epsind]+eps[0]**2)**(3/2.))
		H2eps1 = 1./(8*np.pi*(r2[epsind:]+eps[1]**2)**(3/2.))
		H1 = np.append((r2[:epsind] + 2*eps[0]**2)*H2eps0,(r2[epsind:] + 2*eps[1]**2)*H2eps1)
		H2 = np.append(H2eps0,H2eps1)
		M[ind,0::3] =   H1 + d2[:,0]*H2
		M[ind,1::3] =   d[:,0]*d[:,1]*H2
		M[ind,2::3] =   d[:,0]*d[:,2]*H2
		M[ind+1,0::3] = d[:,1]*d[:,0]*H2
		M[ind+1,1::3] = H1 + d2[:,1]*H2
		M[ind+1,2::3] = d[:,1]*d[:,2]*H2
		M[ind+2,0::3] = d[:,2]*d[:,0]*H2
		M[ind+2,1::3] = d[:,2]*d[:,1]*H2
		M[ind+2,2::3] = H1 + d2[:,2]*H2
	return M
		
def arrayDims(v):
	'''Deal with 1D vs 2D arrays.'''
	if len(v.shape) == 1:
		v = np.array([v])
		P = 1
	else:
		P = v.shape[0]
	return v, P		
	
def specifyForces(U,l,eps,mu):
	'''Specifies the force magnitude for the two point model. U is velocity magnitude, a vector with one magnitude for each organism. l is the prescribed distance between the head and tail points (scalar, the same for all organisms). eps has two elements, a blob parameter for the head and a larger blob parameter for the tail in that order. mu is dynamic viscosity of the fluid. The output is the required force magnitude for the head and tail to maintain speed U in isolation. (The force vectors are equal and opposite at the head and tail.)'''
	F = 4*np.pi*mu*U/(1./eps[0] - 1./np.sqrt(l**2 + eps[1]**2))
	N = len(U)
	F.shape = (N,1)
	return F
	
	
def calcFlow(xu,xh,v,U,l,eps,mu):
	'''Calculates the velocity at a group of points x (Px3 ndarray) in 3D from a group of organisms with heads located at xh (Nx3 ndarray), body orientations v (Nx3 ndarray), and speeds (in isolation) of U (Nx1 ndarray, all positive). All organisms have the same body length l (scalar). mu is fluid dynamic viscosity (scalar). eps has two elements (list or ndarray), a blob parameter for the head first and then a larger one for the tail.'''
	xt = xh - l*v
	xf = np.row_stack([xh,xt])
	F = specifyForces(U,l,eps,mu)
	fh = F*v 
	f = np.row_stack([fh,-fh])
	epsind = xf.shape[0]/2
	u = calcVelMultiplePointsTwoEps(xu,xf,f,eps,mu,epsind)
	return u, xt
	
def scenario1pt():
	pass
	
def scenario3pts():
	l = 1.e-5
	U = np.array([6.e-5,4.5e-5,5.e-5])
	mu = 1.e-3
	eps = [l/2.,3*l/2.]
	v = -1*np.array([[np.cos(np.pi/6),np.sin(np.pi/6),0],[np.cos(5*np.pi/6),np.sin(5*np.pi/6),0],[np.cos(-np.pi/2),np.sin(-np.pi/2),0]])
	xh = -4.e-6*v
	x = np.arange(-1.e-3,1.01e-3,1.e-5)
	X,Y = np.meshgrid(x,x)
	xu = np.column_stack([np.column_stack([X.flatten(),Y.flatten()]),np.zeros(X.flatten().shape)])
	u, xt = calcFlow(xu,xh,v,U,l,eps,mu)
	yu = xu[:,1].reshape((len(x),len(x)))
	xu = xu[:,0].reshape((len(x),len(x)))
	v = u[:,1].reshape((len(x),len(x)))
	u = u[:,0].reshape((len(x),len(x)))
	# plotVels(xu,yu,u,v,xh,xt)
	mat2py.write('/Users/bcummins/colonyvel.mat', {'xu':xu,'yu':yu,'xh':xh,'xt':xt,'u':u,'v':v,'eps':eps})
	return xu,yu,u,v,xh,xt
	
def plotVels(xu,yu,u,v,xh,xt):
	f = plt.figure(1)
	plt.plot(xh[:,0],xh[:,1],'r.')
	plt.plot(xt[:,0],xt[:,1],'r.')
	plt.quiver(xu,yu,u,v) #for use with arrays
	# plt.gca().set_aspect('equal')
	f.canvas.draw()
	plt.savefig('/Users/bcummins/colonyvel.png')

	