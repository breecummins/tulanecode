#!/usr/bin/env python

import numpy as nm

# Slice the matrix. 
def xStack(u):
	'''Periodic BCs.'''
	um1 = nm.row_stack([u[-1,:],u[:-1,:]])
	up1 = nm.row_stack([u[1:,:],u[0,:]])
	return um1, up1

def yStack(u):
	'''Periodic BCs, cell edge values.'''
	um1 = nm.column_stack([u[:,-1],u[:,:-1]])
	up1 = nm.column_stack([u[:,1:],u[:,0]])
	return um1, up1

# Interpolate from nearest 4 points.	
def uInterp(u):
	'''Average nearest 4 u points to get middle v point, making use of periodic BCs. See MAC grid diagram in class docs.'''
	uxm1, uxp1 = xStack(u)
	uym1, uyp1 = yStack(u)
	uxp1ym1, uxp1yp1 = yStack(uxp1)
	return (u+uxp1+uym1+uxp1ym1)/4

def vInterp(v):
	'''Average nearest 4 v points to get middle u point, making use of periodic BCs. See MAC grid diagram in class docs.'''
	vxm1, vxp1 = xStack(v)
	vym1, vyp1 = yStack(v)
	vxm1ym1, vxm1yp1 = yStack(vxm1)
	return (v+vxm1+vyp1+vxm1yp1)/4

# The following are center differences in 2h (e.g. u -> m1, v -> m2)
def xDerivWide(u,h):
	'''First derivative in x using centered differences with periodic boundary conditions in two dimensions. u is a two dimensional array with the ij-th value cooresponding to the (xi,yj) point on the grid (cell edges). h is the grid spacing.'''
	um1, up1 = xStack(u)
	return (up1-um1)/(2*h)	

def yDerivWide(u,h):
	'''First derivative in y using centered differences with periodic boundary conditions in two dimensions. u is a two dimensional array with the ij-th value cooresponding to the (xi,yj) point on the grid (cell edges). h is the grid spacing.'''
	um1, up1 = yStack(u)
	return (up1-um1)/(2*h)	

def Laplacian5pt(u,h):
	'''5 point Laplacian with periodic boundary conditions in two dimensions. u is a two dimensional array with the ij-th value corresponding to the (xi,yj) point on the grid (either cell edges or centers). h is the grid spacing.'''
	uxm1, uxp1 = xStack(u)
	uym1, uyp1 = yStack(u)
	return (uxm1+uxp1+uym1+uyp1-4*u)/h^2

# The following are half-wide center differences (e.g. edges to centers or vice versa).	
def xDerivSlimPlus(u,h):
	'''Returning first derivatives at half points.'''
	um1, up1 = xStack(u)
	return (up1-u)/h	

def xDerivSlimMinus(u,h):
	'''Returning first derivatives at half points.'''
	um1, up1 = xStack(u)
	return (u-um1)/h	

def yDerivSlimPlus(u,h):
	'''Returning first derivatives at half points.'''
	um1, up1 = yStack(u)
	return (up1-u)/h	

def yDerivSlimMinus(u,h):
	'''Returning first derivatives at half points.'''
	um1, up1 = yStack(u)
	return (u-um1)/h	

def divEdge2Cen(u,v,h):
	'''Divergence at cell centers using centered differences with periodic boundary conditions in two dimensions. u,v are two dimensional arrays at the cell edges. h is the grid spacing.'''
	return xDerivSlimPlus(u,h) + yDerivSlimPlus(v,h)

