#!/usr/bin/env python

import numpy as nm

class QuadMethod:
	'''This is the base class for quadrature methods. Composite trapezoid, Simpsons, and Riemann rules were unit tested on 05/06/2010 using: 1) f(x) = 2*pi on [0,3]; 2) f(x) = x*sqrt(2) on [0,3]; and 3) f(x) = sqrt(4 - x^2) on [-2,2].'''
	def __init__(self, quadspacing):
		self.h = quadspacing
		
	def performQuad(self,f):
		'''f is a 1D array containing the function values at x_i (equally spaced points in a 1D domain).'''
		P = f.shape[0]
		coeffs = self.makeCoeffs(P)
		quad = nm.sum(coeffs*f)
		return quad

	
class CompTrapRule(QuadMethod):
	'''Composite trapezoidal rule.'''
	def __init__(self,quadspacing):
		QuadMethod.__init__(self,quadspacing)
		
	def makeCoeffs(self,P):
		'''P is the number of equally-spaced points in the quadrature.'''
		coeffs = self.h*nm.asarray([0.5] + (P-2)*[1] + [0.5])
		return coeffs

class CompSimpRule(QuadMethod):
	'''Composite Simpson\'s rule.'''
	def __init__(self,quadspacing):
		QuadMethod.__init__(self,quadspacing)

	def makeCoeffs(self,P):
		'''P is the number of equally-spaced points in the quadrature.'''
		if not nm.mod(P,2):
			print("Warning: P is even. P must be odd for Simpson's rule.")	
			import sys
			sys.exit()
		coeffs = self.h*nm.asarray([1./3] + ((P-3)/2)*[4./3, 2./3] + [4./3, 1./3])
		return coeffs


class RiemannSum(QuadMethod):
	'''Easy-peasy Riemann sum.'''
	def __init__(self,quadspacing):
		QuadMethod.__init__(self,quadspacing)

	def makeCoeffs(self,P):
		'''P is the number of equally-spaced points in the quadrature.'''
		coeffs = self.h*nm.ones(P)
		return coeffs

