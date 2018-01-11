#!/usr/bin/env python

import numpy as nm
import sys

def solver(func,IC,Tol=1.e-9, NumIts=100000):
	dif = 10*Tol
	its = 0
	x = IC
	while dif > Tol:
		xnew = func(x)
		dif = abs(xnew - x)
		x = xnew
		its+=1
		if its > NumIts:
			print('Maximum number of iterations reached. Last difference was ' + str(dif) + '.')
			sys.exit()
	print(x)
	return x
		
def testcase1(x):
	'''Exception: Can only get zero fixed point; function too steep for 1 or -1.'''
	return x**3

def testcase2(x):
	'''Exception: Cannot get zero fixed point. Approximates x near the fixed point; always have a nonzero value that approximates sin(x) arbitrarily well.'''
	return nm.sin(x)
	
def testcase3(x):
	'''Can only get lower fixed point.'''
	return nm.exp(x) - 2
	
def testcase4(x):
	'''Rearrange test case 3 -- only lower fixed point.'''
	return nm.log(x+2)
	
def testcase5(x):
	return (x**2-1)/3.
	
if __name__ == '__main__':
	solver(testcase5,3.4)

	

