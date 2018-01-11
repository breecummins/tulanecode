#!/usr/bin/env python

import numpy as nm
import sys

def solver(F,a,b,Tol=1.e-9, NumIts=100000):
	lval = F(a)
	rval = F(b)
	if (lval >0 and rval > 0 ) or (lval <0 and rval <0):
		print('Choose a different interval; F(a) and F(b) do not fall on opposite sides of zero.')
		sys.exit()
	h = 0
	its = 0
	while its < NumIts:
		h = (b-a)/2.
		p = a + h
		newval = F(p)
		if newval == 0 or h < Tol:
			print(p)
			return p
		if lval*newval > 0:
			a = p
			lval = newval
		else:
			b = p		
		its+=1
	print('Failed to find root after ',str(NumIts),' iterations.')
		
def testcase1(x):
	return x**3 - x

	
if __name__ == '__main__':
	solver(testcase1,-1.5,-0.5)

	

