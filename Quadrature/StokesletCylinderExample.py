#!/usr/bin/env python

import numpy as nm
import os
import sys

from scipy.sparse.linalg import gmres
import mat2py as mf

import BoundaryElements as be
import BasisFunctions as bf
import QuadratureMethods as qm
import IntegrationKernels as ik


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
	
	#BEM method, both blobs, varying N
	mu = 1.0
	rad = 0.25
	obspts = constructpatches()
	exactu, exactv = exactSoln(obspts,rad,mu)
	M = 50
	hat = bf.HatBasisFunc(M,'t')
	numchords = [2**k for k in range(3,9)]
	for N in numchords:
		numbps = 100
		chordlen = rad*2*nm.sin(nm.pi/N)
		blobparams = nm.linspace((chordlen/M)/100,8*chordlen/M,numbps)
		BCs = nm.asarray(N*[1,0])
		#BEM
		circBEM = be.Circle(N,M,rad)
		circBEM.LinearElements2D()
		trap = qm.CompTrapRule(circBEM.nodespacing)
		#begin iterating over blob params
		u_blob1 = nm.empty((0,))
		v_blob1 = nm.empty((0,))
		u_blob2 = nm.empty((0,))
		v_blob2 = nm.empty((0,))
		for k in range(numbps):
			bp = blobparams[k]
			print('Hat/trap BEM, N = ' + str(N)+ ', Blob = '+str(bp))
			blob1 = ik.Stokeslet2DReg(bp,mu)
			blob2 = ik.Stokeslet2DReg2(bp,mu)
			fmat1 = blob1.makeMatrix2DBEM(circBEM.nodes,circBEM,hat,trap)
			fmat2 = blob2.makeMatrix2DBEM(circBEM.nodes,circBEM,hat,trap)
			forces1 = useGMRES(fmat1,BCs)
			forces2 = useGMRES(fmat2,BCs)
			vmat1 = blob1.makeMatrix2DBEM(obspts,circBEM,hat,trap)
			vmat2 = blob2.makeMatrix2DBEM(obspts,circBEM,hat,trap)
			vel1 = nm.dot(vmat1,forces1)
			vel2 = nm.dot(vmat2,forces2)
			if len(u_blob1) == 0:
				u_blob1 = vel1[0::2]                   	
				v_blob1 = vel1[1::2]                    	
				u_blob2 = vel2[0::2]                   	
				v_blob2 = vel2[1::2]                    	
			else:                                                     	
				u_blob1 = nm.append(u_blob1,vel1[0::2])	
				v_blob1 = nm.append(v_blob1,vel1[1::2])	
				u_blob2 = nm.append(u_blob2,vel2[0::2])	
				v_blob2 = nm.append(v_blob2,vel2[1::2])	
		fname = 'StokesCylTest_hattrapBEM_bothblobs'+ (4-len(str(N)))*'0' + str(N) + '.mat'
		mf.write(fname,{'u_blob1':u_blob1,'v_blob1':v_blob1,'u_blob2':u_blob2,'v_blob2':v_blob2,'exactu':exactu,'exactv':exactv,'obspts':obspts,'blobparams':blobparams,'numchords':numchords,'M':M})
		print('Results in ' + fname )
		
	
		#original method, both blobs, varying N
		mu = 1.0
		rad = 0.25
		obspts = constructpatches()
		exactu, exactv = exactSoln(obspts,rad,mu)
		numchords = [2**k for k in range(3,9)]
		for N in numchords:
			numbps = 50
			chordlen = rad*2*nm.sin(nm.pi/N)
			blobparams = nm.linspace(chordlen/15,chordlen,numbps)
			BCs = nm.asarray(N*[1,0])
			circ = be.Circle(N,1,rad)
			circ.LinearElements2D()
			#begin iterating over blob params
			u_blob1 = nm.empty((0,))
			v_blob1 = nm.empty((0,))
			u_blob2 = nm.empty((0,))
			v_blob2 = nm.empty((0,))
			for k in range(numbps):
				bp = blobparams[k]
				print('Original method, N = ' + str(N)+ ', Blob = '+str(bp))
				blob1 = ik.Stokeslet2DReg(bp,mu)
				blob2 = ik.Stokeslet2DReg2(bp,mu)
				fmat1 = blob1.makeMatrix2D(circ.nodes,circ.nodes)
				fmat2 = blob2.makeMatrix2D(circ.nodes,circ.nodes)
				#fmat = origmeth.makeMatrix2DBEM(circ.nodes,circ,cst,riem)
				forces1 = useGMRES(fmat1,BCs)
				forces2 = useGMRES(fmat2,BCs)
				vmat1 = blob1.makeMatrix2D(obspts,circ.nodes)
				vmat2 = blob2.makeMatrix2D(obspts,circ.nodes)
				# vmat = origmeth.makeMatrix2DBEM(obspts,circ,cst,riem)
				vel1 = nm.dot(vmat1,forces1)
				vel2 = nm.dot(vmat2,forces2)
				if len(u_blob1) == 0:
					u_blob1 = vel1[0::2]                   	
					v_blob1 = vel1[1::2]                    	
					u_blob2 = vel2[0::2]                   	
					v_blob2 = vel2[1::2]                    	
				else:                                                     	
					u_blob1 = nm.append(u_blob1,vel1[0::2])	
					v_blob1 = nm.append(v_blob1,vel1[1::2])	
					u_blob2 = nm.append(u_blob2,vel2[0::2])	
					v_blob2 = nm.append(v_blob2,vel2[1::2])	
			fname = 'StokesCylTest_origmethod_bothblobs'+ (4-len(str(N)))*'0' + str(N) + '.mat'
			mf.write(fname,{'u_blob1':u_blob1,'v_blob1':v_blob1,'u_blob2':u_blob2,'v_blob2':v_blob2,'exactu':exactu,'exactv':exactv,'obspts':obspts,'blobparams':blobparams})
			print('Results in ' + fname )


	# #try to match Ricardo's pictures -- original method and with hat/trap
	# N = 128
	# mu = 1.0
	# rad = 0.25
	# obspts = constructpatches()
	# exactu, exactv = exactSoln(obspts,rad,mu)
	# # circ = be.Circle(N,1,rad)
	# # circ.LinearElements2D()
	# # cst = bf.CstBasisFunc(1,'r')
	# # riem = qm.RiemannSum(circ.nodespacing)
	# # blobparams = nm.arange(1.e-3,0.013,5.e-4)
	# M = 40
	# circBEM = be.Circle(N,M,rad)
	# circBEM.LinearElements2D()
	# hat = bf.HatBasisFunc(M,'t')
	# trap = qm.CompTrapRule(circBEM.nodespacing)
	# blobparams = nm.arange(1.e-5,5.e-4,2.e-6)
	# BCs = nm.asarray(N*[1,0])
	# u = nm.empty((0,))
	# v = nm.empty((0,))
	# for bp in blobparams:
	# 	print('Blob = '+str(bp))
	# 	# origmeth = ik.Stokeslet2DReg2(bp,mu)
	# 	# # fmat = origmeth.makeMatrix2D(circ.nodes,circ.nodes)
	# 	# fmat = origmeth.makeMatrix2DBEM(circ.nodes,circ,cst,riem)
	# 	# forces = useGMRES(fmat,BCs)
	# 	# # vmat = origmeth.makeMatrix2D(obspts,circ.nodes)
	# 	# vmat = origmeth.makeMatrix2DBEM(obspts,circ,cst,riem)
	# 	# vel = nm.dot(vmat,forces)
	# 	# if len(u) == 0:
	# 	# 	u = vel[0::2]                   	
	# 	# 	v = vel[1::2]                    	
	# 	# else:                                                     	
	# 	# 	u = nm.append(u,vel[0::2])	
	# 	# 	v = nm.append(v,vel[1::2])	
	# 	hattrapBEM = ik.Stokeslet2DReg2(bp,mu)
	# 	fmat = hattrapBEM.makeMatrix2DBEM(circBEM.nodes,circBEM,hat,trap)
	# 	forces = useGMRES(fmat,BCs)
	# 	vmat = hattrapBEM.makeMatrix2DBEM(obspts,circBEM,hat,trap)
	# 	vel = nm.dot(vmat,forces)
	# 	if len(u) == 0:
	# 		u = vel[0::2]                   	
	# 		v = vel[1::2]                    	
	# 	else:                                                     	
	# 		u = nm.append(u,vel[0::2])	
	# 		v = nm.append(v,vel[1::2])
	# # fname = 'StokesCylTest_matchRicardopatches_origmethodnewblob'+ (4-len(str(N)))*'0' + str(N) + '.mat'
	# fname = 'StokesCylTest_matchRicardopatches_hattrapnewblob'+ (4-len(str(N)))*'0' + str(N) + '.mat'
	# mf.write(fname,{'u':u,'v':v,'exactu':exactu,'exactv':exactv,'obspts':obspts,'blobparams':blobparams})
	# print('Results in ' + fname )
	# 	
	# 
	# 
		