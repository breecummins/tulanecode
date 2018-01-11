import numpy as np
from scipy.sparse.linalg import gmres

#FIXME I blow up in 2D

def takeStepNeumann1(n,un,C,g):
    '''
    Set ghost point equal to exact solution.
    n is the time step number.
    un is an NxN matrix with un_ij = un(x_i,y_j).
    g is a handle to a formula for the exact values of the function at the ghost points.
    Returns un+1 -- the next step in the time solver.
    '''
    N = un.shape[0]
    M = np.zeros((N**2,N**2))
    RHS = np.zeros((N**2))
    row = 0
    for i in range(N):
        for j in range(N):
            M[row,row] = 1-4*C
            if j < N-1 and j >0 and i < N-1 and i >0:
                M[row,row-1] = C 
                M[row,row+1] = C 
                M[row,(i-1)*j] = C 
                M[row,(i+1)*j] = C 
                RHS[row] = (1+4*C)*un[i,j] - C*(un[i-1,j]+un[i+1,j]+un[i,j-1]+un[i,j+1])
            elif i == 0:
                M[row,row+1] = C 
                M[row,(i+1)*j] = C 
                if j < N-1 and j >0:
                    M[row,row-1] = C 
                    RHS[row] = (1+4*C)*un[i,j] - C*(un[i+1,j]+un[i,j-1]+un[i,j+1]+g(N,i-1,j,n)+g(N,i-1,j,n+1))
                elif j == 0:
                    RHS[row] = (1+4*C)*un[i,j] - C*(un[i+1,j]+un[i,j+1]+g(N,i-1,j,n)+g(N,i-1,j,n+1)+g(N,i,j-1,n)+g(N,i,j-1,n+1))                    
                elif j == N-1:
                    RHS[row] = (1+4*C)*un[i,j] - C*(un[i+1,j]+un[i,j-1]+g(N,i-1,j,n)+g(N,i-1,j,n+1)+g(N,i,j+1,n)+g(N,i,j+1,n+1))                    
            elif i == N-1:
                M[row,row-1] = C 
                M[row,(i-1)*j] = C 
                if j < N-1 and j >0:
                    M[row,row+1] = C 
                    RHS[row] = (1+4*C)*un[i,j] - C*(un[i-1,j]+un[i,j-1]+un[i,j+1]+g(N,i+1,j,n)+g(N,i+1,j,n+1))
                elif j == 0:
                    RHS[row] = (1+4*C)*un[i,j] - C*(un[i-1,j]+un[i,j+1]+g(N,i,j-1,n)+g(N,i,j-1,n+1)+g(N,i+1,j,n)+g(N,i+1,j,n+1))                    
                elif j == N-1:
                    RHS[row] = (1+4*C)*un[i,j] - C*(un[i-1,j]+un[i,j-1]+g(N,i,j+1,n)+g(N,i,j+1,n+1)+g(N,i+1,j,n)+g(N,i+1,j,n+1))                    
            elif j == 0:
                M[row,row+1] = C 
                M[row,(i+1)*j] = C 
                M[row,(i-1)*j] = C 
                RHS[row] = (1+4*C)*un[i,j] -  C*(un[i-1,j]+un[i+1,j]+un[i,j+1]+g(N,i,j-1,n)+g(N,i,j-1,n+1))
            elif j == N-1:
                M[row,row-1] = C 
                M[row,(i+1)*j] = C 
                M[row,(i-1)*j] = C 
                RHS[row] = (1+4*C)*un[i,j] -  C*(un[i-1,j]+un[i+1,j]+un[i,j-1]+g(N,i,j+1,n)+g(N,i,j+1,n+1))
            row += 1
    unp1, info = gmres(M,RHS)
    return unp1.reshape(N,N)
            
        

def takeStepNeumann2(n,un,C,G,h):
    '''
    Set ghost point equal to extrapolation from exact gradient at boundary.
    un is an NxN matrix with un_ij = un(x_i,y_j).
    C is a numerical coefficient.
    G is a handle to a formula for the exact values of the gradient on the boundary.
    h is the spatial spacing.
    '''
    pass


def exactSoln2HeatEqn(N,xind=None,yind=None,tind=None):
    h=0.5/N
    mu = 1.0
    b = 1.0
    if xind == None:
        #return initial condition with b = mu = 1
        X, Y = np.meshgrid(np.arange(h/2,0.5,h), np.arange(h/2,0.5,h))
        return np.exp( -( (X-0.25)**2 + (Y-0.25)**2 )/(4*b) )/b, h, mu
    else:
        x = h/2 + xind*h
        y = h/2 + yind*h
        t = tind*h
        return (1./(b+mu*t))*np.exp( -( (x-0.25)**2 + (y-0.25)**2 )/(4*(b+mu*t)) )


def CrankNicholson(N):
    '''
    Iterate the differential equation.
    '''
    g0, h, mu = exactSoln2HeatEqn(N)
    Tf = 10
    dt = h
    C = - (mu*dt)/(2*h**2)
    un = g0
    savevals = np.zeros( (np.floor(Tf/dt),un.shape[0],un.shape[1]) )
    for n in range(np.floor(Tf/dt).astype('int32')):
        savevals[n,:,:] = un
        un = takeStepNeumann1(n,un,C,exactSoln2HeatEqn)
        print(n)
        print(np.max(np.abs(un)))
    return savevals



































