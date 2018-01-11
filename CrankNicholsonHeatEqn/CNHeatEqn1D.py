import numpy as np
import scipy.sparse.linalg as spl
from scipy.sparse import spdiags
import matplotlib.pyplot as plt

def takeStepNeumann1(n,un,C,g,p):
    '''
    Set ghost point equal to exact solution.
    n is the time step number.
    un is an N vector with un_i = un(x_i).
    g is a handle to a formula for the exact values of the function at the ghost points.
    p are function parameters
    Returns right hand side of matrix equation
    '''
    #build the right hand side of the matrix equation accounting for boundary conditions
    RHS = (1+2*C)*un
    RHS[1:-1] -= C*( un[:-2] + un[2:])
    #overwrite boundary edges
    RHS[0] =  g(-p['h']/2.,(n+1)*p['dt'],p)
    RHS[-1] = g(p['N']*p['h'] + p['h']/2.,(n+1)*p['dt'],p)
    return RHS
            
        

def takeStepNeumann2(n,un,C,G,p):
    '''
    Set ghost point equal to extrapolation from exact gradient at boundary.
    un is an NxN matrix with un_ij = un(x_i,y_j).
    C is a numerical coefficient.
    G is a handle to a formula for the exact values of the gradient on the boundary.
    p are function parameters
    Returns right hand side of matrix equation
    '''
    #build the right hand side of the matrix equation accounting for boundary conditions
    RHS = (1+2*C)*un
    RHS[1:-1] -= C*( un[:-2] + un[2:])
    #overwrite boundary edges
    RHS[0] = - G(0,(n+1)*p['dt'],p)*p['h'] 
    RHS[-1] = G(p['N']*p['h'],(n+1)*p['dt'],p)*p['h'] 
    return RHS


def exactSoln2HeatEqn(x,t,p):
    denom = np.sqrt(p['b']+2*p['mu']*t)
    exponum = -(x-0.25)**2
    expodenom = 2*(p['b']+2*p['mu']*t)
    return np.exp( exponum/expodenom ) / denom


def exactgradSoln2HeatEqn(x,t,p):
    denom = (p['b']+2*p['mu']*t)**(1.5)
    exponum = -(x-0.25)**2
    expodenom = 2*(p['b']+2*p['mu']*t)
    return -(x-0.25)*np.exp( exponum/expodenom ) / denom

def testFunc(x,t,p):
    return 0

def setParams(N,Tf):
    h=0.5/N
    mu = 0.01
    b = 0.01
    dt = h**2
    return {'N':N, 'h':h, 'mu':mu, 'b':b, 'dt':dt, 'Tf':Tf}

def buildMatrix(C,N,flag):
    # build C-N matrix for different Neumann BC methods
    sidep = C*np.ones(N+2)
    sidem = C*np.ones(N+2)
    if flag == 1:
        sidep[1] = 0
        sidem[-2] = 0
    elif flag == 2:
        sidep[1] = -1
        sidem[-2] = -1
    #construct the C-N tridiagonal sparse matrix
    main = (1-2*C)*np.ones(N+2)
    main[0] = 1
    main[-1] = 1
    diagrows = np.array([sidem,main,sidep])
    diags = [-1,0,1]
    M=spdiags(diagrows, diags, N+2, N+2)  
    return M     


def CrankNicholson(N,Tf,flag):

    '''
    Iterate the differential equation.
    N = number of spatial points.
    Tf = total time to run.
    flag = 1 uses exact assignment method
    flag = 2 uses gradient extrapolation method
    '''

    # initial conditions and parameters
    p = setParams(N,Tf)
    x = np.arange(-p['h']/2, (N+1)*p['h'], p['h']) #grid points + 2 ghost points
    un = exactSoln2HeatEqn(x,0,p)
    C = - (p['mu']*p['dt'])/(2*p['h']**2)
    if flag == 1:
        myUpdater = takeStepNeumann1
        myfunc = exactSoln2HeatEqn
    elif flag == 2:
        myUpdater = takeStepNeumann2
        myfunc = exactgradSoln2HeatEqn 

    # make C-N matrix
    M = buildMatrix(C,N,flag)
    Minv = np.linalg.inv(M.todense())

    # iterate
    numsteps = np.floor(Tf/p['dt']).astype('int32')
    savevals = np.zeros( (numsteps,N+2) )
    for n in range(numsteps):
        savevals[n,:] = un
        RHS = myUpdater(n,un.squeeze(),C,myfunc,p)
        #solve the equation for the next time step
        un = np.dot(Minv,RHS)
        un = np.asarray(un)
#        un = np.linalg.solve(M.todense(),RHS) #fastest for N up to 128
#        un = spl.spsolve(M,RHS)  #complains about matrix format, but seems to be right
#        un, info = spl.gmres(M,RHS) # slow and potentially buggy. Don't use!
    return savevals, p

def doMany(Nlist,Tf,flag):
    plt.close('all')
    results=[]
    params=[]
    import time
    for N in Nlist:
        print(N)
        start = time.clock()
        r, p = CrankNicholson(N,Tf,flag)
        results.append( r )
        params.append( p )
        end = time.clock()
        print('Time elapsed is %f' % (end-start) + ' seconds')
    vizAll(Nlist,results,flag,params)
    return Nlist, results, params

def vizAll(Nlist,results,flag,params,back=0):
    l1err, l2err, linferr, l2relerr = vizError(Nlist,results,flag,params)
    plotExact(params[0],0.2)
    for k in range(len(Nlist)):
        plotNumerical(Nlist[k],results[k],flag,params[k],0.2)    
    if back:
        return l1err, l2err, linferr, l2relerr


def vizResults(N,savevals,p):
    l2err = []
    l1err = []
    linferr = []
    l2relerr = []
    x = np.arange(-p['h']/2, (N+1)*p['h'], p['h'])
    for n in range(savevals.shape[0]):
        exactans = exactSoln2HeatEqn(x,n*p['dt'],p)
        l2err.append( np.sqrt( p['h']*( ( (savevals[n,:] - exactans)**2 ).sum() ) ) )
        l1err.append( p['h']*( ( np.abs(savevals[n,:] - exactans) ).sum() ) ) 
        linferr.append( np.max( np.abs(savevals[n,:] - exactans) ) )
        l2relerr.append( np.sqrt( p['h']*( ( ( (savevals[n,:] - exactans)/exactans )**2 ).sum() ) ) )
    plt.figure(1)
    plt.plot(p['h']*np.arange(savevals.shape[0]),l2err)
    plt.figure(2)
    plt.plot(p['h']*np.arange(savevals.shape[0]),l1err)
    plt.figure(3)
    plt.plot(p['h']*np.arange(savevals.shape[0]),linferr)
    plt.figure(4)
    plt.plot(p['h']*np.arange(savevals.shape[0]),l2relerr)
    return l1err, l2err, linferr, l2relerr

def vizSpatialError(N,savevals,n,p):
    x = np.arange(-p['h']/2, (N+1)*p['h'], p['h'])
    exactans = exactSoln2HeatEqn(x,n*p['dt'],p)
    plt.figure()
    plt.plot(x,np.abs(savevals[n,:] - exactans))
 
        
def vizError(Nlist,results,flag,params):
    plt.close('all')
    l2err = []
    l1err = []
    linferr = []
    l2relerr = []
    for k in range(len(Nlist)):
        l1, l2, li, l2r =  vizResults(Nlist[k],results[k],params[k]) 
        l2err.append(l2)
        l1err.append(l1)
        linferr.append(li)
        l2relerr.append(l2r)
    plt.figure(1)
    plt.title('Neumann method %d' % flag + ', $L_2$ error')
    plt.legend([str(N) for N in Nlist])
    plt.figure(2)
    plt.title('Neumann method %d' % flag + ', $L_1$ error')
    plt.legend([str(N) for N in Nlist])
    plt.figure(3)
    plt.title('Neumann method %d' % flag + ', $L_\infty$ error')
    plt.legend([str(N) for N in Nlist])
    plt.figure(4)
    plt.title('Neumann method %d' % flag + ', $L_2$ relative error')
    plt.legend([str(N) for N in Nlist])
    return l1err, l2err, linferr, l2relerr

    
def plotExact(p,tstep):
#    plt.close('all')
    plt.figure()
    timevec=np.arange(0,min(p['Tf'],10),tstep)
    if p['Tf'] > 10:
        rem = p['Tf'] - 10
        timevec2 = np.linspace(10,p['Tf'],20)
        timevec = np.append(timevec,timevec2)
    x = np.linspace(0,0.5,200)
    for t in timevec:
        exactans = exactSoln2HeatEqn(x,t,p)
        plt.plot(x,exactans)
    if len(timevec) <= 5:
        titlestr= str(timevec)
    else:
        titlestr = ('Exact, Times %f - %f' % (timevec[0],timevec[-1]))
    plt.title(titlestr)
        
def plotNumerical(N,results,flag,p,tstep):
#    plt.close('all')
    plt.figure()
    nt = np.round(tstep/p['dt']).astype('int32')    
    tindvec=range(0,min(np.int_(10/p['dt']),results.shape[0]),max(nt,1))
    if p['Tf'] > 10:
        rem = (p['Tf'] - 10)/p['dt']
        skip = np.round(rem/20)
        tindvec2 = np.arange(10/p['dt'],p['Tf']/p['dt'],skip)
        tindvec = np.append(tindvec,tindvec2)
    x = np.arange(-p['h']/2, (N+1)*p['h'], p['h'])
    for tind in tindvec:
        y = results[tind,:]
        plt.plot(x,y)
    plt.axis([0,p['N']*p['h'],0,1./np.sqrt(p['b'])])
    titlestr = ('Neumann method %d, N = %d, Times %f - %f' % (flag,N,0,tind*p['dt']))
    plt.title(titlestr)


if __name__ == '__main__':
    out1 = doMany([32],50,1)
    vizAll(out1[0],out1[1],1,out1[2])




























