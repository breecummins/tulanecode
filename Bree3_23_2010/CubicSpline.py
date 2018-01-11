from Numeric import *
from numpy.oldnumeric.linear_algebra import *
from KSPSolvers import *
###############################################################################
##  Python code to fit a given set of data points with a cubic spline. 
##                                                   John Chrispell (11/5/2008)
##                                                           Modified:   Today
###############################################################################
class Spline_base:
    """ 
    This is the base class for spline routines.  
    
    members:
       cubic
    """
    def __init__(self):
        pass

class cubic_spline(Spline_base):
    """
    Returns the coefficients for a cubic spline of a given set of data points.  
    """	
    def __init__(self,xvec,yvec):
	Spline_base.__init__(self)
	self.xvec = xvec
	self.yvec = yvec
        
    def periodic(self):
        # Initilizeation 
        n_pts = len(self.xvec)
        n     = 4*n_pts
        A     = zeros(n*n, Float); A.shape=(n,n)
        x     = ones(n, Float);    x.shape=(n,1)
        b     = zeros(n, Float);   b.shape=(n,1)
        #      S[i](x) = a[i] + b[i](x) + c[i](x)**2 + d[i](x)**3
        # Loop over the given vector of points 
        for i in range(n_pts): 
            if(self.xvec[(i) % n_pts] > self.xvec[(i-1) % n_pts]): 
                # Match Left end point
                A[4*i  ,(4*i  )% n] = 1.0
                A[4*i  ,(4*i+1)% n] =  self.xvec[(i) % n_pts]
                A[4*i  ,(4*i+2)% n] = (self.xvec[(i) % n_pts])**2.0
                A[4*i  ,(4*i+3)% n] = (self.xvec[(i) % n_pts])**3.0
                # Match Right end point
                A[4*i+1,(4*(i-1)  )% n] = 1.0
                A[4*i+1,(4*(i-1)+1)% n] =  self.xvec[(i) % n_pts]
                A[4*i+1,(4*(i-1)+2)% n] = (self.xvec[(i) % n_pts])**2.0
                A[4*i+1,(4*(i-1)+3)% n] = (self.xvec[(i) % n_pts])**3.0
                # Match First derivative
                A[4*i+2,(4*i+1)% n] =  1.0
                A[4*i+2,(4*i+2)% n] =  2.0*( self.xvec[((i) % n_pts)])
                A[4*i+2,(4*i+3)% n] =  3.0*((self.xvec[((i) % n_pts)])**2.0)
                A[4*i+2,(4*(i-1)+1)% n] = -1.0
                A[4*i+2,(4*(i-1)+2)% n] = -2.0*( self.xvec[((i) % n_pts)])
                A[4*i+2,(4*(i-1)+3)% n] = -3.0*((self.xvec[((i) % n_pts)])**2.0)
                # Match Second derivative
                A[4*i+3,(4*i+2)% n] =  2.0
                A[4*i+3,(4*i+3)% n] =  6.0*(self.xvec[(i) % n_pts])
                A[4*i+3,(4*(i-1)+2)% n] = -2.0
                A[4*i+3,(4*(i-1)+3)% n] = -6.0*(self.xvec[(i) % n_pts])
                # Set up Right hand side
                b[ 4*i   % n] = self.yvec[(i) % n_pts]
                b[(4*i+1)% n] = self.yvec[(i) % n_pts]
                b[(4*i+2)% n] = 0.0
                b[(4*i+3)% n] = 0.0

            else: 
                # Match Left end point
                A[4*i  ,(4*i  )% n] = 1.0
                A[4*i  ,(4*i+1)% n] = (self.xvec[(i) % n_pts]) 
                A[4*i  ,(4*i+2)% n] = (self.xvec[(i) % n_pts])**2.0
                A[4*i  ,(4*i+3)% n] = (self.xvec[(i) % n_pts])**3.0
                # Match Right end point
                A[4*i+1,(4*(i-1)  )% n] = 1.0
                A[4*i+1,(4*(i-1)+1)% n] =  self.xvec[(i) % n_pts]+2.0*pi
                A[4*i+1,(4*(i-1)+2)% n] = (self.xvec[(i) % n_pts]+2.0*pi)**2.0
                A[4*i+1,(4*(i-1)+3)% n] = (self.xvec[(i) % n_pts]+2.0*pi)**3.0
                # Match First derivative
                A[4*i+2,(4*i+1)% n] =  1.0
                A[4*i+2,(4*i+2)% n] =  2.0*( self.xvec[((i) % n_pts)])
                A[4*i+2,(4*i+3)% n] =  3.0*((self.xvec[((i) % n_pts)])**2.0)
                A[4*i+2,(4*(i-1)+1)% n] = -1.0
                A[4*i+2,(4*(i-1)+2)% n] = -2.0*( self.xvec[((i) % n_pts)]+2.0*pi)
                A[4*i+2,(4*(i-1)+3)% n] = -3.0*((self.xvec[((i) % n_pts)]+2.0*pi)**2.0)
                # Match Second derivative
                A[4*i+3,(4*i+2)% n] =  2.0
                A[4*i+3,(4*i+3)% n] =  6.0*(self.xvec[(i) % n_pts])
                A[4*i+3,(4*(i-1)+2)% n] = -2.0
                A[4*i+3,(4*(i-1)+3)% n] = -6.0*(self.xvec[(i) % n_pts]+2.0*pi)
                # Set up Right hand side
                b[ 4*i   % n] = self.yvec[(i) % n_pts]
                b[(4*i+1)% n] = self.yvec[(i) % n_pts]
                b[(4*i+2)% n] = 0.0
                b[(4*i+3)% n] = 0.0


        #============================================================
        print "Solving System for function coefficients with gmres   " 
        print "======================================================"
        System = gmres(Matrix=A,RHS=b,x=x,Tol=1e-14,maxIts=len(b))
        coefficents, error, totalIters = System.solve()

        # Rearrange the spline coefficents effectivly 
        NiceCoefficents = zeros(4*n_pts, Float); NiceCoefficents.shape=(n_pts,4)
        for i in range(n_pts):
            NiceCoefficents[i,0] = coefficents[4*i  ,0]
            NiceCoefficents[i,1] = coefficents[4*i+1,0]
            NiceCoefficents[i,2] = coefficents[4*i+2,0]
            NiceCoefficents[i,3] = coefficents[4*i+3,0]

	return NiceCoefficents

class fast_cubic_spline(Spline_base):
    """
    Returns the coefficients for a cubic spline of a given set of data points.  
    """	
    def __init__(self,xvec,yvec):
	Spline_base.__init__(self)
	self.xvec = xvec
	self.yvec = yvec
        
    def periodic(self):
        # Initilizeation 
        n_pts = len(self.xvec)
        n     = n_pts
        A     = zeros(n*n, Float); A.shape=(n,n)
        x     = ones(n, Float);    x.shape=(n,1)
        b     = zeros(n, Float);   b.shape=(n,1)
        # This is a much faster implementation using a periodic tridagonal.
        # Spline is in the form: 
        #      S[i](x) = a[i] + b[i](x-x[i]) + c[i](x-x[i])**2 + d[i](x-x[i])**3
        # Loop over the given vector of points 
        for i in range(n_pts): 
            if(self.xvec[(i+1) % n_pts] > self.xvec[(i) % n_pts] > self.xvec[(i-1) % n_pts]): 
                hi   = (self.xvec[(i+1)% n_pts] - self.xvec[(i)  % n_pts])
                himo = (self.xvec[(i)  % n_pts] - self.xvec[(i-1)% n_pts])     
                # Set up a matrix for the c[i]'s 
                A[i  ,(i  )% n] = 2.0*(himo + hi)
                A[i  ,(i+1)% n] = hi
                A[i  ,(i-1)% n] = himo
                # Set up Right hand side
                b[i] = ((3.0/hi)*(self.yvec[(i+1) % n_pts] - self.yvec[(i) % n_pts]) 
                      - (3.0/himo)*(self.yvec[(i) % n_pts] - self.yvec[(i-1) % n_pts]))
            
            if(self.xvec[(i-1) % n_pts] > self.xvec[(i+1) % n_pts] > self.xvec[(i) % n_pts]):
                hi   = (self.xvec[(i+1)% n_pts]-self.xvec[(i)% n_pts])
                himo = ((self.xvec[(i)% n_pts]+2.0*pi)-self.xvec[(i-1)% n_pts])
                # Set up a matrix for the c[i]'s 
                A[i  ,(i  )% n] = 2.0*(himo + hi)
                A[i  ,(i+1)% n] = hi
                A[i  ,(i-1)% n] = himo
                # Set up Right hand side
                b[i] = ((3.0/hi)*(self.yvec[(i+1) % n_pts] - self.yvec[(i) % n_pts]) 
                       - (3.0/himo)*(self.yvec[(i) % n_pts] - self.yvec[(i-1) % n_pts]))

            if(self.xvec[(i) % n_pts] > self.xvec[(i-1) % n_pts] > self.xvec[(i+1) % n_pts]):
                hi   = ((self.xvec[(i+1)% n_pts]+2.0*pi)-self.xvec[(i)% n_pts])
                himo = (self.xvec[(i)% n_pts]-self.xvec[(i-1)% n_pts])                
                # Set up a matrix for the c[i]'s 
                A[i  ,(i  )% n] = 2.0*(himo + hi)
                A[i  ,(i+1)% n] = hi
                A[i  ,(i-1)% n] = himo
                # Set up Right hand side
                b[i] = ((3.0/hi)*(self.yvec[(i+1) % n_pts] - self.yvec[(i) % n_pts]) 
                      - (3.0/himo)*(self.yvec[(i) % n_pts] - self.yvec[(i-1) % n_pts]))
        
        #============================================================
        print "Solving System for function coefficients with gmres   " 
        print "======================================================"

        System = gmres(Matrix=A,RHS=b,x=x,Tol=1e-14,maxIts=len(b))
        coefficents, error, totalIters = System.solve()

        # Rearrange the spline coefficents effectivly 
        NiceCoefficents = zeros(4*n_pts, Float); NiceCoefficents.shape=(n_pts,4)
        for i in range(n_pts):            
            if(self.xvec[(i+1) % n_pts] > self.xvec[(i) % n_pts]): 
                hi   = (self.xvec[(i+1)% n_pts]-self.xvec[(i)% n_pts])
            else:
                hi   = (self.xvec[(i+1)% n_pts]+2*pi -self.xvec[(i)% n_pts])
            
            NiceCoefficents[i,0] = self.yvec[i]
            NiceCoefficents[i,1] = ((1.0/hi)*(self.yvec[(i+1) %n_pts] - self.yvec[i]) 
                                     - (hi/3.0)*(2.0*coefficents[i,0] + coefficents[(i+1) % n_pts,0]))
            NiceCoefficents[i,2] = coefficents[i,0]
            NiceCoefficents[i,3] = (1.0/(3.0*hi))*(coefficents[(i+1) % n_pts,0] - coefficents[i,0])

	return NiceCoefficents
