from Numeric import *
from numpy.oldnumeric.linear_algebra import *
###############################################################################
##  Python implementations of the GMRES and bicgstab algorithms. 
##                                                   John Chrispell (11/4/2008)
##                                                           Modified:   Today
###############################################################################
class Solver_base:
    """ 
    This is the base class for some KSP solvers.  
    
    members:
       bcgstab
    """
    def __init__(self):
        pass

class gmres(Solver_base):
    """
    This is a function that solves Ax = b using the gmres routine. 
    """	
    def __init__(self,Matrix,RHS, x, Tol=1e-11, maxIts=50):
		Solver_base.__init__(self)
		self.Matrix = Matrix
		self.RHS = RHS
		self.x = x
	    self.Tol = Tol
		self.maxIts = maxIts

    def backsub(self,U,b):
        n = len(b)
        x = b
        
        j = n-1
        while(j >= 0):
            if(U[j,j] == 0):
                print "Singular Matrix Error!"
            x[j] = b[j]/U[j,j]
            for i in range(j):
                b[i] = b[i] - U[i,j]*x[j]
            j = j-1
        return x
            

    def givapp(self,c,s,vin,k):
        # Apply a sequence of givens rotations. 

        self.k = k 
        if(self.k <= 1):
            self.c = []
            self.c.append(c)
            self.s = [] 
            self.s.append(s)
        else:
            self.c = c
            self.s = s

        self.vrot = vin

        for i in range(self.k):
            w1 = self.c[i]*self.vrot[i] - self.s[i]*self.vrot[i+1]
            w2 = self.s[i]*self.vrot[i] + self.c[i]*self.vrot[i+1]

            self.vrot[i] = w1 
            self.vrot[i+1] = w2

        return self.vrot
        
    def solve(self):
        # Initilizeation 
        n = len(self.RHS)
        h = zeros((self.maxIts+1)*(self.maxIts+1), Float); h.shape=(self.maxIts+1,self.maxIts+1)
        v = zeros((self.maxIts+1)*(n), Float);             v.shape=(n,self.maxIts+1)
        c = zeros((self.maxIts + 1), Float);               c.shape=((self.maxIts+1),1)
        s = zeros((self.maxIts + 1), Float);               s.shape=((self.maxIts+1),1)
        
        if( sqrt(float(matrixmultiply(transpose(self.x),self.x))) != 0.0): 
            r = self.RHS - matrixmultiply(self.Matrix, self.x) 
        else:
            r = self.RHS 
        
        rho = sqrt(float(matrixmultiply(transpose(r),r)))

        g = zeros((self.maxIts + 1),Float); g.shape=(self.maxIts+1,1)
        g[0,0] = rho
        errtol = self.Tol*(sqrt(float(matrixmultiply(transpose(self.RHS),self.RHS))))
        error = []
        error.append(rho)
        totalIters = 0

        if(rho < errtol):
            return (self.x)

        v[:,0] = (1.0/rho)*r[:,0] 
        beta = rho
        k = 0

        # The GMRES Iteration 
        while((rho > errtol) and (k < self.maxIts)):
            k = k + 1
            v[:,k] = matrixmultiply(self.Matrix,v[:,(k-1)])
            normav = sqrt(matrixmultiply(transpose(v[:,k]),v[:,k]))

            # Modified Gram-Schmidt
            for j in range(k):
                h[j,k-1] = matrixmultiply(transpose(v[:,j]),v[:,k])
                v[:,k] = v[:,k] - h[j,k-1]*(v[:,j])

            h[k,k-1] = sqrt(matrixmultiply(transpose(v[:,k]),v[:,k]))
            normav2 = h[k,k-1]

            # Reorthoganlaize (Brown-Hindmarsh condition)
            if( (normav + 0.001*normav2) == normav ):
                for j in range(k):
                    hr = matrixmultiply(transpose(v[:,j]),v[:,k])
                    h[j,k-1] = h[j,k-1] + hr
                    v[:,k] = v[:,k] - hr*v[:,j]
                h[k,k-1] = sqrt(matrixmultiply(transpose(v[:,k]),v[:,k]))
                
            if(h[k,k-1] != 0):
                v[:,k] = v[:,k]/h[k,k-1]

            # Form the new Givens Rotation 
            if(k > 1): 
                h[0:k,k-1] = self.givapp(c=c[0:k-1,0],s=s[0:k-1,0],vin = h[0:k,k-1],k=k-1)

            nu = sqrt( matrixmultiply(transpose(h[k-1:k+1,k-1]),h[k-1:k+1,k-1] ))

            if(nu != 0):
                c[k-1,0] = float(h[k-1,k-1]/nu)
                s[k-1,0] = float(-h[k,k-1]/nu)
                h[k-1,k-1] = float(c[k-1]*h[k-1,k-1] - s[k-1]*h[k,k-1])
                if(k < self.maxIts):
                    h[k,k-1] = 0
                g[k-1:k+1,0] = self.givapp(c=c[k-1,0],s=s[k-1,0],vin=g[k-1:k+1,0],k=1)

            # Update the residual 
            rho = abs(g[k,0])
            error.append(float(rho))

        # Compute x 
        y = self.backsub(U=h[0:k,0:k],b=g[0:k,0])
        totalIters = k
        self.x[:,0] = self.x[:,0] + matrixmultiply(v[0:n,0:k],y)

	return self.x, error, totalIters


class bicgstab(Solver_base):
    """
    This is a function that solves Ax = b using the bicgstab routine. 
    """	
    def __init__(self,Matrix,RHS, x, Tol=1e-11, maxIts=50):
	Solver_base.__init__(self)
	self.Matrix = Matrix
	self.RHS = RHS
	self.x = x
        self.Tol = Tol
	self.maxIts = maxIts

    def norm(self,vec):
        norm = float(sqrt(matrixmultiply(transpose(vec),vec)))
        return norm
    
    def solve(self):

        # Initialize the Algorithm 
        n = len(self.RHS)
        errortol = self.Tol*self.norm(self.RHS)
        rho = zeros((self.maxIts+3), Float);   rho.shape=(self.maxIts+3,1)
        hatr0 = zeros(n, Float); hatr0.shape=(n,1)
        r     = zeros(n, Float); r.shape=(n,1)
        v     = zeros(n, Float); v.shape=(n,1)
        p     = zeros(n, Float); p.shape=(n,1)
        s     = zeros(n, Float); s.shape=(n,1)
        error = [] 
        

        if(self.norm(self.x) != 0.0): 
            r = self.RHS - matrixmultiply(self.Matrix, self.x) 
        else:
            r = self.RHS

        k        = 0
        hatr0[:,0] = r[:,0]
        alpha    = 1.0
        omega    = 1.0 
        rho[0,0] = 1.0
        rho[1,0] = float(matrixmultiply(transpose(hatr0[:,0]),r[:,0]))
        zeta = self.norm(r)
        error.append(zeta)

        # Bi-CGSTAB iteration 
        while((k <= self.maxIts) and (zeta > errortol)): 
            k=k+1
            if(omega == 0.0):
                print "Division by zero in bicgstab Algorithm flag 0"
                return
            beta = float((float(rho[k,0])/float(rho[k-1,0]))*(alpha/omega))
            p[:,0] =  r[:,0] + beta*(p[:,0]-omega*v[:,0])
            v[:,0] =  matrixmultiply(self.Matrix,p[:,0])
            tau = float(matrixmultiply(transpose(hatr0[:,0]),v[:,0]))
            if(tau == 0.0):
                print "Division by zero in bicgstab Algorithm flag 1"
                return
            alpha = float(float(rho[k,0])/tau)
            s[:,0] = r[:,0] - alpha*v[:,0]
            t = matrixmultiply(self.Matrix,s)
            tau = float(matrixmultiply(transpose(t[:,0]),t[:,0]))
            if(tau == 0.0):
                print "Division by zero in bicgstab Algorithm flag 2"
                return
            omega = float((matrixmultiply(transpose(t[:,0]),s[:,0]))/tau) 
            rho[k+1,0] = (-1.0)*omega*(matrixmultiply(transpose(hatr0[:,0]),t[:,0]))
            self.x[:,0] = self.x[:,0] + alpha*p[:,0] + omega*s[:,0]
            r[:,0] = s[:,0] - omega*t[:,0]
            zeta = self.norm(r[:,0])
            error.append(zeta)

        totalIters = k
	return self.x, error, totalIters
