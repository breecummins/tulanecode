      	program SCVT_By_Lloyd
      	include "scvt.m"

C*******************************************************************
C  Read the data sets for initialization                            
C*******************************************************************
      	open(15,file='scvt.in',status='unknown')
        read(15,*) 
      	read(15,*) ntemp
        read(15,*) 
C  Read the maximal number of iterations for Lloyd's method
      	read(15,*) max_iter
        read(15,*)
C  Read the tolerance for Lloyd's method         
	read(15,*) eps
        read(15,*)
C  Should we place one node at the poles ?(0--no -1--south 1--North)
        read(15,*) ip
        close(15)

C  The file "scvt_s.dat" must exist and then read it
       	open(16,file='scvt_s.dat',status='unknown')
C  Read the number of generators
       	read(16,*) n
        if (n.ge.nmax) then
           print *,"The number of generators must be less than ",nmax
           stop
        endif
C  Read the initial coordinates of the generators
       	do node = 1,n
           read(16,*) ntemp,x(node),y(node),z(node)
      	enddo
      	close(16)

        open (35,file='fixedpoints.m',status='unknown')
        read (35,*) ntemp
        close(35)

        eps = eps/(1.0*n)**0.5
        print *,"Number of Generators = ",n
        print *,"Number of Maximum Iterations = ",max_iter
        print *,"Maximum Tolerance = ",eps

C  Fix the south/north pole if necessary
        indx = 0   
        if (ip.ne.0) then
           distmin = 1000000.0
           do node = 1,n
              dist = abs(x(node))+abs(y(node))+abs(z(node)-ip*1.0)           
              if (dist<distmin) then
                 distmin = dist
                 indx = node
              endif
           enddo
           x(indx) = 0.0
           y(indx) = 0.0
           z(indx) = ip*1.0
        endif
      	print *,'Initialization is done!'
        print *,"Fixed node index = ",indx
        print *,"Number of Fixed nodes on interface = ",ntemp

C*******************************************************************
C  Start the Lloyd's algorithm 
C*******************************************************************
        print *,'Start Lloyd iteration ...'
      	ic = 0
      	do iloop = 1,max_iter
           call grid(n,x,y,z,xc,yc,zc,neigh,neisz,vortx,vorsz,ic)
           dmove = 0.0
           do node = 1,n
              ns = neisz(node)
              snw(1:3) = 0.0
              vsnw  = 0.0
              v1(1) = x(node)
              v1(2) = y(node)
              v1(3) = z(node)
              dfw1 = dens_f(v1) 
              do i1 = 1,ns
                 i2 = mod(i1,ns)+1 
                 v2(1) = xc(vortx(node,i1))
                 v2(2) = yc(vortx(node,i1))
                 v2(3) = zc(vortx(node,i1))
                 dfw2 = dens_f(v2)
                 v3(1) = xc(vortx(node,i2))
                 v3(2) = yc(vortx(node,i2))
                 v3(3) = zc(vortx(node,i2))
                 dfw3 = dens_f(v3)
                 T_area = areas(v1,v2,v3)            
                 do j = 1,3 
                    crdwei = v1(j)*dfw1+v2(j)*dfw2+v3(j)*dfw3
                    snw(j) = snw(j)+T_area*crdwei
                 enddo
                 vsnw = vsnw+T_area*(dfw1+dfw2+dfw3)
              enddo  
C  Compute the centroid according to the density function
              snw(1:3) = snw(1:3)/vsnw
              st = sqrt((snw(1))**2+(snw(2))**2+(snw(3))**2)
              snw(1:3) = snw(1:3)/st
              dm = sqrt((x(node)-snw(1))**2+(y(node)-snw(2))**2
     .                  +(z(node)-snw(3))**2)
              if (dm>dmove) dmove = dm
              if (node.ne.indx) then
	         if (node>ntemp) then
                    x(node) = snw(1)
                    y(node) = snw(2)
                    z(node) = snw(3)
                 endif
              endif
           enddo 
           if (dmove.lt.eps) goto 10  
       	enddo


C  Write the final coornidate of nodes into file "scvt_lloyd.dat"
10      print *,"Iter = ",iloop-1," Maximal movement = ",dmove
        if (dmove.gt.eps) print *,"The maixmum iterations are reached!"
	
        open (17,file='scvt_lloyd.dat',status='unknown')
      	write(17,100) n
      	do node = 1,n             
           write(17,200) node,x(node),y(node),z(node)
      	enddo
      	close(17) 
100     format(I9)
200     format(I9,3X,F16.8,3X,F16.8,3X,F16.8)


C  Write the Delaunay triangles into the file "deltri.dat"
C  Write the Voronoi regions into the file "voronoi.dat"
      	ic = 1
      	call grid(n,x,y,z,xc,yc,zc,neigh,neisz,vortx,vorsz,ic)

      	end program SCVT_By_Lloyd
