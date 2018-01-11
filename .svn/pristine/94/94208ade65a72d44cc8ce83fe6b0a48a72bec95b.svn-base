        program SCVT_By_MonteCarlo
        implicit real(a-h,o-z)
        real,allocatable :: nodes_crd(:,:)
	parameter(max_test=10000)

C  Read the number of generators ...
        open (15,file='scvt.in',status='unknown')
        read (15,*)
        print *,'Reading the number of generators ...'
        read (15,*) n
        print *,'Number of Generators =',n
        close(15)
 
C  Allocate the memory space ...
        allocate(nodes_crd(3,n))

C  Set the seed for random number generator
        nrank = 0
        call set_random_seed(nrank)
C  Find the maximum of the density function
        dens_max = get_max_density(max_test)

C  Keep fixed points
        open (35,file='fixedpoints.m',status='unknown')
        read (35,*) ntemp
        print *,"Number of Fixed Generators = ",ntemp
	if (ntemp>0) then
       	  do node = 1,ntemp
           read(35,*) nodes_crd(1,node),nodes_crd(2,node),
     .                nodes_crd(3,node)
      	  enddo
        endif
        close(35)

C  Generate the SCVT using Monte Carlo method
        do node = ntemp+1,n
           call random_generator(nodes_crd(1,node),dens_max)
        enddo
     
C  Write the coordinates of the resulted SCVT
        open (16, file='scvt_mc.dat',status='unknown')
        write(16,100) n
        do node = 1,n            
           write(16,200) node,nodes_crd(1:3,node)
        enddo
        close(16)     
100    	format(I9)
200     format(I9,3X,F16.8,3X,F16.8,3X,F16.8) 

C  Free the memory space...
        deallocate(nodes_crd)

        end program SCVT_By_MonteCarlo


C*******************************************************************
C  Get the maximum of the density function dens_f
C*******************************************************************
        function get_max_density(max_t)
        real    get_max_density,pt_test(3),ff
        integer i,max_t
                                                                                
        get_max_density = 0.0
        do i = 1,max_t
           call unif_random(pt_test)
           ff = dens_f(pt_test)
           if (ff.gt.get_max_density) get_max_density = ff
        enddo
                                                                                
        end function get_max_density
