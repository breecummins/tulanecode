C*******************************************************************
C     Define the density function for the SCVT
C*******************************************************************
        function dens_f(pt)

        real  pt(3),dens_f
        
        dens_f = 1.0
c        if ((pt(3)>0.65) .or. (pt(3)<-0.65)) then
c		dens_f = 2.0
c        endif

        end function dens_f
