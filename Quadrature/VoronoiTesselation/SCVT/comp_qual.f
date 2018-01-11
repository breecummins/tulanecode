        program Comp_Quality
        include "scvt.m"
	parameter(pi=3.1415926)

        open(16,file='nodes.dat',status='unknown')
        read(16,*) n
        if (n.ge.nmax) then
           print *,"The number of generators must be less than ",nmax
           stop
        endif
        do node = 1,n
           read(16,*) ntemp,x(node),y(node),z(node)
        enddo
        close(16)
        print *,"Numebr of Generators = ",n

C  Generate the mesh 
        call grid(n,x,y,z,xc,yc,zc,neigh,neisz,vortx,vorsz,ic)
       
C  Measure the quality 
        dist_max = 0
        dist_min = 10000
        area_max = 0
        area_min = 10000
	angle_min = 10000
	angle_max = 0
        iobangles = 0
        do node = 1,n
           V_area = 0
           ns = neisz(node)
           v1(1) = x(node)
           v1(2) = y(node)
           v1(3) = z(node)
           do j = 1,ns
              ind_neig = neigh(node,j)
              v2(1) = x(ind_neig)
              v2(2) = y(ind_neig)
              v2(3) = z(ind_neig)
              dist = (v1(1)-v2(1))**2+(v1(2)-v2(2))**2+(v1(3)-v2(3))**2
              dist = sqrt(dist)
              if (dist<dist_min) dist_min = dist
              if (dist>dist_max) dist_max = dist
              i1 = j
              i2 = mod(i1,ns)+1
              v2(1) = xc(vortx(node,i1))
              v2(2) = yc(vortx(node,i1))
              v2(3) = zc(vortx(node,i1))
              v3(1) = xc(vortx(node,i2))
              v3(2) = yc(vortx(node,i2))
              v3(3) = zc(vortx(node,i2))
              V_area = V_area+areas(v1,v2,v3)
              i2 = neigh(node,j)
              i3 = neigh(node,mod(j,ns)+1)
              if ((node<i2).and.(node<i3)) then
                 v2(1) = x(i2)
                 v2(2) = y(i2)
                 v2(3) = z(i2)
                 v3(1) = x(i3)
                 v3(2) = y(i3)
                 v3(3) = z(i3)
                 call get_angles(v1,v2,v3,angles)
                 if (angles(1)>angle_max) angle_max = angles(1)
                 if (angles(1)<angle_min) angle_min = angles(1)
                 if (angles(1)>0.5*pi) iobangles = iobangles+1
                 if (angles(2)>angle_max) angle_max = angles(2)
                 if (angles(2)<angle_min) angle_min = angles(2)
                 if (angles(2)>0.5*pi) iobangles = iobangles+1
                 if (angles(3)>angle_max) angle_max = angles(3)
                 if (angles(3)<angle_min) angle_min = angles(3)
                 if (angles(3)>0.5*pi) iobangles = iobangles+1
              endif
           enddo
           if (V_area<area_min) area_min = V_area
           if (V_area>area_max) area_max = V_area
        enddo
	angle_max = angle_max*180/pi
        angle_min = angle_min*180/pi
        print *,'Minimal Angle = ',angle_min
        print *,'Maximal Angle = ',angle_max
        print *,'*** ','Number of Obtuse Triangles = ',iobangles,'***'
        print *,'Maximal Dist = ',dist_max,' Minimal Dist = ',dist_min
        print *,'*** ','Dist Ratio = ',dist_max/dist_min,'***'
        print *,'Maximal Area = ',area_max,' Minimal Area = ',area_min
        print *,'*** ','Area Ratio = ',area_max/area_min,'***'

	end program Comp_Quality
