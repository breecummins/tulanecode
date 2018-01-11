      	subroutine grid(n,x,y,z,xc,yc,zc,neigh,neisz,vortx,vorsz,ic)
      	include "scvt.m"

C  Generate the Delaunay triangules
      	CALL TRMESH (N,X,Y,Z,LIST,LPTR,LEND,LNEW,IWK,IWK(N+1),DS,IER)

      	DO NODE = 1,N   
           LPL = LEND(NODE)
           LP = LPL
           K = 0
10         K = K + 1
           LP = LPTR(LP)
           ND = LIST(LP)
           NEIGH(NODE,K) = ND
           IF (LP .NE. LPL) GO TO 10
           NEISZ(NODE) = K
      	ENDDO  
      	K = 1
      	DO NODE = 1,N
           NS = NEISZ(NODE)
           DO J = 1,NS
              J1 = J
              J2 = MOD(J,NS)+1 
              NV1 = NEIGH(NODE,J1)
              NV2 = NEIGH(NODE,J2)
              IF ((NV1.GT.NODE).AND.(NV2.GT.NODE)) then
                 NTRI(K,1) = NODE
                 NTRI(K,2) = NV1
                 NTRI(K,3) = NV2
                 K = K+1
              ENDIF
           ENDDO
      	ENDDO
        NT = K-1

C  Generate the Voronoi diagrams
      	CALL CRLIST (N,NCOL,X,Y,Z,LIST,LEND, LPTR,LNEW,
     .               LBTRI, LISTC,NB,XC,YC,ZC,RC,IER)
	DO NODE = 1,N
           NV = 0
           LPL = LEND(NODE)
           LP = LPL
20         LP = LPTR(LP)
           KT = LISTC(LP)
           NV = NV + 1
           IWK(NV) = KT
           VORTX(NODE,NV) = KT 
           IF (LP .NE. LPL) GO TO 20
           VORSZ(NODE) = NV
      	ENDDO 

C  Save the Delaunay triangules and Voronoi diagrams if ic==1
        if (ic.eq.1) then
           open(15,file='deltri.dat',status='unknown')
           write(15,*) "Number of nodes:",N
           do node = 1,N
              write(15,200) node,X(node),Y(node),Z(node)
           enddo
           write(15,*) "Number of triangles:",NT
           do i = 1,NT
              write(15,300) i,NTRI(i,1),NTRI(i,2),NTRI(i,3)
           enddo
           close(15)

           open(35,file='points.dat',status='unknown')
           write(35,*) N
           do node = 1,N
              write(35,100) X(node),Y(node),Z(node)
           enddo
           close(35)

           open(45,file='triangles.dat',status='unknown')
           do i = 1,NT
              write(45,300) i,NTRI(i,1),NTRI(i,2),NTRI(i,3)
           enddo
           close(45)

100     FORMAT(3X,F16.8,3X,F16.8,3X,F16.8)
200   	FORMAT(I9,3X,F16.8,3X,F16.8,3X,F16.8)
300   	FORMAT(4I10)


           open(16,file='voronoi.dat',status='unknown')
           write(16,*) "Number of generators:",N
           do node = 1,N
              write(16,200) node,X(node),Y(node),Z(node)
           enddo
           write(16,*) "Neighbor generators:"
           do node = 1,N
              NS = NEISZ(node)
              write(16,400) NS,(NEIGH(node,i),i=1,NS)
           enddo
           write(16,*) "Voronoi vortices:"
           do node = 1,N
              NS = NEISZ(node)
              write(16,400) NS,(VORTX(node,i),i=1,NS)
           enddo
400   	FORMAT(16I9)
500     FORMAT(16I2)
           write(16,*) "Number of vortices:",NT
           do node = 1,NT
              write(16,200) node,XC(node),YC(node),ZC(node)
           enddo   
           close(16)

           open(46,file='edges.dat',status='unknown')
           do node = 1,N
              NS = NEISZ(node)
              do i = 1,NS
                if (NEIGH(node,i)>node) then  
                    write(46,400) node,NEIGH(node,i)
                endif
              enddo
           enddo
	       close(46)
        endif

	end subroutine grid
