%
%
% FIND THE NECESSARY FUNCTIONS
%
%
 syms r d real
 syms a b c real
 syms s real

 %start with a blob
 blob = 2*d^4/(pi*(r^2+d^2)^3);

%GEt G in two steps
 Gtmp = int(r*blob,r)/r + 1/(2*pi*r);
 G=int(Gtmp,r);
 
 %get Bprime
 Bp = int(r*G,r)/r;
 
 %now the H1 and H2
 H1 = Bp/r - G + 0/(8*pi) ;
 H1 = simplify( 8*pi*H1+1 )/(8*pi);
 H2 = simplify( (diff(Bp,r)*r - Bp)/r^3 );
 
 disp(' computed H1 and H2' )
 %%
 
 %Now find the integrals needed
 
 tmp= subs(H1,r,sqrt(a*s^2+b*s+c));
 tmp1 = int(s*tmp,s);
 I1 = subs(tmp1,s,1)-subs(tmp1,s,0);
 
 tmp= subs(H2,r,sqrt(a*s^2+b*s+c));
 tmp1 = int(s*tmp,s); 
 I2 = subs(tmp1,s,1)-subs(tmp1,s,0);
  
 tmp1 = int(s*s*tmp,s);
 I3 = subs(tmp1,s,1)-subs(tmp1,s,0);
 
 tmp1 = int(s*s*s*tmp,s);
 I4 = subs(tmp1,s,1)-subs(tmp1,s,0);
 
 disp(' computed the necessary integrals')
 %%
 
 
 %TEST TEST TEST
 % Constant forces on a circle
 %%%%%%%%%%%%%%%%%%%%%%%%
 N = 4;  %number of points around the circle
 h = 2*pi/N;
 ang = (0 : h : 2*pi-h/2)';
 x = cos(ang);
 y = sin(ang);
 
 %use some fake forces to compute the velocity
 f1 = 4*ones(N,1) * h;
 f2 = zeros(N,1)  * h;
 
 epsi = 0.1;%guess a good value of epsilon
 
 %eval point (xj,yj)
 xj = 1;   yj = 0;
         
 %make the necessary piece of matrix
 FullMatrix = zeros(2,2*N);
 
 
 for k = 1:N  %go around the circle
    matblock = zeros(2,2);
 
    %set the two neighbors of x(k)
    knext = [k-1,k+1]; 
    if k==1, knext = [N,2]; end
    if k==N, knext = [1,N-1]; end
 
    for kcount = 1:2
       km1 = knext(kcount);
 
       dxj = xj-x(km1);  dyj = yj-y(km1);
       dxk = x(k)-x(km1);  dyk = y(k)-y(km1);
 
       a1 = dxk^2 + dyk^2;
       b1 = -2*( dxj*dxk + dyj*dyk );
       c1 = dxj^2+dyj^2;
 
       C =  [ [dxj^2,   dxj*dyj]; ...
              [dxj*dyj, dyj^2] ];

       B = -[ [2*dxj*dxk,       dxj*dyk+dyj*dxk]; ...
              [dxj*dyk+dyj*dxk, 2*dyj*dyk] ];
          
       A =  [ [dxk^2,   dxk*dyk]; ...
              [dxk*dyk, dyk^2] ];
 
       matblock=matblock ...
     + subs(I1,{a,b,c,d},{a1,b1,c1,epsi})*eye(2) ...
     + subs(I2,{a,b,c,d},{a1,b1,c1,epsi})*C ...
     + subs(I3,{a,b,c,d},{a1,b1,c1,epsi})*B ...
     + subs(I4,{a,b,c,d},{a1,b1,c1,epsi})*A;
 
    end
    FullMatrix(:,2*k-1:2*k) = matblock;
 end
 
 FullMatrix

 force = [f1';f2'];  force = force(:);  uj = FullMatrix * force 
 plot(x,y,xj,yj,'r.'),axis equal,hold on
 quiver(xj,yj,uj(1),uj(2),0,'r');pause(.1)
 
 hold off