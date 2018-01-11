
% THIS FILE SHOWS AN ANIMATION OF THE SOLUTION
% OF THE 'HEAT' EQUATION  u_t = u_xx USING
% CRANK-NICOLSON.  THE EQUATION IS SOLVED BETWEEN
% x=0 AND x=1 WITH BC GIVEN BY phiL(t) AND phiR(t).
% INITIAL CONDITIONS ARE GIVEN BY f(x).

clear

M=64; h = 1/(M);  %SPATIAL DISCRETIZATION
x=(h/2:h:1)'; 
xaug = [-h/2;x;1+h/2];
% xL = -h/2;  xR = 1+h/2;
%------------------------
a = 1;  b = 0.01;
uex  =@(x,t) 0.1*          exp(-(x-1/2).^2/(2*(b+2*a*t)))/sqrt(b+2*a*t);
duex =@(x,t) -0.1*(x-1/2).*exp(-(x-1/2).^2/(2*(b+2*a*t)))/((b+2*a*t)^(3/2));
%------------------------
Un = [uex(h/2,0)-h*duex(0,0);uex(x,0);uex(1-h/2,0)+h*duex(1,0)];


%SET TIME STEP
dt = h/10;  %TIME STEP
kappa = a*dt/(2*h^2);

A = zeros(M+2,M+2);
rhs = zeros(M+2,1);
for i=1:M+2
  A(i,i) = 1+2*kappa;
end
for i=1:M+1
  A(i+1,i) = -kappa;
  A(i,i+1) = -kappa;
end
A(1,1) = 1;  A(1,2) = -1; A(1,3:end) = zeros(1,M);
A(M+2,M+2) = 1;  A(M+2,M+1) = -1;  A(M+2,1:M) = zeros(1,M);


figure(1)
Tf = 100.0;  nt = fix(Tf/dt);
tt = (0:dt:(nt-1)*dt);
for j=1:nt
	 tn = (j-1)*dt; tnp1 = j*dt;
	err = Un - uex(xaug,tn); %plot(xaug,err),axis([-h 1+h -2e-3 2e-3]),pause(.1)
	err2(j) = sqrt( h*sum(abs(err.^2)) );
	 %--------------------
	%   if mod(j,1)==0,
	%   xaug = [-h/2;x;1+h/2];
	%   plot(xaug,Un,'b'),hold on
	%   plot(xaug,Un,'r.'),
	%   plot(xaug,uex(xaug,tn),'k--'),hold off
	%   axis([-2*h 1+2*h -1 2]),grid
	%   xlabel('x')
	%   title(['TIME = ',num2str(tn)]),pause(.1)
	%   end
	%----------------
	% SET THE NEW RHS

	 rhs(1) = -h*duex(0,tnp1);
	 for i=2:M+1
	   rhs(i) = Un(i) + kappa*( Un(i-1)-2*Un(i)+Un(i+1) );
	 end
	 rhs(M+2) = h*duex(1,tnp1);
 
	 Un = A\rhs;
end
figure(2),plot(tt,err2,'r'),hold on,grid on
