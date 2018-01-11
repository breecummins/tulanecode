function Bree_ChoanoCode();

s.L =0.25; s.F=1; s.a=0.4; s.delh = 0.5;  s.delt = 0.55; s.critdist = 2*s.delh; s.mu=1; 
s.forces = @(ff,fnew) [ff+[fnew;-fnew;-fnew;fnew]];
%s.forces = @(ff,fnew) [ff+[fnew;-fnew;-0*fnew;0*fnew]];

xh0 = [ 2 0 0 ; 3.5 0 0];
v0 = [ 0 1 0 ; 0 1 0];
% [theta0,phi0,r0] = cart2sph(v0(:,1),v0(:,2),v0(:,3));
% angs = [theta0,phi0].';
theta0 = atan2(v0(:,2),v0(:,1));
Y0 = xh0.';
y0 = [Y0(:);theta0];

% plot(xh0(:,1),xh0(:,2),'b.','MarkerSize',20), hold on
% plot(xt0(:,1),xt0(:,2),'r.','MarkerSize',20)
% keyboard

[t,y] = ode45(@RRForces,[0:20:800],y0,[],s);

d=[];
figure(1)
for k = 1:length(t);
	xh = reshape(y(k,1:6).',3,[]).'; 
	tmpd = sqrt( sum((xh(1,:)-xh(2,:)).^2) );
    d = [d; tmpd];
	theta = y(k,7:8).'; 
	b = [cos(theta),sin(theta),zeros(size(theta))];
	xt = xh-b*s.L;
	ff = [s.F*b;-s.F*b];
    %check distance between heads and add forces
    if (tmpd<s.critdist), 
        fnew = s.a*(1-tmpd/s.critdist)*(xh(1,:)-xh(2,:))/tmpd;
        ff = s.forces(ff,fnew);
    end
	%plot critters
	plot(xh(:,1),xh(:,2),'b.','MarkerSize',20), hold on
 	plot(xt(:,1),xt(:,2),'r.','MarkerSize',20)
 	quiver(xh(:,1),xh(:,2),ff(1:2,1),ff(1:2,2),'r')
 	quiver(xt(:,1),xt(:,2),ff(3:4,1),ff(3:4,2),'r')

	xx = [xh;xt];
 	xmin = min(xx(:,1));  xmax = max(xx(:,1));
 	ymin = min(xx(:,2));  ymax = max(xx(:,2));
 	hold off,axis equal,axis([xmin-4*s.L xmax+4*s.L ymin-4*s.L ymax+4*s.L])
 	grid on,title(['time = ',num2str(t(k))])
 	pause(.02)
end

figure(2), hold on
plot( t.',d,'r' ),xlabel('time'),ylabel('distance')

end %function

function yp = RRForces(t,y,s);
	
	% y = [xh, b], where b is orientation
	% yp = [dxh/dt, db/dt]
	
	Np = 2;
	dh2 = s.delh^2;  dt2 = s.delt^2;
		
	xh = reshape(y(1:6).',3,[]).'; 
	theta = y(7:8); 
	b = [cos(theta),sin(theta),zeros(size(theta))];
	xt = xh-b*s.L;
    xx = [xh;xt]; ff = [s.F*b;-s.F*b];
	
    %check distance between heads and add forces
    tmpd = sqrt( sum((xh(1,:)-xh(2,:)).^2) );
    if (tmpd<s.critdist), 
        fnew = s.a*(1-tmpd/s.critdist)*(xh(1,:)-xh(2,:))/tmpd;
        ff = s.forces(ff,fnew);
    end

    u = zeros(size(xx(:,1)));  v = u;   w = u;

    for k=1 : Np
      dx = xx(:,1)-xx(k,1);
      dy = xx(:,2)-xx(k,2);
      dz = xx(:,3)-xx(k,3);
   
      r2 = dx.^2 + dy.^2 + dz.^2;
      R  = sqrt(r2+dh2);
      % All the relevant functions
      H1 = (1./R + dh2./R.^3)/(8*pi*s.mu);
      H2 = (1./R.^3)/(8*pi*s.mu);
      fdotx = ff(k,1)*dx + ff(k,2)*dy + ff(k,3)*dz;

      u = u + ff(k,1)*H1 + fdotx.*dx.*H2;
      v = v + ff(k,2)*H1 + fdotx.*dy.*H2;
      w = w + ff(k,3)*H1 + fdotx.*dz.*H2;
      
      % now do the tails
      dx = xx(:,1)-xx(Np+k,1);
      dy = xx(:,2)-xx(Np+k,2);
      dz = xx(:,3)-xx(Np+k,3);
   
      r2 = dx.^2 + dy.^2 + dz.^2;
      R  = sqrt(r2+dt2);
      % All the relevant functions
      H1 = (1./R + dt2./R.^3)/(8*pi*s.mu);
      H2 = (1./R.^3)/(8*pi*s.mu);
      fdotx = ff(Np+k,1)*dx + ff(Np+k,2)*dy + ff(Np+k,3)*dz;

      u = u + ff(Np+k,1)*H1 + fdotx.*dx.*H2;
      v = v + ff(Np+k,2)*H1 + fdotx.*dy.*H2;
      w = w + ff(Np+k,3)*H1 + fdotx.*dz.*H2;
    end

	U = [u,v,w];
	uh = U(1:2,1); 	vh = U(1:2,2); 	wh = U(1:2,3);
	ut = U(3:4,1); 	vt = U(3:4,2); 	wt = U(3:4,3);
	dtheta = ( cos(theta) .* (vh-vt) - sin(theta) .* (uh-ut) )/s.L;
	Yp = U(1:2,:).';
	yp = [Yp(:);dtheta];

end %function	

