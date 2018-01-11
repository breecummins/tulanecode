%DifferentBlobs.m

clear

%get sphere discretization
load voronoivertices2400.dat
vv = voronoivertices2400;
rad = 1; %fortran code discretizes sphere of radius one.


% GIVE THE SPHERE A POSITIVE X VELOCITY.
u = zeros(size(vv,1)*size(vv,2),1);
umag = 1;
u(1:3:end) = umag; 


% DEFINE PARAMETERS
tv = vv(2:end,:);
tv(:,1) = tv(:,1) - vv(1,1);
tv(:,2) = tv(:,2) - vv(1,2);
tv(:,3) = tv(:,3) - vv(1,3);
epsilon = min(calcnorm(tv))./2;
mu = 1; %fluid viscosity
blob =1;

%test different spread
epsilon=epsilon./5;

% CALCULATE FORCES
M = constructM(vv,vv,epsilon,size(vv,1),size(vv,1),blob);
f = mu*M\u;
F=zeros(size(vv.'));
F(:)=f;
forces=-F.';

% CALCULATE DRAG
Da = sum(forces)

save ~/rsyncfolder/data/Quadrature/voronoi2400blob1eps5.mat vv rad u epsilon mu blob forces Da


