%plotDrag.m

clear

D1=[];
D1de=[];
D2=[];

load ~/rsyncfolder/data/Quadrature/voronoi0300blob1.mat
D1(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi0300blob1eps5.mat
D1de(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi0300blob2.mat
D2(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi0600blob1.mat
D1(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi0600blob1eps5.mat
D1de(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi0600blob2.mat
D2(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi1200blob1.mat
D1(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi1200blob1eps5.mat
D1de(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi1200blob2.mat
D2(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi2400blob1.mat
D1(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi2400blob1eps5.mat
D1de(end+1,:) = Da;
load ~/rsyncfolder/data/Quadrature/voronoi2400blob2.mat
D2(end+1,:) = Da;

xax=[300,600,1200,2400];
De = [-6*pi*mu*max(u)*rad,0,0];

figure
plot(xax,D1(:,1),'b-')
hold on
plot(xax,D2(:,1),'r-')
plot(xax,D1de(:,1),'k-')
plot(xax,De(1)*ones(size(xax)),'b--')
hold off
title('x direction')
xlabel('# of generators')
ylabel('Total drag')
legend('Blob 1', 'Blob 2','Blob 1, eps/5','Exact')


figure
plot(xax,D1(:,2),'b-')
hold on
plot(xax,D2(:,2),'r-')
plot(xax,D1de(:,2),'k-')
plot(xax,De(2)*ones(size(xax)),'b--')
hold off
title('y direction')
xlabel('# of generators')
ylabel('Total drag')
legend('Blob 1', 'Blob 2','Blob 1, eps/5','Exact')


figure
plot(xax,D1(:,3),'b-')
hold on
plot(xax,D2(:,3),'r-')
plot(xax,D1de(:,3),'k-')
plot(xax,De(3)*ones(size(xax)),'b--')
hold off
title('z direction')
xlabel('# of generators')
ylabel('Total drag')
legend('Blob 1', 'Blob 2','Blob 1, eps/5','Exact')


