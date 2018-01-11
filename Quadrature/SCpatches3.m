%SCpatches3.m

%use this code for blob parameters that change with N.  othrewise use SCpatches2.m

clear 
close all

plotvblob = 1;
plotbestblob = 1;

base1 = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_farnearpatches_bpwithN';
base2 = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_farnearpatches_bpwithN_originalmethod';

l2errBEM = [];
l2errOrig = [];
allbps = [];
allbpsBEM = [];
ctr=0;

for N = 2.^(3:8);
	ctr=ctr+1;
	str = sprintf('%04d',N);
	
	load([base1,str])
	ptsperpatch =size(obspts,1)/2;
	
	l2 = [];
	for k = ptsperpatch:ptsperpatch:length(sqrerr);
		l2(end+1) = sqrt(sum( sqrerr(k-ptsperpatch+1:k) ));  
	end
	l2errBEM(:,:,ctr) = [l2(1:2:end).',l2(2:2:end).'];

	allbpsBEM(:,ctr) = blobparams;

	load([base2,str])
	ptsperpatch =size(obspts,1)/2;
		
	l2_orig = [];
	for k = ptsperpatch:ptsperpatch:length(sqrerr);
		l2_orig(end+1) = sqrt(sum( sqrerr(k-ptsperpatch+1:k) ));
	end
	l2errOrig(:,:,ctr)=[l2_orig(1:2:end).',l2_orig(2:2:end).'];
	
	allbps(:,ctr) = blobparams;
end

N = (2.^(3:8));
chordlen = 2*0.25*sin(pi./N).'; 

if plotvblob == 1;
	clrs = {'k','b','r','g','c','m'};
	set(0,'DefaultAxesFontSize',16)
	set(0,'DefaultLineLineWidth',1.5)

	leg={};
	for k = 1:length(N);
		leg{end+1} = ['N=',int2str(N(k)),',',sprintf('\n'),' clen=',sprintf('%0.3g',chordlen(k))];
	end

	figure(1)
	plot(allbpsBEM,squeeze(l2errBEM(:,1,:)))
	title('Hat/trap BEM, far patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
	figure(2)
	plot(allbpsBEM,squeeze(l2errBEM(:,2,:)))
	title('Hat/trap BEM, near patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
	figure(3)
	plot(allbps,squeeze(l2errOrig(:,1,:)))
	title('Original method, far patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
	figure(4)
	plot(allbps,squeeze(l2errOrig(:,2,:)))
	title('Original method, near patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
end

if plotbestblob == 1;
	bp=[];
	bpBEM=[];
	for k = 1:2;
		err = squeeze(l2errOrig(:,k,:));
		for l = 1:6;
			ind = find(err(:,l)==min(err(:,l)));
			bp(l,k)=allbps(ind,l);
		end
		err = squeeze(l2errBEM(:,k,:));
		for l = 1:6;
			ind = find(err(:,l)==min(err(:,l)));
			bpBEM(l,k)=allbpsBEM(ind,l);
		end
	end

	leg={'orig','0.21clen','hat/trap','0.17clen^{3/2}'};

	figure
	plot(chordlen,bp(:,1),'b.')
	hold on
	plot(chordlen,0.21*chordlen,'r')
	plot(chordlen,bpBEM(:,1),'k.')
	plot(chordlen,0.17*chordlen.^(3/2),'m-')
	title('blob param v chord len, far patch')
	xlabel('chordlen')
	ylabel('blob param')
	legend(leg,'Location','Northwest')

	leg={'orig','0.16clen','hat/trap','0.15clen^{3/2}'};

	figure
	plot(chordlen,bp(:,2),'b.')
	hold on
	plot(chordlen,0.16*chordlen,'r')
	plot(chordlen,bpBEM(:,2),'k.')
	plot(chordlen,0.15*chordlen.^(3/2),'m-')
	title('blob param v chord len, near patch')
	xlabel('chordlen')
	ylabel('blob param')
	legend(leg,'Location','Northwest')
end		


