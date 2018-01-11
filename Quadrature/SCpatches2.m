%SCpatches2.m

clear 
close all

plotvblob = 0;
plotvblobnf = 1;
plotbestblob = 0;
plotconv = 0;

base1 = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_farnearpatches_mostbps';
base2 = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_farnearpatches_mostbps_originalmethod';

l2errBEM = [];
l2errOrig = [];
ctr=0;
for N = 2.^(3:8);
	ctr=ctr+1;
	str = sprintf('%04d',N);
	
	load([base1,str])
	ptsperpatch =size(obspts,1)/2;
	
	l2 = [];
	for k = ptsperpatch:ptsperpatch:length(sqrerr);
		patcherr = sqrerr(k-ptsperpatch+1:k);
		if ~all(size(patcherr) == [1 ,25]);
			size(patcherr)
			error('There is a bug.')
		end
		l2(end+1) = sqrt(sum( patcherr ));  
	end
	l2errBEM(:,:,ctr) = [l2(1:2:end).',l2(2:2:end).'];

	load([base2,str])
	ptsperpatch =size(obspts,1)/2;
		
	l2_orig = [];
	for k = ptsperpatch:ptsperpatch:length(sqrerr);
		patcherr = sqrerr(k-ptsperpatch+1:k);
		if ~all(size(patcherr) == [1 ,25]);
			size(patcherr)
			error('There is a bug.')
		end
		l2_orig(end+1) = sqrt(sum( patcherr ));
	end
	l2errOrig(:,:,ctr)=[l2_orig(1:2:end).',l2_orig(2:2:end).'];
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
	plot(blobparams,squeeze(l2errBEM(:,1,:)))
	title('Hat/trap BEM, far patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
	figure(2)
	plot(blobparams,squeeze(l2errBEM(:,2,:)))
	title('Hat/trap BEM, near patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
	figure(3)
	plot(blobparams,squeeze(l2errOrig(:,1,:)))
	title('Original method, far patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
	figure(4)
	plot(blobparams,squeeze(l2errOrig(:,2,:)))
	title('Original method, near patch')
	xlabel('Blob parameters')
	ylabel('L_2 error')
	legend(leg,'Location','BestOutside')
end

if plotvblobnf == 1;
	set(0,'DefaultAxesFontSize',16)
	set(0,'DefaultLineLineWidth',1.5)

	leg={'Far field','Near field'};
	bind =5;

	figure
	plot(blobparams,squeeze(l2errBEM(:,1,bind)),'b.')
	hold on
	plot(blobparams,squeeze(l2errBEM(:,2,bind)),'r.')
	title(['Hat/trap BEM, N = ', int2str(N(bind))])
	xlabel('Blob parameters')
	ylabel('L_2 error')
	xlim([0,0.014])
	legend(leg,'Location','BestOutside')
	figure
	plot(blobparams,squeeze(l2errBEM(:,1,bind)),'b.')
	hold on
	plot(blobparams,squeeze(l2errBEM(:,2,bind)),'r.')
	title(['Original method, N = ', int2str(N(bind))])
	xlabel('Blob parameters')
	ylabel('L_2 error')
	xlim([0,0.014])
	legend(leg,'Location','BestOutside')
end

if plotbestblob == 1;
	bp=[];
	bpBEM=[];
	for k = 1:2;
		err = squeeze(l2errOrig(:,k,:));
		for l = 1:6;
			ind = find(err(:,l)==min(err(:,l)));
			bp(l,k)=blobparams(ind);
		end
		err = squeeze(l2errBEM(:,k,:));
		for l = 1:6;
			ind = find(err(:,l)==min(err(:,l)));
			bpBEM(l,k)=blobparams(ind);
		end
	end

	leg={'orig','0.21N','hat/trap','0.17N^{3/2}'};

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

	leg={'orig','0.16N','hat/trap','0.15N^{3/2}'};

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


if plotconv == 1;
	% bv = abs(blobparams - 0.2);
	% %ind = find(bv == min(bv));
	ind = 50;
	blobp = blobparams(ind)
	indBEM = 50;
	blobpBEM = blobparams(indBEM)

	%t = { 'Hat/trap BEM, far', 'Hat/trap BEM, near', 'Original method, far', 'Original method, near' };
	t = { 'Original method, far', 'Original method, near' };
	tBEM = { 'Hat/trap BEM, far', 'Hat/trap BEM, near'};

	errmat=[];
	errmatBEM=[];
	coeffs = [];
	coeffsBEM = [];

	fg = linspace(min(N),max(N),50);

	for k = 1:2;
		err = squeeze(l2errOrig(ind,k,:));
		p = polyfit(log(N).',log(err),1);
	
		figure
		loglog(N, err,'.')
		hold on
		loglog(fg,exp(p(2)).*fg.^(p(1)),'r-')
		xlabel('N')
		ylabel('L_2 error')
		title(t{k})
	
		errmat(:,end+1) = err;
		coeffs(:,end+1) = p;

		err = squeeze(l2errBEM(indBEM,k,:));
		p = polyfit(log(N).',log(err),1);
	
		figure
		loglog(N, err,'.')
		hold on
		loglog(fg,exp(p(2)).*fg.^(p(1)),'r-')
		xlabel('N')
		ylabel('L_2 error')
		title(tBEM{k})
	
		errmatBEM(:,end+1) = err;
		coeffsBEM(:,end+1) = p;

	end

	power = abs(coeffs(1,:))
	powerBEM = abs(coeffsBEM(1,:))

	% cmat = [chordlen.^power(1),chordlen.^power(2),chordlen.^power(3),chordlen.^power(4)];
	cmat = [chordlen.^power(1),chordlen.^power(2)];
	ordermat = errmat./cmat;
	errmat;
end




% figure(5)
% ang=0:2*pi/500:2*pi;
% circ = 0.25*[cos(ang).',sin(ang).'];
% plot(circ(:,1),circ(:,2),'r-')
% hold on
% plot(obspts(:,1),obspts(:,2),'b.')
% title('Near and far patches')	
% 
% figure(6)
% plot(N,squeeze(l2err(1,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams(1))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(7)
% plot(N,squeeze(l2err(50,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams(50))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(8)
% plot(N,squeeze(l2err(100,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams(100))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(9)
% plot(N,squeeze(l2err(150,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams(150))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(10)
% plot(N,squeeze(l2err(200,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams(200))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% % figure(6)
% % plot(2.^(3:8),squeeze(l2err(1:50,1,:)).')
% % title('Hat/trap BEM, far patch, bps = 1.e-4 - 4.2e-4')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(7)
% % plot(2.^(3:8),squeeze(l2err(51:100,1,:)).')
% % title('Hat/trap BEM, far patch, bps = 4.36e-4 - 1.84e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(8)
% % plot(2.^(3:8),squeeze(l2err(101:150,1,:)).')
% % title('Hat/trap BEM, far patch, bps = 1.90e-3 - 8.03e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(9)
% % plot(2.^(3:8),squeeze(l2err(151:200,1,:)).')
% % title('Hat/trap BEM, far patch, bps = 8.27e-3 - 3.5e-2')
% % xlabel('N')
% % ylabel('L_2 error')
% % 
% % 
% % 
% % figure(10)
% % plot(2.^(3:8),squeeze(l2err(1:50,2,:)).')
% % title('Hat/trap BEM, near patch, bps = 1.e-4 - 4.2e-4')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(11)
% % plot(2.^(3:8),squeeze(l2err(51:100,2,:)).')
% % title('Hat/trap BEM, near patch, bps = 4.36e-4 - 1.84e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(12)
% % plot(2.^(3:8),squeeze(l2err(101:150,2,:)).')
% % title('Hat/trap BEM, near patch, bps = 1.90e-3 - 8.03e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(13)
% % plot(2.^(3:8),squeeze(l2err(151:200,2,:)).')
% % title('Hat/trap BEM, near patch, bps = 8.27e-3 - 3.5e-2')
% % xlabel('N')
% % ylabel('L_2 error')
% % 
% % 
% % 
% % figure(14)
% % plot(2.^(3:8),squeeze(l2err(1:50,3,:)).')
% % title('Orig method, far patch, bps = 1.e-4 - 4.2e-4')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(15)
% % plot(2.^(3:8),squeeze(l2err(51:100,3,:)).')
% % title('Orig method, far patch, bps = 4.36e-4 - 1.84e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(16)
% % plot(2.^(3:8),squeeze(l2err(101:150,3,:)).')
% % title('Orig method, far patch, bps = 1.90e-3 - 8.03e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(17)
% % plot(2.^(3:8),squeeze(l2err(151:200,3,:)).')
% % title('Orig method, far patch, bps = 8.27e-3 - 3.5e-2')
% % xlabel('N')
% % ylabel('L_2 error')
% % 
% % 
% % 
% % figure(18)
% % plot(2.^(3:8),squeeze(l2err(1:50,4,:)).')
% % title('Orig method, near patch, bps = 1.e-4 - 4.2e-4')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(19)
% % plot(2.^(3:8),squeeze(l2err(51:100,4,:)).')
% % title('Orig method, near patch, bps = 4.36e-4 - 1.84e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(20)
% % plot(2.^(3:8),squeeze(l2err(101:150,4,:)).')
% % title('Orig method, near patch, bps = 1.90e-3 - 8.03e-3')
% % xlabel('N')
% % ylabel('L_2 error')
% % figure(21)
% % plot(2.^(3:8),squeeze(l2err(151:200,4,:)).')
% % title('Orig method, near patch, bps = 8.27e-3 - 3.5e-2')
% % xlabel('N')
% % ylabel('L_2 error')
% % 
% % 
% % 
