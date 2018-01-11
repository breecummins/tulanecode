%SCpatches.m

clear 
close all

load ~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTestExactSoln_farnearpatches
exactu = exactsoln.data(:,3);
exactv = exactsoln.data(:,4);
ptsperpatch = size(exactu,1)/2;

base1 = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_farnearpatches';
base2 = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_farnearpatches_originalmethod';

l2err = [];
ctr=0;
for N = 2.^(3:8);
	ctr=ctr+1;
	str = sprintf('%04d',N);
	
	load([base2,str])
	uverr_orig = u_v_err.data;
	
	load([base1,str])
	uverr = u_v_err.data;
	
	l2 = [];
	l2_orig = [];
	%tog = 1;
	for k = ptsperpatch:ptsperpatch:size(uverr,1);
		% long version commented out; used to check the short version below
		% if tog ==1;
		% 	eu = exactu(1:25);
		% 	ev = exactv(1:25);
		% elseif tog == 0;
		% 	eu = exactu(26:end);
		% 	ev = exactv(26:end);
		% else;
		% 	error('Toggle is failing.')
		% end
		% l2(end+1) = sqrt(sum( (uverr(k-ptsperpatch+1:k,1)-eu).^2 + (uverr(k-ptsperpatch+1:k,2)-ev).^2 ));  
		% l2_orig(end+1) = sqrt(sum( (uverr_orig(k-ptsperpatch+1:k,1)-eu).^2 + (uverr_orig(k-ptsperpatch+1:k,2)-ev).^2 ));
		% tog =mod(tog+1,2);
		l2(end+1) = sqrt(sum( uverr(k-ptsperpatch+1:k,3) ));  %for farnearpatches files, the error column is (u-uexact).^2 + (v-vexact).^2. so to get the L2 norm, sum all of these over the patch and take the square root
		l2_orig(end+1) = sqrt(sum( uverr_orig(k-ptsperpatch+1:k,3) ));
	end
	l2err(:,:,ctr)=[l2(1:2:end).',l2(2:2:end).',l2_orig(1:2:end).',l2_orig(2:2:end).'];
end

bv = abs(blobparams.data - 0.0025);
ind = find(bv == min(bv));
bp = blobparams.data(ind)

t = { 'Hat/trap BEM, far', 'Hat/trap BEM, near', 'Original method, far', 'Original method, near' };

errmat=[];
coeffs = [];

ff=0;

chordlen = 2*sin(pi./(2.^((3+ff):8))).'; %skip N = 8
fg = linspace(min(chordlen),max(chordlen),50);

for k = 1:4;
	err = squeeze(l2err(ind,k,(ff+1):end)); %skip N = 8
	p = polyfit(chordlen,err,2);
	
	figure
	hold on
	plot(chordlen, err,'.')
	plot(fg, polyval(p,fg),'r-')
	xlabel('chord length')
	ylabel('L2 error')
	title(t{k})
	
	errmat(:,end+1) = err;
	coeffs(:,end+1) = p;
end

coeffs

power = 2*(coeffs(1,:) > coeffs(2,:)) + (coeffs(1,:) <= coeffs(2,:))

cmat = [chordlen.^power(1),chordlen.^power(2),chordlen.^power(3),chordlen.^power(4)];
ordermat = errmat./cmat

% clrs = {'k','b','r','g','c','m'};
% set(0,'DefaultAxesFontSize',16)
% set(0,'DefaultLineLineWidth',1.5)
% 
% leg={};
% for k = 1:6;
% 	N = 2.^(k+2);
% 	chordlen = 2*sin(pi/N);
% 	str = sprintf('%0.3g',chordlen);
% 	leg{end+1} = ['N=',int2str(N),',',sprintf('\n'),' clen=',str];
% end
% 
% figure(1)
% plot(blobparams.data,squeeze(l2err(:,1,:)))
% title('Hat/trap BEM, far patch')
% xlabel('Blob parameters')
% ylabel('L_2 error')
% legend(leg,'Location','BestOutside')
% figure(2)
% plot(blobparams.data,squeeze(l2err(:,2,:)))
% title('Hat/trap BEM, near patch')
% xlabel('Blob parameters')
% ylabel('L_2 error')
% legend(leg,'Location','BestOutside')
% figure(3)
% plot(blobparams.data,squeeze(l2err(:,3,:)))
% title('Original method, far patch')
% xlabel('Blob parameters')
% ylabel('L_2 error')
% legend(leg,'Location','BestOutside')
% figure(4)
% plot(blobparams.data,squeeze(l2err(:,4,:)))
% title('Original method, near patch')
% xlabel('Blob parameters')
% ylabel('L_2 error')
% legend(leg,'Location','BestOutside')
% 
% figure(5)
% ang=0:2*pi/500:2*pi;
% circ = 0.25*[cos(ang).',sin(ang).'];
% plot(circ(:,1),circ(:,2),'r-')
% hold on
% plot(obspts.data(:,1),obspts.data(:,2),'b.')
% title('Near and far patches')	
% 
% figure(6)
% plot(2.^(3:8),squeeze(l2err(1,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams.data(1))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(7)
% plot(2.^(3:8),squeeze(l2err(50,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams.data(50))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(8)
% plot(2.^(3:8),squeeze(l2err(100,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams.data(100))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(9)
% plot(2.^(3:8),squeeze(l2err(150,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams.data(150))])
% legend({'Hat/trap, far', 'Hat/trap, near', 'Orig, far', 'Orig, near'})
% xlabel('N')
% ylabel('L_2 error')
% axis([0,256,0,0.2])
% 
% figure(10)
% plot(2.^(3:8),squeeze(l2err(200,1:4,:)).')
% title(['Blob parameter = ',num2str(blobparams.data(200))])
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
