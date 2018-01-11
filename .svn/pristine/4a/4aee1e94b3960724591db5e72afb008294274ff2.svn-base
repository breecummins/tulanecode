%viewvelocity3.m
	
	clear 
	close all
	
	
	patchspacing = 0.02/5;
	
	bestblobfarorig=[];
	bestblobnearorig=[];
	bestblobfarexactint=[];
	bestblobnearexactint=[];
	
	l2errfarorig=[];
	l2errnearorig=[];
	l2errfarexactint=[];
	l2errnearexactint=[];
	
	n=0;
	for N = 2.^(3:8);
		n=n+1;
		str = sprintf('%04d',N);
		load(['~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_exactintegralsANDorig_patches_singleblobfortable_N',str])
	
		numpts = length(exactu);
		for k = 1:length(origbps);		
			jnd = ((k-1)*numpts +1):(k*numpts);
			uoerr = (origu(jnd) - exactu).^2;
			voerr = (origv(jnd) - exactv).^2;
			ueerr = (u(jnd) - exactu).^2;
			veerr = (v(jnd) - exactv).^2;
			l2errfarorig(n,k) =  patchspacing*sqrt( sum( uoerr(1:end/2) + voerr(1:end/2) ) );
			l2errnearorig(n,k) = patchspacing*sqrt( sum( uoerr((end/2+1):end) + voerr((end/2+1):end) ) );
			l2errfarexactint(n,k) =  patchspacing*sqrt( sum( ueerr(1:end/2) + veerr(1:end/2) ) );
			l2errnearexactint(n,k) = patchspacing*sqrt( sum( ueerr((end/2+1):end) + veerr((end/2+1):end) ) );
		end
		
		indfo = find(l2errfarorig(n,:) == min(l2errfarorig(n,:)));
		indno = find(l2errnearorig(n,:) == min(l2errnearorig(n,:)));
		indfe = find(l2errfarexactint(n,:) == min(l2errfarexactint(n,:)));
		indne = find(l2errnearexactint(n,:) == min(l2errnearexactint(n,:)));
		
		bestblobfarorig(end+1)=origbps(indfo(1));
		bestblobnearorig(end+1)=origbps(indno(1));
		bestblobfarexactint(end+1)=eintbps(indfe(1));
		bestblobnearexactint(end+1)=eintbps(indne(1));

		% figure
		% plot(origbps,l2errfarorig,'b.')
		% hold on
		% plot(origbps,l2errnearorig,'r.')
		% title(['Original method, N=',int2str(N)])
		% legend('far field','near field')
		% xlabel('blob parameter')
		% ylabel('patch error')
		% 	
		% figure
		% plot(eintbps,l2errfarexactint,'b.')
		% hold on
		% plot(eintbps,l2errnearexactint,'r.')
		% title(['Exact integral method, N=',int2str(N)])
		% legend('far field','near field')
		% xlabel('blob parameter')
		% ylabel('patch error')
		% 	
end	

N = 2.^(3:8);
chordlen = 0.25*2*sin(pi./N).';

% c={ 'b.', 'r.', 'k.', 'g.' };
% leg = { 'orig, far', 'orig, near', 'exact int, far','exact int, near'};
% str = { 'farorig','nearorig','farexactint','nearexactint' };
% coeffs=[];
% figure
% for k = 1:4;
% 	loglog(chordlen,eval(['bestblob',str{k}]),c{k})
% 	pause(0.5)
% 	hold on
% 	coeffs(end+1,:) = polyfit(log(chordlen),log(eval(['bestblob',str{k}])),1);
% end
% title('Original and exact integral methods')
% xlabel('chord length')
% ylabel('minimal error blob')
% legend(leg,'Location','BestOutside')
% 
% coeffs = [coeffs(:,1),exp(coeffs(:,2))]

%Make table
table = [l2errfarorig, l2errfarorig./chordlen.^(1.1), l2errnearorig, l2errnearorig./chordlen.^(1.1), l2errfarexactint, l2errfarexactint./chordlen.^2, l2errnearexactint, l2errnearexactint./chordlen.^(2.5)]

















	
	
