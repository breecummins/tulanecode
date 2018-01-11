%viewvelocity.m
	
	clear 
	close all
	
	
	patchspacing = 0.02/5;
	
	bestblobf1o=[];
	bestblobn1o=[];
	bestblobf2o=[];
	bestblobn2o=[];
	bestblobf1b=[];
	bestblobn1b=[];
	bestblobf2b=[];
	bestblobn2b=[];
	
	for N = 2.^(3:8);
		str = sprintf('%04d',N);
		load(['~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_origmethod_bothblobs',str])
		u1o=u_blob1;
		v1o=v_blob1;
		u2o=u_blob2;
		v2o=v_blob2;
		bps_orig = blobparams;
		load(['~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_hattrapBEM_bothblobs',str])
		u1b=u_blob1;
		v1b=v_blob1;
		u2b=u_blob2;
		v2b=v_blob2;
		bps_BEM = blobparams;
	
		numpts = length(exactu);
		l2errfar1o = [];
		l2errnear1o = [];
		l2errfar2o = [];
		l2errnear2o = [];
		l2errfar1b = [];
		l2errnear1b = [];
		l2errfar2b = [];
		l2errnear2b = [];
		for k = 1:length(bps_orig);		
			jnd = ((k-1)*numpts +1):(k*numpts);
			u1oerr = (u1o(jnd) - exactu).^2;
			v1oerr = (v1o(jnd) - exactv).^2;
			u2oerr = (u2o(jnd) - exactu).^2;
			v2oerr = (v2o(jnd) - exactv).^2;
			l2errfar1o(end+1) =  patchspacing*sqrt( sum( u1oerr(1:end/2) + v1oerr(1:end/2) ) );
			l2errnear1o(end+1) = patchspacing*sqrt( sum( u1oerr((end/2+1):end) + v1oerr((end/2+1):end) ) );
			l2errfar2o(end+1) =  patchspacing*sqrt( sum( u2oerr(1:end/2) + v2oerr(1:end/2) ) );
			l2errnear2o(end+1) = patchspacing*sqrt( sum( u2oerr((end/2+1):end) + v2oerr((end/2+1):end) ) );
		end
		for k = 1:length(bps_BEM);		
			u1berr = (u1b(jnd) - exactu).^2;
			v1berr = (v1b(jnd) - exactv).^2;
			u2berr = (u2b(jnd) - exactu).^2;
			v2berr = (v2b(jnd) - exactv).^2;
			l2errfar1b(end+1) =  patchspacing*sqrt( sum( u1berr(1:end/2) + v1berr(1:end/2) ) );
			l2errnear1b(end+1) = patchspacing*sqrt( sum( u1berr((end/2+1):end) + v1berr((end/2+1):end) ) );
			l2errfar2b(end+1) =  patchspacing*sqrt( sum( u2berr(1:end/2) + v2berr(1:end/2) ) );
			l2errnear2b(end+1) = patchspacing*sqrt( sum( u2berr((end/2+1):end) + v2berr((end/2+1):end) ) );
		end
		
		indf1o = find(l2errfar1o == min(l2errfar1o));
		indn1o = find(l2errnear1o == min(l2errnear1o));
		indf2o = find(l2errfar2o == min(l2errfar2o));
		indn2o = find(l2errnear2o == min(l2errnear2o));
		indf1b = find(l2errfar1b == min(l2errfar1b));
		indn1b = find(l2errnear1b == min(l2errnear1b));
		indf2b = find(l2errfar2b == min(l2errfar2b));
		indn2b = find(l2errnear2b == min(l2errnear2b));
		
		bestblobf1o(end+1)=bps_orig(indf1o(1));
		bestblobn1o(end+1)=bps_orig(indn1o(1));
		bestblobf2o(end+1)=bps_orig(indf2o(1));
		bestblobn2o(end+1)=bps_orig(indn2o(1));
		bestblobf1b(end+1)=bps_BEM(indf1b(1));
		bestblobn1b(end+1)=bps_BEM(indn1b(1));
		bestblobf2b(end+1)=bps_BEM(indf2b(1));
		bestblobn2b(end+1)=bps_BEM(indn2b(1));

		figure
		plot(bps_orig,l2errfar1o,'b.')
		hold on
		plot(bps_orig,l2errnear1o,'r.')
		title('Original method, blob 1')
		legend('far field','near field')
		xlabel('blob parameter')
		ylabel('patch error')
			
		figure
		plot(bps_orig,l2errfar2o,'b.')
		hold on
		plot(bps_orig,l2errnear2o,'r.')
		title('Original method, blob 2')
		legend('far field','near field')
		xlabel('blob parameter')
		ylabel('patch error')
			
		figure
		plot(bps_BEM,l2errfar1b,'b.')
		hold on
		plot(bps_BEM,l2errnear1b,'r.')
		title('Hat/trap BEM, blob 1')
		legend('far field','near field')
		xlabel('blob parameter')
		ylabel('patch error')
			
		figure
		plot(bps_BEM,l2errfar2b,'b.')
		hold on
		plot(bps_BEM,l2errnear2b,'r.')
		title('Hat/trap BEM, blob 2')
		legend('far field','near field')
		xlabel('blob parameter')
		ylabel('patch error')
end	

N = 2.^(3:8);
chordlen = 0.25*2*sin(pi./N);
M = 50;

c={ 'b.', 'r.', 'k.', 'g.' };
leg = { 'blob 1, far', 'blob 1, near,' 'blob 2, far', 'blob 2, near' };
str = { 'f1o','n1o','f2o','n2o' };
origcoeffs=[];
figure
for k = 1:4;
	loglog(chordlen,eval(['bestblob',str{k}]),c{k})
	hold on
	origcoeffs(end+1,:) = polyfit(log(chordlen),log(eval(['bestblob',str{k}])),1);
end
title('Original method')
xlabel('chord length')
ylabel('minimal error blob')
legend(leg,'Location','BestOutside')

str = { 'f1b','n1b','f2b','n2b' };
BEMcoeffs=[];
figure
for k = 1:4;
	loglog(chordlen,eval(['bestblob',str{k}]),c{k})
	hold on
	BEMcoeffs(end+1,:) = polyfit(log(chordlen),log(eval(['bestblob',str{k}])),1);
end
title('Hat/trap BEM')
xlabel('chord length')
ylabel('minimal error blob')
legend(leg,'Location','BestOutside')

origcoeffs = [origcoeffs(:,1),exp(origcoeffs(:,2))]
BEMcoeffs = [BEMcoeffs(:,1),exp(BEMcoeffs(:,2))]





















	
	
