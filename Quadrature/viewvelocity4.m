%viewvelocity4.m
	
	clear 
	close all
	
	basefname = '~/rsyncfolder/data/Quadrature/stokes2Dcylinder/StokesCylTest_exactintegralsANDorig_2points_singleblob4poster_varyingeps_N';
	singleblob = 1;
	Nvec = 2.^(4:7);
	chordlen = 0.25*2*sin(pi./Nvec);
	
	err_farp_origmeth=[];
	err_near_origmeth=[];
	err_farp_exactint=[];
	err_near_exactint=[];
	
	allorigbps = [];
	alleintbps = [];
	
	n=0;
	for N = Nvec;
		n=n+1;
		str = sprintf('%04d',N);
		load([basefname,str])
		
		allorigbps(n,:) = origbps;
		alleintbps(n,:) = eintbps;
	
		uoerr = (origu - repmat(exactu,1,length(origbps))).^2;
		voerr = (origv - repmat(exactv,1,length(origbps))).^2;
		ueerr = (u - repmat(exactu,1,length(eintbps))).^2;
		veerr = (v - repmat(exactv,1,length(eintbps))).^2;
		err_farp_origmeth(n,:) =  sqrt( uoerr(1:2:end) + voerr(1:2:end)  );
		err_near_origmeth(n,:) =  sqrt( uoerr(2:2:end) + voerr(2:2:end)  );
		err_farp_exactint(n,:) =  sqrt( ueerr(1:2:end) + veerr(1:2:end)  );
		err_near_exactint(n,:) =  sqrt( ueerr(2:2:end) + veerr(2:2:end)  );
	end	

if singleblob == 0;
	% figure
	% plot(allorigbps.',err_near_origmeth.')
	% title('Trap, near')
	% 
	% figure
	% plot(allorigbps.',err_farp_origmeth.')
	% title('Trap, far')

	figure
	plot(alleintbps.',err_near_exactint.')
	title('BEM, near')
	
	figure
	plot(alleintbps.',err_farp_exactint.')
	title('BEM, far')

	% me_farporig = min(err_farp_origmeth,[],2);
	% me_nearorig = min(err_near_origmeth,[],2);
	me_farpeint = min(err_farp_exactint,[],2);
	me_neareint = min(err_near_exactint,[],2);
	
	% figure
	% loglog(chordlen,me_farporig,'b.')
	% xlabel('Chord length (=0.25*2*sin(\pi /N))')
	% ylabel('minimal error over all blobs tested')
	% title('Far, orig')
	% p_fo = polyfit(log(chordlen.'),log(me_farporig),1);
	% hold on
	% loglog(chordlen,exp(p_fo(2))*chordlen.^(p_fo(1)),'r')
	% hold off
	% legend('minimal error',['slope ',num2str(p_fo(1))])
	% 
	% 
	% figure
	% loglog(chordlen,me_nearorig,'b.')
	% xlabel('Chord length (=0.25*2*sin(\pi /N))')
	% ylabel('minimal error over all blobs tested')
	% title('Near, orig')
	% p_no = polyfit(log(chordlen.'),log(me_nearorig),1);
	% hold on
	% loglog(chordlen,exp(p_no(2))*chordlen.^(p_no(1)),'r')
	% hold off
	% legend('minimal error',['slope ',num2str(p_no(1))])

	figure
	loglog(chordlen,me_farpeint,'b.')
	xlabel('Chord length (=0.25*2*sin(\pi /N))')
	ylabel('minimal error over all blobs tested')
	title('Far, exact int')
	p_fe = polyfit(log(chordlen.'),log(me_farpeint),1);
	hold on
	loglog(chordlen,exp(p_fe(2))*chordlen.^(p_fe(1)),'r')
	hold off
	legend('minimal error',['slope ',num2str(p_fe(1))])
	
	
	figure
	loglog(chordlen,me_neareint,'b.')
	xlabel('Chord length (=0.25*2*sin(\pi /N))')
	ylabel('minimal error over all blobs tested')
	title('Near, exact int')
	p_ne = polyfit(log(chordlen.'),log(me_neareint),1);
	hold on
	loglog(chordlen,exp(p_ne(2))*chordlen.^(p_ne(1)),'r')
	hold off
	legend('minimal error',['slope ',num2str(p_ne(1))])
	
	% bb_farporig = [];
	% bb_nearorig = [];
	bb_farpeint = [];
	bb_neareint = [];

	for k = 1:length(Nvec);
		% ind = find(err_farp_origmeth(k,:) == me_farporig(k));
		% bb_farporig(k) = allorigbps(k,ind);
		% ind = find(err_near_origmeth(k,:) == me_nearorig(k));
		% bb_nearorig(k) = allorigbps(k,ind);
		ind = find(err_farp_exactint(k,:) == me_farpeint(k));
		bb_farpeint(k) = alleintbps(k,ind);
		ind = find(err_near_exactint(k,:) == me_neareint(k));
		bb_neareint(k) = alleintbps(k,ind);
	end

	% figure
	% loglog(chordlen,bb_farporig,'b.')
	% xlabel('Chord length (=0.25*2*sin(\pi /N))')
	% ylabel('best blob')
	% title('Far, orig')
	% p_fo = polyfit(log(chordlen),log(bb_farporig),1);
	% hold on
	% loglog(chordlen,exp(p_fo(2))*chordlen.^(p_fo(1)),'r')
	% hold off
	% legend('best blobs',['slope ',num2str(p_fo(1))])
	% 
	% figure
	% loglog(chordlen,bb_nearorig,'b.')
	% xlabel('Chord length (=0.25*2*sin(\pi /N))')
	% ylabel('best blob')
	% title('Near, orig')
	% p_no = polyfit(log(chordlen),log(bb_nearorig),1);
	% hold on
	% loglog(chordlen,exp(p_no(2))*chordlen.^(p_no(1)),'r')
	% hold off
	% legend('best blobs',['slope ',num2str(p_no(1))])
	
	figure
	loglog(chordlen,bb_farpeint,'b.')
	xlabel('Chord length (=0.25*2*sin(\pi /N))')
	ylabel('best blob')
	title('Far, exact int')
	p_fe = polyfit(log(chordlen),log(bb_farpeint),1);
	hold on
	loglog(chordlen,exp(p_fe(2))*chordlen.^(p_fe(1)),'r')
	hold off
	legend('best blobs',['slope ',num2str(p_fe(1))])
	
	figure
	loglog(chordlen,bb_neareint,'b.')
	xlabel('Chord length (=0.25*2*sin(\pi /N))')
	ylabel('best blob')
	title('Near, exact int')
	p_ne = polyfit(log(chordlen),log(bb_neareint),1);
	hold on
	loglog(chordlen,exp(p_ne(2))*chordlen.^(p_ne(1)),'r')
	hold off
	legend('best blobs',['slope ',num2str(p_ne(1))])
	
	

else;
	
	chordlen=chordlen.';
	%Make table
	table = [err_farp_origmeth, err_farp_origmeth./chordlen.^(1), err_near_origmeth, err_near_origmeth./chordlen.^(1), err_farp_exactint, err_farp_exactint./chordlen.^(2), err_near_exactint, err_near_exactint./chordlen.^(2)]
end

















	
	
