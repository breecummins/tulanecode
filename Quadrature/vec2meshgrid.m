function varargout=vec2meshgrid(A);
	
	%Take an N by M array A and split it up into M arrays that are each sqrt(N) square matrices. Note that this will fail if N ~= n^2 for all n.
	% outputs will probably look like either X,Y,U,V or X,Y,Z,U,V,W.
	
	if size(A,1) == 1;
		A = A.';
	end
	
	n = sqrt(size(A,1));
	M = size(A,2);
	
	for k = 1:M;
		temp = [];
		for j = 1:n:n^2;
			temp(:,end+1) = A(j:j+n-1,k);
		end
		varargout(k) = {temp};
	end