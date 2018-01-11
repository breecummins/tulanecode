function v = calcnorm(M)

%This function calculates the 2-norms of real-valued row vectors stacked in a matrix. The output is a column vector with the same number of rows as the input.

if ~all(all(imag(M) == 0));
	error('Inputs must be real-valued!')
end
v = sqrt( sum(M.^2, 2) );