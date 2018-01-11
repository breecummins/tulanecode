function A = MatrixStokeslet(x,y,d)

%
% x = evaluation points
% y = location of the Stokeslets
% d = delta
% 
%
d2 = d^2;

nx = length(x(:,1));
ny = length(y(:,1));

A = zeros(3*nx,3*ny);

for i = 1 : ny

      columnid = 3*(i-1);
      columnid1 = columnid+1;
      columnid2 = columnid+2;
      columnid3 = columnid+3;

      dx1 = x(:,1)-y(i,1);
      dx2 = x(:,2)-y(i,2);
      dx3 = x(:,3)-y(i,3);
      r2  = dx1.^2 + dx2.^2 + dx3.^2;
      R  = sqrt(r2+d2);
      % BEALE BLOB
%       H1 = 1/2./R + 3/2*d2*d2./R.^5;
%       H2 = 1/2./R.^3 + 3/2*d2./R.^5;
      % REGULAR BLOB
      H1 = 1/2./R + 1/2*d2./R.^3;
      H2 = 1/2./R.^3;

      
      A(1:3:end,columnid1) = H1 + dx1.*dx1.*H2;
      A(1:3:end,columnid2) = dx1.*dx2.*H2 ;
      A(1:3:end,columnid3) = dx1.*dx3.*H2 ;

      A(2:3:end,columnid1) = dx2.*dx1.*H2 ;
      A(2:3:end,columnid2) = H1 + dx2.*dx2.*H2 ;
      A(2:3:end,columnid3) = dx2.*dx3.*H2 ;

      A(3:3:end,columnid1) = dx3.*dx1.*H2 ;
      A(3:3:end,columnid2) = dx3.*dx2.*H2 ;
      A(3:3:end,columnid3) = H1 + dx3.*dx3.*H2 ;
end
A = A/(4*pi);
