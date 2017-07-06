function wvec=autoquad3d(klocs_d1,klocs_d2,klocs_d3,delta_d1,delta_d2,delta_d3,test,p)

% The ith weight w_i = 1 / sum sinc^2(delta_d1 * (klocs_d1_i-klocs_d1_j)) *
%                           j
%                              sinc^2(delta_d2 * (klocs_d2_i-klocs_d2_j)) *
%
%                              sinc^2(delta_d3 * (klocs_d3_i-klocs_d3_j))
%
% where sinc(x)=sin(pi * x) / (pi * x)
% If test=1, will test result for correctness with allowable (L2) error up to p
% If test=0, p is not used

ifl=1;
wvec=1./sincsq3d(ifl,klocs_d1*delta_d1,klocs_d2*delta_d2,klocs_d3*delta_d3,ones(size(klocs_d1)),1e-16);

if test==1
    test3d(klocs_d1,klocs_d2,klocs_d3,delta_d1,delta_d2,delta_d3,wvec,p);
end

function test3d(klocs_d1,klocs_d2,klocs_d3,delta_d1,delta_d2,delta_d3,wvec,err_thresh)
% Check taken from sincsq3d

[a1,b1]=ndgrid(klocs_d1*delta_d1,klocs_d1*delta_d1);
[a2,b2]=ndgrid(klocs_d2*delta_d2,klocs_d2*delta_d2);
[a3,b3]=ndgrid(klocs_d3*delta_d3,klocs_d3*delta_d3);
x=sin(pi*(a1-b1))./(pi*(a1-b1));
y=sin(pi*(a2-b2))./(pi*(a2-b2));
z=sin(pi*(a3-b3))./(pi*(a3-b3));
x(arrayfun(@isnan,x))=1;
y(arrayfun(@isnan,y))=1;
z(arrayfun(@isnan,z))=1;
sincmat=x.*y.*z;
sincmat=sincmat.^2;
correct_fast=sum(sincmat,2); % column vector
correct_fast=1./correct_fast;

err=sqrt(sum((correct_fast-wvec).*(correct_fast-wvec)));
if err>err_thresh
    fprintf("Error: autoquad3d calculation failed with error %g instead of requested %g\n",err,err_thresh);
end