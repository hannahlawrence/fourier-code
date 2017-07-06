function wvec=autoquad2d(klocs_d1,klocs_d2,delta_d1,delta_d2,test,p)

% The ith weight w_i = 1 / sum sinc^2(delta_d1 * (klocs_d1_i-klocs_d1_j)) *
%                           j
%                              sinc^2(delta_d2 * (klocs_d2_i-klocs_d2_j))
%
% where sinc(x)=sin(pi * x) / (pi * x)
% If test=1, will test result for correctness with allowable (L2) error up to p
% If test=0, p is not used

ifl=1;
wvec=1./sincsq2d(ifl,klocs_d1*delta_d1,klocs_d2*delta_d2,ones(size(klocs_d1))',1e-16);

if test==1
    test2d(klocs_d1,klocs_d2,delta_d1,delta_d2,wvec,p);
end

function test2d(klocs_d1,klocs_d2,delta_d1,delta_d2,wvec,err_thresh)

% Check taken from sincsq2d
[a1,b1]=ndgrid(klocs_d1*delta_d1,klocs_d1*delta_d1);
[a2,b2]=ndgrid(klocs_d2*delta_d2,klocs_d2*delta_d2);
x=sin(pi*(a1-b1))./(pi*(a1-b1));
y=sin(pi*(a2-b2))./(pi*(a2-b2));
x(arrayfun(@isnan,x))=1;
y(arrayfun(@isnan,y))=1;
sincmat=x.*y;
sincmat=sincmat.^2;
correct_fast=sum(sincmat,2); % column vector
correct_fast=1./correct_fast;
err=sqrt(sum((correct_fast-wvec).*(correct_fast-wvec)));
if err>err_thresh
    fprintf("Error: autoquad2d calculation failed with error %g instead of requested %g\n",err,err_thresh);
end
