function wvec=autoquad2d(klocs_d1,klocs_d2,delta_d1,delta_d2)

% The ith weight w_i = 1 / sum sinc^2(delta_d1 * (klocs_d1_i-klocs_d1_j)) *
%                           j
%                              sinc^2(delta_d2 * (klocs_d2_i-klocs_d2_j))
%
% where sinc(x)=sin(pi * x) / (pi * x)

ifl=1;
wvec=1./sincsq2d(ifl,klocs_d1*delta_d1,klocs_d2*delta_d2,ones(size(klocs_d1)),1e-16);

% Check correctness (can comment out next line for speed)
test2d(klocs_d1,klocs_d2,delta_d1,delta_d2,wvec,1e-7);

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
if max(abs(correct_fast-wvec))>err_thresh
    fprintf("Error: autoquad2d calculation failed\n");
end
