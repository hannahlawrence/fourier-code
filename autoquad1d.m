function wvec=autoquad1d(klocs,delta,test,p)

% The ith weight w_i = 1 / sum sinc^2(delta * (klocs_i-klocs_j))
%                           j
% where sinc(x)=sin(pi * x) / (pi * x)
% If test=1, will test result for correctness with allowable (L2) error up to p
% If test=0, p is not used

ifl=1;
wvec=1./sincsq1d(ifl,delta*klocs,ones(size(klocs)),1e-16);

if test==1
    test1d(klocs,delta,wvec,p);
end

function test1d(klocs,delta,wvec,err_thresh)
% Check taken from sincsq1d

[a,b]=ndgrid(delta*klocs,delta*klocs);
sincmat=(sin(pi*(a-b))./(pi*(a-b))).^2;
sincmat(arrayfun(@isnan,sincmat))=1;
correct_wtrans=sum(sincmat,2);
correct_wtrans=1./correct_wtrans;
err=sqrt(sum((correct_wtrans-wvec).*(correct_wtrans-wvec)));
if err>err_thresh
    fprintf("Error: autoquad2d calculation failed with error %g instead of requested %g\n",err,err_thresh);
end
