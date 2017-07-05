function wvec=autoquad1d(klocs,delta)

% The ith weight w_i = 1 / sum sinc^2(delta * (klocs_i-klocs_j))
%                           j
% where sinc(x)=sin(pi * x) / (pi * x)

ifl=1;
wvec=1./sincsq1d(ifl,delta*klocs,ones(size(klocs)),1e-16);

% Check correctness (can comment out next line for speed)
test1d(klocs,delta,wvec,1e-7);

function test1d(klocs,delta,wvec,err_thresh)
% Check taken from sincsq1d

[a,b]=ndgrid(delta*klocs,delta*klocs);
sincmat=(sin(pi*(a-b))./(pi*(a-b))).^2;
sincmat(arrayfun(@isnan,sincmat))=1;
correct_wtrans=sum(sincmat,2);
correct_wtrans=1./correct_wtrans;
if max(abs(correct_wtrans-wvec))>err_thresh
    fprintf("Error: autoquad1d calculation failed\n");
end
