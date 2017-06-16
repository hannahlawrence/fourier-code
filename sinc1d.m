function wtrans= sinc1d(ifl,klocs,q,rsamp)

if(nargin<1), test_sinc1d; return; end

%  
% wtrans(j) = sum sinc(klocs(k)-klocs(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs = (real) sample locations
% q = sample strengths
% rsamp = oversampling parameter in discretization of sinc quadrature

% ensure vector inputs are column vectors
q=q(:);
klocs=klocs(:);

rkmax=max(bsxfun(@max,zeros(size(klocs)),abs(klocs)));

if ifl==1
    rkmax=pi*rkmax;
    klocs=pi*klocs;
end

nx=ceil(rsamp*round(rkmax+3)); % following fortran code

[xx,ww]=lgwt(nx,-1,1);
h_at_xx=finufft1d3(klocs,q,-1,1e-15,xx);

% factor of 1/2 required for change of integration bounds
wtrans=0.5*finufft1d3(xx,h_at_xx.*ww,1,1e-15,klocs);

% only return real values
wtrans=real(wtrans);

function test_sinc1d

numtrials=20;
n=100;
resultmat=zeros(n,numtrials);

for t=1:numtrials
    klocs=(-pi)+(pi*2*rand(n,1));
    q=.1+rand(1,n)*5;
    ifl=1;
    my_wtrans= sinc1d(ifl,klocs,q,2);
    correct_wtrans=slowsinc1d(ifl,klocs,q);
    if(size(correct_wtrans,1)~=size(my_wtrans,1))
        correct_wtrans=correct_wtrans.';
    end
    resultmat(:,t)=correct_wtrans-my_wtrans;
    fprintf("Error: %g\n", mean(abs(my_wtrans-correct_wtrans)));
end

function correct_wtrans=slowsinc1d(ifl,klocs,q) 
    [a,b]=ndgrid(klocs,klocs);
    if ifl==1
        sincmat=sin(pi*(a-b))./(pi*(a-b));
    else
        sincmat=sin(a-b)./(a-b);
    end
    sincmat(arrayfun(@isnan,sincmat))=1;
    correct_wtrans=sum(repmat(q,length(q),1).*sincmat,2);
