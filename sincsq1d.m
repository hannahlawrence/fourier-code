function wtrans= sincsq1d(ifl,klocs,q,rsamp)

if(nargin<1), test_sincsq1d; return; end

%  
% wtrans(j) = sum sinc^2(klocs(k)-klocs(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs = (real) sample locations
% q = sample strengths

rkmax=max(bsxfun(@max,zeros(size(klocs)),abs(klocs)));

if ifl==1
    rkmax=pi*rkmax;
    klocs=pi*klocs;
end

nx=rsamp*round(rkmax+3);
[xx_pre,ww_pre]=lgwt(nx,-1,1);

% quadrature points and weights to cover [-2,2]
xx=vertcat((xx_pre-1),(xx_pre+1));
ww=vertcat((ww_pre-1),(ww_pre+1));

h_at_xx=finufft1d3(klocs,q,-1,1e-15,xx);

trianglevec=zeros(size(ww));
for i=1:(2*nx)
    val=abs(xx(i));
    trianglevec(i)=2-val;
end

wtrans=.25*finufft1d3(xx,h_at_xx.*ww.*trianglevec,1,1e-15,klocs);
wtrans=real(wtrans);

function test_sincsq1d
numtrials=5;
n=10;
for numt=1:numtrials
    klocs=rand(n,1);
    q=rand(1,n); 
    ifl=1;

    [a,b]=ndgrid(klocs,klocs);
    repmat(q,length(q),1);
    if ifl==1
        sincmat=(sin(pi*(a-b))./(pi*(a-b))).^2;
    else
        sincmat=(sin(a-b)./(a-b)).^2;
    end
    sincmat(arrayfun(@isnan,sincmat))=1;
    correct_wtrans=sum(repmat(q,length(q),1).*sincmat,2);
    my_wtrans= sincsq1d(ifl,klocs,q,2);
    fprintf("Error: %g\n", mean(abs(my_wtrans-correct_wtrans)));
end




