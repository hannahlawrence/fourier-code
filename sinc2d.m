function wtrans=sinc2d(ifl,klocs_d1,klocs_d2,q,rsamp)

if(nargin<1), test_sinc2d; return; end

%  
% wtrans(j) = sum sinc(klocs_d1(k)-klocs_d1(j)) * sinc(klocs_d2(k)-klocs_d2(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs_d1 = (real) sample locations in dimension 1
% klocs_d2 = (real) sample locations in dimension 2
% q = sample strengths

rkmaxx=max(bsxfun(@max,zeros(size(klocs_d1)),abs(klocs_d1)));
rkmaxy=max(bsxfun(@max,zeros(size(klocs_d2)),abs(klocs_d2)));

if ifl==1
    rkmaxx=pi*rkmaxx;
    rkmaxy=pi*rkmaxy;
    klocs_d1=pi*klocs_d1;
    klocs_d2=pi*klocs_d2;
end

nx=ceil(rsamp*round(rkmaxx+3));
ny=ceil(rsamp*round(rkmaxy+3));

[xx,wwx]=lgwt(nx,-1,1);
[yy,wwy]=lgwt(ny,-1,1);

% create 2D grid of points and corresponding weights
[c,d]=ndgrid(xx,yy);
allxx=c(:);
allyy=d(:);
[e,f]=ndgrid(wwx,wwy);
allww=e(:).*f(:);

h_at_xxyy=finufft2d3(klocs_d1,klocs_d2,q,-1,1e-15,allxx,allyy);

wtrans=0.25*finufft2d3(allxx,allyy,h_at_xxyy.*allww,1,1e-15,klocs_d1,klocs_d2);
wtrans=real(wtrans);

function test_sinc2d
n=100;
numtrials=5;
for t=1:numtrials
    klocs_d1=-pi+(2*pi*rand(n,1));
    klocs_d2=-pi+(2*pi*rand(n,1));
    q=rand(1,n);
    ifl=0;
    correct=slowsinc2d(ifl,klocs_d1,klocs_d2,q);
    myresult=sinc2d(ifl,klocs_d1,klocs_d2,q,2);
    fprintf("Error: %g\n",mean(abs(correct-myresult)));
end

function correct = slowsinc2d(ifl,klocs_d1,klocs_d2,q)
    [a1,b1]=ndgrid(klocs_d1,klocs_d1);
    [a2,b2]=ndgrid(klocs_d2,klocs_d2);
    if ifl==1
        x=sin(pi*(a1-b1))./(pi*(a1-b1));
        y=sin(pi*(a2-b2))./(pi*(a2-b2));
    else
        x=sin(a1-b1)./(a1-b1);
        y=sin(a2-b2)./(a2-b2);
    end
    x(arrayfun(@isnan,x))=1;
    y(arrayfun(@isnan,y))=1;
    sincmat=x.*y;
    correct=sum(repmat(q,length(klocs_d1),1).*sincmat,2);
