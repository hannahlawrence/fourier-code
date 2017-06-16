function wtrans=sinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q,rsamp)

if(nargin<1), test_sinc3d; return; end

%  
% wtrans(j) = sum sinc(klocs_d1(k)-klocs_d1(j)) * sinc(klocs_d2(k)-klocs_d2(j)) * sinc(klocs_d3(k)-klocs_d3(j)) * q(j)
%              k
%
% ifl = sinc convention
%   0: sinc(x) = sin(x)/x
%   1: sinc(x)=sin(pi*x)/(pi*x)
% klocs_d1 = (real) sample locations in dimension 1
% klocs_d2 = (real) sample locations in dimension 2
% klocs_d3 = (real) sample locations in dimension 3
% q = sample strengths

rkmaxx=max(bsxfun(@max,zeros(size(klocs_d1)),abs(klocs_d1)));
rkmaxy=max(bsxfun(@max,zeros(size(klocs_d2)),abs(klocs_d2)));
rkmaxz=max(bsxfun(@max,zeros(size(klocs_d3)),abs(klocs_d3)));

if ifl==1
    rkmaxx=pi*rkmaxx;
    rkmaxy=pi*rkmaxy;
    rkmaxz=pi*rkmaxz;
    klocs_d1=pi*klocs_d1;
    klocs_d2=pi*klocs_d2;
    klocs_d3=pi*klocs_d3;
end
nx=ceil(rsamp*round(rkmaxx+3)); 
ny=ceil(rsamp*round(rkmaxy+3));
nz=ceil(rsamp*round(rkmaxz+3));

[xx,wwx]=lgwt(nx,-1,1);
[yy,wwy]=lgwt(ny,-1,1);
[zz,wwz]=lgwt(nz,-1,1);

[c,d,e]=ndgrid(xx,yy,zz);
allxx=c(:);
allyy=d(:);
allzz=e(:);
[f,g,h]=ndgrid(wwx,wwy,wwz);
allww=f(:).*g(:).*h(:);

h_at_xxyyzz=finufft3d3(klocs_d1,klocs_d2,klocs_d3,q,-1,1e-15,allxx,allyy,allzz);

wtrans=(1/8)*finufft3d3(allxx,allyy,allzz,h_at_xxyyzz.*allww,1,1e-15,klocs_d1,klocs_d2,klocs_d3);
wtrans=real(wtrans);

function test_sinc3d
n=10;
numtrials=5;
for t=1:numtrials
    klocs_d1=-pi+(2*pi*rand(n,1));
    klocs_d2=-pi+(2*pi*rand(n,1));
    klocs_d3=-pi+(2*pi*rand(n,1));
    q=rand(1,n)*3; 
    ifl=0;
    correct=slowsinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q);
    myresult=sinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q,2);
    fprintf("Error: %g\n",mean(abs(correct-myresult)));
end

function correct=slowsinc3d(ifl,klocs_d1,klocs_d2,klocs_d3,q)
    [a1,b1]=ndgrid(klocs_d1,klocs_d1);
    [a2,b2]=ndgrid(klocs_d2,klocs_d2);
    [a3,b3]=ndgrid(klocs_d3,klocs_d3);
    if ifl==1
        x=sin(pi*(a1-b1))./(pi*(a1-b1));
        y=sin(pi*(a2-b2))./(pi*(a2-b2));
        z=sin(pi*(a3-b3))./(pi*(a3-b3));
    else
        x=sin(a1-b1)./(a1-b1);
        y=sin(a2-b2)./(a2-b2);
        z=sin(a3-b3)./(a3-b3);
    end
    x(arrayfun(@isnan,x))=1;
    y(arrayfun(@isnan,y))=1;
    z(arrayfun(@isnan,z))=1;
    sincmat=x.*y.*z;
    correct=sum(repmat(q,length(klocs_d1),1).*sincmat,2);
