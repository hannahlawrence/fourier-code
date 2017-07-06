function [image] = recon2d(klocs_d1,klocs_d2,kdata,unifpts_d1,unifpts_d2,method)

if(nargin<1), test_recon2d; return; end

    function Aout = FtranspF(ps)
        fwd=finufft2d3(2*pi*unifpts_d1,2*pi*unifpts_d2,ps,-1,1e-15,klocs_d1,klocs_d2); % should be type 2, possibly
        if length(find(isnan(fwd)))>0
            fprintf("fwd is nan on input:\n");
            display(ps);
        end
        Aout=finufft2d3(2*pi*klocs_d1,2*pi*klocs_d2,fwd,1,1e-15,unifpts_d1,unifpts_d2);
        if length(find(isnan(fwd)))>0
            fprintf("Aout is nan on input:\n");
            display(ps);
        end
    end

    function out=F(ps)
        out=finufft1d3(2*pi*unifpts,ps,-1,1e-15,klocs);
    end
    function out=preconR(y_bar)
        t1=sqwmat*y_bar; %might have to check dimensions: make sure y a column vector
        if length(find(isnan(t1)))>0
            fprintf("t1 is nan on input:\n");
            display(y_bar);
        end
        t2=finufft2d3(2*pi*klocs_d1,2*pi*klocs_d2,t1,1,1e-15,unifpts_d1,unifpts_d2);
        t3=finufft2d3(2*pi*unifpts_d1,2*pi*unifpts_d1,t2,-1,1e-15,klocs_d1,klocs_d2);
        out=sqwmat*t3;
    end
        
if strcmp(method,'direct')
    unknownfactor=1;
    numpts=floor(sqrt(length(unifpts_d1))); % assume same #x y for testing
    unifspace=(max(unifpts_d1)-min(unifpts_d1))/numpts;
    %unknownfactor=1/(numpts*unifspace);
    weights=autoquad2d(klocs_d1,klocs_d2,1/unknownfactor,1/unknownfactor,1,1e-7); %delta?
    image=(1/length(klocs_d1))*finufft2d3(2*pi*klocs_d1,2*pi*klocs_d2,kdata.*weights,1,1e-15,unifpts_d1,unifpts_d2); %2 pi?
else
    if strcmp(method,'manual') 
        % actually: can't do conjugate-gradient because not square matrix
        % ------------------------------
        % just apply F* to kdata to see
        image=(1/length(kdata))*finufft2d3(2*pi*klocs_d1,2*pi*klocs_d2,kdata,1,1e-15,unifpts_d1,unifpts_d2);
    else
        if strcmp(method,'PCG-left')
            b=finufft2d3(2*pi*klocs_d1,2*pi*klocs_d2,kdata,1,1e-15,unifpts_d1,unifpts_d2); 
            image=myCG(@FtranspF,b,zeros(size(b)),50,1e-5);
        else
            if strcmp(method,'PCG-right')
                unknownfactor=1;
                weights=autoquad2d(klocs_d1,klocs_d2,unknownfactor,unknownfactor,1,1e-7);
                wmat=diag(weights);
                sqwmat=sqrt(wmat);
                b=sqrt(weights).*kdata;
                if b~=(sqwmat*kdata)
                    fprintf("strange\n");
                end
                if size(b,2)~=1
                    fprintf("dim issue\n");
                end
                y_bar=pcg(@preconR,b); 
                y=(1./sqrt(weights)).*y_bar;
                image=finufft2d3(2*pi*klocs_d1,2*pi*klocs_d2,weights.*y,1,1e-15,unifpts_d1,unifpts_d2);
            end
        end
    end
end
end

function test_recon2d

numunif=20;
unif_lb=-2;
unif_ub=2;
numNU=1000;
NU_lb=-3;
NU_ub=3;

unifpts=fillmesh1d(numunif,unif_lb,unif_ub);
allunifx=zeros(numunif^2,1);
allunify=zeros(numunif^2,1);
c=1;
for i=1:numunif
    for j=1:numunif
        allunifx(c)=unifpts(i);
        allunify(c)=unifpts(j);
        c=c+1;
    end
end
unifspace=((max(unifpts)-min(unifpts))/(length(unifpts)-1))
f_true=sin(exp(allunifx))+sin(exp(allunify)); 
klocs_d1=NU_lb+(NU_ub-NU_lb)*rand(numNU,1);
klocs_d2=NU_lb+(NU_ub-NU_lb)*rand(numNU,1);
kdata=finufft2d3(2*pi*allunifx,2*pi*allunify,f_true,-1,1e-15,klocs_d1,klocs_d2); 
recon_direct=recon2d(klocs_d1,klocs_d2,kdata,allunifx,allunify,'direct');
recon_manual=recon2d(klocs_d1,klocs_d2,kdata,allunifx,allunify,'manual');
recon_PCG_left=recon2d(klocs_d1,klocs_d2,kdata,allunifx,allunify,'PCG-left');
recon_PCG_right=recon2d(klocs_d1,klocs_d2,kdata,allunifx,allunify,'PCG-right');
numpts=length(allunifx);
plot(1:numpts,f_true,'k',1:numpts,real(recon_direct),'r',1:numpts,real(recon_manual),'b',1:numpts,real(recon_PCG_left),'g',1:numpts,real(recon_PCG_right),'m');
xlabel('Black: True. Red: Direct. Blue: Manual. Green: PCG-left. Magenta: PCG-right.');
end