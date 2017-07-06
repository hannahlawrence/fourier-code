function [image] = recon1d(klocs,kdata,unifpts,method)

if(nargin<1), test_recon1d; return; end

    function Aout = FtranspF(ps)
        fwd=finufft1d3(2*pi*unifpts,ps,-1,1e-15,klocs);
        if length(find(isnan(fwd)))>0
            fprintf("fwd is nan on input:\n");
            display(ps);
        end
        Aout=finufft1d3(2*pi*klocs,fwd,1,1e-15,unifpts);
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
        t2=finufft1d3(2*pi*klocs,t1,1,1e-15,unifpts);
        t3=finufft1d3(2*pi*unifpts,t2,-1,1e-15,klocs);
        out=sqwmat*t3;
    end
    function out=M(a)
        image=sinc1d(1,klocs,a,1e-15);
        
    end
        
if strcmp(method,'direct')
    %unknownfactor=(2*pi)/()
    unknownfactor=1;
    space=(max(unifpts)-min(unifpts)) / (length(unifpts)-1);
    unknownfactor=1/(space*length(unifpts));
    Fbound=2*max(unifpts);
    weights=1./(Fbound*sincsq1d(1,Fbound*klocs,ones(size(klocs)),1e-15));
    image=finufft1d3(2*pi*klocs,kdata.*weights,1,1e-15,unifpts); %2 pi?
else
    if strcmp(method,'manual')
        % can't do conjugate-gradient because not square matrix!
        % just to see
        image=(1/length(kdata))*finufft1d3(2*pi*klocs,kdata,1,1e-15,unifpts); %hugely too big!!!!!!!
        
    else
        if strcmp(method,'PCG-left')
            b=finufft1d3(2*pi*klocs,kdata,1,1e-15,unifpts);
            image=myCG(@FtranspF,b,zeros(size(b)),20,1e-5);
        else
            if strcmp(method,'PCG-right')
                unknownfactor=1;
                weights=autoquad1d(klocs,unknownfactor,1,1e-7);
                wmat=diag(weights);
                sqwmat=sqrt(wmat);
                b=sqrt(weights).*kdata;
                if b~=(sqwmat*kdata)
                    fprintf("err\n");
                end
                if size(b,2)~=1
                    fprintf("dim issue\n");
                end
                y_bar=myCG(@preconR,b,zeros(size(b)),20,1e-5);
                y=(1./sqrt(weights)).*y_bar;
                image=finufft1d3(2*pi*klocs,weights.*y,1,1e-15,unifpts);
            else
                if strcmp(method,'sinc')
                    a=myCG(@M,kdata,30,1e-5); %could add additional precon using weights??
                    image=finufft1d3(2*pi*klocs,a,1,1e-15,unifpts);
                end
            end
            end
        end
    end
end

function test_recon1d

numunif=200;
unif_lb=-.5; % --> F=1
unif_ub=.5;
numNU=400;
NU_lb=-3;
NU_ub=3;

unifpts=fillmesh1d(numunif,unif_lb,unif_ub);
unifspace=((max(unifpts)-min(unifpts))/(length(unifpts)-1));
f_true=exp(sin(5*unifpts));
klocs=NU_lb+6*rand(numNU,1);
kdata=finufft1d3(2*pi*unifpts,f_true,-1,1e-15,klocs); 
recon_direct=recon1d(klocs,kdata,unifpts,'direct');
recon_manual=recon1d(klocs,kdata,unifpts,'manual'); % is not working for some reason:
recon_PCG_left=recon1d(klocs,kdata,unifpts,'PCG-left');
recon_PCG_right=recon1d(klocs,kdata,unifpts,'PCG-right');
plot(unifpts,f_true,'k',unifpts,real(recon_direct),'r',unifpts,real(recon_manual),'b',unifpts,real(recon_PCG_left),'g',unifpts,real(recon_PCG_right),'m');
xlabel('Black: True. Red: Direct. Blue: Manual. Green: PCG-left. Magenta: PCG-right.');
end