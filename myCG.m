function [finalx,rvec]=myCG(A,b,x_init,numiter,tol)

if(nargin<1), test_myCG; return; end

% b should be a column vector
if size(b,2)~=1
    b=b.';
end
numiter=min(numiter, size(b,1)); %not sure if this is a good idea or not...

count=0;
if isa(A,'function_handle')
    FH=1;
else
    FH=0;
end

x=x_init; %initial guess
if FH
    r=b-A(x);
else
    r=b-A*x;
end
p=r;
rsold=r.'*r;
rvec=rsold;
besttol=sqrt(sum(r.^2));
pbest=x_init;
while sqrt(sum(r.^2))>tol && count<numiter
    if FH
        Ap=A(p);
    else
        Ap=A*p;
    end
    alpha=rsold/(p.'*Ap);
    x=x+alpha*p;
    r=r-alpha*Ap;
    rsnew=r.'*r;
    p=r+(rsnew/rsold)*p;
    rsold=rsnew;
    rvec=vertcat(rvec,rsold);
    if sqrt(sum(r.^2))<besttol
        besttol=sqrt(sum(r.^2));
        pbest=x;
    end
    count=count+1;
end
if count<numiter
    finalx=x;
else %exited b/c hit max numiter
    finalx=pbest;
end

end

function test_myCG
n=10;
tol=1e-5;
maxit=20;
A=rand(n); % A needs to be symmetric!
A=.5*(A+A'); % All entries of A <= 1
A=A+n*eye(n); % A is diagonally dominant
% A = diagonally dominant with real + diagonal entries! = positive definite
b=10*rand(n,1);
[x,rvec]=myCG(@sampleA,b,zeros(size(b)),maxit,tol)
display(sqrt(sum((sampleA(x)-b).^2)));
matlabx=pcg(@sampleA,b,tol,maxit)
display(sqrt(sum((A*matlabx-b).^2))); %note: still returns same answer if #iter same
end

function vec=sampleA(x)
A=magic(length(x))/100;
A=.5*(A+A');
m=max(max(A));
A=A+m*eye(length(x)); % A is positive definite and symmetric
vec=A*x;
end
