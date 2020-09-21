function [x_zero,fval,exitflag,value,iter] = mybrent(a,b,crit,maxiter,tau_s,Parameters,Grids,flags,Vguess)
% -------------------------------------------------------------------------
% AIM: This function is a rootfinding algorithm based on Brent's method. It
% finds x_zero such that F(x_zero)=0
%
% INPUTS:    "fun"    is a function handle that returns F(x)
%            "a,b"    are the endpoints of the interval
%
% OUTPUT:    "x_zero" is the zero of the function "fun"
%            "fval"   is the function evaluated at x_zero
%
% -------------------------------------------------------------------------
% Copyright ? 2019 by Alessandro Di Nola. All rights
% reserved.
% -------------------------------------------------------------------------

%% Brent

savea = a;
saveb = b;

%% Lower endpoint
if isempty(Vguess)
    [fa,~,~,~,~,value,~,~] = PE_interpolation1(a,tau_s,Parameters,Grids,flags,[]);
else
    [fa,~,~,~,~,value,~,~] = PE_interpolation1(a,tau_s,Parameters,Grids,flags,Vguess);
end
Vguess = value.V1; %save Vguess at first iteration
%fa=fun(a);

%% Higher endpoint
[fb,~,~,~,~,value,~,~] = PE_interpolation1(b,tau_s,Parameters,Grids,flags,Vguess);
Vguess = value.V1;
%fb=fun(b);

if fa*fb>0
    error('f(a) und f(b) have the same sign!');
end

c=a; fc=fa;   %At the beginning, c = a

c=a; fc=fa; d=b-a; e=d;

%% Main block
iter=0;
%maxiter=1000;

while iter<=maxiter
    iter=iter+1;
    
    if fb*fc>0
        c=a; fc=fa; d=b-a; e=d;
    end
    
    if abs(fc)<abs(fb)
        a=b; b=c; c=a;
        fa=fb; fb=fc; fc=fa;
    end
    
    tol=2*eps*abs(b)+crit; m=(c-b)/2; %Tolerance
    
    if (abs(m)>tol) && (abs(fb)>0) %Proceed with iterations
        
        if (abs(e)<tol) || (abs(fa)<=abs(fb))
            d=m; e=m;
        else
            s=fb/fa;
            if a==c
                p=2*m*s; q=1-s;
            else
                q=fa/fc; r=fb/fc;
                p=s*(2*m*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            end
            if p>0
                q=-q;
            else
                p=-p;
            end
            s=e; e=d;
            if ( 2*p<3*m*q-abs(tol*q) ) && (p<abs(s*q/2))
                d=p/q;
            else
                d=m; e=m;
            end
        end
        
        a=b; fa=fb;
        
        if abs(d)>tol
            b=b+d;
        else
            if m>0
                b=b+tol;
            else
                b=b-tol;
            end
        end
    else
        break;
    end
    %New function evalutation
    
    [fb,~,~,~,~,value,~,~] = PE_interpolation1(b,tau_s,Parameters,Grids,flags,Vguess);
    %fb=fun(b);
    Vguess = value.V1; %update Vguess
    
end

exitflag = 1;


%% Check convergence and update flags
if iter>=maxiter && (abs(m)>tol) && (abs(fb)>0)
    warning('Brent reached maximum number of iterations!')
    exitflag = 0;
else
    disp(['Brent stopped after ',num2str(iter),' iterations'])
    disp('Zero Found In Interval')
    disp([savea, saveb])
end

x_zero=b;
fval = fb;


end %END FUNCTION

