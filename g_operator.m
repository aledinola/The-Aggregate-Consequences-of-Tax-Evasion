function [Gop,fail] = g_operator(policy,Parameters,Grids,flags,UC,w0,r0)
%{
INPUTS:
    - Policy functions in Benchmark economy, on finer grid for assets
    - Parameters: Matlab structure with parameters
    - Grids,flags: other Matlab structures
    - UC: function handle for consumption piece of U(C,L) 
    - w0, r0: prices (needed by external function trans_operator.m)
    
OUTPUTS:
    - Gop is defined as G^B(x) in Ale's memo. It is the denominator of the
    CEV formula.
    
EXTERNAL FUNCTIONS:
    - This function calls the external function "prob_audit"
    pk = prob_audit([],[],[],kpol_bench(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
    - "trans.operator.m" which gives out Transition T(a,eps,theta,a',eps',theta'), big 6-D

NOTES:
    %Note: all policy functions are defined on (agrid_dist,negrid,ntgrid). To
    %be pedantic, they should all have the suffix "det"
    %Note: all policy functions are from the benchmark (B) economy
    %Note: If I use transition function T I don't need policies for assets.
    %However, T is a 6-D array so maybe is better to avoid it and use directly
    %the assets policy functions.

    %Reference: Alessandro's welfare memo

%Updated by Alessandro Di Nola on September 18, 2020.

%}

%Unpack 
beta              = Parameters.beta;
pflag             = Parameters.pflag;
pn_1              = Parameters.pn_1;
pn_2              = Parameters.pn_2;
pn_3              = Parameters.pn_3;

occpol    = policy.occpoldet;
cpolse0   = policy.cpolse0det;
cpolse1   = policy.cpolse1det;
cpolw     = policy.cpolwdet;
policycap = policy.policycapdet;

[nagrid_dist,negrid,ntgrid] = size(occpol);
sizeC = [nagrid_dist,negrid,ntgrid];

%Check consumption policy functions
%Consumption choice has to be >0
validateattributes(cpolw, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'positive', 'size', sizeC})
validateattributes(cpolse0, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'positive', 'size', sizeC})
validateattributes(cpolse1, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'positive', 'size', sizeC})

Gop     = zeros(nagrid_dist,negrid,ntgrid);
%Gop = (1/(1-beta))*UC(cpolw);

Gop_new = zeros(nagrid_dist,negrid,ntgrid);

%Compute T
T = trans_operator(Parameters,Grids,policy,flags,w0,r0);

if 0==1
    
    sumrowsT = zeros(nagrid_dist,negrid,ntgrid);
    for t=1:ntgrid % current theta
        for j=1:negrid % current eps
            for i=1:nagrid_dist % current a (assets)
                temp = squeeze(T(i,j,t,:,:,:));
                sumrowsT(i,j,t) = sum(temp(:));
                
            end
        end
    end
end

iter    = 0;
tol     = 1e-4;
dist    = tol+1;
itermax = 1000;
fail    = 0;
verbose = 1;

tic
while dist>tol && iter<=itermax
    
    iter = iter+1;
    if verbose==1; fprintf('Iter: %d,',iter); end
    
    for t=1:ntgrid
        for j=1:negrid
            for i=1:nagrid_dist
                
                %First compute GopNext
                %Tx_to_xp =  squeeze(T(i,j,t,:,:,:));
                %GopNext = sum(sum(sum(Tx_to_xp.*Gop)));
                GopNext = 0;
                for tp=1:ntgrid
                    for jp=1:negrid
                        for ip=1:nagrid_dist
                            GopNext = GopNext + T(i,j,t,ip,jp,tp)*Gop(ip,jp,tp);
                        end
                    end
                end

                %Add up all terms (eq.12 of Ale's memo)
                pk = prob_audit([],[],[],policycap(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
                Gop_new(i,j,t)= (1-occpol(i,j,t))*(UC(cpolw(i,j,t)))+... %workers
                    occpol(i,j,t)*((1-pk)*UC(cpolse0(i,j,t)) +pk*UC(cpolse1(i,j,t))  )+... %entre, 0=nd,1=d
                    beta*GopNext;
                
            end
        end
    end
    
    %Compute distance
    dist = max(abs(Gop_new(:)-Gop(:)));
    if verbose==1; fprintf(' Dist: %f \n',dist); end
    
    %Update
    Gop = Gop_new;
    
end %end while
toc

if iter>=itermax
    warning('G operator did not converge!')
    fail = 1;
else
    fprintf('\n')
    disp('G operator converged succesfully!')
    fprintf('\n')
end

end %END FUNCTION "georgi_operator"

