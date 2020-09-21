function [ED,agg,mu_work,mu_se] = aggregates(kn,r0,w0,tau_s,policy,mu,Parameters,Grids)

%--------------------------- LEGEND --------------------------------------%
%{
%DESCRIPTION:
 Given policy function and distribution, compute some aggregate variables
 and excess demand (for asset market clearing)

%INPUTS:
 kn:          capital/labor implied by Cobb-Douglas conditions
 r0,w0,tau_s: interest rate, wage, tax function parameter
 policy:      matlab structure with policy functions
 mu:          distribution
 Parameters:  structure with parameters
 Grids:       structure with grids

%OUTPUTS
 ED:      excess demand (capital demand minus capital supply)
 agg:     structure with aggregate variables
 mu_work: distribution conditional on occpol=work
 mu_se:   distribution conditional on occpol=se

Notes: Updated by Alessandro Di Nola on September 18, 2020

%}

Parameters.tau_s = tau_s; 

%% Unpack structures

delta         = Parameters.delta;
vi            = Parameters.vi;
gamma         = Parameters.gamma;
pn_1          = Parameters.pn_1 ;
pn_2          = Parameters.pn_2 ;
pn_3          = Parameters.pn_3 ;
s             = Parameters.s;
pflag         = Parameters.pflag;
nagrid_dist   = Parameters.nagrid_dist; 
negrid        = Parameters.negrid;
ntgrid        = Parameters.ntgrid;
agrid_dist    = Grids.agrid_dist;
eps_grid      = Grids.eps_grid;
theta         = Grids.theta;
occpoldet     = policy.occpoldet;
policycapdet  = policy.policycapdet;
lpolwdet      = policy.lpolwdet;
policyndet    = policy.policyndet;
lepoldet      = policy.lepoldet;
policyphidet  = policy.policyphidet;

%% Distribution of assets of workers and self-employed

% MU is defined on the finer grid for assets (nagrid_dist)

%Compute marginal pdf of assets by integrating out (eps,theta), i.e. mu(a) = sum_z(mu(a,z))
%mu_a = sum(sum(mu,2),3); % dimension: nagrid_dist,1

% Given mu, compute the conditional distributions mu_work and mu_entre 
% for workers and self-employed, respectively
mu_se = zeros(nagrid_dist,negrid,ntgrid); 
mu_work = zeros(nagrid_dist,negrid,ntgrid); 

for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            if ( occpoldet(i,j,t)==1 )
                mu_se(i,j,t) = mu(i,j,t);
            else
                mu_work(i,j,t) = mu(i,j,t);
                
            end
        end
    end
end

%Normalize so they sum to one
mu_se = mu_se/sum(mu_se(:));
mu_work = mu_work/sum(mu_work(:));


%% Compute aggregate variables

% net capital supply
% K_C + K_E = aggregate savings
% k_supply = aggregate savings - K_E
% In equilibrium k_supply = K_C

k_supply = 0;
for t=1: ntgrid
    for j=1:negrid
        for i = 1:nagrid_dist
            k_supply = k_supply + agrid_dist(i)*mu(i,j,t);
            %k_supply = k_supply + agrid_dist(i)*mu(i,j,t) - occpoldet(i,j,t)*policycapdet(i,j,t)*mu(i,j,t);
            if occpoldet(i,j,t)==1 % this is redundant, we can multiply by occpoldet next in the next line
                k_supply = k_supply - policycapdet(i,j,t)*mu(i,j,t);
            end
        end
    end
end

% net labor supply 
n_supply = 0;
for t=1: ntgrid 
    for j=1:negrid
        for i = 1:nagrid_dist
            % sum only if the guy is a worker, ie. OCC = 0
            n_supply = n_supply + (1-occpoldet(i,j,t))*eps_grid(j)*lpolwdet(i,j,t)*mu(i,j,t)...
                -occpoldet(i,j,t)*policyndet(i,j,t)*mu(i,j,t);
        end
    end
end

% capital-labor ratio implied by market clearing
kn_implied = k_supply/n_supply;

if kn_implied<0
    disp('kn_implied is negative')
    keyboard
end

% Put excess demand in relative terms as (kn - kn_implied)/kn
%ED = kn - kn_implied; % kn: K/N from firm FOCs
ED = (kn - kn_implied)/kn;

%% Government budget constraint
% Here we do the following: 
% (1) we compute aggregate revenues (i.e. taxes plus fines collected by Gov)
% (2) we compute aggregate revenues separatly from workers and entrep
% (3) we compute tau_s needed to balance the gov budget and we call it tau_s_implied
% NB: (3) is still to be implemented (there is a likely mistake, check!!!)

taxes_w = 0;
taxes_e = 0;

for j=1:negrid
    for t=1: ntgrid
        for i = 1:nagrid_dist
            % Compute pk
            
            kopt   = policycapdet(i,j,t);
            nopt   = policyndet(i,j,t);
            phiopt = policyphidet(i,j,t);
            leopt  = lepoldet(i,j,t);
            profitopt = busincfun(theta(t),kopt,leopt,nopt,r0,w0,gamma,vi,delta);% scalar 
            %profitopt_pos = max(profitopt,0)
            p_k = prob_audit(phiopt,profitopt,theta(t),kopt,nopt,pflag,pn_1,pn_2,pn_3,leopt,gamma,vi);
            
            inc_work = w0*eps_grid(j)*lpolwdet(i,j,t)+ r0*agrid_dist(i);
            % Income entre includes also financial income 
            % bus_inc_entre is "pi"
            % inc_entre is "y_e = pi+ra"
            bus_inc_entre = busincfun(theta(t),kopt,leopt,nopt,r0,w0,gamma,vi,delta); 
            %inc_entre     = bus_inc_entre+r0*agrid_dist(i);
            % can evade only business income (baseline)
            tax_entre_paid = tax_entre(bus_inc_entre*(1-phiopt)+r0*agrid_dist(i),Parameters);
            taxes_w = taxes_w + (1-occpoldet(i,j,t))*mu(i,j,t)*tax_work(inc_work,Parameters);
            taxes_e = taxes_e + occpoldet(i,j,t)* mu(i,j,t)*(tax_entre_paid  +  ...
                p_k*s*( tax_entre(bus_inc_entre +r0*agrid_dist(i),Parameters )-tax_entre_paid )); % SE
        end
    end
end

% taxes_e = Taxes paid by entre + fine revenues
taxes_total = taxes_w + taxes_e;

%% Some aggregate moments: output, K/Y, tfp and gov_spending/output
% For other targets, see targets_compute_clean.m

%% Pack some results

agg.k_supply    = k_supply;
agg.n_supply    = n_supply;
agg.taxes_total = taxes_total;
agg.taxes_w     = taxes_w;
agg.taxes_e     = taxes_e;

end %END FUNCTION "aggregates.m"
