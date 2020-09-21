function [ModelResults,modelTargetsS] = targets_compute_clean(r0,w0,tau_s,scale,Parameters,Grids,value,distrib,policy,agg,flags,ED)
%{
DESCRIPTION:
 In this function we compute model targets. After computing the targets, 
 we export them in txt file "targets_model_manual.txt" and we return them 
 to the calling program in an M-structure.
---------------------------------------------------------------------------
INPUTS:
 r0,w0: interest rate and wage
 tau_s,scale: tax function parameters
 Parameters,Grids: structures with model primitives
 value,distrib,policy,agg: structures with equilibrium objects
 flags: structure with some flags
 ED: excess demand (a scalar)
---------------------------------------------------------------------------
OUTPUTS:
 ModelResults: M-structure with model results
 modelTargetsS: M-structure with model targets
---------------------------------------------------------------------------
 AUXILIARY FILES:
 This file calls the following functions:
 + fun_exit_rate
    share_entre_by_wealth
    share_entre_by_income
    income_gap_all
    income_gap_only_entre
    taxes_gap_only_entre
    taxes_gap_all
    leverage_fun
    prob_businc
 + 
---------------------------------------------------------------------------
 NOTES:
 %% Updated by Alessandro Di Nola on September 16, 2020. Folder: RED.
%}

% Validate inputs
if ~isstruct(Parameters)
    error('Input "Parameters" must be a structure!')
end
if ~isstruct(Grids)
    error('Input "Grids" must be a structure!')
end
if ~isstruct(value)
    error('Input "value" must be a structure!')
end
if ~isstruct(distrib)
    error('Input "distrib" must be a structure!')
end
if ~isstruct(policy)
    error('Input "policy" must be a structure!')
end
if ~isstruct(agg)
    error('Input "agg" must be a structure!')
end
if ~isstruct(flags)
    error('Input "flags" must be a structure!')
end

%% Unpack structures
% The structure "policy" is visible in the host function, "f_solve_model"

Parameters.tau_s    = tau_s;
Parameters.scale    = scale;
Parameters.scale_se = 1;%scale_se;

%Parameters
alpha          = Parameters.alpha;
vi             = Parameters.vi;
delta          = Parameters.delta;  
lambda         = Parameters.lambda;  
pn_1           = Parameters.pn_1;
pn_2           = Parameters.pn_2;
gamma          = Parameters.gamma;
cc0            = Parameters.cc0;
cc1            = Parameters.cc1;
cc2            = Parameters.cc2;
nagrid_dist    = Parameters.nagrid_dist;
ntgrid         = Parameters.ntgrid;
negrid         = Parameters.negrid;
nngrid         = Parameters.nngrid; 

%Grids
agrid          = Grids.agrid;
agrid_dist     = Grids.agrid_dist;
eps_grid       = Grids.eps_grid;
theta          = Grids.theta;
ngrid          = Grids.ngrid;


%Value functions
Vw  = value.Vw;
Vse = value.Vse;
V1  = value.V1;

%Distributions
mu      = distrib.mu;
mu_work = distrib.mu_work;
mu_se   = distrib.mu_se;

%Policy functions
occpoldet     = policy.occpoldet;
policycapdet  = policy.policycapdet;
policycap     = policy.policycap;
lpolwdet      = policy.lpolwdet;
lepoldet      = policy.lepoldet;
policyndet    = policy.policyndet;
policyphidet  = policy.policyphidet;
apolse1det    = policy.apolse1det;
apolse0det    = policy.apolse0det;
apolwdet      = policy.apolwdet;
% Aggregates
k_supply      = agg.k_supply;
n_supply      = agg.n_supply;  
taxes_total   = agg.taxes_total;
taxes_w       = agg.taxes_w;
taxes_e       = agg.taxes_e;

%% Fraction of entrepreneurs    
share_entre = sum(sum(sum(occpoldet.*mu))); 

%% Share of entre income

% Share of entre income over total income
taxable_income_work1         = zeros(nagrid_dist,1);
taxable_income_entre1        = zeros(nagrid_dist,1);    %TRUE taxable income
taxable_income_entre1_decl   = zeros(nagrid_dist,1);    %DECLARED income 
for j=1:negrid
    for t=1: ntgrid
        % work: added endog labor supply
        taxable_income_work1 = taxable_income_work1  + (1-occpoldet(:,j,t)).*(w0*repmat(eps_grid(j),nagrid_dist,1).*lpolwdet(:,j,t)+ r0*agrid_dist).*mu(:,j,t);
        businc = busincfun(repmat(theta(t),nagrid_dist,1),policycapdet(:,j,t),lepoldet(:,j,t),policyndet(:,j,t),r0,w0,gamma,vi,delta);
        taxable_income_entre1 = taxable_income_entre1 +  occpoldet(:,j,t).*(businc+r0*agrid_dist).*mu(:,j,t);
        taxable_income_entre1_decl = taxable_income_entre1_decl+occpoldet(:,j,t).*((1-policyphidet(:,j,t)).*businc+r0*agrid_dist).*mu(:,j,t);
    end
end

taxable_income_work        = sum(taxable_income_work1);
taxable_income_entre       = sum(taxable_income_entre1);
taxable_income_entre_decl  = sum(taxable_income_entre1_decl);

share_inc_entre       = (taxable_income_entre/(taxable_income_entre+taxable_income_work)); % share of entre income to total income
share_inc_entre_decl  = (taxable_income_entre_decl/(taxable_income_entre_decl+taxable_income_work)); % share of declared entre income to total income 

%% Exit rate (overall) exit_ew

if flags.lowMemory==0
    %----------------- FROM ENTRE to WORK ------------------------------------%
    % fun_exit_rate calls the function trans_operator
    exit_ew = fun_exit_rate(Parameters,Grids,policy,flags,w0,r0,mu);
elseif flags.lowMemory==1
    exit_ew = nan;
end


%% T7 - Capital used by entre OUT

% Vectorize the distribution of agents over (a,eps,theta)
mu_vec      = mu(:); % [nagrid_dist*negrid*ntgrid,1]
mu_work_vec = mu_work(:);
mu_se_vec   = mu_se(:);

% Compute useful arrays
wealth      = repmat(agrid_dist,[1 negrid ntgrid]); % array [nagrid_dist,negrid,ntgrid]
wealth_vec  = wealth(:); % array [nagrid_dist*negrid*ntgrid,1]
wealth_se   = zeros(nagrid_dist,negrid, ntgrid);
wealth_work = zeros(nagrid_dist,negrid, ntgrid);

for j=1:negrid
    for t=1: ntgrid
        for i = 1:nagrid_dist
            wealth_se(i,j,t)   = agrid_dist(i);
            wealth_work(i,j,t) = agrid_dist(i);
        end
    end
end 

wealth_se_vec    = wealth_se(:);
wealth_work_vec  = wealth_work(:);
wealth_tot       = wealth_vec'*mu_vec;

%% Assets owned by entre
%assets held by entre
assets_e1 = zeros(nagrid_dist,1);
assets_w1 = zeros(nagrid_dist,1); %assets held by workers

for j=1:negrid
    for t=1:ntgrid
        assets_e1  = assets_e1 + occpoldet(:,j,t).*agrid_dist.*mu(:,j,t);
        assets_w1  = assets_w1 + (1-occpoldet(:,j,t)).*agrid_dist.*mu(:,j,t);
    end
end

assets_e = sum(assets_e1);
% Fraction of total assets held by ENTRE
assets_e_to_work = (assets_e/wealth_tot); 

%% Ratio of median assets (entre/worker)

% Compute quantiles of wealth
med_assets_work   = quantili(wealth_vec,mu_work_vec,0.50);
med_assets_se     = quantili(wealth_vec,mu_se_vec,0.50);
wealth_med_ratio_E_W = med_assets_se/med_assets_work;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          EQUILIBRIUM DECISION RULES                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute decision rules as policy_function(a,e,theta)

% Decision rule for WORK INCOME
%pregov_incw1 = w0*repmat(eps_grid',nagrid_dist,1,ntgrid)+r0*repmat(agrid_dist',1,negrid,ntgrid); 
pregov_incw = zeros(nagrid_dist,negrid,ntgrid);
for j=1:negrid
    for t=1:ntgrid
        % ENDO LABOR
        pregov_incw(:,j,t) = w0*eps_grid(j)*lpolwdet(:,j,t)+ r0*agrid_dist;
    end
end
% Decision rule for ENTRE INCOME
pregov_incse        = zeros(nagrid_dist,negrid,ntgrid);
pregov_incse_bus    = zeros(nagrid_dist,negrid,ntgrid); % NOT include asset income
pregov_incse_bus_decl = zeros(nagrid_dist,negrid,ntgrid);
pregov_incse_decl   = zeros(nagrid_dist,negrid,ntgrid); % declared SE income (PSID)
for j=1:negrid
    for t=1:ntgrid
        % *true* SE income w/out financial assets (i.e. business income)
        pregov_incse_bus(:,j,t) = busincfun(repmat(theta(t),nagrid_dist,1),...
            policycapdet(:,j,t),lepoldet(:,j,t),policyndet(:,j,t),r0,w0,gamma,vi,delta);
        % *declared* SE income w/out financial assets (i.e. business income)
        pregov_incse_bus_decl(:,j,t) = (1-policyphidet(:,j,t)).*pregov_incse_bus(:,j,t);
        % *true* SE income with financial assets
        pregov_incse(:,j,t)      = pregov_incse_bus(:,j,t)+r0*agrid_dist-cost_evasion(policyphidet(:,j,t),cc0,cc1,cc2);
        % *declared* SE income with financial assets
        pregov_incse_decl(:,j,t) = (1-policyphidet(:,j,t)).*pregov_incse_bus(:,j,t)+r0*agrid_dist-cost_evasion(policyphidet(:,j,t),cc0,cc1,cc2);
    end
end

% Decision rule for INCOME (all)
pregov_inc      = zeros(nagrid_dist,negrid,ntgrid);
pregov_inc_decl = zeros(nagrid_dist,negrid,ntgrid);
for j=1:negrid
    for t=1: ntgrid
        pregov_inc(:,j,t)      = (1-occpoldet(:,j,t)).*pregov_incw(:,j,t)+occpoldet(:,j,t).*pregov_incse(:,j,t);
        pregov_inc_decl(:,j,t) = (1-occpoldet(:,j,t)).*pregov_incw(:,j,t)+occpoldet(:,j,t).*pregov_incse_decl(:,j,t);
    end
end

% Decision rule for TAX GAP
% UNPAID TAXES AND TRUE TAX LIABILITY
% TAX GAP = unpaid taxes/true tax liability
unpaid_taxes   = zeros(nagrid_dist,negrid,ntgrid);
paid_taxes     = zeros(nagrid_dist,negrid,ntgrid);
true_tax       = zeros(nagrid_dist,negrid,ntgrid);
tax_gap_single = zeros(nagrid_dist,negrid,ntgrid);

for i=1:nagrid_dist
    for j=1:negrid
        for t=1: ntgrid
            if     occpoldet(i,j,t)==1 % SE
                   true_tax(i,j,t)     = tax_entre(pregov_incse_bus(i,j,t) + r0*agrid_dist(i),Parameters);
                   paid_taxes(i,j,t)   = tax_entre((1-policyphidet(i,j,t))*pregov_incse_bus(i,j,t) + r0*agrid_dist(i),Parameters);
                   unpaid_taxes(i,j,t) = true_tax(i,j,t)-paid_taxes(i,j,t);
            elseif occpoldet(i,j,t)==0 % WORK
                   true_tax(i,j,t)     = tax_work(pregov_incw(i,j,t),Parameters);
                   paid_taxes(i,j,t)   = tax_work(pregov_incw(i,j,t),Parameters);
                   unpaid_taxes(i,j,t) = true_tax(i,j,t)-paid_taxes(i,j,t);% For workers unpaid taxes=0 by construction
            end
                   tax_gap_single(i,j,t) = unpaid_taxes(i,j,t)/abs(true_tax(i,j,t));
        end
    end
end

%% Vectorize equilibrium decision rules

% Vectorize stuff columnwise
% From 3-D [nagrid_dist,negrid,ntgrid] to 1-D [nagrid_dist*negrid*ntgrid,1]
pregov_incse_vec      = pregov_incse(:);
pregov_incse_decl_vec = pregov_incse_decl(:);
%pregov_incse_bus_vec  = pregov_incse_bus(:);
%pregov_incse_bus_decl_vec = pregov_incse_bus_decl(:);
pregov_inc_vec        = pregov_inc(:);

unpaid_taxes_vec      = unpaid_taxes(:);
true_tax_vec          = true_tax(:);
tax_gap_single_vec    = tax_gap_single(:);

lpolwdet_vec          = lpolwdet(:);
%lepoldet_vec          = lepoldet(:);
policyndet_vec        = policyndet(:);
apolse1det_vec        = apolse1det(:);
apolse0det_vec        = apolse0det(:);

temp_input = zeros(nagrid_dist,negrid,ntgrid);
pkmat      = zeros(nagrid_dist,negrid,ntgrid);
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            temp_input(i,j,t) = policyphidet(i,j,t)*pregov_incse_bus(i,j,t); %phi*pi
            pkmat(i,j,t) =  1./(1 + pn_1*exp(-pn_2*temp_input(i,j,t)));
       end
    end
end

%pkmat                 = repmat(pk,[1 negrid ntgrid]); % array [nagrid_dist,negrid,ntgrid]
pkmat_vec             = pkmat(:);

% Compute mean/median ratios for wealth and income
wealth_median              = quantili(wealth_vec,mu_vec,0.50);
wealth_mean                = sum(wealth_vec.*mu_vec);
wealth_mean_to_median      = wealth_mean/wealth_median;

% Only for SE
wealth_median_se           = quantili(wealth_vec,mu_se_vec,0.50);
wealth_mean_se             = sum(wealth_se_vec.*mu_se_vec);
wealth_mean_to_median_se   = wealth_mean_se/wealth_median_se;

% Only for WORK
wealth_median_work         = quantili(wealth_vec,mu_work_vec,0.50);
wealth_mean_work           = sum(wealth_work_vec.*mu_work_vec);
wealth_mean_to_median_work = wealth_mean_work/wealth_median_work;

% For ALL
inc_median                 = quantili(pregov_inc_vec,mu_vec,0.50);
inc_mean                   = sum(pregov_inc_vec.*mu_vec);
inc_mean_to_median         = inc_mean/inc_median;

% For ENTRE true
incse_median                 = quantili(pregov_incse_vec,mu_se_vec,0.50);
incse_mean                   = sum(pregov_incse_vec.*mu_se_vec);
incse_mean_to_median         = incse_mean/incse_median; %used in txt_export

% For ENTRE declared
incse_median_decl                 = quantili(pregov_incse_decl_vec,mu_se_vec,0.50);
incse_mean_decl                   = sum(pregov_incse_decl_vec.*mu_se_vec);
incse_mean_to_median_decl         = incse_mean_decl/incse_median_decl;

%% Share of entre by wealth
% share_entre is the share in the whole population
%[share_entre_wealth_quint] = share_entre_by_wealth(policy,distrib,wealth,share_entre);

%% .. and by income
[share_entre_inc_quint] = share_entre_by_income(policy,distrib,pregov_inc);

%% Average prob. of auditing for each business income quintile
%We compute it using two different definitions of SE income
%[audit_businc]      = prob_businc(pregov_incse_bus_vec,mu_se_vec,policycapdet,pn_1,pn_2,pn_3,pflag);
%[audit_businc_decl] = prob_businc(pregov_incse_decl_vec,mu_se_vec,policycapdet,pn_1,pn_2,pn_3,pflag);
%[audit_businc_decl] = prob_businc(pregov_incse_bus_decl_vec,mu_se_vec,policycapdet,pn_1,pn_2,pn_3,pflag);

%% Tax evasion statistics (see Slemrod 2010)

% Call file to compute tax gap: aggregate, by deciles/quintiles of true income
%--------------------------------------------%
% We have 4 possible options:                %
% 1 - income gap ALL                         %
% 2 - income gap ONLY ENTRE                  %
% 3 - tax gap ALL                            %
% 4 - tax gap ONLY ENTRE                     %
%--------------------------------------------%
% consider whole population
[inc_gap_all,taxev_quint_inc] = income_gap_all(Parameters,policy,distrib,...
    pregov_incse_bus,pregov_incse,pregov_incw,pregov_inc);        
% consider only self-employed
[inc_gap_entre] = income_gap_only_entre(Parameters,....
    policy,distrib,pregov_incse,pregov_incse_bus);
[tax_gap_entre] = taxes_gap_only_entre(Parameters,policy,distrib,Grids,r0,pregov_incse_vec,pregov_incse_bus,tax_gap_single,unpaid_taxes_vec,true_tax_vec,tax_gap_single_vec);  % consider only self-employed
[tax_gap_all] = taxes_gap_all(Parameters,policy,distrib,Grids,r0,pregov_inc_vec,pregov_incse_bus,pregov_incw,tax_gap_single,unpaid_taxes_vec,true_tax_vec);         % consider whole population


%% Labor targets, revised version RED

% Average hours worked by workers. Policy function is lpolwdet
ave_n_work = sum(lpolwdet_vec.*mu_work_vec);

% Average hours worked by SE
%ave_n_entre_num = sum(lepoldet_vec.*mu_se_vec);

%Total labor supply of self-employed
le_entre   = 0;
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            le_entre = le_entre+occpoldet(i,j,t)*lepoldet(i,j,t)*mu(i,j,t);
        end
    end
end
%Can be written more compactly as 
%le_entre   = sum(lepoldet_vec.*mu_se_vec)*share_entre;

% Share of Hiring Self-employed (data: share of unincorp SE with paid
% employees, BLS)
if nngrid>1
    dist1 = 0.5*(ngrid(2)-ngrid(1)); % ngrid(1)=0 by construction
    share_hiring_se = sum((policyndet_vec>(ngrid(1)+dist1)).*mu_se_vec);
else
    share_hiring_se = 0;
end
% Share of workers hired by SE 
n_corp = n_supply;%((1-alpha)/w0)^(1/alpha)*k_supply; % labor demand from Corp.sec
n_entre  = sum(policyndet_vec.*mu_se_vec)*share_entre;
%n_entre should be exactly equal to n_entre1
%n_entre1 = sum(sum(sum(occpoldet.*policyndet.*mu)));
share_workInSe = n_entre/(n_entre+n_corp);

% Distribution of hiring self-employed
% Rescaling parameter to translate eff_units to the number of hired workers
% 
% E(ability|worker) = 1 and we target ave_n_work to be 1/3, hence  ksi_n = ave_n_work*cond_mean_eps_work = 1/3;
% We use it to translate grid points from data to the model
% n_bins = (0, 2.5/3, 7/3, 14.5/3, 30/3) in model units of efficiency labor

% Calculate distribution of hirings

scale_hiring = 3; %I use 2, Anna uses 3

unc_n_firm_size_dist = zeros(nngrid,1);

A = repmat(ngrid,[1 length(policyndet_vec)]);
%policyndet_vec_help = repmat(policyndet_vec,[1 nngrid]);
%[minValue,closestIndex] = min(abs(A-policyndet_vec_help'));
[~,closestIndex] = min(abs(A-policyndet_vec'));
for i= 1:nngrid
    unc_n_firm_size_dist(i) = sum(mu_se_vec(closestIndex==i));
end

cond_firm_size_dist = zeros(4,1);

norm_need = sum(unc_n_firm_size_dist(2:end)); %hiring labor firms - mass

%these two numbers should be the same
%disp([share_hiring_se, norm_need])

if nngrid < 6
        
    if norm_need > 0
        
        cond_firm_size_dist = 100 * unc_n_firm_size_dist(2:end) /  norm_need;
        
    end
    
else
    
    bins = [0, 5, 10, 20, ngrid(end)*scale_hiring+1];
    
    for i = 2: length(bins)
    
     ind_ngrid_bodies = ngrid*scale_hiring<=bins(i)&ngrid*scale_hiring>bins(i-1);
     
     if norm_need > 0
    
     cond_firm_size_dist(i-1)  = 100  * sum(unc_n_firm_size_dist(ind_ngrid_bodies)) /  norm_need;
     
     end
     
    end
    
end

%% Capital in both sectors

k_entre = sum(sum(sum(occpoldet.*policycapdet.*mu)));
k_corp  = k_supply; %for illustrative purpose :) we use k_supply but it is k_corp
% Aggregate capital in the economy
capital = k_entre + k_supply;

%% Output in both sectors
%  Aggregate output in entrep sector
output_se_mat = zeros(nagrid_dist,negrid,ntgrid);
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            k = policycapdet(i,j,t);
            n = policyndet(i,j,t);
            l = lepoldet(i,j,t);
            output_se_mat(i,j,t) = theta(t).*prodfun(k,l,n,gamma,vi);
        end
    end
end

output_se = sum(sum(sum(occpoldet.*output_se_mat.*mu)));

% Aggregate output in the corporate sector
output_corp = k_supply^alpha*n_supply^(1-alpha); %n_supply or n_corp

% Aggregate output in the economy
output = output_se+output_corp;

%% K/Y ratio
K_to_Y = capital/output;

taxes_total_toY = taxes_total / output; %taxes_total calculated in aggregates.m
taxes_total_toY_p = taxes_total_toY*100;

%% Net worth by percentiles of asset holdings 

%% Sort things by net worth:
% [income_sorted, ind_income_sorted] = sort(pregov_incse_vec);
% mu_sorted = mu_se_vec(ind_income_sorted);

[wealth_se_vec_sorted, ind_wealth_sorted] = sort(wealth_se_vec); % maybe use apol?
mu_se_sorted                              = mu_se_vec(ind_wealth_sorted);


% Sort policy functions for a' and pk by net worth

apolse1det_vec_sorted   = apolse1det_vec(ind_wealth_sorted);
apolse0det_vec_sorted   = apolse0det_vec(ind_wealth_sorted);
apol_vec_sorted         = pkmat_vec(ind_wealth_sorted).*apolse1det_vec_sorted+(1-pkmat_vec(ind_wealth_sorted)).*apolse0det_vec_sorted;

% Cumulative distribution
mu_se_sorted_norm            = mu_se_sorted / sum(mu_se_sorted);
mu_se_sorted_norm_cumulative = cumsum(mu_se_sorted_norm);

%% Deciles

[~, ind_wealth_perc10] = min(abs(mu_se_sorted_norm_cumulative-0.10));
[~, ind_wealth_perc20] = min(abs(mu_se_sorted_norm_cumulative-0.20));
[~, ind_wealth_perc30] = min(abs(mu_se_sorted_norm_cumulative-0.30));
[~, ind_wealth_perc40] = min(abs(mu_se_sorted_norm_cumulative-0.40));
[~, ind_wealth_perc50] = min(abs(mu_se_sorted_norm_cumulative-0.50));
[~, ind_wealth_perc60] = min(abs(mu_se_sorted_norm_cumulative-0.60));
[~, ind_wealth_perc70] = min(abs(mu_se_sorted_norm_cumulative-0.70));
[~, ind_wealth_perc80] = min(abs(mu_se_sorted_norm_cumulative-0.80));
[~, ind_wealth_perc90] = min(abs(mu_se_sorted_norm_cumulative-0.90));

%Decile Indices for income:
dec_wealth{1}  =  1:ind_wealth_perc10;
dec_wealth{2}  =  ind_wealth_perc10+1:ind_wealth_perc20;
dec_wealth{3}  =  ind_wealth_perc20+1:ind_wealth_perc30;
dec_wealth{4}  =  ind_wealth_perc30+1:ind_wealth_perc40;
dec_wealth{5}  =  ind_wealth_perc40+1:ind_wealth_perc50;
dec_wealth{6}  =  ind_wealth_perc50+1:ind_wealth_perc60;
dec_wealth{7}  =  ind_wealth_perc60+1:ind_wealth_perc70;
dec_wealth{8}  =  ind_wealth_perc70+1:ind_wealth_perc80;
dec_wealth{9}  =  ind_wealth_perc80+1:ind_wealth_perc90;
dec_wealth{10} = ind_wealth_perc90+1:length(mu_se_sorted);

apol_dec_wealth          = zeros(10,1);
apol_dec_wealth_norm     = zeros(10,1); %normalized by E(a|occ=e)
apol_dec_wealth_bis      = zeros(10,1);

for i=1:10
    apol_dec_wealth_bis(i)  = sum( mu_se_sorted(dec_wealth{i}).*apol_vec_sorted(dec_wealth{i}) )/...
                              sum( mu_se_sorted(dec_wealth{i}));
    
    apol_dec_wealth(i)      = sum( mu_se_sorted(dec_wealth{i}).*wealth_se_vec_sorted(dec_wealth{i}) )/...
                              sum( mu_se_sorted(dec_wealth{i}));
    
    apol_dec_wealth_norm(i) = apol_dec_wealth(i)/wealth_mean_se;
end

%% Lorenz Curve and Gini of WEALTH and INCOME

% Lorenz for wealth ALL
[pw, ~, ~, gini_wealth]          = lrzcurve(mu_vec,wealth_vec,'g');
% Lorenz for wealth SE 
[pw_se, ~, ~, gini_wealth_se]    = lrzcurve(mu_se_vec,wealth_se_vec,'g');
% Lorenz for wealth W 
[pw_work, ~, ~, gini_wealth_work] = lrzcurve(mu_work_vec,wealth_work_vec,'g');
% Lorenz for Income SE true
[pw_incse, ~, ~, gini_incse]      = lrzcurve(mu_se_vec,pregov_incse_vec,'b'); 
% Lorenz for Income SE declared
[pw_incse_decl, ~, ~, gini_incse_decl] = lrzcurve(mu_se_vec,pregov_incse_decl(:),'b'); 
% Lorenz for Income ALL, true
[pw_inc, ~, ~, gini_inc_true]          = lrzcurve(mu_vec,pregov_inc_vec,'b');%use for targets
% Lorenz for Income ALL, declared
[~, ~, ~, gini_inc_decl]               = lrzcurve(mu_vec,pregov_inc_decl(:),'b');%use for targets

%% Compute specific percentiles of the wealth/income distribution
% ALL AGENTS
% Fraction of total wealth held by...
% Bottom 40% /% Top 20% /% Top 10% /% Top 1%
a_bot40 = interp1q(pw(:,1),pw(:,2),0.40);
a_bot80 = interp1q(pw(:,1),pw(:,2),0.80);
a_top20 = 1-a_bot80;
a_bot90 = interp1q(pw(:,1),pw(:,2),0.90);
a_top10 = 1-a_bot90;
a_bot99 = interp1q(pw(:,1),pw(:,2),0.99);
a_top1  = 1-a_bot99;

wealth_perc = [a_bot40 a_top20 a_top10 a_top1];

% SELF-EMPLOYED ONLY
% Fraction of total wealth held by...
% Bottom 40%/% Top 20% /% Top 10% /% Top 1%
a_bot40_se = interp1q(pw_se(:,1),pw_se(:,2),0.40);
a_bot80_se = interp1q(pw_se(:,1),pw_se(:,2),0.80);
a_top20_se = 1-a_bot80_se;
a_bot90_se = interp1q(pw_se(:,1),pw_se(:,2),0.90);
a_top10_se = 1-a_bot90_se;
a_bot99_se = interp1q(pw_se(:,1),pw_se(:,2),0.99);
a_top1_se  = 1-a_bot99_se;

wealth_perc_se = [a_bot40_se a_top20_se a_top10_se a_top1_se];

% WORKERS ONLY
% Fraction of total wealth held by...
% Bottom 40% / % Top 20% /% Top 10% /% Top 1%
a_bot40_work = interp1q(pw_work(:,1),pw_work(:,2),0.40);
a_bot80_work = interp1q(pw_work(:,1),pw_work(:,2),0.80);
a_top20_work = 1-a_bot80_work;
a_bot90_work = interp1q(pw_work(:,1),pw_work(:,2),0.90);
a_top10_work = 1-a_bot90_work;
a_bot99_work = interp1q(pw_work(:,1),pw_work(:,2),0.99);
a_top1_work  = 1-a_bot99_work;

wealth_perc_work = [a_bot40_work a_top20_work a_top10_work a_top1_work];

% FRACTION OF CAPITAL HELD BY

% Fraction of total income held by...
% Bottom 40%/ % Top 20% /% Top 10% /% Top 1%
I_bot40 = interp1q(pw_inc(:,1),pw_inc(:,2),0.40);
I_bot80 = interp1q(pw_inc(:,1),pw_inc(:,2),0.80);
I_top20 = 1-I_bot80;
I_bot90 = interp1q(pw_inc(:,1),pw_inc(:,2),0.90);
I_top10 = 1-I_bot90;
I_bot99 = interp1q(pw_inc(:,1),pw_inc(:,2),0.99);
I_top1  = 1-I_bot99;

inc_perc = [I_bot40 I_top20 I_top10 I_top1];

% SELF-EMPLOYED ONLY -- true income distribution
Ise_bot40 = interp1q(pw_incse(:,1),pw_incse(:,2),0.40);
Ise_bot80 = interp1q(pw_incse(:,1),pw_incse(:,2),0.80);
Ise_top20 = 1-Ise_bot80;
Ise_bot90 = interp1q(pw_incse(:,1),pw_incse(:,2),0.90);
Ise_top10 = 1-Ise_bot90;
Ise_bot99 = interp1q(pw_incse(:,1),pw_incse(:,2),0.99);
Ise_top1  = 1-Ise_bot99;

incse_perc   = [Ise_bot40 Ise_top20 Ise_top10 Ise_top1];
incse_perc_p = 100*incse_perc;

% SELF-EMPLOYED ONLY -- declared income distribution
Ise_bot40_decl = interp1q(pw_incse_decl(:,1),pw_incse_decl(:,2),0.40);
Ise_bot80_decl = interp1q(pw_incse_decl(:,1),pw_incse_decl(:,2),0.80);
Ise_top20_decl = 1-Ise_bot80_decl;
Ise_bot90_decl = interp1q(pw_incse_decl(:,1),pw_incse_decl(:,2),0.90);
Ise_top10_decl = 1-Ise_bot90_decl;
Ise_bot99_decl = interp1q(pw_incse_decl(:,1),pw_incse_decl(:,2),0.99);
Ise_top1_decl  = 1-Ise_bot99_decl;

incse_perc_decl = [Ise_bot40_decl Ise_top20_decl Ise_top10_decl Ise_top1_decl];
incse_perc_decl_p = 100*incse_perc_decl; %used by txt_export


%% Firm size distribution

% Average capital given occpol=1
ave_k_entre = 0;
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            ave_k_entre = ave_k_entre+policycapdet(i,j,t)*mu_se(i,j,t);
        end
    end
end

n_cut = 0;%0.5*(ngrid(2)-ngrid(1)); % ngrid(1)=0 by construction
ave_n_entre_num = 0;
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            if policyndet(i,j,t)>=n_cut
                ave_n_entre_num = ave_n_entre_num+policyndet(i,j,t)*mu_se(i,j,t);
            end
        end
    end
end
ave_n_entre_den = 0; %this should be equal to the share of hiring SE
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            if policyndet(i,j,t)>=n_cut
                ave_n_entre_den = ave_n_entre_den+mu_se(i,j,t);
            end
        end
    end
end

ave_n_entre = (ave_n_entre_num/ave_n_entre_den)+1;
%First bin starts at 1 (not 0)

%% Conditional distribution (pdf) of theta

prob_th_unc       = reshape(sum(sum(mu,1),2),[ntgrid,1]);
%temp              = Ptheta^1000;
% disp('Consistency check:')
% disp([prob_th_unc,prob_th_unc_check])

prob_th_cond = zeros(ntgrid,1);
for t=1:ntgrid
    prob_th_cond(t) = sum(sum(occpoldet(:,:,t).*mu(:,:,t)));
end

% "cond_mean_theta" Conditional average of theta (given occpol=1)
sum_prob_th_cond = sum(prob_th_cond); % check this is % of entre
%unc_mean_theta   = theta'*prob_th_unc;
cond_mean_theta  = 0;
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            cond_mean_theta = cond_mean_theta+ theta(t)*((occpoldet(i,j,t)*mu(i,j,t))/sum_prob_th_cond);
        end
    end
end

%% Welfare
 
%interpolate value functions
Vsedet = zeros(nagrid_dist,negrid,ntgrid);
Vwdet = zeros(nagrid_dist,negrid,ntgrid);
V1det = zeros(nagrid_dist,negrid,ntgrid);

for j=1:negrid
    for t=1:ntgrid
        Vsedet(:,j,t) = interp1(agrid, Vse(:,j,t), agrid_dist);
        Vwdet(:,j,t)  = interp1(agrid, Vw(:,j,t), agrid_dist);
        V1det(:,j,t)  = interp1(agrid, V1(:,j,t), agrid_dist); %V1 is max(Vse, Vw)
    end
end
 
%Welfare of entrepreneurs
We=0;
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            We = We + Vsedet(i,j,t)*mu_se(i,j,t);
            % mu_se = occ*mu/(sum(occ*mu))
        end
    end
end
 
%Welfare of workers
Ww=0;
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            Ww = Ww + Vwdet(i,j,t)*mu_work(i,j,t);
        end
    end
end
 
%Total welfare (average welfare) 
Wtot=We*share_entre + Ww*(1-share_entre); % 

%%
 %Calculate the fraction of borrowing constrained individuals
 
frac_bor=0; %those who have very little assets
 
 for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            if apolwdet(i,j,t) < 1e-3
                frac_bor = frac_bor + mu_se(i,j,t);
            end
        end
    end
 end
 
%fraction of SE who are borrowing constrained 
disp_leverage = 0; %flag 0-1
[frac_bor1,leverage_percentiles,leverage_mean] = leverage_fun(policycapdet,...
    agrid_dist,lambda,mu_se,nagrid_dist,negrid,ntgrid,disp_leverage,agrid,policycap); 
 
% Convert some model targets from fraction to percentage
share_entre_p      = share_entre*100; %share of entre in percent
exit_ew_p          = exit_ew*100;     %exit rate in percent 
tax_gap_entre_p    = tax_gap_entre*100; %tax gap in percent
tax_gap_all_p      = tax_gap_all*100;
inc_gap_all_p      = inc_gap_all*100;
inc_gap_entre_p    = inc_gap_entre*100;
share_inc_entre_p  = share_inc_entre_decl*100;
share_inc_entre_true_p = share_inc_entre*100;
assets_e_to_work_p = assets_e_to_work*100;
wealth_perc_p      = wealth_perc*100;
wealth_perc_se_p   = wealth_perc_se*100;
wealth_perc_work_p = wealth_perc_work*100;
inc_perc_p         = inc_perc*100;
gini_wealth_p      = gini_wealth*100; %gini wealth total in percent
gini_wealth_se_p   = gini_wealth_se*100;   %gini wealth se in percent
gini_wealth_work_p = gini_wealth_work*100;
gini_inc_decl_p    = gini_inc_decl*100;    %gini declared income in percent
gini_inc_true_p    = gini_inc_true*100;    %gini declared income in percent
gini_incse_p       = gini_incse*100;
gini_incse_decl_p  = gini_incse_decl*100;
taxev_quint_inc_p  = taxev_quint_inc*100;
%ave_n_work_p       = ave_n_work*100;
share_hiring_se_p  = share_hiring_se*100;
%share_workInSe_p   = share_workInSe*100;
share_entre_quint_p = share_entre_inc_quint*100; 

%%
%% 
% Show for manual run
% For income and wealth: Gini - Mean/Median - Bottom 40% - Top 20% - Top 10% - Top 1%
model_wealth_all    =   [gini_wealth_p, wealth_mean_to_median, wealth_perc_p(1), wealth_perc_p(2), wealth_perc_p(3), wealth_perc_p(4)]';
model_wealth_entre  =   [gini_wealth_se_p, wealth_mean_to_median_se, wealth_perc_se_p(1), wealth_perc_se_p(2), wealth_perc_se_p(3), wealth_perc_se_p(4)]';
model_wealth_work   =   [gini_wealth_work_p, wealth_mean_to_median_work, wealth_perc_work_p(1), wealth_perc_work_p(2), wealth_perc_work_p(3), wealth_perc_work_p(4)]';
model_income_all    =   [gini_inc_decl_p, inc_mean_to_median, inc_perc_p(1), inc_perc_p(2), inc_perc_p(3), inc_perc_p(4)]';




%-------------------------------------------------------------------------%
% Model targets % order follows memo: Revision: Model and Calibration,
% Table 2, C:\Users\tkhir\Dropbox\ProjectMariia\revision\model_revision_v3
model_targets = [ave_n_work; K_to_Y; share_inc_entre_p; ...
                 share_hiring_se_p; ...
                 exit_ew_p; share_entre_quint_p; share_entre_p; assets_e_to_work_p; ... 
                 taxev_quint_inc_p; inc_gap_all_p;...
                 taxes_total_toY_p;...
                 cond_firm_size_dist];  
             
%% Collect model targets into structure "ModStat"

ModStat = v2struct(share_entre_p,cond_mean_theta,ave_k_entre,ave_n_entre,...
k_entre,k_supply,capital,output,output_corp,output_se,r0,w0,taxes_total_toY_p,...
taxes_total,taxes_w,taxes_e,Wtot,We,Ww,frac_bor,frac_bor1,model_targets,...
share_inc_entre_true_p,wealth_med_ratio_E_W,inc_gap_entre_p,tax_gap_entre_p,...
tax_gap_all_p,ED,leverage_mean,model_wealth_all,model_wealth_entre,model_wealth_work,...
gini_inc_true_p,gini_inc_decl_p,model_income_all,gini_incse_p,incse_mean_to_median,...
incse_perc_p,gini_incse_decl_p,incse_mean_to_median_decl,incse_perc_decl_p);

%% Export in txt files

txt_export(flags,Parameters,ModStat);

%% Collect some moments into structure "modelTargetsS"

modelTargetsS = v2struct(ED,ave_n_work,K_to_Y,share_inc_entre_true_p,share_hiring_se_p,...
    exit_ew_p,gini_incse_p,share_entre_p,taxev_quint_inc_p,inc_gap_all_p,...
    taxes_total_toY_p,leverage_mean,cond_firm_size_dist);

%% Collect some moments into structure "ModelResults" (similar to "Table results in excel")
ModelResults.share_entre_p       = share_entre_p;   %share of SE in pop.
ModelResults.cond_mean_theta     = cond_mean_theta; %E(theta|occpol=SE)
ModelResults.ave_k_entre         = ave_k_entre;     %k_entre/share of SE
ModelResults.ave_n_entre         = ave_n_entre;     %given n>0
ModelResults.k_entre             = k_entre;         %integral of k(x) given occpol=SE
ModelResults.n_entre             = n_entre;         %integral of n(x) given occpol=SE
ModelResults.le_entre            = le_entre;        %integral of l^SE(x) given occpol=SE
ModelResults.output_se           = output_se;       %integral of y(x) given occpol=SE
ModelResults.frac_bor_business   = frac_bor1;       %share of SE who face binding borrowing constraint
%-------------------------------------------------------------------------%
ModelResults.k_corp              = k_corp;          %K^C
ModelResults.n_corp              = n_corp;          %N^C
ModelResults.output_corp         = output_corp;     %Y^C
%-------------------------------------------------------------------------%
ModelResults.share_hiring_se_p   = share_hiring_se_p; %Share of hiring SE
ModelResults.share_workInSe      = share_workInSe;    %Share of workers hired by SE
ModelResults.ave_n_work          = ave_n_work; %Avg. hours worked
%-------------------------------------------------------------------------%
ModelResults.r0                  = r0;
ModelResults.w0                  = w0;
%-------------------------------------------------------------------------%
ModelResults.taxes_total_toY_p   = taxes_total_toY_p;  %total taxes/total output
ModelResults.gini_wealth_p       = gini_wealth_p; % Gini wealth ALL\
ModelResults.gini_wealth_se_p    = gini_wealth_se_p; % Gini wealth SE
ModelResults.gini_wealth_work_p  = gini_wealth_work_p; % Gini wealth WORK
%-------------------------------------------------------------------------%
ModelResults.cond_firm_size_dist = cond_firm_size_dist; %conditional firm size distrib. (n>0)
ModelResults.prob_th_cond        = prob_th_cond; %distrib. of theta conditional on occpol=SE
ModelResults.prob_th_unc         = prob_th_unc; %unconditional distrib. of theta
ModelResults.adist               = sum(sum(distrib.mu,2),3);
%-------------------------------------------------------------------------%
ModelResults.output               = output; %output_se+output_corp
ModelResults.K_to_Y               = K_to_Y;
ModelResults.frac_bor1            = frac_bor1;
ModelResults.leverage_percentiles = leverage_percentiles;
ModelResults.leverage_mean        = leverage_mean;
%-------------------------------------------------------------------------%


end %END FUNCTION "targets_compute_clean"
