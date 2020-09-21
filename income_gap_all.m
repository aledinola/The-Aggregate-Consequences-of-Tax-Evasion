function [inc_gap_all,taxev_quint_inc,share_entre_all_quint] = ...
    income_gap_all(Parameters,policy,distrib,pregov_incse_bus,pregov_incse,pregov_incw,pregov_inc)

% This function "INCOME GAP ALL" is called by targets_compute_clean.m
% We follow Johns and Slemrod (2010) as closely as possible

% INPUTS
% Parameters,policy,distrib: structures with equilibrium objects
% pregov_incse_bus: decision rule for business income [nagrid_dist,negrid,ntgrid]
% pregov_incse: decision rule for business income plus r*a [nagrid_dist,negrid,ntgrid]
% pregov_inc: decision rule for income all, considering occpol

% OUTPUTS: 
% inc_gap_all
% taxev_quint_inc
% share_entre_all_quint
% % taxev_dec_inc1 NOT USED
%-------------------------------------------------------------------------%

occpoldet    = policy.occpoldet;
policyphidet = policy.policyphidet;
mu           = distrib.mu;


%% Income Gap
% Share of misreported income out of total true income (declared+undeclared) 
% For all agents. Wrt income gap entre, only the denominator changes
unrep_inc_entre = 0;
true_inc_all    = 0;
true_inc_w      = 0;
true_inc_entre  = 0;
for i=1:Parameters.nagrid_dist
    for j=1:Parameters.negrid
        for t=1:Parameters.ntgrid
            unrep_inc_entre = unrep_inc_entre + occpoldet(i,j,t)*policyphidet(i,j,t)*pregov_incse_bus(i,j,t)*mu(i,j,t);
            true_inc_all = true_inc_all + abs(pregov_inc(i,j,t))*mu(i,j,t);
            true_inc_w = true_inc_w + (1-occpoldet(i,j,t))*abs(pregov_incw(i,j,t))*mu(i,j,t);
            true_inc_entre = true_inc_entre + occpoldet(i,j,t)*abs(pregov_incse(i,j,t))*mu(i,j,t);
            % We take abs value of each term in denominator as suggested by
            % Johns and Slemrod (2010, footnote 11 pag.403)
        end 
    end
end
% Income gap = NMP for income, see Johns and Slemrod (2010, Table 2, column 1)
inc_gap_all = unrep_inc_entre/true_inc_all;
inc_gap_check = unrep_inc_entre/(true_inc_w+true_inc_entre);

if abs(inc_gap_all-inc_gap_check)>1e-5
    warning('Check income_gap_all')
end

%% Vectorize stuff columnwise
% From 3-D [nagrid_dist,negrid,ntgrid] to 1-D [nagrid_dist*negrid*ntgrid,1]
pregov_incse_bus_vec  = pregov_incse_bus(:);
pregov_inc_vec        = pregov_inc(:);
mu_vec                = mu(:);
policyphidet_vec      = policyphidet(:);
occpoldet_vec         = occpoldet(:);

%% Sort things by income:
% Consider self-employed & workers (ALL)
[income_sorted, ind_income_sorted] = sort(pregov_inc_vec);
mu_sorted = mu_vec(ind_income_sorted);
pregov_incse_bus_vec_sorted = pregov_incse_bus_vec(ind_income_sorted); 
% We sort self-employed by their total income but then we consider only
% their business income for tax evasion

%taxevasion_vec_sorted = taxevasion_vec(ind_income_sorted);
policyphidet_vec_sorted = policyphidet_vec(ind_income_sorted);
occpoldet_vec_sorted = occpoldet_vec(ind_income_sorted);

% Cumulative distribution
mu_sorted_norm = mu_sorted / sum(mu_sorted);
mu_sorted_norm_cumulative = cumsum(mu_sorted_norm);


%% Quintiles

[~, ind_inc_perc20] = min(abs(mu_sorted_norm_cumulative-0.20));
[~, ind_inc_perc40] = min(abs(mu_sorted_norm_cumulative-0.40));
[~, ind_inc_perc60] = min(abs(mu_sorted_norm_cumulative-0.60));
[~, ind_inc_perc80] = min(abs(mu_sorted_norm_cumulative-0.80));

%Quintile Indices for income:
quint_inc{1} =  1:ind_inc_perc20;
quint_inc{2} =  ind_inc_perc20+1:ind_inc_perc40;
quint_inc{3} =  ind_inc_perc40+1:ind_inc_perc60;
quint_inc{4} =  ind_inc_perc60+1:ind_inc_perc80;
quint_inc{5} =  ind_inc_perc80+1:length(mu_sorted);

% Income gap by quintile of self-employed income
taxev_quint_inc       = zeros(5,1);
weight5               = zeros(5,1);
share_entre_all_quint = zeros(5,1);

for i=1:5
    taxev_quint_inc(i) = sum( mu_sorted(quint_inc{i}).*occpoldet_vec_sorted(quint_inc{i}).*policyphidet_vec_sorted(quint_inc{i}).*pregov_incse_bus_vec_sorted(quint_inc{i}) )/...
        sum( mu_sorted(quint_inc{i}).*income_sorted(quint_inc{i}) );
    weight5(i) = sum( mu_sorted(quint_inc{i}).*income_sorted(quint_inc{i}) )/sum(mu_sorted.*income_sorted);
    share_entre_all_quint(i) = sum( mu_sorted(quint_inc{i}).*occpoldet_vec_sorted(quint_inc{i}) )/sum(mu_sorted(quint_inc{i}));
  
end

% check 
sum(weight5); % this should be exactly 1
sum(taxev_quint_inc.*weight5); % this should be exactly inc_gap_all
%junk = sum(share_entre_all_quint*0.2); % this should give exactly share_entre

end %END FUNCTION
