function [inc_gap_entre,taxev_quint_inc_entre] = income_gap_only_entre(Parameters,....
    policy,distrib,pregov_incse,pregov_incse_bus)

% This function "INCOME GAP ONLY ENTRE" is called by targets_compute.m
% We follow Johns and Slemrod (2010) as closely as possible

% INPUTS
% Parameters,policy,distrib: structures with equilibrium objects
% pregov_incse_bus: decision rule for business income [nagrid_dist,negrid,ntgrid]
% pregov_incse: decision rule for business income plus r*a [nagrid_dist,negrid,ntgrid]
% pregov_inc: decision rule for income all, considering occpol

% OUTPUTS: 
% inc_gap_entre
% taxev_quint_inc1
%-------------------------------------------------------------------------%

occpoldet    = policy.occpoldet;
policyphidet = policy.policyphidet;
mu           = distrib.mu;
mu_se        = distrib.mu_se;


%% Tax Gap
% Share of misreported income out of total true income (declared+undeclared)
% for self-employed only
unrep_inc_entre = 0;
true_inc_entre = 0;
for i=1:Parameters.nagrid_dist
    for j=1:Parameters.negrid
        for t=1:Parameters.ntgrid
            unrep_inc_entre = unrep_inc_entre + occpoldet(i,j,t)*policyphidet(i,j,t)*pregov_incse_bus(i,j,t)*mu(i,j,t);
            true_inc_entre  = true_inc_entre + abs(occpoldet(i,j,t)*pregov_incse(i,j,t))*mu(i,j,t);
            % We take abs value of each term in denominator as suggested by
            % Johns and Slemrod (2010, footnote 11 pag.403)
        end 
    end
end

% Income gap = NMP for income, see Johns and Slemrod (2010, Table 2, column 1)
inc_gap_entre = unrep_inc_entre/true_inc_entre;

%% Vectorize stuff columnwise
% From 3-D [nagrid_dist,negrid,ntgrid] to 1-D [nagrid_dist*negrid*ntgrid,1]

mu_se_vec            = mu_se(:);
pregov_incse_vec     = pregov_incse(:);
pregov_incse_bus_vec = pregov_incse_bus(:);
policyphidet_vec     = policyphidet(:);


%% Sort things by income:
% Consider ONLY self-employed
[income_sorted, ind_income_sorted] = sort(pregov_incse_vec);
mu_sorted = mu_se_vec(ind_income_sorted);
pregov_incse_bus_vec_sorted = pregov_incse_bus_vec(ind_income_sorted); 
% We sort self-employed by their total income but then we consider only
% their business income for tax evasion


%taxevasion_vec_sorted = taxevasion_vec(ind_income_sorted);
policyphidet_vec_sorted = policyphidet_vec(ind_income_sorted);

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
taxev_quint_inc_entre = zeros(5,1);
weight5         = zeros(5,1);

for i=1:5
    taxev_quint_inc_entre(i) = sum( mu_sorted(quint_inc{i}).*policyphidet_vec_sorted(quint_inc{i}).*pregov_incse_bus_vec_sorted(quint_inc{i}) )/...
        sum( mu_sorted(quint_inc{i}).*income_sorted(quint_inc{i}) );
    weight5(i) = sum( mu_sorted(quint_inc{i}).*income_sorted(quint_inc{i}) )/sum(mu_sorted.*income_sorted);
  
end

% check (law of total probabilities)
sum(weight5); % this should be exactly 1!!!
sum(taxev_quint_inc_entre.*weight5); % this should be exactly inc_gap_entre!!!

end %end function "income_gap_only_entre"
