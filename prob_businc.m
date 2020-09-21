function [audit_quint_inc_entre] = prob_businc(pregov_incse_bus_vec,mu_se_vec,policycapdet,pn_1,pn_2,pn_3,pflag)

%{
Aim: compute expected prob. of auditing by quintiles of business income
INPUTS
  pregov_incse_bus_vec: decision rule for business income, dim: [nagrid_dist*negrid*ntgrid,1]
  mu_se_vec: distribution SE, dim: [nagrid_dist*negrid*ntgrid,1]
  policycapdet: decision rule for capital, dim: [nagrid_dist*negrid*ntgrid,1]
OUTPUTS
  audit_quint_inc_entre: 5 elements vector with E(pk|quintile)

EXTERNAL
This function also calls the function 
+ pk = prob_audit(phi,profit,theta,k,n,pflag,pn_1,pn_2,pn_3,le,gamma,vi)

NOTES:
  Updated by Alessandro Di Nola on September 16, 2020.
%}


%% Sort things by income:
% pregov_incse_bus
% Consider ONLY self-employed. SURE?
[~, ind_income_sorted_se] = sort(pregov_incse_bus_vec);
mu_sorted = mu_se_vec(ind_income_sorted_se);
%pregov_incse_bus_vec_sorted = pregov_incse_bus_vec(ind_income_sorted_se); 

%policycapdet_vec = policycapdet(:);
%policycapdet_vec_sorted = policycapdet_vec(ind_income_sorted_se);

policycapdet_vec = policycapdet(:);
policycapdet_vec_sorted = policycapdet_vec(ind_income_sorted_se);


% Cumulative distribution
mu_sorted_norm = mu_sorted / sum(mu_sorted);
mu_sorted_norm_cumulative = cumsum(mu_sorted_norm);

% Deciles
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

% Prob of auditing by quintile of self-employed income
audit_quint_inc_entre  = zeros(5,1);
weight5                = zeros(5,1);

for i=1:5
    weight5(i) = sum( mu_sorted(quint_inc{i}));
    %pk = prob_audit(phi,profit,theta,k,n,pflag,pn_1,pn_2,pn_3,le,gamma,vi)
    pk_fun = prob_audit([],[],[],policycapdet_vec_sorted(quint_inc{i}),[],pflag,pn_1,pn_2,pn_3,[],[],[] );
    audit_quint_inc_entre(i) = sum( pk_fun.*mu_sorted(quint_inc{i}))/weight5(i);
end

%disp(weight5)

end %END function
