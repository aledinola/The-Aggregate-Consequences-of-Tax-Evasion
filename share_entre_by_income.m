function [share_entre_inc_quint,check] = share_entre_by_income(policy,distrib,pregov_inc)

%{
DESCRIPTION: This function computes the share of entre by income quintiles
---------------------------------------------------------------------------
INPUTS:
policy,distrib: structures with equilibrium objects
pregov_inc: decision rule for income all, 
            considering occpol [nagrid_dist,negrid,ntgrid]
share_entre: share of SE in the whole population
---------------------------------------------------------------------------

OUTPUTS: 
share_entre_inc_quint:    share entre by income quintiles
check
---------------------------------------------------------------------------

NOTES: 
What is pregov_inc?
pregov_inc(:,j,t)=(1-occpoldet(i,j,t))*pregov_incw(i,j,t)+occpoldet(i,j,t)*pregov_incse(i,j,t);

** Updated by Alessandro Di Nola on September 16, 2020.
%}

%Vectorize stuff columnwise
% From 3-D [nagrid_dist,negrid,ntgrid] to 1-D [nagrid_dist*negrid*ntgrid,1]
pregov_inc_vec        = pregov_inc(:);
mu_vec                = distrib.mu(:);
occpoldet_vec         = policy.occpoldet(:);

%Sort things by income:
% Consider self-employed & workers (ALL)
[~, ind_income_sorted] = sort(pregov_inc_vec);
mu_sorted              = mu_vec(ind_income_sorted);
occpoldet_vec_sorted   = occpoldet_vec(ind_income_sorted);

% Cumulative distribution
mu_sorted_norm = mu_sorted / sum(mu_sorted);
mu_sorted_norm_cumulative = cumsum(mu_sorted_norm);


[~, ind_perc20] = min(abs(mu_sorted_norm_cumulative-0.20));
[~, ind_perc40] = min(abs(mu_sorted_norm_cumulative-0.40));
[~, ind_perc60] = min(abs(mu_sorted_norm_cumulative-0.60));
[~, ind_perc80] = min(abs(mu_sorted_norm_cumulative-0.80));


%Quintile Indices for total wealth:
quint_inc{1}    =  1:ind_perc20;
quint_inc{2}    =  ind_perc20+1:ind_perc40;
quint_inc{3}    =  ind_perc40+1:ind_perc60;
quint_inc{4}    =  ind_perc60+1:ind_perc80;
quint_inc{5}    =  ind_perc80+1:length(mu_sorted);

share_entre_inc_quint      = zeros(5,1);

for i=1:5
    share_entre_inc_quint(i) = sum( mu_sorted(quint_inc{i}).*occpoldet_vec_sorted(quint_inc{i}) )...
        /sum(mu_sorted(quint_inc{i}));    
end

check = sum(share_entre_inc_quint*0.2); % this should give exactly share_entre

end %END FUNCTION

