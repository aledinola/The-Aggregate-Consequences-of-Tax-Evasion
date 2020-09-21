function [share_entre_wealth_quint,check] = share_entre_by_wealth(policy,distrib,wealth,share_entre)
%{
DESCRIPTION: This function computes the share of entre by wealth quintiles
-------------------------------------------------------------------------%
INPUTS
policy,distrib: structures with equilibrium objects
wealth:         array with "a", dim [nagrid_dist,negrid,ntgrid]
share_entre:    share of SE in the whole population

OUTPUTS: 
share_entre_wealth_quint: share entre by wealth quintiles

NOTES:
** Updated by Alessandro Di Nola on September 16, 2020.
-------------------------------------------------------------------------%
%}

%Unpack structure
occpoldet     = policy.occpoldet;
mu            = distrib.mu;
occpoldet_vec = occpoldet(:);
mu_vec        = mu(:);
wealth_vec    = wealth(:);

%Sort things by wealth:
[~, ind_wealth_sorted] = sort(wealth_vec);
mu_sorted              = mu_vec(ind_wealth_sorted);
occpoldet_vec_sorted   = occpoldet_vec(ind_wealth_sorted);

% Cumulative distribution
mu_sorted_norm            = mu_sorted / sum(mu_sorted);
mu_sorted_norm_cumulative = cumsum(mu_sorted_norm);


[~, ind_wealth_perc20] = min(abs(mu_sorted_norm_cumulative-0.20));
[~, ind_wealth_perc40] = min(abs(mu_sorted_norm_cumulative-0.40));
[~, ind_wealth_perc60] = min(abs(mu_sorted_norm_cumulative-0.60));
[~, ind_wealth_perc80] = min(abs(mu_sorted_norm_cumulative-0.80));


%Quintile Indices for total wealth:
quint_wealth{1}    =  1:ind_wealth_perc20;
quint_wealth{2}    =  ind_wealth_perc20+1:ind_wealth_perc40;
quint_wealth{3}    =  ind_wealth_perc40+1:ind_wealth_perc60;
quint_wealth{4}    =  ind_wealth_perc60+1:ind_wealth_perc80;
quint_wealth{5}    =  ind_wealth_perc80+1:length(mu_sorted);

share_entre_wealth_quint      = zeros(5,1);

for i=1:5
    share_entre_wealth_quint(i) = sum( mu_sorted(quint_wealth{i}).*occpoldet_vec_sorted(quint_wealth{i}) )...
        /sum(mu_sorted(quint_wealth{i}));
end

check = sum(share_entre_wealth_quint*0.2); % this should give exactly share_entre

if abs(check-share_entre)>1e-5
    warning('Please check file share_entre_by')
end

end %END FUNCTION "share_entre_by_wealth"

