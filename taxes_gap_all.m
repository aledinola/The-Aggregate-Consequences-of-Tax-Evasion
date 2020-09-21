function [tax_gap_all,taxev_dec_inc1,taxev_quint_inc1,tax_gap_ave] = taxes_gap_all(Parameters,policy,distrib,grids,r0,...
    pregov_inc_vec,pregov_incse_bus,pregov_incw,tax_gap_single,unpaid_taxes_vec,true_tax_vec)

%{
 Here we compute the tax gap, total and over the distribution of income.
 For the definition of the tax gap, see Johns and Slemrod (2011). We follow 
 this article quite closely. See Online Appendix

 INPUTS:
  TBA
    
 OUTPUTS:
 tax_gap_all:      overall tax gap
 taxev_dec_inc1:   tax gap by income deciles
 taxev_quint_inc1: tax gap by income quintiles

 NOTES:
 %% Updated by Alessandro Di Nola on September 16, 2020. Folder: RED.
%}
    
% Validate inputs
if ~isstruct(Parameters)
    error('Input "Parameters" must be a structure!')
end
if ~isstruct(policy)
    error('Input "policy" must be a structure!')
end
if ~isstruct(distrib)
    error('Input "distrib" must be a structure!')
end
if ~isstruct(grids)
    error('Input "grids" must be a structure!')
end

%Unpack
occpoldet    = policy.occpoldet;
occpoldet_vec = occpoldet(:);
policyphidet = policy.policyphidet;
agrid_dist   = grids.agrid_dist;
mu           = distrib.mu;
mu_vec       = mu(:);
tax_gap_single_vec = tax_gap_single(:);

%% Tax Gap
% Share of unpaid/misreported taxes out of total taxes (paid+unpaid) 
unpaid_tax_entre = 0;
total_tax_entre = 0;
total_tax_work = 0;
for i=1:Parameters.nagrid_dist
    for j=1:Parameters.negrid
        for t=1:Parameters.ntgrid
            % true_tax
            tobepaid = tax_entre(pregov_incse_bus(i,j,t)+r0*agrid_dist(i),Parameters);
            paid = tax_entre((1-policyphidet(i,j,t))*pregov_incse_bus(i,j,t)+r0*agrid_dist(i),Parameters);
            unpaid_tax_entre = unpaid_tax_entre + occpoldet(i,j,t)*(tobepaid-paid)*mu(i,j,t);
            total_tax_entre = total_tax_entre + occpoldet(i,j,t)*tax_entre(pregov_incse_bus(i,j,t)+r0*agrid_dist(i),Parameters)*mu(i,j,t);
            total_tax_work = total_tax_work + (1-occpoldet(i,j,t))*tax_work(pregov_incw(i,j,t),Parameters)*mu(i,j,t);
        end 
    end
end
% Slemrod's definition: tax gap = ratio of total unpaid/total true liab
tax_gap_all = unpaid_tax_entre/(total_tax_entre+total_tax_work);

% NOT USED Second definition: tax gap is average of individual agents' tax gaps
tax_gap_ave = sum(tax_gap_single(:).*mu_vec); % including also workers

clear tobepaid paid unpaid_tax_entre total_tax_entre total_tax_work

%% Sort things by income:
% [true_inc_vec_sorted, ind_income_sorted] = sort(pregov_incse_vec);
% mu_sorted = mu_se_vec(ind_income_sorted);

[~, ind_income_sorted] = sort(pregov_inc_vec);
mu_sorted = mu_vec(ind_income_sorted);

%taxevasion_vec_sorted = taxevasion_vec(ind_income_sorted); hidden income
unpaid_taxes_vec_sorted   = unpaid_taxes_vec(ind_income_sorted); 
true_tax_vec_sorted       = true_tax_vec(ind_income_sorted);
tax_gap_single_vec_sorted = tax_gap_single_vec(ind_income_sorted);
%policyphidet_vec_sorted   = policyphidet_vec(ind_income_sorted);
occpoldet_vec_sorted      = occpoldet_vec(ind_income_sorted);

% Cumulative distribution
mu_sorted_norm = mu_sorted / sum(mu_sorted);
mu_sorted_norm_cumulative = cumsum(mu_sorted_norm);


%% Quintiles

%[dev_perc10, ind_inc_perc10] = min(abs(mu_sorted_norm_cumulative-0.10));
[~, ind_inc_perc20] = min(abs(mu_sorted_norm_cumulative-0.20));
%[dev_perc30, ind_inc_perc30] = min(abs(mu_sorted_norm_cumulative-0.30));
[~, ind_inc_perc40] = min(abs(mu_sorted_norm_cumulative-0.40));
%[dev_perc50, ind_inc_perc50] = min(abs(mu_sorted_norm_cumulative-0.50));
[~, ind_inc_perc60] = min(abs(mu_sorted_norm_cumulative-0.60));
%[dev_perc70, ind_inc_perc70] = min(abs(mu_sorted_norm_cumulative-0.70));
[~, ind_inc_perc80] = min(abs(mu_sorted_norm_cumulative-0.80));
%[dev_perc90, ind_inc_perc90] = min(abs(mu_sorted_norm_cumulative-0.90));

%Quintile Indices for income:
quint_inc{1} =  1:ind_inc_perc20;
quint_inc{2} =  ind_inc_perc20+1:ind_inc_perc40;
quint_inc{3} =  ind_inc_perc40+1:ind_inc_perc60;
quint_inc{4} =  ind_inc_perc60+1:ind_inc_perc80;
quint_inc{5} =  ind_inc_perc80+1:length(mu_sorted);

% Tax gap by quintile of total income
taxev_quint_inc1  = zeros(5,1);
taxev_quint_inc   = zeros(5,1);
share_entre_quint = zeros(5,1);
true_tax_quint    = zeros(5,1); % true tax quintile i / total true taxes

for i=1:5
  % WE HAVE TO USE METHOD 1!!!
  % Method 1: ratio of total unpaid taxes in decile i / total true taxes
  taxev_quint_inc1(i) = sum( mu_sorted(quint_inc{i}).*unpaid_taxes_vec_sorted(quint_inc{i}) )/sum( mu_sorted(quint_inc{i}).*true_tax_vec_sorted(quint_inc{i}) );
  % Method 2: Conditional average of tax gap (unpaid/due) for all agents in quintile i 
  taxev_quint_inc(i) = sum( mu_sorted(quint_inc{i}).*tax_gap_single_vec_sorted(quint_inc{i}) )/sum(mu_sorted(quint_inc{i}));
  share_entre_quint(i) = sum( mu_sorted(quint_inc{i}).*occpoldet_vec_sorted(quint_inc{i}) )/ sum(mu_sorted(quint_inc{i})) ;
  true_tax_quint(i) = sum( mu_sorted(quint_inc{i}).*true_tax_vec_sorted(quint_inc{i}) )/sum( mu_sorted.*true_tax_vec_sorted );

end

% check (law of total probabilities)
sum(true_tax_quint); % this has to be 1
sum(share_entre_quint*0.2); % this should give exactly (share_entre|occ=entre)
sum(taxev_quint_inc1.*true_tax_quint); % % this should be exactly tax_gap_all
sum(taxev_quint_inc*0.2); % this should be exactly tax_gap_ave

end %END FUNCTION FILE





