function [cev_total,cev_total_se,cev_total_work,cev_se,cev_work] = welfare_computations_ce1...
    (policy_bench,value_bench,distrib_bench,policy_reform,value_reform,distrib_reform,Parameters,Grids,flags,w0,r0)

%% HERE WE FOLLOW CE1 DEFINITION - micro based. 
%It is our preferred definition used for the paper

%{
INPUTS
- policy_bench,value_bench,distrib_bench: structures from benchmark economy
- policy_reform,value_reform,distrib_reform: structures from reform 
                                             (counterfactual) economy
- Parameters: structure with parameters
- Grids: structure with grids
- flags: structure with flags
- Parameters,Grids,flags are of course the *same* in both economies

OUTPUTS
- cev_total: CEV whole population, averaging CE1(x). See eq.18 memo
- cev_total_se: CEV conditional on j=self-employed
- cev_total_work: CEV conditional on j=work
- cev_se: CEV by deciles of wealth, for j=self-employed
- cev_work: CEV by deciles of wealth, for j=workers

NOTES:
- sigma is sigma1 in the paper, i=1,..,10 (wealth deciles) and j=SE,W
- Policy functions are defined on agrid_dist to be compatible with the
distribution, so they all have _det

EXTERNAL FUNCTIONS:
- g_operator.m
- quantili.m

%Updated by Alessandro Di Nola on September 16, 2020.
%}

%% Unpack objects

%Paramerters
nagrid_dist   = Parameters.nagrid_dist;
negrid        = Parameters.negrid;
ntgrid        = Parameters.ntgrid;
%beta          = Parameters.beta;
sigma         = Parameters.sigma;
% pflag         = Parameters.pflag;
% pn_1          = Parameters.pn_1;
% pn_2          = Parameters.pn_2;
% pn_3          = Parameters.pn_3;

%Grids
agrid        = Grids.agrid;
agrid_dist   = Grids.agrid_dist;
wealth_vec   = Grids.wealth_vec;


%Benchmark economy
occpol_bench     = policy_bench.occpoldet;
% cpolw_bench      = policy_bench.cpolwdet;
% cpolse1_bench    = policy_bench.cpolse1det;
% cpolse0_bench    = policy_bench.cpolse0det;
% kpol_bench_bench = policy_bench.policycapdet;
Vse_bench        = interp1(agrid,value_bench.Vse,agrid_dist);
Vw_bench         = interp1(agrid,value_bench.Vw,agrid_dist);
V_bench          = max(Vse_bench,Vw_bench);
mu_bench         = distrib_bench.mu;
mu_work_bench    = distrib_bench.mu_work;
mu_se_bench      = distrib_bench.mu_se;

%Reform economy
% occpol_reform    = policy_reform.occpoldet;
% cpolw_reform     = policy_reform.cpolwdet;
% cpolse1_reform   = policy_reform.cpolse1det;
% cpolse0_reform   = policy_reform.cpolse0det;

Vse_reform     = interp1(agrid,value_reform.Vse,agrid_dist);
Vw_reform      = interp1(agrid,value_reform.Vw,agrid_dist);
V_reform       = max(Vse_reform,Vw_reform);
% mu_reform      = distrib_reform.mu;
% mu_work_reform = distrib_reform.mu_work;
% mu_se_reform   = distrib_reform.mu_se;

% % cev totals, single number
% cev_total = 0;

%Utility from consumption (w/out labor part!)
UC = @(c) (c.^(1-sigma)/(1-sigma));

%% CE1 single number (all,workers,self-employed)

%First we need to compute the expected discounted utility from consumption
%in the benchmark economy (see eq.12 memo)
%METHOD G. This is the recursive method explained in Section 5.3.1 of
%my welfare memo: only policies in (B) econ are used

disp('Start iteration on G operator.. please be patient!')
[Georgi,fail_georgi] = g_operator(policy_bench,Parameters,Grids,flags,UC,w0,r0);
if fail_georgi==1; error('Smth wrong with G operator'); end
%G operator is G^(x) in the welfare memo

%Second, compute CE1 according to eq.17 memo, for each x
% CE1 is a 3-D array
CE1       = (((V_reform-V_bench)./Georgi)+1).^(1/(1-sigma))-1;

%Compute average as in eq.18 of memo
cev_total      = sum(sum(sum(CE1.*mu_bench)));
%Conditional average by occupational status
cev_total_work = sum(sum(sum(CE1.*mu_work_bench)));
cev_total_se   = sum(sum(sum(CE1.*mu_se_bench)));

% if Parameters.endo_ls_work==0 && Parameters.endo_ls_entre==0
%     %Standard formula for CE1 with CRRA w/out labor
%     cev_total = (EV_reform/EV_bench)^(1/(1-sigma))-1;
% else
%     %Formula for CE2 with separable u(c,l), see Bruggeman (2019)
%     cev_total = ( ((EV_reform-EV_bench)/E0) + 1)^(1/(1-sigma))-1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Deciles

% wealth deciles
q = 0.1:0.1:1;
assets_all_bench = quantili(wealth_vec, mu_bench(:),q); 
assets_all_bench(10) = max(agrid_dist);


%% SELF-EMPLOYED, CE1 by wealth deciles
% See eq.19 of my memo
%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ce1_se_numerator   = zeros(10,1);
ce1_se_denominator = zeros(10,1);

%First interval do separetely
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            % if you are in the first decile and a self-employed
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==1) 
                ce1_se_numerator(1) = ce1_se_numerator(1) + CE1(i,j,t)*mu_bench(i,j,t);
                ce1_se_denominator(1) = ce1_se_denominator(1) + mu_bench(i,j,t);
            end
        end
    end
end

%Other intervals
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            for id=2:10
                % if you are in decile id and a self-employed
                if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==1)
                    ce1_se_numerator(id) = ce1_se_numerator(id) + CE1(i,j,t)*mu_bench(i,j,t);
                    ce1_se_denominator(id) = ce1_se_denominator(id) + mu_bench(i,j,t);
                end
                
            end
        end
    end
end

%check_dist_se = XXX; %TBC

%Write in a vector %10 deciles for self-employed
cev_se = ce1_se_numerator./ce1_se_denominator;

%% WORKERS, CE1 by wealth deciles
% See eq.19 of my memo
%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ce1_work_numerator   = zeros(10,1);
ce1_work_denominator = zeros(10,1);

% First interval do separetely
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            % if you are in the first decile and a worker
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==0) 
                ce1_work_numerator(1) = ce1_work_numerator(1) + CE1(i,j,t)*mu_bench(i,j,t);
                ce1_work_denominator(1) = ce1_work_denominator(1) + mu_bench(i,j,t);
            end
        end
    end
end

%Other intervals
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            for id=2:10
                % if you are in decile id and a worker
                if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==0)
                    ce1_work_numerator(id) = ce1_work_numerator(id) + CE1(i,j,t)*mu_bench(i,j,t);
                    ce1_work_denominator(id) = ce1_work_denominator(id) + mu_bench(i,j,t);
                end
                
            end
        end
    end
end

%check_dist_work = YYY; %TBC

%Write in a vector %10 deciles for workers
cev_work = ce1_work_numerator./ce1_work_denominator;


% %%% check
% 
% value_work = sum(W_B_w.*sum_w_bench)
% value_entre = sum(W_B_se.*sum_se_bench)
% 
% chec = 0.135*value_entre+(1-0.135)*value_work

end %end function


