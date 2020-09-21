function [cev_total,cev_total_se,cev_total_work,cev_se,cev_w] = welfare_computations_ce2...
    (policy_bench,value_bench,distrib_bench,policy_reform,value_reform,distrib_reform,Parameters,Grids,flags,w0,r0)

%Updated by Ale on 15Jan2019

%% HERE WE FOLLOW CE2 DEFINITION (Lucas aggregate number, see my memo)

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
sigma is sigma1 in the paper, i=1,..,10 (wealth deciles) and j=SE,W

EXTERNAL FUNCTIONS:
- g_operator.m
- quantili.m

%}

%% Unpack objects

%Paramerters
nagrid_dist   = Parameters.nagrid_dist;
negrid        = Parameters.negrid;
ntgrid        = Parameters.ntgrid;
beta          = Parameters.beta;
sigma         = Parameters.sigma;
pflag         = Parameters.pflag;
pn_1          = Parameters.pn_1;
pn_2          = Parameters.pn_2;
pn_3          = Parameters.pn_3;

%Grids
agrid        = Grids.agrid;
agrid_dist   = Grids.agrid_dist;
wealth_vec_bench = Grids.wealth_vec;


%Benchmark economy
occpol_bench  = policy_bench.occpoldet;
cpolw         = policy_bench.cpolwdet;
cpolse1       = policy_bench.cpolse1det;
cpolse0       = policy_bench.cpolse0det;
kpol_bench    = policy_bench.policycapdet;
Vse_bench     = interp1(agrid,value_bench.Vse,agrid_dist);
Vw_bench      = interp1(agrid,value_bench.Vw,agrid_dist);
V_bench       = max(Vse_bench,Vw_bench);
mu_bench      = distrib_bench.mu;
mu_work_bench = distrib_bench.mu_work;
mu_se_bench   = distrib_bench.mu_se;

%Reform economy
occpol_reform  = policy_reform.occpoldet;
Vse_reform     = interp1(agrid,value_reform.Vse,agrid_dist);
Vw_reform      = interp1(agrid,value_reform.Vw,agrid_dist);
V_reform       = max(Vse_reform,Vw_reform);
mu_reform      = distrib_reform.mu;
mu_work_reform = distrib_reform.mu_work;
mu_se_reform   = distrib_reform.mu_se;


% % cev totals, single number
% cev_total = 0;

% consumption equivalent variation by deciles of wealth, vectors (10,1)
cev_se = zeros(10,1); 
cev_w  = zeros(10,1);

%Utility from consumption
UC = @(c) (c.^(1-sigma)/(1-sigma));


% Calculations

%% CEV single number
EV_bench = 0;

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            EV_bench = EV_bench + (occpol_bench(i,j,t)*Vse_bench(i,j,t) + (1-occpol_bench(i,j,t))*Vw_bench(i,j,t))*mu_bench(i,j,t);
        end
    end
end

EV_reform = 0;

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            EV_reform = EV_reform + (occpol_reform(i,j,t)*Vse_reform(i,j,t) + (1-occpol_reform(i,j,t))*Vw_reform(i,j,t))*mu_reform(i,j,t);
        end
    end
end

%We call E0 the term at the denominator of Bruggeman(2019, eq.C3). See also
%our memo on welfare analysis.

E0temp = 0;
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            pk = prob_audit([],[],[],kpol_bench(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
            E0temp = E0temp + ((1-occpol_bench(i,j,t))*UC(cpolw(i,j,t))+...
                occpol_bench(i,j,t)*(pk*UC(cpolse1(i,j,t))+(1-pk)*UC(cpolse0(i,j,t))))*mu_bench(i,j,t); 
				% we use consumption from the benchmark
        end
    end
end

E0 = E0temp/(1-beta); %This is the denominator in eq.(2) of Ale's memo

if Parameters.endo_ls_work==0 && Parameters.endo_ls_entre==0
    cev_total = (EV_reform/EV_bench)^(1/(1-sigma))-1;
else
    cev_total = ( ((EV_reform-EV_bench)/E0) + 1)^(1/(1-sigma))-1;
end
fprintf('CEV total = %f \n', cev_total);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EV_bench_work = 0;

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            EV_bench_work = EV_bench_work + Vw_bench(i,j,t).*mu_work_bench(i,j,t);
        end
    end
end

EV_reform_work = 0;

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            EV_reform_work = EV_reform_work + Vw_reform(i,j,t).*mu_work_reform(i,j,t);
        end
    end
end

cev_total_work = (EV_reform_work/EV_bench_work)^(1/(1-sigma))-1;
disp(cev_total_work)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EV_bench_se = 0;

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            EV_bench_se = EV_bench_se + Vse_bench(i,j,t).*mu_se_bench(i,j,t);
        end
    end
end

EV_reform_se = 0;

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            EV_reform_se = EV_reform_se + Vse_reform(i,j,t).*mu_se_reform(i,j,t);
        end
    end
end

cev_total_se = (EV_reform_se/EV_bench_se)^(1/(1-sigma))-1;
disp(cev_total_se)

%% Deciles

% wealth deciles
q = 0.1:0.1:0.9;
assets_all_bench = quantili(wealth_vec_bench, mu_bench(:),q); 
assets_all_bench = [assets_all_bench;
                    agrid_dist(end)];


%% Compute W_B by deciles for SE as in Ale's memo

W_B_se = zeros(10,1);
ind_se = zeros(nagrid_dist, negrid, ntgrid);

% First interval do separetely
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==1) % if you are in the first decile and a self-employed
                ind_se(i,j,t) = 1;
                pk = prob_audit([],[],[],kpol_bench(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
                W_B_se(1) = W_B_se(1) + ind_se(i,j,t)*((1-occpol_bench(i,j,t))*UC(cpolw(i,j,t))+...
                    occpol_bench(i,j,t)*(pk*UC(cpolse1(i,j,t))+(1-pk)*UC(cpolse0(i,j,t)))) *mu_bench(i,j,t);
            end
            
        end
    end
end


for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            for id=2:10
                
                if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==1) % if you are in the first decile and a self-employed
                    ind_se(i,j,t) = 1;
                    pk = prob_audit([],[],[],kpol_bench(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
                    W_B_se(id) = W_B_se(id) + ind_se(i,j,t)*((1-occpol_bench(i,j,t))*UC(cpolw(i,j,t))+...
                        occpol_bench(i,j,t)*(pk*UC(cpolse1(i,j,t))+(1-pk)*UC(cpolse0(i,j,t)))) *mu_bench(i,j,t);
                end
                
            end
        end
    end
end

W_B_se = W_B_se/(1-beta);

%% Compute W_B by deciles for Workers as in Ale's memo

W_B_w  = zeros(10,1);
ind_w  = zeros(nagrid_dist, negrid, ntgrid);

% First interval do separetely
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==0) % if you are in the first decile and a self-employed
                ind_w(i,j,t) = 1;
                pk = prob_audit([],[],[],kpol_bench(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
                W_B_w(1) = W_B_w(1) + ind_w(i,j,t)*((1-occpol_bench(i,j,t))*UC(cpolw(i,j,t))+...
                    occpol_bench(i,j,t)*(pk*UC(cpolse1(i,j,t))+(1-pk)*UC(cpolse0(i,j,t)))) *mu_bench(i,j,t);
            end
            
        end
    end
end


for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            for id=2:10
                
                if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==0) % if you are in the first decile and a self-employed
                    ind_w(i,j,t) = 1;
                    pk = prob_audit([],[],[],kpol_bench(i,j,t),[],pflag,pn_1,pn_2,pn_3,[],[],[]);
                    W_B_w(id) = W_B_w(id) + ind_w(i,j,t)*((1-occpol_bench(i,j,t))*UC(cpolw(i,j,t))+...
                        occpol_bench(i,j,t)*(pk*UC(cpolse1(i,j,t))+(1-pk)*UC(cpolse0(i,j,t)))) *mu_bench(i,j,t);
                end
                
            end
        end
    end
end

W_B_w = W_B_w/(1-beta);


%%
%%%%%%%%%%%%%%%%%%%%%% SELF-EMPLOYED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

agrid_dist_states = repmat(agrid_dist, 1 , negrid, ntgrid);

welfare_se_bench  = zeros(10,1);
sum_se_bench      = zeros(10,1);
welfare_se_reform = zeros(10,1);
sum_se_reform     = zeros(10,1);

ind_se            = zeros(nagrid_dist, negrid, ntgrid);
 

% First interval do separetely
% welfare_se_bench is the numerator of EV^B
% sum_se_bench is the denominator of EV^B
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==1) % if you are in the first decile and a self-employed
                ind_se(i,j,t) = 1;
                welfare_se_bench(1) = welfare_se_bench(1) + ind_se(i,j,t).*(occpol_bench(i,j,t)*Vse_bench(i,j,t)...
                    + (1-occpol_bench(i,j,t))*Vw_bench(i,j,t))*mu_bench(i,j,t);
                sum_se_bench(1)     = sum_se_bench(1) + ind_se(i,j,t).*mu_bench(i,j,t);
            end
            
        end
    end
end


for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            for id=2:10
                
            if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==1) % if you are in the first decile and a self-employed
                ind_se(i,j,t) = 1;
                welfare_se_bench(id) = welfare_se_bench(id) + ind_se(i,j,t).*(occpol_bench(i,j,t)*Vse_bench(i,j,t)...
                    + (1-occpol_bench(i,j,t))*Vw_bench(i,j,t))*mu_bench(i,j,t);
                sum_se_bench(id)     = sum_se_bench(id) + ind_se(i,j,t).*mu_bench(i,j,t);
            end
            
            end
        end
    end
end

check_dist_se = sum(sum_se_bench) %share of self-employed

% first interval do separetely
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==1) % if you are in the first decile and a self-employed
                ind_se(i,j,t) = 1;
                welfare_se_reform(1) = welfare_se_reform(1) + ind_se(i,j,t).*(occpol_reform(i,j,t)*Vse_reform(i,j,t)...
                    + (1-occpol_reform(i,j,t))*Vw_reform(i,j,t))*mu_reform(i,j,t);
                sum_se_reform(1)     = sum_se_reform(1) + ind_se(i,j,t).*mu_reform(i,j,t);
            end
            
        end
    end
end

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            for id=2:10
            if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==1) % if you are in the first decile and a self-employed
                ind_se(i,j,t) = 1;
                welfare_se_reform(id) = welfare_se_reform(id) + ind_se(i,j,t).*(occpol_reform(i,j,t)*Vse_reform(i,j,t)...
                    + (1-occpol_reform(i,j,t))*Vw_reform(i,j,t))*mu_reform(i,j,t);
                sum_se_reform(id)     = sum_se_reform(id) + ind_se(i,j,t).*mu_reform(i,j,t);
            end
            end
        end
    end
end

check_dist_se1 = sum(sum_se_reform) %share of self-employed

% Write in a vector %10 deciles for self-employed
if Parameters.endo_ls_work==0 && Parameters.endo_ls_entre==0
    %Old formula with exogenous labor
    for i = 1:10 % 10 deciles
        cev_se(i) = ( (welfare_se_reform(i)/sum_se_reform(i) ) / (welfare_se_bench(i)/sum_se_bench(i)))^(1/(1-sigma))-1;
    end
else
    %New formula with endogenous labor
    for i = 1:10 % 10 deciles
        EV_C_i1 = welfare_se_reform(i)/sum_se_reform(i); %for each i, given j=1 (i.e. SE)
        EV_B_i1 = welfare_se_bench(i)/sum_se_bench(i);
        %Compute W_B !!!!!!!!!!!!!!!!!!!!!!!!!!
        cev_se(i) = (((EV_C_i1-EV_B_i1)/(W_B_se(i)/sum_se_bench(i))  ) + 1)^(1/(1-sigma))-1;
    end
end


%%
%%%%%%%%%%%%%%%%% WORKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


welfare_w_bench  = zeros(10,1);
sum_w_bench      = zeros(10,1);
welfare_w_reform = zeros(10,1);
sum_w_reform     = zeros(10,1);

ind_w            = zeros(nagrid_dist, negrid, ntgrid);
 

% First interval do separetely
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==0) % if you are in the first decile and a self-employed
                ind_w(i,j,t) = 1;
                welfare_w_bench(1) = welfare_w_bench(1) + ind_w(i,j,t).*(occpol_bench(i,j,t)*Vse_bench(i,j,t)...
                    + (1-occpol_bench(i,j,t))*Vw_bench(i,j,t))*mu_bench(i,j,t);
                sum_w_bench(1)     = sum_w_bench(1) + ind_w(i,j,t).*mu_bench(i,j,t);
            end
            
        end
    end
end


for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            for id=2:10
                
            if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==0) % if you are in the first decile and a worker
                ind_w(i,j,t) = 1;
                welfare_w_bench(id) = welfare_w_bench(id) + ind_w(i,j,t).*(occpol_bench(i,j,t)*Vse_bench(i,j,t)...
                    + (1-occpol_bench(i,j,t))*Vw_bench(i,j,t))*mu_bench(i,j,t);
                sum_w_bench(id)     = sum_w_bench(id) + ind_w(i,j,t).*mu_bench(i,j,t);
            end
            
            end
        end
    end
end

check_dist_w = sum(sum_w_bench) %share of self-employed

% first interval do separetely
for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            
            if (agrid_dist(i) <= assets_all_bench(1) && occpol_bench(i,j,t)==0) % if you are in the first decile and a self-employed
                ind_w(i,j,t) = 1;
                welfare_w_reform(1) = welfare_w_reform(1) + ind_w(i,j,t).*(occpol_reform(i,j,t)*Vse_reform(i,j,t)...
                    + (1-occpol_reform(i,j,t))*Vw_reform(i,j,t))*mu_reform(i,j,t);
                sum_w_reform(1)     = sum_w_reform(1) + ind_w(i,j,t).*mu_reform(i,j,t);
            end
            
        end
    end
end

for i=1:nagrid_dist
    for j=1:negrid
        for t=1:ntgrid
            for id=2:10
            if (agrid_dist(i) > assets_all_bench(id-1) && agrid_dist(i) <= assets_all_bench(id) && occpol_bench(i,j,t)==0) % if you are in the first decile and a self-employed
                ind_w(i,j,t) = 1;
                welfare_w_reform(id) = welfare_w_reform(id) + ind_w(i,j,t).*(occpol_reform(i,j,t)*Vse_reform(i,j,t)...
                    + (1-occpol_reform(i,j,t))*Vw_reform(i,j,t))*mu_reform(i,j,t);
                sum_w_reform(id)     = sum_w_reform(id) + ind_w(i,j,t).*mu_reform(i,j,t);
            end
            end
        end
    end
end

check_dist_w1 = sum(sum_w_reform) %share of self-employed

% Write in a vector %10 deciles for workers
if Parameters.endo_ls_work==0 && Parameters.endo_ls_entre==0
    %Old formula with exogenous labor
    for i = 1:10 % 10 deciles
        cev_w(i) = ( (welfare_w_reform(i)/sum_w_reform(i) ) / (welfare_w_bench(i)/sum_w_bench(i)))^(1/(1-sigma))-1;
    end
else
    %New formula with endogenous labor
    for i = 1:10 % 10 deciles
        EV_C_i0 = welfare_w_reform(i)/sum_w_reform(i); %for each i, given j=0 (i.e. W)
        EV_B_i0 = welfare_w_bench(i)/sum_w_bench(i);
        %Compute W_B !!!!!!!!!!!!!!!!!!!!!!!!!!
        cev_w(i) = (((EV_C_i0-EV_B_i0)/(W_B_w(i)/sum_w_bench(i))  ) + 1)^(1/(1-sigma))-1;
    end
end


% %%% check
% 
% value_work = sum(W_B_w.*sum_w_bench)
% value_entre = sum(W_B_se.*sum_se_bench)
% 
% chec = 0.135*value_entre+(1-0.135)*value_work

end %end function


