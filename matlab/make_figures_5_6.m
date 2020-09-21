%% This script generates Figures (5-6), The Impact of the Fine on Aggregate Outcomes and Welfare
%{
External files called by this script:
    + plots_no_redistrib
    + plots_lump_sum
    + plots_tax_cut
    + plots_tax_cut_se
Notes: Updated by Alesandro Di Nola on September 17, 2020.
%}

clear;clc;close all

%% Specify useful paths

ResultsDir = 'results\mat\'; %folder where .mat files are stored
SaveDir = 'results\Figures\'; %Specify here where you want to save the figures
disp(['Saving figures in subfolder: ' SaveDir])

%% Replication flags
% To replicate Figure 5 (The Impact of the Fine on Aggregate Outcomes)
% No fiscal redistribution
% Top panel:    set method='PE' & experiment='no_redistrib'
% Bottom panel: set method='GE' & experiment='no_redistrib'
%-------------------------------------------------------------------------%
% To replicate Figure 6 (The Impact of the Fine on Welfare). Analysis in
% general equilibrium with (3 ways of) balancing the government budget.
% Panels (a)-(b): set method='GE' & experiment='lump_sum'
% Panels (c)-(d): set method='GE' & experiment='tax_cut'
% Panels (e)-(f): set method='GE' & experiment='tax_cut_se'
%-------------------------------------------------------------------------%

method        = 'PE'; %options: 'PE' or 'GE'

experiment    = 'no_redistrib';      %Options: 
                                 %'no_redistrib'
                                 %'lump_sum'
                                 %'tax_cut'
                                 %'tax_cut_se'

interpolation = 'YES'; %Options are: 'YES','NO'
axis_tight    = 0; %Options are: 0 or 1. If 1, activate "axis tight" in plots


%% Options for plotting

ldw     = 2; %line width
fontw   = 13; %font size for xlabel, legend, etc.
fileExt = '-depsc'; % '-depsc', File extension: png or eps
ts      = 0; %0/1 display titles

%% Load results

if strcmp(experiment,'no_redistrib')
    
    if strcmp(method,'PE')
        disp('FINE EXPERIMENT: Load PARTIAL equilibrium results, no redistr.')
        load([ResultsDir 'results_fine_pe_no_red_0406.mat'])
        load([ResultsDir 'fine_cev_PE_no_redistrib.mat'])
    elseif strcmp(method,'GE')
        disp('FINE EXPERIMENT: Load GENERAL equilibrium results, no redistr.')
        load([ResultsDir 'results_fine_ge_no_red_0406.mat'])
        load([ResultsDir 'fine_cev_GE_no_redistrib.mat'])
    else
        error('Selected method does not exist!')
    end
    
elseif strcmp(experiment,'lump_sum')
    if strcmp(method,'GE')
        disp('FINE EXPERIMENT: Load GENERAL equilibrium results,lump-sum')
        load([ResultsDir 'results_fine_ge_lump_sum_0406.mat'])
        load([ResultsDir 'fine_cev_GE_lump_sum.mat'])
    else
        error('Selected method does not exist!')
    end
    
elseif strcmp(experiment,'tax_cut')
    
    if strcmp(method,'GE')
        disp('FINE EXPERIMENT: Load GENERAL equilibrium results, tax_cut')
        load([ResultsDir 'results_fine_ge_tax_cut_0406.mat'])
        load([ResultsDir 'fine_cev_GE_tax_cut.mat'])
    else
        error('Selected method does not exist!')
    end
elseif strcmp(experiment,'tax_cut_se')
    if strcmp(method,'GE')
        disp('FINE EXPERIMENT: Load GENERAL equilibrium results, tax_cut_se')
        load([ResultsDir 'results_fine_ge_tax_cut_se_0406.mat'])
        load([ResultsDir 'fine_cev_GE_tax_cut_se.mat'])
    else
        error('Selected method does not exist!')
    end
    
else
    error('Selected experiment does not exist!')
end

%% Some checks

for s_c =1:length(s_vec)
    if ResultsFine(s_c).fail>0
        fprintf('Something went wrong for s = %f \n',s_vec(s_c));
        fprintf('FAIL CHECK = %d \n',ResultsFine(s_c).fail);
    end
end

%Set max value s to plot
s_max = 10;
indmax = find(s_vec<=s_max, 1, 'last' );
if isempty(indmax)
    error('s_max for plotting is out of range!')
end

%% Prepare data
%Income gap all as a function of s
inc_gap_vec = nan(length(s_vec),1);
for s_c =1:length(s_vec)
    inc_gap_vec(s_c) = ResultsFine(s_c).ModelResults.inc_gap_all_p;
end

% Tax revenues: taxes_total,taxes_e,taxes_w
taxes_all_vec = nan(length(s_vec),1);
taxes_e_vec   = nan(length(s_vec),1);
taxes_w_vec   = nan(length(s_vec),1);
for s_c =1:length(s_vec)
    taxes_all_vec(s_c) = ResultsFine(s_c).ModelResults.taxes_total;
    taxes_e_vec(s_c)   = ResultsFine(s_c).ModelResults.taxes_e;
    taxes_w_vec(s_c)   = ResultsFine(s_c).ModelResults.taxes_w;
end

% share_entre_p,cond_mean_theta
share_entre_p_vec   = nan(length(s_vec),1);
cond_mean_theta_vec = nan(length(s_vec),1);
for s_c =1:length(s_vec)
    share_entre_p_vec(s_c)   = ResultsFine(s_c).ModelResults.share_entre_p;
    cond_mean_theta_vec(s_c) = ResultsFine(s_c).ModelResults.cond_mean_theta;
end

% capital,output (both total)
capital_vec   = nan(length(s_vec),1); %k_entre+k_corp
output_vec    = nan(length(s_vec),1); %output_corp + output_se
labor_vec     = nan(length(s_vec),1); %n_entre+n_corp
for s_c =1:length(s_vec)
    capital_vec(s_c)= ResultsFine(s_c).ModelResults.capital;
    output_vec(s_c) = ResultsFine(s_c).ModelResults.output;
    labor_vec(s_c) = ResultsFine(s_c).ModelResults.n_entre+...
        ResultsFine(s_c).ModelResults.n_corp;
end

%Welfare (as integral of V(.)). (We compute the CEVs in another file).
welfare_all_vec   = nan(length(s_vec),1); 
welfare_entre_vec = nan(length(s_vec),1); 
welfare_work_vec  = nan(length(s_vec),1); 

for s_c =1:length(s_vec)
    welfare_all_vec(s_c)   = ResultsFine(s_c).ModelResults.welfare_all;
    welfare_entre_vec(s_c) = ResultsFine(s_c).ModelResults.welfare_entre;
    welfare_work_vec(s_c)  = ResultsFine(s_c).ModelResults.welfare_work;
end

%ave_k_entre,ave_n_entre
ave_k_entre_vec    = nan(length(s_vec),1); 
ave_n_entre_vec    = nan(length(s_vec),1); 
small_firms_vec    = nan(length(s_vec),1); %maybe plot NON-norm share of firms
large_firms_vec    = nan(length(s_vec),1); %maybe plot NON-norm share of firms

for s_c =1:length(s_vec)
    ave_k_entre_vec(s_c)= ResultsFine(s_c).ModelResults.ave_k_entre;
    ave_n_entre_vec(s_c) = ResultsFine(s_c).ModelResults.ave_n_entre;
    small_firms_vec(s_c) = ResultsFine(s_c).ModelResults.cond_firm_size_dist(1);
    large_firms_vec(s_c) = sum(ResultsFine(s_c).ModelResults.cond_firm_size_dist(2:end));
end


%Data for Figure on (r,w) only GE
if strcmp(method,'GE')
    r_vec   = nan(length(s_vec),1);
    w_vec    = nan(length(s_vec),1);
    ED_vec    = nan(length(s_vec),1);
    
    for s_c =1:length(s_vec)
        r_vec(s_c)= ResultsFine(s_c).ModelResults.r0;
        w_vec(s_c)= ResultsFine(s_c).ModelResults.w0;
        ED_vec(s_c) = ResultsFine(s_c).ED;
    end
end

%% Normalization wrt benchmark: divide by value at 1.75 and *100

load([ResultsDir 'taxevasion_ge.mat'],'Parameters')
indx = find(s_vec==Parameters.s);
if isempty(indx)
    error('s=1.75 (benchmark) not found in vector s_vec!')
end
clear Parameters

if strcmp(method,'GE')
    r_vec_norm = (r_vec/r_vec(indx))*100;
    w_vec_norm = (w_vec/w_vec(indx))*100;
end

% taxes_all_vec_norm = (taxes_all_vec/taxes_all_vec(indx))*100;
% taxes_e_vec_norm   = (taxes_e_vec/taxes_e_vec(indx))*100;
% taxes_w_vec_norm   = (taxes_w_vec/taxes_w_vec(indx))*100;

strlist = {'taxes_all_vec','taxes_e_vec','taxes_w_vec','share_entre_p_vec',...
    'cond_mean_theta_vec','capital_vec','output_vec','labor_vec','ave_k_entre_vec',...
    'ave_n_entre_vec','small_firms_vec','large_firms_vec'};

for strct = 1:numel(strlist)
   eval([char(strlist(strct)) '_norm =' char(strlist(strct)) '/' char(strlist(strct)) '(' num2str(indx) ')*100']);
end

welfare_all_vec_norm = (((welfare_all_vec)-welfare_all_vec(indx))/abs(welfare_all_vec(indx)))*100+100;
welfare_entre_vec_norm = (((welfare_entre_vec)-welfare_entre_vec(indx))/abs(welfare_entre_vec(indx)))*100+100;
welfare_work_vec_norm = (((welfare_work_vec)-welfare_work_vec(indx))/abs(welfare_work_vec(indx)))*100+100;

%% Data for CEV Figure
% s_vec, cev_total_p, cev_total_se_p, cev_total_w_p are VECTORS
% indx is s.t. s_vec(indx)=1.75

cev_total_p(indx)    = 0.0;
cev_total_se_p(indx) = 0.0;
cev_total_w_p(indx)  = 0.0;


%% Smooth out results

smooth = 0.3; %Smoothing parameter for HP filter
if strcmp(method,'PE') && strcmp(experiment,'tax_cut')
    knots = [1 2 3 4 6 7 8 9 10 11]';
elseif strcmp(method,'GE') && strcmp(experiment,'tax_cut')
    knots = [1 3 4 5 6 7 9 10 11]';
elseif strcmp(method,'PE') && strcmp(experiment,'tax_cut_se')
    knots = [1 2 3 4 5 6 8 9 10 11]';
elseif strcmp(method,'GE') && strcmp(experiment,'tax_cut_se')
    knots = [1 3 4 6 7 8 9 10 11]';
elseif strcmp(method,'PE') && strcmp(experiment,'lump_sum')
    knots = [1 3 4 6 9 10 11]';
elseif strcmp(method,'GE') && strcmp(experiment,'lump_sum')
    knots = [1 3 4 6 9 10 11]';
else
    knots = (1:length(s_vec))';
end
%knots for interpolation 
interp_method = 'pchip'; %'spline','pchip','makima'

if strcmp(interpolation,'YES')
     inc_gap_vec                  = hpfilter(inc_gap_vec,smooth);
    inc_gap_vec              = interp1(s_vec(knots),inc_gap_vec(knots),s_vec,interp_method,'extrap');
    taxes_all_vec_norm       = interp1(s_vec(knots),taxes_all_vec_norm(knots),s_vec,interp_method,'extrap');
    taxes_e_vec_norm         = interp1(s_vec(knots),taxes_e_vec_norm(knots),s_vec,interp_method,'extrap');
    taxes_w_vec_norm         = interp1(s_vec(knots),taxes_w_vec_norm(knots),s_vec,interp_method,'extrap');
    %!
     share_entre_p_vec_norm   = interp1(s_vec(knots),share_entre_p_vec_norm(knots),s_vec,interp_method,'extrap');
     cond_mean_theta_vec_norm = interp1(s_vec(knots),cond_mean_theta_vec_norm(knots),s_vec,interp_method,'extrap');
     capital_vec_norm         = interp1(s_vec(knots),capital_vec_norm(knots),s_vec,interp_method,'extrap');
     output_vec_norm          = interp1(s_vec(knots),output_vec_norm(knots),s_vec,interp_method,'extrap');
    ave_k_entre_vec_norm     = interp1(s_vec(knots),ave_k_entre_vec_norm(knots),s_vec,interp_method,'extrap');
    ave_n_entre_vec_norm     = interp1(s_vec(knots),ave_n_entre_vec_norm(knots),s_vec,interp_method,'extrap');
    welfare_all_vec_norm     = interp1(s_vec(knots),welfare_all_vec_norm(knots),s_vec,interp_method,'extrap');
    %!
    welfare_entre_vec_norm   = interp1(s_vec(knots),welfare_entre_vec_norm(knots),s_vec,interp_method,'extrap');
    welfare_work_vec_norm    = interp1(s_vec(knots),welfare_work_vec_norm(knots),s_vec,interp_method,'extrap');
    %!
    cev_total_p              = interp1(s_vec(knots),cev_total_p(knots),s_vec,interp_method,'extrap');
    cev_total_se_p           = interp1(s_vec(knots),cev_total_se_p(knots),s_vec,interp_method,'extrap');
    cev_total_w_p            = interp1(s_vec(knots),cev_total_w_p(knots),s_vec,interp_method,'extrap');
    %!
    small_firms_vec          = interp1(s_vec(knots),small_firms_vec(knots),s_vec,interp_method,'extrap');
    large_firms_vec          = interp1(s_vec(knots),large_firms_vec(knots),s_vec,interp_method,'extrap');
    if strcmp(method,'GE')
        r_vec                    = interp1(s_vec(knots),r_vec(knots),s_vec,interp_method,'extrap');
        w_vec                    = interp1(s_vec(knots),w_vec(knots),s_vec,interp_method,'extrap');
    end
end

% Normalize income gap after interpolation
inc_gap_vec_norm = (inc_gap_vec/inc_gap_vec(indx))*100;

%% Plots

fig_ct = 0;

if strcmp(experiment,'no_redistrib')
    %Here do plots for the case without fiscal redistribution
    
    plots_no_redistrib
    
elseif strcmp(experiment,'lump_sum')
    %Here do plots for the case with fiscal redistrib. via lump-sum
    %transfers
    
    plots_lump_sum
    
elseif strcmp(experiment,'tax_cut')
    %Here do plots for the case with fiscal redistrib. via tax cut on all
    
    plots_tax_cut
    
elseif strcmp(experiment,'tax_cut_se')
    %Here do plots for the case with fiscal redistrib. via tax cut for
    %self-employed only
    
    plots_tax_cut_se
    
else
    warning('Selected experiment does not exist!')
    
end
