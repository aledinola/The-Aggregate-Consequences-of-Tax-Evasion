%------------------ LEGEND -----------------------------------------------%

% Main program to replicate model results in 
% "The Aggregate Consequences of Tax Evasion" forthcoming Review of Economic
% Dynamics 
% by Alessandro Di Nola, Georgi Kocharkov, Almuth Scholl and Anna-Mariia Tkhir
% 
% Notes: Updated on September 17, 2020. 
       
%-------------------------------------------------------------------------%

%{
Options to compute benchmark economy and several counterfactuals and
robustness checks:

- If do_GE = 0, the code computes a *partial* equilibrium
  for given vector of parameters "guess". Specifically, it fixes the
  interest rate (at 4 percent) and does not clear the asset market.

- If do_GE = 1, the code computes an equilibrium
  for given vector of parameters "guess". Specifically, it fixes beta at
  the calibrated value and find the interest rate r that clears the market.

- Counterfactual with perfect tax enforcement: set no_evasion = 1.

- Replication instructions: See comments below and README file provided in 
  the replication folder. Additional details can be found at 
  https://github.com/aledinola/The-Aggregate-Consequences-of-Tax-Evasion

- To replicate "Tax Evasion Benchmark", column 1 in Table 7, set
  do_GE=0,do_gov_loop=0 and no_evasion=0 and do_lambda_exp=0.
  To replicate counterfactual with "Perfect Tax Enforcement", column 2 in
  Table 7, set do_GE=1,do_gov_loop=0 and no_evasion=1
  To replicate counterfactual with "Perfect Tax Enforcement" and lump-sum 
  redistribution, set do_GE=1,do_gov_loop=1,no_evasion=1,balance='lump_sum'.
  To replicate counterfactual with "Perfect Tax Enforcement" and cut tax all 
  redistribution, set do_GE=1,do_gov_loop=1,no_evasion=1,balance='tax_cut'.
  To replicate counterfactual with "Perfect Tax Enforcement" and cut tax 
  for self-employed only, set do_GE=1,do_gov_loop=1,no_evasion=1 and 
  balance='tax_cut_se'.

- In all the above economies, set decomp_case = 0. If decompositions are
  desired, set decomp_case=1:6, depending on which column of Table (8) the
  user wants to replicate.

- Lambda experiments (Section 5.5, Tax evasion and credit constraints)
  For lambda=1.2, set do_lambda_exp=1, do_GE=0, do_gov_loop=0, no_evasion=0  
  For lambda=1.8, set do_lambda_exp=2, do_GE=0, do_gov_loop=0, no_evasion=0 

- Auxiliary files called by this script:
  + function f_solve_model_withGOV.m, which in turn calls:
            + fun_parameters.m
            + PE_interpolation1.m
            + targets_compute_clean.m
            + fun_estimation.m

%}

clear; close all;clc;format long g;

addpath(genpath('utilities'));

%% Computational Parameters

do_GE           = 0; % 0 = PE, 1 = GE
do_gov_loop     = 0; % 0 = no gov loop, 1 = gov loop
no_evasion      = 0; % if = 1 then an economy W/OUT tax evasion is computed (where \phi is restricted to be zero)
do_lambda_exp   = 0; %Options: 0,1,2.

par_fortran     = 1; % set it = 1 if OpenMP for Fortran VFI
display_iter    = 1; % 0=NOT display any fortran messages,1=low verbosity, 2=high verbosity
display_howard  = 0; % 1 = display iterations Howard improv. step
display_mu      = 0; % 1 = display iterations for distribution
flags.Verbose   = 0; % 1 = display intermediate calibration results
debug_mode      = 0; % set = 1 to perform checks (it runs slower)
plot_figures    = 0; % set = 1 to plot figures
flags.lowMemory = 0; % set = 1 if you don't have much RAM on your system
n_cores         = 8; % number of cores
flags.doSave    = 0; % 1=save final results in .mat file (after targets_compute_clean)
flags.loadV     = 1; % 1=use pre-exisiting V as initial guess for VFI
flags.saveV     = 0; % 1=save new V as initial guess
do_estimation   = 0;    %NOT USED but present in the code, leave it ==0
est_algo        = 'ga'; %NOT USED but present in the code

% Flags for decompositions
% Optionns: {0:6}. See Table (8) in the paper. If 0, no decompositions
% Options from 1 to 6 replicate columns 1-6 in Table (8).
flags.decomp_case = 0; 
% 1 - No tax evasion, GE
% 2 - Fixed: o(x),k(x),n(x) ==> subsidy channel
% 3 - Fixed: o(x) ==> subsidy channel + detection channel
% 4 - Fixed: k(x) and n(x) ==> subsidy channel + selection channel
% 5 - All decisions are endogenous, prices still fixed at (1)
% 6 - Benchmark with tax evasion, GE
% Note: In {2:5}, prices are fixed at 1.

% Some useful paths
main_folder     = fullfile(pwd,'');
exe_folder      = fullfile(main_folder,'exe');
results_folder  = fullfile(main_folder,'results');

% Pack into M-structures
flags.par_fortran    = par_fortran;
flags.display_iter   = display_iter;
flags.main_folder    = main_folder;
flags.exe_folder     = exe_folder;
flags.results_folder = results_folder;
flags.do_estimation  = do_estimation;
flags.plot_figures   = plot_figures;
flags.debug_mode     = debug_mode;
flags.do_GE          = do_GE;
flags.do_gov_loop    = do_gov_loop;
flags.n_cores        = n_cores;
flags.display_howard = display_howard;
flags.display_mu     = display_mu;
flags.display_iter    = display_iter;
Parameters.no_evasion = no_evasion;

if do_estimation~=0
    error('do_estimation must be set equal to 0')
end

%% First, set exogenous parameters (not part of internal calibration)
Parameters.sigma    = 2;     % risk aversion (called \sigma_1 in the paper)
Parameters.sigma2   = 1.67;  % labor supply elasticity (inverse of Frisch elast.)
Parameters.alpha    = 0.38;  % share of capital - corporate sector
Parameters.s        = 1.75;  % tax evasion fine
Parameters.rho      = 0.89;  % persistence labor productivity shock (called \rho_{\varepsilon} in the paper)
Parameters.sigmaeps = 0.21;  % standard dev. labor productivity shock
Parameters.cc1      = 0;     % NOT USED but present in the code, hence==0
Parameters.cc2      = 0;     % NOT USED but present in the code, hence==0
Parameters.pn_3     = 0;     % NOT USED but present in the code, hence==0

%% Set bounds for internal parameters
%[First number is lower bound, second number is the upper bound]
%See also Table 2 in the paper
bounds.beta          = [0.90, 0.96]; % Discount factor
bounds.psi           = [0.5, 0.9];   % Disutility hours worked
bounds.delta         = [0.0, 10];    % Capital depreciation
bounds.vi            = [0.5, 0.9];   % Span of control
bounds.gamma         = [0.1, 1.0];   % Capital share, self-employed
bounds.rho_theta     = [0.1, 0.99];  % persistence of ability shock
bounds.sigma_theta   = [0.1, 2.0];   % stdev of ability shock
bounds.uncmean_theta = [-2, -0.5];   % unconditional mean of ability shock
bounds.pn_1          = [10, 10000];  % parameter prob detection: intercept
bounds.pn_2          = [0.005, 0.5]; % parameter prob detection: slope
bounds.cc0           = [0.0, 0.3];   % fixed cost of evasion (\kappa in the paper)
bounds.ksi_tax       = [4.0, 15.0];  % Rescale for tax function
bounds.lambda        = [1.2, 1.8];   % Leverage ratio collateral constraint


%% Set parameter names IN THE SAME ORDER AS THEY APPEAR in bounds and guess
pnames = {'beta';
          'psi';
          'delta';
          'vi';
          'gamma';
          'rho_theta';
          'sigma_theta';
          'uncmean_theta';
          'pn_1';
          'pn_2';
          'cc0';
          'ksi_tax';
          'lambda'};

%% Manual runs and experiments

if Parameters.no_evasion==1
    disp('Counterfactual economy with perfect tax enforcement')
elseif Parameters.no_evasion==0
    disp('Benchmark economy with tax evasion')
end

if do_lambda_exp == 0
    disp('MANUAL RUN: parameters are fixed at their calibrated values')
    beta1 = 0.945332328222695; 
    guess =  [beta1   %beta
              0.83    %psi
              0.11    %delta
              0.74    %vi
              0.73    %gamma
              0.952   %rho_theta
              0.67    %sigma_theta
              -1.12   %mu_theta
              2250    %p1
              0.35    %p2
              0.1315  %c0 (or kappa)
              11.0    %ksi
              1.50];  %lambda
   
    [obj_smm1,model_targets,policy,r,ED,ModelResults,agg] = f_solve_model_withGOV(guess, Parameters, bounds, flags,pnames);
    
elseif do_lambda_exp ==1
    disp('Lambda Experiment: lambda=1.2')
    % lambda = 1.2
    beta1 = 0.9448074;
    guess=  [beta1 0.83 0.11 0.74 0.73 0.952 0.67 -1.1 2250 0.35 0.137 11.0 1.2]';  
    [obj_smm1,model_targets,policy,r,ED,ModelResults,agg] = f_solve_model_withGOV(guess, Parameters, bounds, flags,pnames);
    
elseif do_lambda_exp ==2
    disp('Lambda Experiment: lambda=1.8')
    % lambda = 1.8
    beta1 = 0.9458858;
    guess=  [beta1 0.83 0.11 0.74 0.73 0.952 0.67 -1.133 1750 0.35 0.127 11.0 1.8]';
    [obj_smm1,model_targets,policy,r,ED,ModelResults,agg] = f_solve_model_withGOV(guess, Parameters, bounds, flags,pnames);
else
    error('do_lambda_exp must be either 0,1,2')
end 

%% Display some model results on the screen
% For more detailed output, see M-files make_tables_*

fprintf('\n')
disp('Some results from the model:')
fprintf('\n')
fprintf('K/Y: %f \n',ModelResults.K_to_Y)
fprintf('Interest rate: %f \n',ModelResults.r0)
fprintf('Wage: %f \n',ModelResults.w0)
fprintf('Share of self-employed: %f \n',ModelResults.share_entre_p)
fprintf('E(theta|SE): %f \n',ModelResults.cond_mean_theta)
fprintf('E(k|SE): %f \n',ModelResults.ave_k_entre)
fprintf('Taxes total: %f \n',agg.taxes_total)
