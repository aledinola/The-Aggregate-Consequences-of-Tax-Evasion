function [Parameters,Grids] = fun_parameters(Parameters,flags)

%{
This function assigns some additional parameters to the M-structure
Parameters, generates discretized grids and stores them in the
M-structure Grids

INPUTS:
- Parameters: structure with parameters
- flags:      structure with some flags
 
OUTPUTS:
- Parameters: structure with parameters
- Grids:      structure with grids for discretized assets,shocks, etc.

Notes: Updated by Alesandro Di Nola on 14 September 2020
%}

%% Some additional flags

Parameters.labor_hiring  = 1; %0/1 labor hirings by entre
Parameters.endo_ls_work  = 1; %0/1 exogenous/endogenous labor supply workers
Parameters.endo_ls_entre = 1; %0/1 exogenous/endogenous labor supply entre
minLwork = 0.0001; maxLwork = 2.9999; %total available time for workers
minEwork = 0.0001; maxEwork = 2.9999; %total available time for entre

%% Computational Parameters

Parameters.tolV        = 10^-4; % Tolerance for value function (it can be looser than tol_dist)
Parameters.tol_dist    = 10^-6; % Tolerance for the distribution (please set it lower than 1e-5)
Parameters.maxitV      = 20;    % Maximum number of iterations for the value function (with howard algo, you never need more than 15-20 iterations)
Parameters.maxit_dist  = 10000; % Maximum number of iterations for the distribution mu
grid_method            = 2;    % if = 1 equispaced grid, if = 2 then grid is finer for smaller values of a

% We use these flags for decompositions. If decomp_case=0,these flags
% remain at default values. O/w they are reset accordingly.
Parameters.fix_occpol  = 0;    % 0 is free, 1 is fixed. o(x), occupational choice
Parameters.fix_kpol    = 0;    % 0 is free, 1 is fixed. k(x), capital demand
Parameters.fix_npol    = 0;    % 0 is free, 1 is fixed. n(x), labor demand
Parameters.fix_lpol    = 0;    % 0 is free, 1 is fixed. l(x), labor supply
Parameters.fix_apol_nd = 0;    % NOT USED but present in the code, leave it ==0

switch flags.decomp_case
% See Table(8), Sec 5.3. In cases {2:5}, prices are fixed at (1)    
    case 1 %No tax evasion GE
        disp('(1) No tax evasion GE')
        Parameters.no_evasion = 1; %1=WITHOUT tax evasion,0=WITH tax evasion
        Parameters.fix_occpol = 0;
        Parameters.fix_kpol   = 0;
        Parameters.fix_npol   = 0;
        Parameters.fix_lpol   = 0;
        
    case 2 %Fixed: o(x),k(x),n(x) ==> subsidy channel
        disp('(2) Fixed: o(x),k(x),n(x) ==> subsidy channel')
        Parameters.no_evasion = 0;
        Parameters.fix_occpol = 1;
        Parameters.fix_kpol   = 1;
        Parameters.fix_npol   = 1;
        Parameters.fix_lpol   = 0;
        
    case 3 %Fixed: o(x) ==> subsidy channel + detection channel
        disp('(3) Fixed: o(x) ==> subsidy channel + detection channel')
        Parameters.no_evasion = 0;
        Parameters.fix_occpol = 1;
        Parameters.fix_kpol   = 0;
        Parameters.fix_npol   = 0;
        Parameters.fix_lpol   = 0;
        
        
    case 4 %Fixed: k(x) and n(x) ==> subsidy channel + selection channel
        disp('(4) Fixed: k(x) and n(x) ==> subsidy channel + selection channel')
        Parameters.no_evasion = 0;
        Parameters.fix_occpol = 0;
        Parameters.fix_kpol   = 1;
        Parameters.fix_npol   = 1;
        Parameters.fix_lpol   = 0;
        
    case 5 %All decisions are endogenous, prices still fixed at (1)
        disp('(5) All decisions are endogenous, prices still fixed at (1)')
        Parameters.no_evasion = 0;
        Parameters.fix_occpol = 0;
        Parameters.fix_kpol   = 0;
        Parameters.fix_npol   = 0;
        Parameters.fix_lpol   = 0;
    case 6 %Benchmark WITH tax evasion, GE (r = 4%)
        disp('(6) Benchmark WITH tax evasion, GE')
        Parameters.no_evasion = 0;
        Parameters.fix_occpol = 0;
        Parameters.fix_kpol   = 0;
        Parameters.fix_npol   = 0;
        Parameters.fix_lpol   = 0;
    otherwise
        disp('NO decompositions, standard analysis')
        
end
disp('-----------------------------------------------------------')

%functional form of p()
pk_flag = 1; % 1 = pk (prob of being caught) logistic, 2=pk constant, 3= pk step function 0-1, 4= convex/concave
             % ONLY logistic is operational in the current version

%argument of p():
pflag = 4; %1=misrep income
           %2=f(theta,k,n)
           %3=f(k,n), new try
           %4=k
           %5=n

% Legend:
%mu_method 1 = 'finer'  ! Distribution is computed on finer grid, size: nagrid_fine with lotteries a la Rios Rull
%mu_method 2 = 'finer2'  ! as in 'finer' but WITHOUT lotteries
Parameters.mu_method         = 1;
Parameters.evade_only_profit = 1; % NOT USED but DO NOT CHANGE

%% Fortran controls
Parameters.fortran_flag  = 1;  % 0 = MATLAB; 1 = FORTRAN i.e. part of the code is executed in Fortran
                               % Note: We do not provide the MATLAB
                               % implementation, so fortran_flag must be
                               % set to 1
                               
                               
% The flags below switch on/off some speed-up tricks.                                
Parameters.howard_flag   = 1;  % 1 = do Howard policy improvement algorithm
Parameters.accmax        = 100; % num. of accelerations for Howard algo
Parameters.vfi_tricks    = 1;  % 1=yes, 0=no tricks
Parameters.concavity     = 1;  %concavity


%local search upward, suggested value: 0.1
%local search down, suggested value: 0.02
%small grid for n, set tricks_n_down/up = 1
Parameters.tricks_k_down = 0.02; %(- *100 percentage around which the FORTRAN/MATLAB executable will search for k(a) based on k(a(-1))
Parameters.tricks_k_up   = 0.2; %(+ *100 percentage around which the FORTRAN/MATLAB executable will search for k(a) based on k(a(-1))
Parameters.tricks_n_down = 0.02; %(- *100 percentage around which the FORTRAN/MATLAB executable will search for n(a) based on n(a(-1))
Parameters.tricks_n_up   = 0.2;%0.1; %(- *100 percentage around which the FORTRAN/MATLAB executable will search for n(a) based on n(a(-1))

%NOT USED but DO NOT CHANGE
%We imposed monotonicity and concavity on a'(a) ==> no need to do local
%search on a' anymore
Parameters.tricks_a_up   = 0.1; %NOT USED
Parameters.tricks_a_down = 0.1; %NOT USED
Parameters.tricks_eval_k = 1.0; %NOT USED


%% Income tax function

% Flag for tax function
% Legend: 1;    Heathcote, Storesletten and Violante (2016)
%         2;    Gouveia and Strauss (1994) tax function
%         3;    simple proportional tax rate
Parameters.taxfunc = 2; 

% (1) HSV tax function
Parameters.lambda_entre  = 0.9; % 10% at mean income
Parameters.lambda_work   = 0.88;
Parameters.tau_entre     = 0.087;
Parameters.tau_work      = 0.09;
% (2) Gouveia and Strauss tax function
% (we borrowed estimates for b,p from Cagetti & DeNardi 2009)
% s_w and s_e are determined in equilibrium to clear gbc
% scale_a2 is chosen in the GOV loop
Parameters.p_work  = 0.76; % a1 in the paper
Parameters.b_work  = 0.32;
Parameters.p_entre = 1.40;
Parameters.b_entre = 0.26;
Parameters.s_work = (Parameters.ksi_tax^Parameters.p_work)*0.22; %a2_w = 0.22 (De Nardi)
Parameters.s_entre = (Parameters.ksi_tax^Parameters.p_entre)*0.44; %a2_e = 0.44 (De Nardi)

% (3) Proportional tax function
Parameters.taxrate_entre = 0.1;
Parameters.taxrate_work  = 0.1;

%% Grid for assets: we use 3 grids of different size

nagrid      = 30;%40;  %number of asset grid points
nagrid_fine = nagrid*30;%nagrid*30; % choice set for a' (next-period assets)
nagrid_dist = nagrid*10;%nagrid*30; % grid for distribution
% #grid for a < #grid for a' <= #grid for mu(a)
% nagrid      < nagrid_dist  <= nagrid_fine

abar = 1e-6; % Lower bound for asset holdings. Please do not change 
amin = abar; % min value of assets
amax = 250;  % max value of assets, should not bind

if  grid_method == 1     % grid when equally spaced
    agrid = linspace(amin,amax,nagrid)';
    agrid_fine = linspace(amin,amax,nagrid_fine)';
    agrid_dist = linspace(amin,amax,nagrid_dist)';
elseif grid_method == 2  % grid when logarithmically spaced
    agridtemp = linspace(0,log(amax+1-amin),nagrid)';
    agrid = exp(agridtemp)-1+amin;
    agridtemp_fine =linspace(0,log(amax+1-amin),nagrid_fine)';
    agrid_fine = exp(agridtemp_fine)-1+amin;
    agridtemp_dist =linspace(0,log(amax+1-amin),nagrid_dist)';
    agrid_dist = exp(agridtemp_dist)-1+amin;
end

%plot(agrid,agrid,'o')

%% Labor productivity grid (variable \eps in the paper)
% Unconditional mean for eps normalized to one
% AR1 process, discretize using Tauchen and Hussey (1991), (see Kitao 2008)

siglog       = Parameters.sigmaeps^2/(1-Parameters.rho^2);
negrid       = 7;         % number of epsilon (working ability) grid points
uncmean_eps1 = -siglog/2; % unconditional mean normalized to unity

% eps_grid  - nodes for labor productivity, (negrid,1) array
% Peps      - transition probability for epsilon, (negrid,negrid) array
[eps_grid1,Peps] = tauchenhussey(negrid,uncmean_eps1,Parameters.rho,Parameters.sigmaeps, Parameters.sigmaeps);
eps_grid         = exp(eps_grid1);

%Compute unconditional distrib. of epsilon shock
Petemp   = Peps^5000;
eps_dist = Petemp(1,:);
Parameters.uncmean_eps = eps_dist*eps_grid;

%% Enterpreneurial productivity grid (variable \theta in the paper)

ntgrid = 12;    % number of grid points for theta

% parameters of theta process
% rho_theta
% Parameters.sigmaeps_theta
% Parameters.uncmean_theta

% Unconditional mean: Parameters.uncmean_theta,
%                     is a calibrated/estim parameter
% See Cagetti and De Nardi (2006)

% AR1 process, discretize using Tauchen and Hussey (1991) or Rouwenworst

%[theta,Ptheta] = tauchenhussey(ntgrid,uncmean,rho_theta,Parameters.sigmaeps_theta, Parameters.sigmaeps_theta);
[Ptheta_rouwen, theta] = rouwen(Parameters.rho_theta,Parameters.uncmean_theta,Parameters.sigmaeps_theta, ntgrid);
Ptheta = Ptheta_rouwen'; %we want rows summing to one, i.e. sum(Ptheta,2)=1 
theta  = exp(theta);

%% Capital grid (variable k in the paper)

nkap = round(nagrid_dist/1);   % number of grid points for capital
kmin = 1e-6;                   % lower bound for capital 
kmax = Parameters.lambda*amax; % upper bound for capital

if  grid_method == 1       % grid when equally spaced
    kgrid = linspace(kmin, kmax, nkap)';
elseif grid_method == 2    % grid when logarithmically spaced
    kgridtemp = linspace(0,log(kmax+1-kmin),nkap)';
    kgrid = exp(kgridtemp)-1+kmin;
end

%% Evasion grid (variable \phi in the paper)

if Parameters.no_evasion==0
    nphi    = 2;%15; % number of grid points for evasion
    phigrid = linspace(0,1,nphi)'; % grid for tax evasion, always between zero and one
elseif Parameters.no_evasion==1
    nphi    = 1; % number of grid points for evasion
    phigrid = linspace(0,0,nphi)'; % grid for tax evasion, always between zero and one
end

%% Labor supply grid workers (variable \ell in the paper)

switch Parameters.endo_ls_work
    case 1
        nlgrid = 20;
        lgrid  = linspace(minLwork,maxLwork,nlgrid)';
        
    case 0 %Exogenous LS: W supply inelastically 1 unit of time
        nlgrid = 1;  %they are forced to "choose" only one point
        lgrid  = maxLwork;
end


%% Labor supply grid entre (variable \ell in the paper)

switch Parameters.endo_ls_entre
    case 1
        nlegrid = 20;
        legrid  = linspace(minEwork,maxEwork,nlegrid)';
        
    case 0 %Entre supply inelastically le (can be 0.33 or 1)
        nlegrid = 1;     %they are forced to "choose" only one point
        legrid  = maxEwork;
end

%% Labor hiring (employees hired by entrep, variable n in the paper)

switch Parameters.labor_hiring
    
    case 1 %YES, entre can hire workers
        
        Ngrid_method = 3; %see below
        nngrid = 100;
        nmin   = 0;
        nmax   = kmax/3;
        curv   = 1.5; %curv=1, equal spacing, curv>1 more points near the lower bound
        if  Ngrid_method == 1       % grid when equally spaced
            ngrid = linspace(nmin, nmax, nngrid)';
        elseif Ngrid_method == 2    % grid when logarithmically spaced
            ngridtemp = linspace(0,log(nmax+1-nmin),nngrid)';
            ngrid = exp(ngridtemp)-1+nmin;
        elseif Ngrid_method == 3
            ngrid = zeros(nngrid,1);
            for i=1:nngrid
                ngrid(i) = nmin + ((i-1)/(nngrid-1))^curv*(nmax-nmin);
                
            end
        end
        
    case 0  %NO
        nngrid = 1;
        ngrid  = 0;
end

%% Probability of being caught
% Now defined as a function file

% Display some messages
% Summarize dimensions of arrays
if flags.Verbose == 1
    
    disp('--------------------------------------------------------------------')
    disp(' NUMBER OF GRID POINTS ')
    disp('--------------------------------------------------------------------')
    disp('Number of grid points for current assets')
    disp(nagrid)
    disp('Number of grid points for next-period assets')
    disp(nagrid_fine)
    disp('Number of grid points for distribution')
    disp(nagrid_dist)
    disp('Number of grid points for capital')
    disp(nkap)
    disp('Number of grid points for labor hiring')
    disp(nngrid)
    disp('Number of grid points for evasion')
    disp(nphi)
    disp('Number of grid points for labor supply W')
    disp(nlgrid)
    disp('Number of grid points for labor supply SE')
    disp(nlegrid)
    disp('Number of grid points for epsilon')
    disp(negrid)
    disp('Number of grid points for theta')
    disp(ntgrid)
    
    if Parameters.labor_hiring == 0
        disp('No labor hiring in SE prodfun')
    elseif Parameters.labor_hiring == 1
        disp('labor hiring in SE prodfun')
    end
    if Parameters.endo_ls_work == 0
        disp('Exogenous labor supply W')
    elseif Parameters.endo_ls_work == 1
        disp('ENDOgenous labor supply W')
    end
    if Parameters.endo_ls_entre == 0
        disp('Exogenous labor supply ENTRE')
    elseif Parameters.endo_ls_entre == 1
        disp('ENDOgenous labor supply ENTRE')
    end
    
    if Parameters.no_evasion==1
        disp('Simple economy WITHOUT tax evasion')
    elseif Parameters.no_evasion==0
        disp('Benchmark economy with TAX EVASION')
        if pk_flag==1
            disp('p(.) logistic')
        elseif pk_flag==2
            disp('p(.) = p constant')
        elseif pk_flag==3
            disp('p(.) stepwise')
        elseif pk_flag==4
            disp('p(.) convex/concave')
        end
        if pflag==1
            disp('p(.) argument: misrep income, option 5 in memo')
        elseif pflag==2
            disp('p(.) argument: f(theta,k,n), option 6 in memo')
        elseif pflag==3
            disp('p(.) argument: f(k,n), new try')
        elseif pflag==4
            disp('p(.) argument: k, new try')
        end
    end
    if Parameters.vfi_tricks==1
        disp('VFI tricks on k,n: ON')
    else
        disp('VFI tricks on k,n: OFF')
    end
    if Parameters.concavity==1
        disp('VFI using concavity: ON')
    else
        disp('VFI using concavity: OFF')
    end
    if Parameters.howard_flag==1
        disp('Howard acceleration: ON')
    else
        disp('Howard acceleration: OFF')
    end
    disp('--------------------------------------------------------------------')
end %disp if verbose is on

%% Collect some stuff into a structure

% Pack all grids (and Markov chains) into structure "Grids"
Grids.agrid      = agrid;
Grids.agrid_fine = agrid_fine;
Grids.agrid_dist = agrid_dist;
Grids.eps_grid   = eps_grid;
Grids.Peps       = Peps;
Grids.theta      = theta;
Grids.Ptheta     = Ptheta;
Grids.kgrid      = kgrid;
Grids.phigrid    = phigrid;
Grids.lgrid      = lgrid;
Grids.ngrid      = ngrid;
Grids.legrid     = legrid;

% Pack grid dimensions into structure "Parameters"
Parameters.nagrid      = nagrid;
Parameters.nagrid_fine = nagrid_fine;
Parameters.nagrid_dist = nagrid_dist;
Parameters.negrid      = negrid;
Parameters.ntgrid      = ntgrid;
Parameters.nkap        = nkap;
Parameters.nphi        = nphi;
Parameters.nlgrid      = nlgrid;
Parameters.nngrid      = nngrid;
Parameters.pflag       = pflag;
Parameters.nlegrid     = nlegrid;
Parameters.pk_flag     = pk_flag;
Parameters.minLwork    = minLwork;
Parameters.minEwork    = minEwork;
Parameters.maxLwork    = maxLwork;
Parameters.maxEwork    = maxEwork;
Parameters.nmin        = nmin;
Parameters.kmin        = kmin;
Parameters.abar        = abar;

end %END FUNCTION

