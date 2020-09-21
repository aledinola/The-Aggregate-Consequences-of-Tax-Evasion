function writefiles(r0,w0,tau_s,V0,Parameters,Grids,mypath,flags)

%--------------------------- LEGEND --------------------------------------%
%{
DESCRIPTION:
 This functions writes variables to txt files to be ported from MATLAB to FORTRAN

INPUTS:
 r0: interest rate
 w0: wage
 tau_s: tax parameter to clear gov bc
 V0: initial guess for value function
 Parameters: structure containing parameters
 Grids: structure containing grids
 mypath: path for the folder where Fortran executable is stored
 flags: structure containing flags

OUTPUTS:
 ED: excess demand (it has to be the first output!)
 w0: wage implied by Cobb Douglas focs
 agg: structure with aggregate variables (K/N, total taxes, etc)
 distrib: structure with distributions
 policy: structure with policy functions
 value: structure with value functions
 exitflag_vfi,exitflag_mu: exit flags for VFI and MU. If >0, failed

Notes: Updated by Alessandro Di Nola on September 15, 2020

%}


if ~isstruct(Parameters)
    error('writefiles: input "Parameters" must be a structure!')
end
if ~isstruct(Grids)
    error('writefiles: input "Grids" must be a structure!')
end
if ~isstruct(flags)
    error('writefiles: input "flags" must be a structure!')
end

%Change path
cd(mypath)

%Real parameters
real_params(1) = Parameters.beta;
real_params(2) = Parameters.delta;
real_params(3) = Parameters.tolV;
real_params(4) = r0;
real_params(5) = w0;
real_params(6) = Parameters.lambda_work;
real_params(7) = Parameters.tau_work;
real_params(8) = tau_s;
real_params(9) = Parameters.lambda_entre;
real_params(10) = Parameters.tau_entre;
real_params(11) = Parameters.s;
real_params(12) = Parameters.lambda;
real_params(13) = Parameters.sigma;
real_params(14) = Parameters.tol_dist;
real_params(15) = Parameters.vi;
real_params(16) = Parameters.tricks_a_down;
real_params(17) = Parameters.tricks_a_up;
real_params(18) = Parameters.tricks_k_down;
real_params(19) = Parameters.tricks_k_up;
real_params(20) = Parameters.tricks_eval_k;
real_params(21) = Parameters.uncmean_eps;
real_params(22) = Parameters.b_work;
real_params(23) = Parameters.p_work;
real_params(24) = Parameters.s_work;
real_params(25) = Parameters.b_entre;
real_params(26) = Parameters.p_entre;
real_params(27) = Parameters.s_entre;
real_params(28) = Parameters.taxrate_work;
real_params(29) = Parameters.taxrate_entre;
real_params(30) = Parameters.sigma2;
real_params(31) = Parameters.psi;
real_params(32) = Parameters.gamma;
real_params(33) = Parameters.pn_1;
real_params(34) = Parameters.pn_2;
real_params(35) = Parameters.tricks_n_up;
real_params(36) = Parameters.tricks_n_down;
real_params(37) = Parameters.pn_3;
real_params(38) = Parameters.cc0;
real_params(39) = Parameters.cc1;
real_params(40) = Parameters.cc2;


real_params = real_params';
dlmwrite(fullfile(mypath,'real_params.txt'), real_params, 'delimiter', '\t', 'precision', 16);

%Integer parameters
int_params(1) = Parameters.nagrid;
int_params(2) = Parameters.nagrid_fine;
int_params(3) = Parameters.nagrid_dist;
int_params(4) = Parameters.negrid;
int_params(5) = Parameters.ntgrid;
int_params(6) = Parameters.maxitV;
int_params(7) = Parameters.evade_only_profit;
int_params(8) = Parameters.nphi;
int_params(9) = Parameters.nkap;
int_params(10) = Parameters.mu_method;
int_params(11) = Parameters.vfi_tricks;
int_params(12) = Parameters.howard_flag;
int_params(13) = Parameters.accmax;
int_params(14) = Parameters.taxfunc;
int_params(15) = Parameters.maxit_dist;
int_params(16) = Parameters.fix_occpol;
int_params(17) = Parameters.fix_kpol;
int_params(18) = Parameters.fix_apol_nd;
int_params(19) = Parameters.nlgrid;
int_params(20) = Parameters.nngrid;
int_params(21) = Parameters.concavity;
int_params(22) = Parameters.pflag;
int_params(23) = Parameters.nlegrid;
int_params(24) = flags.display_howard;
int_params(25) = flags.display_mu;
int_params(26) = flags.display_iter;
int_params(27) = flags.par_fortran;

int_params = int_params';
dlmwrite(fullfile(mypath,'int_params.txt'), int_params, 'delimiter', '\t', 'precision', 16);

%Dimensions
dimensions = [length(real_params);length(int_params)]; %the dimensions
dlmwrite(fullfile(mypath,'dimensions.txt'), dimensions, 'delimiter', '\t', 'precision', 16);

%Other objects
dlmwrite(fullfile(mypath,'V0.txt'), V0(:), 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'agrid.txt'), Grids.agrid, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'agrid_fine.txt'), Grids.agrid_fine, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'agrid_dist.txt'), Grids.agrid_dist, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'eps_grid.txt'), Grids.eps_grid, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'theta.txt'), Grids.theta', 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'phigrid.txt'), Grids.phigrid, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'kgrid.txt'), Grids.kgrid, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'Peps.txt'), Grids.Peps(:), 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'Ptheta.txt'), Grids.Ptheta(:), 'delimiter', '\t', 'precision', 16);
%dlmwrite(fullfile(mypath,'pk.txt'), pk, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'lgrid.txt'), Grids.lgrid, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'ngrid.txt'), Grids.ngrid, 'delimiter', '\t', 'precision', 16);
dlmwrite(fullfile(mypath,'legrid.txt'), Grids.legrid, 'delimiter', '\t', 'precision', 16);

% If we do decompositions we need to use some polfuns from no tax
% evasion economy (with prices from no tax evasion economy)
if Parameters.fix_occpol==1
    % load notax GE
    load([flags.results_folder,'\notaxevasion_ge.mat'],'policy')
    occpol_fixed = policy.occpol;
    clear policy
    % For table 2 decompositions we keep _note even if the objects are
    % from tax evasion economy
    dlmwrite(fullfile(mypath,'occpol_fixed_note.txt'), occpol_fixed(:), 'delimiter', '\t', 'precision', 16);
end

if Parameters.fix_kpol==1
    % load notax GE
    load([flags.results_folder,'\notaxevasion_ge.mat'],'policy')
    policycap_fixed = policy.policycap; %real
    Icap_fixed      = policy.Icap; %integer
    clear policy
    % For table 2 decompositions we keep _note even if the objects are
    % from tax evasion economy
    dlmwrite(fullfile(mypath,'policycap_fixed_note.txt'), policycap_fixed(:), 'delimiter', '\t', 'precision', 16);
    dlmwrite(fullfile(mypath,'Icap_fixed_note.txt'), Icap_fixed(:), 'delimiter', '\t', 'precision', 16);
end

if Parameters.fix_npol==1
    load([flags.results_folder,'\notaxevasion_ge.mat'],'policy')
    policyn_fixed = policy.policyn; %real
    Inpol_fixed   = policy.Inpol;   %integer
    clear policy
    % For table 2 decompositions we keep _note even if the objects are
    % from tax evasion economy
    dlmwrite(fullfile(mypath,'policyn_fixed_note.txt'), policyn_fixed(:), 'delimiter', '\t', 'precision', 16);
    dlmwrite(fullfile(mypath,'Inpol_fixed_note.txt'), Inpol_fixed(:), 'delimiter', '\t', 'precision', 16);

end

if Parameters.fix_lpol==1
    load([flags.results_folder,'\notaxevasion_ge.mat'],'policy')
    %Integer policies for labor supply, W and SE
    decisl_fixed   = policy.decisl;   
    lepolind_fixed = policy.lepolind;

    clear policy
    % For table 2 decompositions we keep _note even if the objects are
    % from tax evasion economy
    dlmwrite(fullfile(mypath,'decisl_fixed_note.txt'), decisl_fixed(:), 'delimiter', '\t', 'precision', 16);
    dlmwrite(fullfile(mypath,'lepolind_fixed_note.txt'), lepolind_fixed(:), 'delimiter', '\t', 'precision', 16);

end

end %END FUNCTION "writefiles"
