%{
DESCRIPTION:
  This script computes consumption equivalent variation (CEV) for different 
  scenarios. Loads previously computed economies from results\mat.
  Use this program to reproduce results shown in Figure 4 and Table 9.

EXTERNAL FILES:
    +welfare_computations_ce1.m
        +g_operator.m
    +welfare_computations_ce2.m
        +g_operator.m

** Updated by Alessandro Di Nola on September 18, 2020.

%}

clear;clc

addpath(genpath('utilities'));

%% Set up useful paths for folders

%Specify here where results (.mat files) are saved
ResultsDir = 'results\mat\';
%Specify here where you want to save the welfare_results/figures/tables 
SaveDir = 'results\Figures\';

%% Set some flags for welfare analysis

% For example, to reproduce Figure(4a) of the paper, set 
% experiment = 'no_redistrib' .and. cev_definition = 'CE1'
% The program will compute the consumption equivalent variation of
% eliminating tax evasion w/out imposing fiscal neutrality

% To reproduce Figure(4b) of the paper, set 
% experiment = 'lump_sum' .and. cev_definition = 'CE1'
% The program will compute the CEV of eliminating tax evasion with lump-sum 
% transfers. And so on...

% Note: the *partial equilibrium* economies (suffix _pe) are not used in the
% published version of the paper.

do_welfare_plots = 1;
experiment     = 'no_redistrib';   %Options:
                                    %'no_redistrib_pe'
                                    %'no_redistrib'
                                    %'lump_sum_pe'
                                    %'lump_sum'
                                    %'cut_tax_all_pe'
                                    %'cut_tax_all'
                                    %'cut_tax_se_pe'
                                    %'cut_tax_se'
cev_definition = 'CE1'; %Options:
                        %'CE1' (micro), preferred and used for the paper
                        %'CE2' (macro)

%% Load benchmark economy
%%% HERE load benchmark economy %%%
load([ResultsDir 'taxevasion_ge.mat'],'policy','value','distrib','Parameters','Grids','flags','w0','r0')

policy_bench  = policy;
value_bench   = value;
distrib_bench = distrib;

%Add smth on the fly - TO BE REMOVED LATER
Grids.wealth      = repmat(Grids.agrid_dist,[1 Parameters.negrid Parameters.ntgrid]); % array [nagrid_dist,negrid,ntgrid]
Grids.wealth_vec  = Grids.wealth(:); % array [nagrid_dist*negrid*ntgrid,1]

clear policy value distrib

%% Load counterfactual (reform) economy
%%% HERE load a counter-factual economy you are interested in %%%
if strcmp(experiment,'no_redistrib')
    load([ResultsDir 'notaxevasion_ge.mat'],'policy','value','distrib')
elseif strcmp(experiment,'no_redistrib_pe')
    load([ResultsDir 'notaxevasion_pe.mat'],'policy','value','distrib')
elseif strcmp(experiment,'lump_sum') %Makes more sense to use GE
    load([ResultsDir 'notaxevasion_ge_lump_sum.mat'],'policy','value','distrib')
elseif strcmp(experiment,'lump_sum_pe')
    load([ResultsDir 'notaxevasion_pe_lump_sum.mat'],'policy','value','distrib')
elseif strcmp(experiment,'cut_tax_all')
    load([ResultsDir 'notaxevasion_ge_cut_tax_all.mat'],'policy','value','distrib')
elseif strcmp(experiment,'cut_tax_all_pe')
    load([ResultsDir 'notaxevasion_pe_cut_tax_all.mat'],'policy','value','distrib')
elseif strcmp(experiment,'cut_tax_se')
    load([ResultsDir 'notaxevasion_ge_cut_tax_se.mat'],'policy','value','distrib')
elseif strcmp(experiment,'cut_tax_se_pe')
    load([ResultsDir 'notaxevasion_pe_cut_tax_se.mat'],'policy','value','distrib')
else
    error('Invalid option')
end

policy_reform  = policy;
value_reform   = value;
distrib_reform = distrib;

clear policy value distrib

%% Main computation is outsourced to external functions
disp(['Comparing baseline with ', experiment])
switch cev_definition
    case 'CE1' %individual measure (Guvenen et al. 2019)
        [cev_total,cev_total_se,cev_total_w,cev_se,cev_w] = welfare_computations_ce1...
            (policy_bench,value_bench,distrib_bench,policy_reform,value_reform,distrib_reform,Parameters,Grids,flags,w0,r0);
        
    case 'CE2' %Lucas 1987
        [cev_total,cev_total_se,cev_total_w,cev_se,cev_w] = welfare_computations_ce2...
            (policy_bench,value_bench,distrib_bench,policy_reform,value_reform,distrib_reform,Parameters,Grids,flags,w0,r0);
        
    otherwise
        error('Selected option is NOT available!')
end

%% Save results 

save([SaveDir,'welfare_results_',experiment]);

%% Display CEV's on Matlab screen
%Convert to percentages
cev_total_p    = cev_total*100;
cev_total_se_p = cev_total_se*100;
cev_total_w_p  = cev_total_w*100;
cev_se_p       = cev_se*100;
cev_w_p        = cev_w*100;

fprintf('\n')
disp(['Comparing baseline with ', experiment])
disp('TABLE 9')
disp(['CEV definition: ', cev_definition])
fprintf('Total CEV(%%): %f \n',cev_total_p)
fprintf('Total CEV Self-employed(%%): %f \n',cev_total_se_p)
fprintf('Total CEV Workers(%%): %f \n',cev_total_w_p)
fprintf('\n')
disp('FIGURE 4')
fprintf('CEV(%%) Self-employed: \n')
for i=1:length(cev_se_p)
    fprintf(['Decile' num2str(i) ': %f \n'],cev_se_p(i))
end
fprintf('\n')
fprintf('CEV(%%) Workers: \n')
for i=1:length(cev_se_p)
    fprintf(['Decile' num2str(i) ': %f \n'],cev_w_p(i))
end

%% Plots

FS1 = 12;
fileExt = '-depsc'; % '-dpng','-depsc', File extension: png or eps

%%% NOTE THAT THE FIGURES ARE SAVED WITH THE SAME NAME THEY HAVE ON OUR
%%% OLD DRAFT. To reproduce Figure 4, see also the file
%%% welfare_new_figure4.m

if ~exist('cev_definition','var')
    error('You MUST specify cev_definition as either CE1 or CE2 (see welfare memo)!')
else
    cev_fig = ['_' cev_definition];
end

if do_welfare_plots == 1
    if strcmp(experiment,'no_redistrib')
        figure(1)
        bar([cev_se_p,cev_w_p],1,'grouped');
        legend('Self-employed','Workers','location','northwest','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_nodistr_c' cev_fig],fileExt)
        
    elseif strcmp(experiment,'no_redistrib_pe')
        figure(1)
        bar([cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','northwest','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_nodistr_pe_c' cev_fig],fileExt)
        
    elseif strcmp(experiment,'lump_sum')
        figure(1)
        bar( [cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','best','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_lump_c' cev_fig],fileExt)
    elseif strcmp(experiment,'lump_sum_pe')
        figure(1)
        bar( [cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','best','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_lump_pe_c' cev_fig],fileExt)
    elseif strcmp(experiment,'cut_tax_all')
        figure(1)
        bar( [cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','northwest','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_tax_all_c' cev_fig],fileExt)
    elseif strcmp(experiment,'cut_tax_all_pe')
        figure(1)
        bar( [cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','northwest','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_tax_all_pe_c' cev_fig],fileExt)
    elseif strcmp(experiment,'cut_tax_se')
        figure(1)
        bar( [cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','northwest','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        print([SaveDir 'welfare_hist_tax_entre_c' cev_fig],fileExt)
    elseif strcmp(experiment,'cut_tax_se_pe')
        figure(1)
        bar( [cev_se_p, cev_w_p],1,'grouped' );
        legend('Self-employed','Workers','location','southeast','FontSize',FS1)
        %mycolor = [0 0 1 ; 1 0 0 ];
        %colormap(mycolor)
        %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
        xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
        ylabel('CEV(%)')
        %xlim([0 20])
        ylim([-10 10])
        set(gca,'FontSize',FS1);
        %welfare_hist_tax_entre_pe_c_CE1
        print([SaveDir 'welfare_hist_tax_entre_pe_c' cev_fig],fileExt)
    else
        error('experiment: Invalid option')
    end
end %END flag do plots
