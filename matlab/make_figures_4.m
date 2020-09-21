%% Figure (4): Tax Evasion and Welfare Across Wealth and Occupation
% Note: updated by Alessandro Di Nola on 15 September 2020

clear;clc;close all

%% Set up useful paths for folders

%Specify here where results (.mat files) are saved
ResultsDir = 'results\mat\';
%Specify here where you want to save the welfare_results/figures/tables
SaveDir = 'results\Figures\';
disp(['Saving figures in subfolder: ' SaveDir])

%% Choose experiment
experiment     = 'cut_tax_se' ;    %Options:
                                    %'no_redistrib'
                                    %'lump_sum'
                                    %'cut_tax_all'
                                    %'cut_tax_se'

%See Ale memo for explanation
disp(['Plot results of experiment: ' experiment])

%% Load welfare results in consumption equivalent units

if strcmp(experiment,'no_redistrib')
    load([ResultsDir 'welfare_results_no_redistrib.mat'],'cev_total','cev_total_se','cev_total_w','cev_se','cev_w')
elseif strcmp(experiment,'lump_sum')
    load([ResultsDir 'welfare_results_lump_sum.mat'],'cev_total','cev_total_se','cev_total_w','cev_se','cev_w')
elseif strcmp(experiment,'cut_tax_all')
    load([ResultsDir 'welfare_results_cut_tax_all.mat'],'cev_total','cev_total_se','cev_total_w','cev_se','cev_w')
elseif strcmp(experiment,'cut_tax_se')
    load([ResultsDir 'welfare_results_cut_tax_se.mat'],'cev_total','cev_total_se','cev_total_w','cev_se','cev_w')
else
    error('Invalid option')
end

%% Convert to percentages

%These are (10*1) vectors
cev_se_p = cev_se*100;
cev_w_p  = cev_w*100;

%% Draw Figure (4)

FS1     = 12; %fontsize of legend. Open issue: if latex interpreter is used, 
              %then fontsize and FontName are turned off!!!
fileExt = '-depsc'; % '-dpng','-depsc', File extension: png or eps

if strcmp(experiment,'no_redistrib')
    %% Figure (4a)
    figure(1)
    b = bar([cev_se_p,cev_w_p],1,'grouped');
    b(1).FaceColor = [0 0 0]; %black
    b(2).FaceColor = [0.8 0.8 0.8]; %grey
    legend('Self-employed','Workers','location','northwest','FontSize',FS1)
    %mycolor = [0 0 1 ; 1 0 0 ];
    %colormap(mycolor)
    %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
    %xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
    xlabel('Deciles of wealth','interpreter','latex')
    %ylabel('Consumption equivalent units , (\%)','interpreter','latex')
    ylabel('Consumption equivalent units (\%)','interpreter','latex')
    %xlim([0 20])
    ylim([-10 10])
    set(gca,'FontSize',FS1);
    print([SaveDir 'fig4a'],fileExt)
    
elseif strcmp(experiment,'lump_sum')
    %% Figure (4b)
    figure(1)
    b = bar([cev_se_p,cev_w_p],1,'grouped');
    b(1).FaceColor = [0 0 0];
    b(2).FaceColor = [0.8 0.8 0.8];
    legend('Self-employed','Workers','location','best','FontSize',FS1)
    %mycolor = [0 0 1 ; 1 0 0 ];
    %colormap(mycolor)
    %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
    %xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
    xlabel('Deciles of wealth','interpreter','latex')
    %ylabel('CEV(%)')
    ylabel('Consumption equivalent units (\%)','interpreter','latex')
    %xlim([0 20])
    ylim([-10 10])
    set(gca,'FontSize',FS1);
    print([SaveDir 'fig4b'],fileExt)
    
elseif strcmp(experiment,'cut_tax_all')
    %% Figure (4c)
    figure(1)
    b = bar([cev_se_p,cev_w_p],1,'grouped');
    b(1).FaceColor = [0 0 0];
    b(2).FaceColor = [0.8 0.8 0.8];
    legend('Self-employed','Workers','location','northwest','FontSize',FS1)
    %mycolor = [0 0 1 ; 1 0 0 ];
    %colormap(mycolor)
    %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
    %xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
    xlabel('Deciles of wealth','interpreter','latex')
    %ylabel('CEV(%)')
    ylabel('Consumption equivalent units (\%)','interpreter','latex')
    %xlim([0 20])
    ylim([-10 10])
    set(gca,'FontSize',FS1);
    print([SaveDir 'fig4c'],fileExt)

elseif strcmp(experiment,'cut_tax_se')
    %% Figure (4d)
    figure(1)
    b = bar([cev_se_p,cev_w_p],1,'grouped');
    b(1).FaceColor = [0 0 0];
    b(2).FaceColor = [0.8 0.8 0.8];
    legend('Self-employed','Workers','location','best','FontSize',FS1)
    %ylabel('CEV, in %','FontSize',FS1,'FontName','Times New Roman')
    %xlabel('Deciles of wealth','FontSize',FS1,'FontName','Times New Roman')
    xlabel('Deciles of wealth','interpreter','latex')
    %ylabel('CEV(%)')
    ylabel('Consumption equivalent units (\%)','interpreter','latex')
    %xlim([0 20])
    ylim([-10 10])
    set(gca,'FontSize',FS1);
    print([SaveDir 'fig4d'],fileExt)
    
else
    error('experiment: Invalid option')
end








