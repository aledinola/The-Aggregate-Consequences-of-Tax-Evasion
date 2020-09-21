%% This script generates Figures (1-2-3) and Figure (7) in Appendix B1
%Notes: Updated by Alesandro Di Nola on September 15, 2020. 

clear;clc;close all

%% Set up some useful paths
ResultsDir = 'results\mat\'; %folder where .mat files are stored
SaveDir = 'results\Figures\'; %Specify here where you want to save the figures
disp(['Saving figures in subfolder: ' SaveDir])

%% Load results from benchmark economy (tax evasion with GE)

load([ResultsDir 'taxevasion_ge.mat'])

%% Unpack structures

nagrid = Parameters.nagrid;
negrid = Parameters.negrid;
ntgrid = Parameters.ntgrid;


%% Options for plotting

ldw     = 2; %line width
fontw   = 12; %font size for xlabel, legend, etc.
fileExt = '-depsc'; % '-depsc', File extension: png or eps
perc    = 100/3; %convert [0,3] into [0% 100%] for hours worked


%% Figure 1: Probability of auditing

pk = prob_audit([],[],[],kgrid,[],pflag,pn_1,pn_2,pn_3,[],[],[]);

klim = find(kgrid>50, 1 ); %this is the index on kgrid

figure(1)
plot(kgrid(1:klim),pk(1:klim),'b','linewidth',ldw)
xlabel('Capital, k')
ylabel('Prob. of auditing')
ylim([0 1])
axis tight
print([SaveDir,'fig1'],fileExt)

%% Figure (2a): Plots "Tax evasion by total income" 


figure(2)
plot(1:5,taxev_quint_inc_p,'b-o','linewidth',ldw)
hold on
plot(1:5,DataTargets.taxev_quint_inc_p,'r--o','linewidth',ldw)
legend('Model','Data','location','northwest')
xlabel('Quintiles of Total Income')
ylabel('Misreported Percentage')
grid on
xticks(1:5)
ylim([0,30])
print([SaveDir,'fig2a'],fileExt)
hold off

%% %% Figure (2b): Plots "Share of self-employed by total income"

figure(3)
plot(1:5,share_entre_quint_p,'b-o','linewidth',ldw)
hold on
plot(1:5,DataTargets.share_entre_quint_p,'r--o','linewidth',ldw)
legend('Model','Data','location','northwest')
xlabel('Quintiles of Total Income')
ylabel('Share of Self-employed')
grid on
xticks(1:5)
ylim([0,40])
print([SaveDir,'fig2b'],fileExt)
hold off

%% Figure (2c): Plots "Expected Probability of auditing by self-employed income"
%expected pk by business income
 
up = 40;

figure(5)
plot(1:5,audit_businc_p,'b-o','linewidth',ldw)
xlabel('Quintiles of Business Income \pi')
ylabel('Prob. of auditing')
grid on
xticks(1:5)
ylim([0,up])
print([SaveDir,'fig2c'],fileExt)

%% Load data for comparison: Tax Evasion vs Perfect Tax Enforcement 

%Load tax evasion
load([ResultsDir 'taxevasion_ge.mat'])
polTE     = policy; %matlab structure with policy funct 
distribTE = distrib;
ModResTE  = ModelResults;

clear policy distrib ModelResults

%Load no tax evasion general equilibrium
load([ResultsDir 'notaxevasion_ge.mat'])
polNoTE     = policy; %matlab structure with policy funct 
distribNoTE = distrib;
ModResNoTE  = ModelResults;

clear policy distrib ModelResults

epsilon_given  = 6;%round(negrid/2);
a_ub = 27; %upper bound for plot


%% Figure (3a): Occupational choice

abar_occpol_te   = abar_fun(polTE.occpol,nagrid,negrid,ntgrid,agrid,epsilon_given);
abar_occpol_note = abar_fun(polNoTE.occpol,nagrid,negrid,ntgrid,agrid,epsilon_given);

figure(8) % occ choice, in the (ability,wealth) space
plot(theta,log(abar_occpol_te),'b',theta,log(abar_occpol_note),'r--','linewidth',ldw)
hold on
grid on
xlabel('Business ability, \theta','FontSize',fontw)
ylabel('Assets (log), a','FontSize',fontw)
text(theta(10),2.5,'Self-employed','FontSize',fontw)
text(theta(3),-0.5,'Workers','FontSize',fontw)
legend('Tax Evasion', 'Perfect Enforcement', 'FontSize',fontw)
xlim([min(theta), max(theta)])
hold off
print([SaveDir 'fig3a'],fileExt)


%% Plotting options for theta,eps and assets upper bound

%fix shock eps at median value
eps_plot                 = min(round(negrid/2),length(eps_grid)); 
%upper bound for assets to plot
assets_upper_bound_value = 40; 
assets_upper_bound_index = find(agrid>=assets_upper_bound_value,1);
al_k                     = min(assets_upper_bound_index,length(agrid)); 

%% Figure(3b): Savings of Self-Employed, apol_te_note_cut_15

theta_ind4b = 9; %theta index to plot

figure(11)
plot(agrid_dist,agrid_dist,'k--','linewidth',ldw)
hold on
plot(agrid_dist,polTE.apolse0det(:,round(negrid/2),theta_ind4b),'b','linewidth',ldw)
hold on
plot(agrid_dist,polNoTE.apolse0det(:,round(negrid/2),theta_ind4b),'r--','linewidth',ldw)
hold on 
legend('45 line','Tax evasion','Perfect Enforcement','Location','best'), grid on
xlabel('Assets, a'),ylabel('Assets, a'''), axis tight
xlim([0 15])
print([SaveDir 'fig3b'],fileExt)
hold off

%% Figure (3c): plot Self-Employed Business Capital

theta_ind_k = 11; %11 theta index to plot !!!!

figure(12)
plot(agrid(1:al_k),lambda*agrid(1:al_k),'k--','linewidth',ldw)
hold on
plot(agrid(1:al_k),polTE.policycap(1:al_k,eps_plot,theta_ind_k),'b','linewidth',ldw)
hold on
plot(agrid(1:al_k),polNoTE.policycap(1:al_k,eps_plot,theta_ind_k ),'r--','linewidth',ldw)
legend('\lambdaa ','Tax Evasion', 'Perfect Enforcement','Location','best','FontSize', fontw)
grid on
xlabel('Assets, a','FontSize', fontw),ylabel('Capital, k','FontSize', fontw), axis tight
%ylim([0 50])
hold off
print([SaveDir 'fig3c'],fileExt)


%% Figure (3d): plot Hired Labor n(x)

al_n      = al_k; %upper bound for assets to plot, specific to npol
theta_ind = theta_ind_k; %theta index to plot !!!!

figure(13)
plot(agrid(1:al_n),polTE.policyn(1:al_n,eps_plot,theta_ind),'b','linewidth',ldw)
hold on
plot(agrid(1:al_n),polNoTE.policyn(1:al_n,eps_plot,theta_ind ),'r--','linewidth',ldw)
legend('Tax Evasion', 'Perfect Enforcement','Location','best','FontSize', fontw)
grid on
xlabel('Assets, a','FontSize', fontw),ylabel('Hired labor, n','FontSize', fontw)
axis tight
%ylim([0 50])
hold off
print([SaveDir 'fig3d'],fileExt)  

%% Figure (3e): plot Workers' Labor Supply

al        = al_k; %upper bound for assets to plot
theta_ind = 12; %theta index to plot !!!!

figure(16)
plot(agrid(1:al),perc*polTE.lpolw(1:al,eps_plot,theta_ind),'b','linewidth',ldw)
hold on
plot(agrid(1:al),perc*polNoTE.lpolw(1:al,eps_plot,theta_ind ),'r--','linewidth',ldw)
legend('Tax Evasion', 'Perfect Enforcement','Location','best','FontSize', fontw), grid on
xlabel('Assets, a','FontSize', fontw),ylabel('Hours worked (%)','FontSize', fontw), axis tight
ylim([perc*0, perc*3])
print([SaveDir 'fig3e'],fileExt)  
hold off

%% Figure (7), Appendix B.1: kpol with unconstrained capital

al_k1 = min(28,length(agrid)); %upper bound for assets to plot, specific to kpol
theta_ind_k = 11; %theta index to plot !!!!

figure(22)
plot(agrid(1:al_k1),lambda*agrid(1:al_k1),'k--','linewidth',ldw)
hold on
plot(agrid(1:al_k1),polTE.policycap(1:al_k1,eps_plot,theta_ind_k),'b','linewidth',ldw)
hold on
plot(agrid(1:al_k1),polNoTE.policycap(1:al_k1,eps_plot,theta_ind_k ),'r--','linewidth',ldw)
legend('\lambdaa ','Tax Evasion', 'Perfect Enforcement','Location','best','FontSize', fontw)
grid on
xlabel('Assets, a','FontSize', fontw),ylabel('Capital, k','FontSize', fontw), axis tight
%ylim([0 50])
hold off
print([SaveDir 'fig7_appendix'],fileExt)


%% Figure (3f): misreporting rate, policy_te_SE

% lth = 2;  %index for 'low theta'
% mth = 9;  %index for 'mid. theta'
% hth = 12; %index for 'high theta'
% al_f = al_k;%19; %upper bound for assets to plot, specific to policyphi
% eps_fix = 7;
% 
% figure(14)
% plot(agrid(1:al_f), 100*polTE.policyphi(1:al_f,eps_fix,lth),':','linewidth',ldw)
% hold on
% plot(agrid(1:al_f), 100*polTE.policyphi(1:al_f,eps_fix,mth),'r--','linewidth',ldw)
% hold on
% plot(agrid(1:al_f), 100*polTE.policyphi(1:al_f,eps_fix,hth),'b','linewidth',ldw)
% legend('high \theta','low \theta','mid. \theta', 'FontSize',fontw)
% xlabel('Assets, a','FontSize',fontw)
% ylabel('Misreporting rate (%)','FontSize',fontw)
% axis tight
% %xlim([0 25])
% hold off
% %print([SaveDir 'policy_te_SE'],fileExt)

%--------------- SUBFUNCTIONS --------------------------------------------%
function [abar_occpol] = abar_fun(occpol,nagrid,negrid,ntgrid,agrid,epsilon_given)
%{
DESCRIPTION:
This function computes the wealth threshold a* such that if a>=a*,
occupational choice = 1, 0 otherwise. The wealth threshold is typically a
decreasing function of entrepreneurial ability theta.
%}

abarind        = zeros(ntgrid,1);
abar_occpol    = zeros(ntgrid,1);
for it=1:ntgrid
    aux = find(occpol(:,epsilon_given,it)==1);
    if isempty(aux)
        abar_occpol_ind(it) = nagrid;
        abar_occpol(it) = agrid(nagrid);
    else
        abar_occpol_ind(it) = min(aux);
        abar_occpol(it) = agrid(abar_occpol_ind(it));
    end
end

end %END FUNCTION "abar_fun"
