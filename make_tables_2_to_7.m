%% This script generates Tables 2-3-4-5-6-7 in the paper

%{
The data files to generate these tables are stored in the folder 
"results\mat". The tables are then saved in "results\Tables"
Notes: Updated by Alessandro Di Nola on 14 September 2020
%}

clear;clc;close all

%% Set up paths for results and save folders
ResultsDir = 'results\mat\'; %folder where .mat files are stored
                                % MUST end with \
SaveDir = 'results\Tables\'; %Specify here where you want to save the tex files
                             % MUST end with \
disp(['Writing latex tables in subfolder: ' SaveDir])

%% Load results from benchmark economy (tax evasion with GE)
load([ResultsDir 'taxevasion_ge.mat']) %this .mat file is saved at the end 
                                       %of "f_solve_model"

%Switch to generate a complete LaTex document or just a table to import in draft:
makeCompleteLatexDocument = 1; %either 0 or 1

%% Table 2: Internally Calibrated Parameters

fid = fopen([SaveDir 'table2.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Internally Calibrated Parameters} \n');
fprintf(fid,'\\label{table:parameters_inside} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{ccll} \\hline \n');
fprintf(fid, 'Parameter & Description & Value & Target \\\\ \n');
fprintf(fid, '\\hline \\hline \n');
fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,'$\\beta$ & Discount factor & %8.4f & 4\\%% Interest rate \\\\ \n',Parameters.beta);
fprintf(fid,'$\\psi$ & Disutility from working & %8.3f & Hours worked \\\\ \n',Parameters.psi);

fprintf(fid, '\\textit{\\underline{Production}} &  & &  \\\\  \n');
fprintf(fid,'$\\delta$ & Capital depreciation & %8.3f & Capital-output ratio \\\\ \n',Parameters.delta);
fprintf(fid,'$\\nu$ & Span of control & %8.3f & Firm Size \\\\ \n', Parameters.vi);
fprintf(fid,'$\\gamma$ & Capital share, self-employed & %8.3f & Share of hiring, self-employed \\\\ \n', Parameters.gamma);
fprintf(fid,'$\\lambda$ & Leverage ratio & %8.3f & Leverage of self-employed \\\\ \n', Parameters.lambda);

fprintf(fid, '\\textit{\\underline{Self-employed ability}} &  & &  \\\\  \n');
fprintf(fid,'$\\rho_{\\theta}$   & Persistence        & %8.3f & Exit rate, SE to W \\\\ \n' ,Parameters.rho_theta);
fprintf(fid,'$\\sigma_{\\theta}$ & Standard deviation & %8.3f & Gini income, SE \\\\ \n', Parameters.sigmaeps_theta);
fprintf(fid,'$\\mu_{\\theta}     $ & Unconditional mean & %8.3f & Share, SE \\\\ [0.5ex] \n'  , Parameters.uncmean_theta);

fprintf(fid, '\\textit{\\underline{Tax evasion detection}} &  & &  \\\\  \n');
fprintf(fid,'$\\kappa$ & Cost ot tax evasion & %8.4f & Misreporting rate \\\\ \n', Parameters.cc0);
fprintf(fid,'$p_{1}$ & Parameter of $p(k)$   & %8.3f & Tax evasion by income \\\\ \n', Parameters.pn_1);
fprintf(fid,'$p_{2}$ & Parameter of $p(k)$   & %8.3f & Tax evasion by income \\\\ \n', Parameters.pn_2);
%fprintf(fid,'$p_{3}$ & Parameter of $p(k)$  & %8.3f & Tax evasion by income \\\\ \n', Parameters.pn_3);

%fprintf(fid,'$\\kappa_{1}$ & Parameter of $\\kappa(\\phi)$  & %8.3f & Tax evasion by income \\\\ \n', Parameters.cc1);
%fprintf(fid,'$\\kappa_{2}$ & Parameter of $\\kappa(\\phi)$  & %8.3f & Tax evasion by income \\\\ \n', Parameters.cc2);
fprintf(fid, '\\textit{\\underline{Tax functions rescale}} &  & &  \\\\  \n');
fprintf(fid,'$\\chi$  & Rescaling parameter       & %8.3f & Tax revenues as share of GDP \\\\ \n', Parameters.ksi_tax);

fprintf(fid,'\\hline \\hline \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);


%% Table 3: Basic Model Statistics

maxLwork = Parameters.maxLwork;

fid = fopen([SaveDir 'table3.tex'], 'wt');
if makeCompleteLatexDocument==1
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[tbp] \n');
fprintf(fid,'\\caption{Basic Model Statistics \\label{tab:model_fit}} \n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\begin{tabular}{lcc} \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'      & Data  & Model \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\hline \n');
fprintf(fid, '\\textit{\\underline{Targeted Statistics}} &       &      \\\\  \n');
fprintf(fid,'Interest rate (\\%%)                  & %8.3f & %8.3f \\\\ \n',4 ,r0*100 );
fprintf(fid,'Capital-output ratio                  & %8.3f & %8.3f \\\\ \n',DataTargets.K_to_Y ,K_to_Y);
fprintf(fid,'Hours worked (\\%%)                   & %8.3f & %8.3f \\\\ \n',33 ,(ave_n_work/maxLwork)*100 );
fprintf(fid,'Share of self-employed (\\%%)         & %8.3f & %8.3f \\\\ \n',DataTargets.share_entre_p,share_entre_p );
fprintf(fid,'Mean leverage of self-employed        & %8.3f & %8.3f \\\\ \n',100*0.289,100*leverage_mean);
fprintf(fid,'Exit rate, self-employed  (\\%%)      & %8.3f & %8.3f \\\\ \n',DataTargets.exit_ew_p,exit_ew_p );
fprintf(fid,'Misreporting rate   (\\%%)            & %8.3f & %8.3f \\\\ \n',DataTargets.inc_gap_all_p,inc_gap_all_p );
fprintf(fid,'Tax revenues/GDP    (\\%%)            & %8.3f & %8.3f \\\\ \n',DataTargets.taxes_total_toY_p,taxes_total_toY_p);

fprintf(fid, '\\textit{\\underline{Non-Targeted Statistics}} &       &      \\\\  \n');
fprintf(fid,'Share of income, self-employed (\\%%)           & %8.3f & %8.3f \\\\ \n',DataTargets.share_inc_entre_p ,share_inc_entre_true_p );
fprintf(fid,'Median wealth ratio, SE to workers              & %8.3f & %8.3f \\\\ \n',4.02,wealth_med_ratio_E_W );
fprintf(fid,'Share of credit-constrained self-employed       & %8.3f & %8.3f \\\\ \n',22.80,100*frac_bor1 );
%fprintf(fid,'Share of assets, self-employed (\\%%) & %8.3f & %8.3f \\\\ \n',DataTargets.assets_e_to_work_p,assets_e_to_work_p );

fprintf(fid,'\\hline \n');
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\end{center} \n');
	
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

%% Table 4bis: Self-employment Firm Size Distribution

fid = fopen([SaveDir 'table4bis.tex'], 'wt');
if makeCompleteLatexDocument==1
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[tbp] \n');
fprintf(fid,'\\caption{Firm Size Distribution: Data and Model \\label{tab:firm_size_employees}} \n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\begin{tabular}{lcc} \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'Moments                   & Data  & Model \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'Share of hiring self-employed (\\%%) & %8.3f & %8.3f \\\\ \n',20.3 ,share_hiring_se_p );
fprintf(fid,'\\hline \n');
fprintf(fid,'1-4 employees (\\%%)    & %8.3f & %8.3f \\\\ \n',75.891,cond_firm_size_dist(1));
fprintf(fid,'5-9 employees           & %8.3f & %8.3f \\\\ \n',14.7  ,cond_firm_size_dist(2));
fprintf(fid,'10-19 employees         & %8.3f & %8.3f \\\\ \n',5.517 ,cond_firm_size_dist(3));
fprintf(fid,'More than 20 employees  & %8.3f & %8.3f \\\\ \n',3.891 ,cond_firm_size_dist(4));
fprintf(fid,'\\hline \n');
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\end{center} \n');
	
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

%% Table 5: Distribution of Self-employed Income

data_income_se = [43.9
1.43
15.702
50.696
36.097
9.776];

fid=fopen([SaveDir 'table5.tex'],'w');
if makeCompleteLatexDocument==1
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');
end
fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Distribution of Self-employed Income}');
fprintf(fid,'\n \\label{tab:income_se}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\begin{tabular}{lcccccc}');
fprintf(fid,'\n \\hline \n');
fprintf(fid,'\n & Gini & Mean/Median & Bottom 40 & Top 20 & Top 10 & Top 1  \\\\');
fprintf(fid,'\n \\hline \n');
fprintf(fid,'\n \\hline \n');
%fprintf(fid,'\n \\textit{\\underline{Income}}  &  &  &  &  & &  \\\\ ');
fprintf(fid,'\n  Model &  $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',[gini_incse_p; incse_mean_to_median; incse_perc_p(1:4)']);
fprintf(fid,'\n  US Data &  $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',data_income_se);
fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

%% Table 6: Wealth Distribution

fid=fopen([SaveDir 'table6.tex'],'w');
if makeCompleteLatexDocument==1
fprintf(fid,'\n \\documentclass[12pt]{article}');
fprintf(fid,'\n \\begin{document}');
end
fprintf(fid,'\n \\begin{table}[h]');
fprintf(fid,'\n \\caption{Wealth Distribution}');
fprintf(fid,'\n \\label{tab:wealth_inc_all}');
fprintf(fid,'\n \\begin{center}');
fprintf(fid,'\n \\begin{tabular}{lcccccc}');
fprintf(fid,'\n \\hline \n');
fprintf(fid,'\n & Gini & Mean/Median & Bottom 40 & Top 20 & Top 10 & Top 1  \\\\');
fprintf(fid,'\n \\hline \n');
fprintf(fid,'\n \\hline \n');
fprintf(fid,'\n \\textit{\\underline{Wealth}}  &  &  &  &  & &  \\\\ ');
fprintf(fid,'\n  Model &  $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',model_wealth_all(1:6));
fprintf(fid,'\n  US Data &  $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ & $%5.2f$ \\\\ ',data_wealth_all);
fprintf(fid,'\n \\hline');
fprintf(fid,'\n \\end{tabular}');
fprintf(fid,'\n \\end{center}');
fprintf(fid,'\n \\end{table}');

if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);


%% Table 7: Aggregate Effects of Tax Evasion 

xx = NaN;

load([ResultsDir 'taxevasion_ge.mat'])
ModResTE = ModelResults; %matlab structure with results from TE economy
                         %Similar to "Table results" in excel

load([ResultsDir 'notaxevasion_pe.mat'])
ModResNoTEpe = ModelResults; %matlab structure with results from No TE economy
                             %prices fixed at the benchmark economy


load([ResultsDir 'notaxevasion_ge.mat'])
ModResNoTE = ModelResults; %matlab structure with results from No TE economy
                           % General equilibrium (r,w clear ED)



fid = fopen([SaveDir 'table7.tex'], 'wt');
if makeCompleteLatexDocument==1
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Aggregate Effects of Tax Evasion} \n');
fprintf(fid,'\\label{table:aggregates} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{l c c c c} \\hline \n');
fprintf(fid, '& Tax Evasion  & Perfect Tax      & Perfect Tax &  Change  \\\\ \n');
fprintf(fid, '& Benchmark  & Enforcement (PE) & Enforcement &   (\\%%) \\\\\ \n');
fprintf(fid, '\\hline \\hline \n');

fprintf(fid, '\\textit{\\underline{Sector of self-employment}} &  & & &  \\\\  \n');
fprintf(fid,'Share of self-employed (\\%%) & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.share_entre_p,ModResNoTEpe.share_entre_p,ModResNoTE.share_entre_p,xx);
fprintf(fid,'$E(\\theta|E)$ & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModResTE.cond_mean_theta,ModResNoTEpe.cond_mean_theta,ModResNoTE.cond_mean_theta,xx);
fprintf(fid,'$E(k|E)$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.ave_k_entre,ModResNoTEpe.ave_k_entre,ModResNoTE.ave_k_entre,xx);
fprintf(fid,'$E(n|E)$       & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.ave_n_entre,ModResNoTEpe.ave_n_entre,ModResNoTE.ave_n_entre,xx);
fprintf(fid,'$K^{E}$        & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.k_entre,ModResNoTEpe.k_entre,ModResNoTE.k_entre,xx);
fprintf(fid,'$N^{E}$        & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.n_entre,ModResNoTEpe.n_entre,ModResNoTE.n_entre,xx); %!!!!!!!!!!
fprintf(fid,'$Y^{E}$        & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.output_se,ModResNoTEpe.output_se,ModResNoTE.output_se,xx);
fprintf(fid,'Share of credit-constrained SE (\\%%) & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*ModResTE.frac_bor1,100*ModResNoTEpe.frac_bor1,100*ModResNoTE.frac_bor1,xx);
%fprintf(fid,'(\\%%) SE borr. constr. & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',0.2084*100,0.2547*100,0.2442*100,xx);


fprintf(fid, '\\textit{\\underline{Corporate sector}} &  & & & \\\\  \n');
fprintf(fid,'$K^{C}$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.k_corp,ModResNoTEpe.k_corp,ModResNoTE.k_corp,xx);
fprintf(fid,'$N^C $    & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.n_corp,ModResNoTEpe.n_corp,ModResNoTE.n_corp,xx);
fprintf(fid,'$Y^{C}$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.output_corp,ModResNoTEpe.output_corp,ModResNoTE.output_corp,xx);

fprintf(fid, '\\textit{\\underline{Labor}} &  & & & \\\\  \n');
fprintf(fid,'Share of hiring SE (\\%%)  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.share_hiring_se_p,ModResNoTEpe.share_hiring_se_p,ModResNoTE.share_hiring_se_p,xx);
%fprintf(fid,'Share of labor hired by SE (\\%%)  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.share_workInSe*100,ModResNoTEpe.share_workInSe*100,ModResNoTE.share_workInSe*100,xx);
fprintf(fid,'Hours worked    & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',100*(ModResTE.ave_n_work/maxLwork),100*(ModResNoTEpe.ave_n_work/maxLwork),100*(ModResNoTE.ave_n_work/maxLwork),xx);

fprintf(fid, '\\textit{\\underline{Prices}} &  & & & \\\\  \n');
fprintf(fid,'$r (\\%%)$   & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.r0*100,ModResNoTEpe.r0*100,ModResNoTE.r0*100,xx);
fprintf(fid,'$w$    & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.w0,ModResNoTEpe.w0,ModResNoTE.w0,xx);

fprintf(fid, '\\textit{\\underline{Tax revenues}} &  & & & \\\\  \n');
fprintf(fid,'$T/Y$ (\\%%)  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.taxes_total_toY_p,ModResNoTEpe.taxes_total_toY_p,ModResNoTE.taxes_total_toY_p,xx);

fprintf(fid, '\\textit{\\underline{Wealth inequality}} &  & & & \\\\  \n');
fprintf(fid,'Gini, all  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.gini_wealth_p,ModResNoTEpe.gini_wealth_p,ModResNoTE.gini_wealth_p,xx);
fprintf(fid,'Gini, self-employed  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.gini_wealth_se_p,ModResNoTEpe.gini_wealth_se_p,ModResNoTE.gini_wealth_se_p,xx);
fprintf(fid,'Gini, workers  & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n',ModResTE.gini_wealth_work_p,ModResNoTEpe.gini_wealth_work_p,ModResNoTE.gini_wealth_work_p,xx);
fprintf(fid,'\\hline  \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);



%% Table 4: Self-employment Firm Size Distribution
%table4bis includes also counterfectual economy

fid = fopen([SaveDir 'table4.tex'], 'wt');
if makeCompleteLatexDocument==1
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[tbp] \n');
fprintf(fid,'\\caption{Self-Employed Firm Size Distribution \\label{tab:firm_size_employees}} \n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\begin{tabular}{lccc} \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'                  & Data  & Tax Evasion & Perfect Tax  \\\\ \n');
fprintf(fid,'                  &       & Benchmark   & Enforcement  \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'Share of hiring self-employed (\\%%) & %8.2f & %8.2f & %8.2f  \\\\ \n',20.3 ,ModResTE.share_hiring_se_p,ModResNoTE.share_hiring_se_p );
fprintf(fid,'\\hline \n');
fprintf(fid,'1-4 employees (\\%%)    & %8.2f & %8.2f & %8.2f  \\\\ \n',75.891,ModResTE.cond_firm_size_dist(1),ModResNoTE.cond_firm_size_dist(1));
fprintf(fid,'5-9 employees           & %8.2f & %8.2f & %8.2f  \\\\ \n',14.7  ,ModResTE.cond_firm_size_dist(2),ModResNoTE.cond_firm_size_dist(2));
fprintf(fid,'10-19 employees         & %8.2f & %8.2f & %8.2f  \\\\ \n',5.517 ,ModResTE.cond_firm_size_dist(3),ModResNoTE.cond_firm_size_dist(3));
fprintf(fid,'More than 20 Employees  & %8.2f & %8.2f & %8.2f  \\\\ \n',3.891 ,ModResTE.cond_firm_size_dist(4),ModResNoTE.cond_firm_size_dist(4));
fprintf(fid,'\\hline \n');
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\end{center} \n');
	
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);
