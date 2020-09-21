%% This script generates table(s) related to decompositions

clear;clc;close all

%Set path for directories
ResultsDir = 'results\'; %folder where .mat files are stored
SaveDir    = 'results\'; %Specify here where you want to save the tex files

disp(['Writing latex tables in subfolder: ' SaveDir])

%Switch to generate a complete LaTex document or just a table to import in draft:
makeCompleteLatexDocument = 1; %either 0 or 1

%% Load (1:6) results into structure array
ModRes(6) = struct(); %Initialize a structure _array_

%ModelResults.cond_mean_theta
%ModelResults.cond_firm_size_dist

%(1) no tax evasion GE
load([ResultsDir 'decomp_col_1.mat'],'ModelResults','inc_gap_all_p')
ModRes(1).ModelResults = ModelResults; 
ModRes(1).inc_gap_all_p = inc_gap_all_p;

%(2) fixed o(x),k(x),n(x),prices
load([ResultsDir 'decomp_col_2.mat'],'ModelResults','inc_gap_all_p')
ModRes(2).ModelResults = ModelResults; 
ModRes(2).inc_gap_all_p = inc_gap_all_p;


%(3) fixed o(x),prices
load([ResultsDir 'decomp_col_3.mat'],'ModelResults','inc_gap_all_p')
ModRes(3).ModelResults = ModelResults; %matlab structure with results 
ModRes(3).inc_gap_all_p = inc_gap_all_p;


%(4) fixed k(x),n(x),prices
load([ResultsDir 'decomp_col_4.mat'],'ModelResults','inc_gap_all_p')
ModRes(4).ModelResults = ModelResults; %matlab structure with results 
ModRes(4).inc_gap_all_p = inc_gap_all_p;


%(5) fixed prices
load([ResultsDir 'decomp_col_5.mat'],'ModelResults','inc_gap_all_p')
ModRes(5).ModelResults = ModelResults; 
ModRes(5).inc_gap_all_p = inc_gap_all_p;


%(6) Benchmark (with taxev) GE
load([ResultsDir 'decomp_col_6.mat'],'ModelResults','inc_gap_all_p')
ModRes(6).ModelResults = ModelResults; 
ModRes(6).inc_gap_all_p = inc_gap_all_p;


%% Create Table

fid = fopen([SaveDir 'table6_decomp.tex'], 'wt');
if makeCompleteLatexDocument
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\usepackage{rotating} \n');
fprintf(fid,'\\usepackage{floatrow} \n');

fprintf(fid,'\\begin{document} \n');
end

fprintf(fid,'\\begin{sidewaystable}[htbp] \\centering \n');
fprintf(fid,'\\caption{Decomposition of Aggregate Effects\\label{tab:decomp}} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{l c c c c c c}  \n');
%fprintf(fid, '& Tax Evasion  & Perfect Tax      &  Change  \\\\ \n');
%fprintf(fid, '& (Benchmark)  & Enforcement &  (\\%%) \\\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid, '\\hline \n');
fprintf(fid,'	&  & \\multicolumn{4}{c}{Tax Evasion Economies}\\\\ \n');
fprintf(fid,'\\cline{3-7} \n');
fprintf(fid,'&  & \\multicolumn{4}{c}{} & \\\\ \n');
fprintf(fid,'& Perfect Tax  & \\multicolumn{4}{c}{Counterfactual Experiments} & Benchmark\\\\ \n');
fprintf(fid,'& Enforcement  & \\multicolumn{4}{c}{Partial Equilibrium} & General Equilibrium \\\\ \n');
fprintf(fid,'\\cline{3-7} \n');
fprintf(fid,'& (1) & (2) & (3) & (4) & (5)&(6) \\\\ \n');
fprintf(fid,'\\hline \\hline \n');
		
fprintf(fid,'\\emph{\\underline{Setup}} &  &  &  &  & & \\\\ \n');
fprintf(fid,'Fixed prices from (1)    & - & $r,w$ & $r,w$ & $r,w$ & $r,w$& -\\\\ \n');
fprintf(fid,'Fixed decisions from (1) & - & $o(x), k(x), n(x)$, & $o(x)$  & $k(x), n(x)$  &-    & -\\\\ \n');
fprintf(fid,'Operational channels & - & \\textit{Subsidy}     & \\textit{Subsidy}      & \\textit{Subsidy}      &\\textit{Subsidy+Detection}& \\textit{Subsidy+Detection}\\\\ \n');
fprintf(fid,'		&   &            &\\textit{+Detection}    &\\textit{+Selection}  &\\textit{+Selection}      & \\textit{+Selection+Prices}  \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid, '\\textit{\\underline{Outcomes}} &  & & & &  &\\\\  \n');
fprintf(fid,'Share of self-employed (\\%%) & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.share_entre_p,ModRes(2).ModelResults.share_entre_p,ModRes(3).ModelResults.share_entre_p,ModRes(4).ModelResults.share_entre_p,ModRes(5).ModelResults.share_entre_p,ModRes(6).ModelResults.share_entre_p);
fprintf(fid,'$E(\\theta|E)$ & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.cond_mean_theta,ModRes(2).ModelResults.cond_mean_theta,ModRes(3).ModelResults.cond_mean_theta,ModRes(4).ModelResults.cond_mean_theta,ModRes(5).ModelResults.cond_mean_theta,ModRes(6).ModelResults.cond_mean_theta);
fprintf(fid,'$E(k|E)$       & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.ave_k_entre,ModRes(2).ModelResults.ave_k_entre,ModRes(3).ModelResults.ave_k_entre,ModRes(4).ModelResults.ave_k_entre,ModRes(5).ModelResults.ave_k_entre,ModRes(6).ModelResults.ave_k_entre);
fprintf(fid,'$K^{E}$        & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.k_entre,ModRes(2).ModelResults.k_entre,ModRes(3).ModelResults.k_entre,ModRes(4).ModelResults.k_entre,ModRes(5).ModelResults.k_entre,ModRes(6).ModelResults.k_entre);
fprintf(fid,'$K^{C}$        & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.k_corp,ModRes(2).ModelResults.k_corp,ModRes(3).ModelResults.k_corp,ModRes(4).ModelResults.k_corp,ModRes(5).ModelResults.k_corp,ModRes(6).ModelResults.k_corp);
fprintf(fid,'$Y^{E}$        & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.output_se,ModRes(2).ModelResults.output_se,ModRes(3).ModelResults.output_se,ModRes(4).ModelResults.output_se,ModRes(5).ModelResults.output_se,ModRes(6).ModelResults.output_se);
fprintf(fid,'$Y^{C}$        & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.output_corp,ModRes(2).ModelResults.output_corp,ModRes(3).ModelResults.output_corp,ModRes(4).ModelResults.output_corp,ModRes(5).ModelResults.output_corp,ModRes(6).ModelResults.output_corp);
fprintf(fid,'Small firms (\\%%) & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).ModelResults.cond_firm_size_dist(1),ModRes(2).ModelResults.cond_firm_size_dist(1),ModRes(3).ModelResults.cond_firm_size_dist(1),ModRes(4).ModelResults.cond_firm_size_dist(1),ModRes(5).ModelResults.cond_firm_size_dist(1),ModRes(6).ModelResults.cond_firm_size_dist(1));
fprintf(fid,'Misreporting rate (\\%%) & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',ModRes(1).inc_gap_all_p,ModRes(2).inc_gap_all_p,ModRes(3).inc_gap_all_p,ModRes(4).inc_gap_all_p,ModRes(5).inc_gap_all_p,ModRes(6).inc_gap_all_p);


fprintf(fid,'\\hline  \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\floatfoot{Note: The outcome Small firms (\\%%) refers to the share of firms with [1-4] employees, denoted as Firm\\_size\\_hiring\\_b2 in the excel file.} \n');
fprintf(fid,'\\end{sidewaystable} \n');

if makeCompleteLatexDocument; fprintf(fid,'\\end{document} \n'); end
fclose(fid);

