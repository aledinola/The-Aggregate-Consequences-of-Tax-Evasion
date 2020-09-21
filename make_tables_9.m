%% This script generates Table 9, Tax Evasion and Aggregate Welfare
%Notes: Updated by Alesandro Di Nola on 15 September 2020. 

clear;clc;close all

%%%% HERE Specify here where results (.mat files) are saved %%%
ResultsDir = 'results\mat\';
%%%% HERE Specify here where you want to save the figures/tables %%%
SaveDir = 'results\Tables\';

%% Load mat files
    %Switch to generate a complete LaTex document or just a table to import in draft:
    makeCompleteLatexDocument = 1; %either 0 or 1
    xx = NaN;
    
    %BENCHMARK
    load([ResultsDir 'taxevasion_ge.mat'],'ModelResults')
    M1 = ModelResults; %matlab structure with results from TE economy
    %M1 is a struct with results for column (1) of Table (7)
    
    
    %COUNTERFACTUAL NO REDISTRIB, column (2)
    load([ResultsDir 'notaxevasion_ge.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_no_redistrib.mat'],'cev_total','cev_total_se','cev_total_w')
    M2 = ModelResults; %matlab structure with results from No TE economy
    C2 = [cev_total;cev_total_se;cev_total_w]*100;
    
    %COUNTERFACTUAL LUMP SUM REDISTRIB, Partial Equilibrium, column (3)
    load([ResultsDir 'notaxevasion_pe_lump_sum.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_lump_sum_pe.mat'],'cev_total','cev_total_se','cev_total_w')
    M3 = ModelResults;
    C3 = [cev_total;cev_total_se;cev_total_w]*100;
    
    %COUNTERFACTUAL LUMP SUM REDISTRIB, General Equilibrium, column (4)
    load([ResultsDir 'notaxevasion_ge_lump_sum.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_lump_sum.mat'],'cev_total','cev_total_se','cev_total_w')
    M4 = ModelResults;
    C4 = [cev_total;cev_total_se;cev_total_w]*100;
    
    %COUNTERFACTUAL TAX CUT ALL, Partial Equilibrium, column (5)
    load([ResultsDir 'notaxevasion_pe_cut_tax_all.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_cut_tax_all_pe.mat'],'cev_total','cev_total_se','cev_total_w')
    M5 = ModelResults;
    C5 = [cev_total;cev_total_se;cev_total_w]*100;
    
    %COUNTERFACTUAL TAX CUT ALL, General Equilibrium, column (6)
    load([ResultsDir 'notaxevasion_ge_cut_tax_all.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_cut_tax_all.mat'],'cev_total','cev_total_se','cev_total_w')
    M6 = ModelResults;
    C6 = [cev_total;cev_total_se;cev_total_w]*100;
    
    %COUNTERFACTUAL TAX CUT SE, Partial Equilibrium, column (7)
    load([ResultsDir 'notaxevasion_pe_cut_tax_se.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_cut_tax_se_pe.mat'],'cev_total','cev_total_se','cev_total_w')
    M7 = ModelResults;
    C7 = [cev_total;cev_total_se;cev_total_w]*100;
    
    %COUNTERFACTUAL TAX CUT SE, General Equilibrium, column (8)
    load([ResultsDir 'notaxevasion_ge_cut_tax_se.mat'],'ModelResults')
    load([ResultsDir 'welfare_results_cut_tax_se.mat'],'cev_total','cev_total_se','cev_total_w')
    M8 = ModelResults;
    C8 = [cev_total;cev_total_se;cev_total_w]*100;
    
%% Table 7, section 5.3, old draft
fid = fopen([SaveDir 'table9.tex'], 'wt');
if makeCompleteLatexDocument
    fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
    fprintf(fid,'\\usepackage{rotating} \n');
    fprintf(fid,'\\usepackage{floatrow} \n');
    fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{sidewaystable}[htbp] \\centering \n');
fprintf(fid,'\\caption{Tax Evasion and Aggregate Welfare} \n');
fprintf(fid,'\\label{table:aggregates} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{l c c c c c c c c} \\hline \n');
fprintf(fid, '& Tax Evasion  & No Redistr  & Lump sum PE & Lump sum GE &  Tax Cut PE & Tax Cut GE & Tax Cut PE & Tax Cut GE \\\\ \n');
fprintf(fid, '& Benchmark    &             & All         & All         &  All        & All        & SE      & SE            \\\\\ \n');
fprintf(fid, '& (1)          &     (2)     & (3)         & (4)         & (5)         & (6)        & (7)     & (8)           \\\\\ \n');
fprintf(fid, '\\hline \\hline \n');
                                                                                                         % Bench(1)       %NoRedistr (2)   %Lump-sum,PE(3)  %Lump-sum,GE(4)  %TaxCutAll,PE(5) %TaxCutAll,GE(6) %TaxCutSE,PE(7)  %TaxCutSE,GE(8)            
fprintf(fid,'Share of self-employed (\\%%) & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',M1.share_entre_p,M2.share_entre_p,M3.share_entre_p,M4.share_entre_p,M5.share_entre_p,M6.share_entre_p,M7.share_entre_p,M8.share_entre_p); 
fprintf(fid,'$Y$                & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',M1.output,       M2.output,       M3.output,       M4.output,       M5.output,       M6.output,       M7.output,       M8.output);
fprintf(fid,'$r(\\%%)$          & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',M1.r0*100,       M2.r0*100,       M3.r0*100,       M4.r0*100,       M5.r0*100,       M6.r0*100,       M7.r0*100,       M8.r0*100);
fprintf(fid,'$w$                & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',M1.w0,           M2.w0,           M3.w0,           M4.w0,           M5.w0,           M6.w0,           M7.w0,           M8.w0);
fprintf(fid,'$CEV total(\\%%)$  & %s    & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n','N.A.',          C2(1),           C3(1),           C4(1),           C5(1),           C6(1),           C7(1),           C8(1)); 
fprintf(fid,'$CEV SE(\\%%)$     & %s    & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n','N.A.',          C2(2),           C3(2),           C4(2),           C5(2),           C6(2),           C7(2),           C8(2)); 
fprintf(fid,'$CEV Work(\\%%)$   & %s    & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n','N.A.',          C2(3),           C3(3),           C4(3),           C5(3),           C6(3),           C7(3),           C8(3)); 

    
fprintf(fid,'\\hline  \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\end{sidewaystable} \n');
if makeCompleteLatexDocument; fprintf(fid,'\\end{document} \n'); end
fclose(fid);
