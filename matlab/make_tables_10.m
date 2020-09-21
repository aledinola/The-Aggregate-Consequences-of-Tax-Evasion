
clear;clc;close all


%% Set up paths for results and save folders
ResultsDir = 'results\lambda_mat\'; %folder where .mat files are stored
                               
SaveDir = 'results\Tables\'; %Specify here where you want to save the tex files
                             
disp(['Writing latex tables in subfolder: ' SaveDir])

%% Load results 
% Column 1
load([ResultsDir 'taxevasion_ge_lambda_1_2.mat']) 
M1 = ModelResults;
load([ResultsDir 'notaxevasion_pe_lambda_1_2.mat']) 
MN1 = ModelResults; 
maxLwork = Parameters.maxLwork;
load([ResultsDir 'welfare_results_no_redistrib_pe_1_2.mat']) 
W1 = [cev_total, cev_total_se, cev_total_w]';
W1 = W1*100; % persent
load ([ResultsDir 'welfare_results_lump_sum_pe_1_2.mat']) 
WL1 = [cev_total, cev_total_se, cev_total_w]';
WL1 = WL1*100;

% Column 2
load([ResultsDir 'taxevasion_ge.mat']) 
M2 = ModelResults;
load([ResultsDir 'notaxevasion_pe.mat']) 
MN2 = ModelResults; 
load([ResultsDir 'welfare_results_no_redistrib_pe.mat']) 
W2 = [cev_total, cev_total_se, cev_total_w]';
W2 = W2*100;
load ([ResultsDir 'welfare_results_lump_sum_pe.mat']) 
WL2 = [cev_total, cev_total_se, cev_total_w]';
WL2 = WL2*100;

% Column 3
load([ResultsDir 'taxevasion_ge_lambda_1_8.mat']) 
M3 = ModelResults;
load([ResultsDir 'notaxevasion_pe_lambda_1_8.mat']) 
MN3 = ModelResults; 
load([ResultsDir 'welfare_results_no_redistrib_pe_1_8.mat']) 
W3 = [cev_total, cev_total_se, cev_total_w]';
W3 = W3*100;
load ([ResultsDir 'welfare_results_lump_sum_pe_1_8.mat']) 
WL3 = [cev_total, cev_total_se, cev_total_w]';
WL3 = WL3*100;

% Column 4
load([ResultsDir 'notaxevasion_ge_lambda_1_2.mat']) 
MN4 = ResultsFine.ModelResults; 
load([ResultsDir 'welfare_results_no_redistrib_1_2.mat']) 
W4 = [cev_total, cev_total_se, cev_total_w]';
W4 = W4*100;
load ([ResultsDir 'welfare_results_lump_sum_1_2.mat']) 
WL4 = [cev_total, cev_total_se, cev_total_w]';
WL4 = WL4*100;

% Column 5
load([ResultsDir 'notaxevasion_ge.mat']) 
MN5 = ModelResults; 
load([ResultsDir 'welfare_results_no_redistrib.mat']) 
W5 = [cev_total, cev_total_se, cev_total_w]';
W5 = W5*100;
load ([ResultsDir 'welfare_results_lump_sum.mat']) 
WL5 = [cev_total, cev_total_se, cev_total_w]';
WL5 = WL5*100;

% Column 6
load([ResultsDir 'notaxevasion_ge_lambda_1_8.mat']) 
MN6 = ResultsFine.ModelResults; 
load([ResultsDir 'welfare_results_no_redistrib_1_8.mat']) 
W6 = [cev_total, cev_total_se, cev_total_w]';
W6 = W6*100;
load ([ResultsDir 'welfare_results_lump_sum_1_8.mat']) 
WL6 = [cev_total, cev_total_se, cev_total_w]';
WL6 = WL6*100;

%Switch to generate a complete LaTex document or just a table to import in draft:
makeCompleteLatexDocument = 0; %either 0 or 1

%% Create Differences
% column 1
M1_vec = [M1.k_entre, M1.output_se, M1.k_corp, M1.output_corp, M1.taxes_total_toY_p,...
    M1.share_entre_p, M1.cond_mean_theta, M1.ave_k_entre, M1.cond_firm_size_dist(1),...
    M1.r0*100, M1.w0]';    
MN1_vec = [MN1.k_entre, MN1.output_se, MN1.k_corp, MN1.output_corp, MN1.taxes_total_toY_p,...
    MN1.share_entre_p, MN1.cond_mean_theta, MN1.ave_k_entre, MN1.cond_firm_size_dist(1),...
    MN1.r0*100, MN1.w0]';
D1  = ((MN1_vec-M1_vec)./M1_vec)*100;
PP1 = MN1_vec-M1_vec;
% column 2
M2_vec = [M2.k_entre, M2.output_se, M2.k_corp, M2.output_corp, M2.taxes_total_toY_p,...
    M2.share_entre_p, M2.cond_mean_theta, M2.ave_k_entre, M2.cond_firm_size_dist(1),...
    M2.r0*100, M2.w0]';    
MN2_vec = [MN2.k_entre, MN2.output_se, MN2.k_corp, MN2.output_corp, MN2.taxes_total_toY_p,...
    MN2.share_entre_p, MN2.cond_mean_theta, MN2.ave_k_entre, MN2.cond_firm_size_dist(1),...
    MN2.r0*100, MN2.w0]';
D2  = ((MN2_vec-M2_vec)./M2_vec)*100;
PP2 = MN2_vec-M2_vec;
% column 3
M3_vec = [M3.k_entre, M3.output_se, M3.k_corp, M3.output_corp, M3.taxes_total_toY_p,...
    M3.share_entre_p, M3.cond_mean_theta, M3.ave_k_entre, M3.cond_firm_size_dist(1),...
    M3.r0*100, M3.w0]';    
MN3_vec = [MN3.k_entre, MN3.output_se, MN3.k_corp, MN3.output_corp, MN3.taxes_total_toY_p,...
    MN3.share_entre_p, MN3.cond_mean_theta, MN3.ave_k_entre, MN3.cond_firm_size_dist(1),...
    MN3.r0*100, MN3.w0]';
D3  = ((MN3_vec-M3_vec)./M3_vec)*100;
PP3 = MN3_vec-M3_vec;
% column 4
MN4_vec = [MN4.k_entre, MN4.output_se, MN4.k_corp, MN4.output_corp, MN4.taxes_total_toY_p,...
    MN4.share_entre_p, MN4.cond_mean_theta, MN4.ave_k_entre, MN4.cond_firm_size_dist(1),...
    MN4.r0*100, MN4.w0]';
D4  = ((MN4_vec-M1_vec)./M1_vec)*100;
PP4 = MN4_vec-M1_vec;
% column 5
MN5_vec = [MN5.k_entre, MN5.output_se, MN5.k_corp, MN5.output_corp, MN5.taxes_total_toY_p,...
    MN5.share_entre_p, MN5.cond_mean_theta, MN5.ave_k_entre, MN5.cond_firm_size_dist(1),...
    MN5.r0*100, MN5.w0]';
D5  = ((MN5_vec-M2_vec)./M2_vec)*100;
PP5 = MN5_vec-M2_vec;
% column 6
MN6_vec = [MN6.k_entre, MN6.output_se, MN6.k_corp, MN6.output_corp, MN6.taxes_total_toY_p,...
    MN6.share_entre_p, MN6.cond_mean_theta, MN6.ave_k_entre, MN6.cond_firm_size_dist(1),...
    MN6.r0*100, MN6.w0]';
D6  = ((MN6_vec-M3_vec)./M3_vec)*100;
PP6 = MN6_vec-M3_vec; 

%% Table 10, section 5.5
fid = fopen([SaveDir 'table10.tex'], 'wt');
if makeCompleteLatexDocument==1
    fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
    fprintf(fid,'\\usepackage{rotating} \n');
    fprintf(fid,'\\usepackage{floatrow} \n');
    fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{sidewaystable}[htbp] \\centering \n');
fprintf(fid,'\\caption{Tax Evasion and the Borrowing Limit} \n');
fprintf(fid,'\\label{table: alevcreditconstr} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{l c c c c c c} \\hline \n');
fprintf(fid, '& Variable & lambda_1.2  & lambda_1.5 & lambda_1.8 & lambda_1.2 & lambda_1.5 &lambda_1.8  \\\\ \n');
fprintf(fid, '\\hline \\hline \n');
fprintf(fid,'$K^{E}$    & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    D1(1), D2(1), D3(1), D4(1), D5(1), D6(1));
                                                                                                                    
fprintf(fid,'$Y^{E}$    & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    D1(2), D2(2), D3(2), D4(2), D5(2), D6(2));

fprintf(fid,'$K^{C}$    & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    D1(3), D2(3), D3(3), D4(3), D5(3), D6(3));

fprintf(fid,'$Y^{C}$    & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    D1(4), D2(4), D3(4), D4(4), D5(4), D6(4));

fprintf(fid,'$T/Y$    & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    PP1(5), PP2(5), PP3(5), PP4(5), PP5(5), PP6(5));
fprintf(fid,'Share of self-employed (\\%%)& %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    PP1(6), PP2(6), PP3(6), PP4(6), PP5(6), PP6(6));
fprintf(fid,'$E(\\theta|E)$ & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n', ...
    D1(7), D2(7), D3(7), D4(7), D5(7), D6(7));
fprintf(fid,'$E(k|E)$       & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', ...
    D1(8), D2(8), D3(8), D4(8), D5(8), D6(8));
fprintf(fid,'Share of small businesses (\\%%)& %8.3f & %8.3f %8.3f & %8.3f%8.3f & %8.3f \\\\ \n', ...
    PP1(9), PP2(9), PP3(9), PP4(9), PP5(9), PP6(9));
fprintf(fid,'$r(\\%%)$          & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    PP1(10), PP2(10), PP3(10), PP4(10), PP5(10), PP6(10));
fprintf(fid,'$w$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    D1(11), D2(11), D3(11), D4(11), D5(11), D6(11));  
fprintf(fid,'$cev$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    W1(1), W2(1), W3(1), W4(1), W5(1), W6(1)); 
fprintf(fid,'$cev_se$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    W1(2), W2(2), W3(2), W4(2), W5(2), W6(2));
fprintf(fid,'$cev_w$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    W1(3), W2(3), W3(3), W4(3), W5(3), W6(3));
fprintf(fid,'$cev$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    WL1(1), WL2(1), WL3(1), WL4(1), WL5(1), WL6(1)); 
fprintf(fid,'$cev_se$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    WL1(2), WL2(2), WL3(2), WL4(2), WL5(2), WL6(2));
fprintf(fid,'$cev_w$                & %8.3f & %8.3f & %8.3f & %8.3f& %8.3f & %8.3f \\\\ \n', ...
    WL1(3), WL2(3), WL3(3), WL4(3), WL5(3), WL6(3));
   
fprintf(fid,'\\hline  \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\end{sidewaystable} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);
    


