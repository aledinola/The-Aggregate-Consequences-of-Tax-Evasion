%Updated on 14Sep2020 by Anna

clear;clc

%%%% HERE Specify here where results (.mat files) are saved %%%
ResultsDir = 'results\lambda_mat\';
%%%% HERE Specify here where you want to save the figures/tables %%%
SaveDir = 'results\Tables\';



%% Load mat files
    %Switch to generate a complete LaTex document or just a table to import in draft:
    makeCompleteLatexDocument = 0; %either 0 or 1
    xx = NaN;
    
    %BENCHMARK
    load([ResultsDir 'taxevasion_ge.mat'],'ModelResults')
    M1 = ModelResults; %matlab structure with results from TE economy
    
    
    %PERFECT TAX ENFORCEMENT lambda=1.65, output_all fixed
    load([ResultsDir 'notaxevasion_ge_output_lambda.mat'],'ModelResults')
    M2 = ModelResults; %matlab structure
    
    
    %PERFECT TAX ENFORCEMENT lambda=1.64, output_se fixed
    load([ResultsDir 'notaxevasion_ge_lambda_out_se.mat'],'ModelResults')
    M3 = ModelResults;
    
    
    %% Create percentage differences
M1_vec = [M1.output
    M1.output_se;
    M1.output_corp;
    M1.k_entre;
    M1.k_corp;
    M1.share_entre_p;
    M1.cond_mean_theta;
    M1.ave_k_entre;
    M1.cond_firm_size_dist(1);
    M1.r0*100;
    M1.w0;
];    
    
M2_vec = [M2.output
    M2.output_se;
    M2.output_corp;
    M2.k_entre;
    M2.k_corp;
    M2.share_entre_p;
    M2.cond_mean_theta;
    M2.ave_k_entre;
    M2.cond_firm_size_dist(1);
    M2.r0*100;
    M2.w0;
]; 

  M3_vec = [M3.output
    M3.output_se;
    M3.output_corp;
    M3.k_entre;
    M3.k_corp;
    M3.share_entre_p;
    M3.cond_mean_theta;
    M3.ave_k_entre;
    M3.cond_firm_size_dist(1);
    M3.r0*100;
    M3.w0;
];   
 
% perc. change column 1
 Perc_change_C1 = ((M2_vec-M1_vec)./M1_vec)*100;  
% perc. change column 2
 Perc_change_C2 = ((M3_vec-M1_vec)./M1_vec)*100;  
 
 
    
    
%% Table 11, section 5.5
fid = fopen([SaveDir 'table11.tex'], 'wt');
if makeCompleteLatexDocument==1
    fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
    fprintf(fid,'\\usepackage{rotating} \n');
    fprintf(fid,'\\usepackage{floatrow} \n');
    fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{sidewaystable}[htbp] \\centering \n');
fprintf(fid,'\\caption{Tax Evasion and the Alleviation of Credit Constraints} \n');
fprintf(fid,'\\label{table: alevcreditconstr} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{l c c } \\hline \n');
fprintf(fid, '& Variable & lambda_1.65  & lambda_1.64  \\\\ \n');
fprintf(fid, '\\hline \\hline \n');
                                                                                                                    
fprintf(fid,'$Y$        & %8.3f & %8.3f  \\\\ \n', 0.00 , Perc_change_C2(1));
fprintf(fid,'$Y^{E}$    & %8.3f & %8.3f  \\\\ \n', Perc_change_C1(2), 0.00 );
fprintf(fid,'$Y^{C}$    & %8.3f & %8.3f  \\\\ \n', Perc_change_C1(3), Perc_change_C2(3));
fprintf(fid,'$K^{E}$    & %8.3f & %8.3f  \\\\ \n', Perc_change_C1(4), Perc_change_C2(4));
fprintf(fid,'$K^{C}$    & %8.3f & %8.3f  \\\\ \n', Perc_change_C1(5), Perc_change_C2(5));
fprintf(fid,'Share of self-employed (\\%%)& %8.3f & %8.3f  \\\\ \n', M2_vec(6)- M1_vec(6), M3_vec(6)- M1_vec(6));
fprintf(fid,'$E(\\theta|E)$ & %8.3f & %8.3f \\\\ \n', Perc_change_C1(7), Perc_change_C2(7));
fprintf(fid,'$E(k|E)$       & %8.3f & %8.3f \\\\ \n', Perc_change_C1(8), Perc_change_C2(8));
fprintf(fid,'Share of small businesses (\\%%)& %8.3f & %8.3f  \\\\ \n', M2_vec(9)- M1_vec(9), M3_vec(9)- M1_vec(9));
fprintf(fid,'$r(\\%%)$          & %8.3f & %8.3f  \\\\ \n', M2_vec(10)- M1_vec(10), M3_vec(10)- M1_vec(10));
fprintf(fid,'$w$                & %8.3f & %8.3f  \\\\ \n', Perc_change_C1(11), Perc_change_C2(11));           
   
fprintf(fid,'\\hline  \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\end{sidewaystable} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);
    






