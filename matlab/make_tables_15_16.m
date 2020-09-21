
clear;clc;close all


%% Set up paths for results and save folders
ResultsDir = 'results\lambda_mat\'; %folder where .mat files are stored
                               
SaveDir = 'results\Tables\'; %Specify here where you want to save the tex files
                             
disp(['Writing latex tables in subfolder: ' SaveDir])

%% Load results from benchmark economy (tax evasion with GE)
load([ResultsDir 'taxevasion_ge_lambda_1_2.mat']) 
M1 = ModelResults; 
maxLwork = Parameters.maxLwork;
M1_vec = [K_to_Y, (ave_n_work/maxLwork)*100, exit_ew_p, inc_gap_all_p, share_inc_entre_true_p, wealth_med_ratio_E_W, frac_bor1]';  
P1 = Parameters;


load([ResultsDir 'taxevasion_ge.mat']) 
M2 = ModelResults;
M2_vec = [K_to_Y, (ave_n_work/maxLwork)*100, exit_ew_p, inc_gap_all_p, share_inc_entre_true_p, wealth_med_ratio_E_W, frac_bor1]';  
P2 = Parameters;


load([ResultsDir 'taxevasion_ge_lambda_1_8.mat']) 
M3 = ModelResults;
M3_vec = [K_to_Y, (ave_n_work/maxLwork)*100, exit_ew_p, inc_gap_all_p, share_inc_entre_true_p, wealth_med_ratio_E_W, frac_bor1]';  
P3 = Parameters;

%Switch to generate a complete LaTex document or just a table to import in draft:
makeCompleteLatexDocument = 0; %either 0 or 1

%% Table 15: Internally Calibrated Parameters: Varying lambda

fid = fopen([SaveDir 'table15.tex'], 'wt');
if (makeCompleteLatexDocument==1)
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[htbp] \\centering \n');
fprintf(fid,'\\caption{Internally Calibrated Parameters} \n');
fprintf(fid,'\\label{table:parameters_inside} \n');
%fprintf(fid,'\\small \n'); % optional
% Begin tabular environment
fprintf(fid, '\\begin{tabular}{lccc} \\hline \n');
fprintf(fid, 'Parameter & Description & lambda=1.2  & lambda=1.8 \\\\ \n');
fprintf(fid, '\\hline \\hline \n');
fprintf(fid, '\\textit{\\underline{Preferences}} &  & &  \\\\  \n');
fprintf(fid,'$\\beta$ & Discount factor & %8.4f & %8.4f \\\\ \n',P1.beta, P3.beta);
fprintf(fid,'$\\psi$ & Disutility from working & %8.3f & %8.3f\\\\ \n',P1.psi, P3.psi);

fprintf(fid, '\\textit{\\underline{Production}} &  & &  \\\\  \n');
fprintf(fid,'$\\delta$ & Capital depreciation & %8.3f & %8.3f\\\\ \n',P1.delta, P3.delta);
fprintf(fid,'$\\nu$ & Span of control & %8.3f & %8.3f\\\\ \n', P1.vi, P3.vi);
fprintf(fid,'$\\gamma$ & Capital share in SE sector & %8.3f & %8.3f\\\\ \n', P1.gamma, P3.gamma);


fprintf(fid, '\\textit{\\underline{Self-employed ability}} &  & &  \\\\  \n');
fprintf(fid,'$\\rho_{\\theta}$   & Persistence        & %8.3f & %8.3f\\\\ \n' ,P1.rho_theta, P3.rho_theta);
fprintf(fid,'$\\sigma_{\\theta}$ & Standard deviation & %8.3f & %8.3f \\\\ \n', P1.sigmaeps_theta, P3.sigmaeps_theta);
fprintf(fid,'$\\mu_{\\theta}     $ & Unconditional mean & %8.3f & %8.3f \\\\ [0.5ex] \n'  , P1.uncmean_theta, P3.uncmean_theta);

fprintf(fid, '\\textit{\\underline{Tax evasion detection}} &  & &  \\\\  \n');
fprintf(fid,'$\\kappa$ & Cost of tax evasion  & %8.3f & %8.3f \\\\ \n', P1.cc0, P3.cc0);
fprintf(fid,'$p_{1}$ & Parameter of $p(k)$  & %8.3f & %8.3f \\\\ \n', P1.pn_1, P3.pn_1);
fprintf(fid,'$p_{2}$ & Parameter of $p(k)$  & %8.3f & %8.3f \\\\ \n', P1.pn_2, P3.pn_2);

fprintf(fid, '\\textit{\\underline{Tax functions rescale}} &  & &  \\\\  \n');
fprintf(fid,'$\\chi$  & Rescaling parameter       & %8.3f & %8.3f \\\\ \n', P1.ksi_tax, P3.ksi_tax);

fprintf(fid,'\\hline \\hline \n \\end{tabular} \n');
% End tabular environment
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);


%% Table 16: Basic Model Statistics: Varying lambda


fid = fopen([SaveDir 'table16.tex'], 'wt');
if makeCompleteLatexDocument==1
fprintf(fid,'\\documentclass[a4paper,10pt]{article} \n');
fprintf(fid,'\\begin{document} \n');
end
fprintf(fid,'\\begin{table}[tbp] \n');
fprintf(fid,'\\caption{Basic Model Statistics \\label{tab:model_fit_lambda}} \n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\begin{tabular}{lcccc} \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'Moments  & Data  & lambda=1.2 & lambda=1.5 & lambda=1.8  \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\hline \n');
fprintf(fid, '\\textit{\\underline{Targeted Statistics}} &  & & &  \\\\  \n');
fprintf(fid,'Interest rate (\\%%)                  & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', 4 , M1.r0*100 , M2.r0*100 , M3.r0*100);
fprintf(fid,'Capital-output ratio                  & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',DataTargets.K_to_Y , M1_vec(1), M2_vec(1), M3_vec(1) );
fprintf(fid,'Hours worked (\\%%)                   & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',33 , M1_vec(2), M2_vec(2), M3_vec(2));
fprintf(fid,'Share of self-employed (\\%%)         & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',DataTargets.share_entre_p, M1.share_entre_p , M2.share_entre_p, M3.share_entre_p);
fprintf(fid,'Exit rate, self-employed  (\\%%)      & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',DataTargets.exit_ew_p, M1_vec(3), M2_vec(3), M3_vec(3) );
fprintf(fid,'Misreporting rate   (\\%%)            & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',DataTargets.inc_gap_all_p, M1_vec(4), M2_vec(4), M3_vec(4) );
fprintf(fid,'Tax revenues/GDP    (\\%%)            & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n',DataTargets.taxes_total_toY_p, M1.taxes_total_toY_p, M2.taxes_total_toY_p, M3.taxes_total_toY_p);
fprintf(fid, '\\textit{\\underline{Non-targeted Statistics}} & & & &  \\\\  \n');
fprintf(fid,'Mean leverage SE                      & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n',0.289,M1.leverage_mean, M2.leverage_mean, M3.leverage_mean);
fprintf(fid,'Share of income, self-employed (\\%%) & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n',DataTargets.share_inc_entre_p, M1_vec(5), M2_vec(5), M3_vec(5) );
fprintf(fid,'Median wealth ratio, SE to W          & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n',4.02, M1_vec(6), M2_vec(6), M3_vec(6) );
fprintf(fid,'Share credit-constrained, SE           & %8.3f & %8.3f & %8.3f & %8.3f\\\\ \n',22.8, M1_vec(7), M2_vec(7), M3_vec(7) );

fprintf(fid,'\\hline \n');
fprintf(fid,'\\end{tabular} \n');
fprintf(fid,'\\end{center} \n');
	
fprintf(fid,'\\end{table} \n');
if makeCompleteLatexDocument==1; fprintf(fid,'\\end{document} \n'); end
fclose(fid);


