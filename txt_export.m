function [] = txt_export(flags,Parameters,ModStat)
%{
 In this function we export model targets in txt files

INPUTS:
 flags: structure with flags
 Parameters: structure with model parameters
 ModStat: structure with model targets ==> txt file
    
OUTPUTS:
 N.A.
 
NOTES:
 %% Updated by Alessandro Di Nola on September 16, 2020. Folder: RED.
%}

% Validate inputs
if ~isstruct(flags)
    error('Input "flags" must be a structure!')
end

if ~isstruct(Parameters)
    error('Input "Parameters" must be a structure!')
end

if ~isstruct(ModStat)
    error('Input "ModStat" must be a structure!')
end

results_folder = flags.results_folder;

fid=fopen([results_folder,'\targets_model_manual.txt'],'wt');  % overwrite

%Flag
fprintf(fid,'\n'); %blank rows (align with excel file)
fprintf(fid,'\n');
fprintf(fid,'Flags \n');
fprintf(fid,'do_GE                        %d \n',flags.do_GE); 
fprintf(fid,'no_evasion                   %d \n',Parameters.no_evasion);
fprintf(fid,'pflag                        %d \n',Parameters.pflag);
fprintf(fid,'vfi_tricks                   %d \n',Parameters.vfi_tricks);
fprintf(fid,'concavity                    %d \n',Parameters.concavity);
fprintf(fid,'labor_hiring                 %d \n',Parameters.labor_hiring);
fprintf(fid,'endo_ls_work                 %d \n',Parameters.endo_ls_work);
fprintf(fid,'endo_ls_entre                %d \n',Parameters.endo_ls_entre);
fprintf(fid,'Num_points_legrid            %d \n',Parameters.nlegrid);
fprintf(fid,'Num_points_nngrid            %d \n',Parameters.nngrid);

%External Parameters
fprintf(fid,'\n');
fprintf(fid,'External_Parameters \n');

fprintf(fid,'sigma                      %5.6f \n', Parameters.sigma);  % risk aversion
fprintf(fid,'sigma2                     %5.6f \n', Parameters.sigma2); % inverse of Frish Elasticity
fprintf(fid,'alpha                      %5.6f \n', Parameters.alpha);  % share of capital
fprintf(fid,'lambda                     %5.6f \n', Parameters.lambda); % borrowing constraint parameter for self-employed
fprintf(fid,'s                          %5.6f \n', Parameters.s);      %fine 

fprintf(fid,'rho                        %5.6f \n',Parameters.rho); 
fprintf(fid,'sigmaeps                   %5.6f \n',Parameters.sigmaeps);

if Parameters.taxfunc==1 % HSV tax function
    fprintf(fid,'tau_work                   %5.6f \n',Parameters.tau_work); 
    fprintf(fid,'lambda_work                %5.6f \n',Parameters.lambda_work); 
    fprintf(fid,'tau_entre                  %5.6f \n',Parameters.tau_entre); 
    fprintf(fid,'lambda_entre               %5.6f \n',Parameters.lambda_entre); 
    fprintf(fid,'s_work                     %5.6f \n',nan); 
    fprintf(fid,'s_entre                    %5.6f \n',nan); 
elseif Parameters.taxfunc==2 % Gouveia and Strauss
    fprintf(fid,'p_work                     %5.6f \n',Parameters.p_work); 
    fprintf(fid,'b_work                     %5.6f \n',Parameters.b_work); 
    fprintf(fid,'p_entre                    %5.6f \n',Parameters.p_entre); 
    fprintf(fid,'b_entre                    %5.6f \n',Parameters.b_entre); 
    fprintf(fid,'s_work                     %5.6f \n',Parameters.s_work); 
    fprintf(fid,'s_entre                    %5.6f \n',Parameters.s_entre); 
elseif Parameters.taxfunc==3 % Proportional tax function
    fprintf(fid,'tau_work                   %5.6f \n',Parameters.tau_work); 
    fprintf(fid,'lambda_work                %5.6f \n',Parameters.lambda_work); 
    fprintf(fid,'tau_entre                  %5.6f \n',Parameters.tau_entre); 
    fprintf(fid,'lambda_entre               %5.6f \n',Parameters.lambda_entre); 
    fprintf(fid,'s_work                     %5.6f \n',nan); 
    fprintf(fid,'s_entre                    %5.6f \n',nan); 
end

fprintf(fid,'\n');
fprintf(fid,'Estimated_Parameters \n'); %order follows Table 2 in memo: Revision: Model and Calibration

% Estimated Parameters
fprintf(fid,'beta                       %5.15f \n', Parameters.beta);   % discount factor
fprintf(fid,'psi                        %5.6f \n', Parameters.psi); % disutility of work
fprintf(fid,'delta                      %5.6f \n', Parameters.delta); % depreciation
fprintf(fid,'vi                         %5.6f \n', Parameters.vi); % span of control (production function) 
fprintf(fid,'gamma                      %5.6f \n', Parameters.gamma); % self-empl. labor share
fprintf(fid,'rho_theta                  %5.6f \n', Parameters.rho_theta); % 
fprintf(fid,'sigmaeps_theta             %5.6f \n', Parameters.sigmaeps_theta); %
fprintf(fid,'uncmean_theta              %5.6f \n', Parameters.uncmean_theta); %
if Parameters.pk_flag == 1 % logistic
    fprintf(fid,'pn_1                       %5.6f \n', Parameters.pn_1); %
    fprintf(fid,'pn_2                       %5.6f \n', Parameters.pn_2); %
    fprintf(fid,'pn_3                       %5.6f \n', Parameters.pn_3); %
elseif Parameters.pk_flag==2  % fixed
    fprintf(fid,'p_fixed                    %5.6f \n', Parameters.p_fixed); %

    fprintf(fid,'p_fixed                    %5.6f \n', Parameters.p_fixed); %

elseif Parameters.pk_flag == 3  % step function
    fprintf(fid,'klim                       %5.6f \n', klim); %
    fprintf(fid,'klim                       %5.6f \n', klim); %
end
fprintf(fid,'cc0                   %5.6f \n', Parameters.cc0); %
fprintf(fid,'cc1                   %5.6f \n', Parameters.cc1); %
fprintf(fid,'cc2                   %5.6f \n', Parameters.cc2); %
fprintf(fid,'ksi_tax               %5.6f \n', Parameters.ksi_tax); %
fprintf(fid,'tau_s                 %5.16f \n',Parameters.tau_s);
fprintf(fid,'scale_a2              %5.16f \n',Parameters.scale);
fprintf(fid,'scale_a2_entre        %5.16f \n',Parameters.scale_se);

% Table results
fprintf(fid,'\n');
fprintf(fid,'Table_Results \n');
fprintf(fid,'share_entre                %5.6f \n',ModStat.share_entre_p);
fprintf(fid,'cond_mean_theta            %5.6f \n',ModStat.cond_mean_theta);
fprintf(fid,'ave_k_entre                %5.6f \n',ModStat.ave_k_entre);
fprintf(fid,'ave_n_entre                %5.6f \n',ModStat.ave_n_entre);
fprintf(fid,'k_entre                    %5.6f \n',ModStat.k_entre);
fprintf(fid,'k_corp                     %5.6f \n',ModStat.k_supply);
fprintf(fid,'k_tot                      %5.6f \n',ModStat.capital);
%fprintf(fid,'TFP                        %5.4f \n',0); %tfp);
fprintf(fid,'Y                          %5.6f \n',ModStat.output);
fprintf(fid,'y_corp                     %5.6f \n',ModStat.output_corp);
fprintf(fid,'y_entre                    %5.6f \n',ModStat.output_se);
fprintf(fid,'r                          %5.15f \n',ModStat.r0);
fprintf(fid,'w                          %5.6f \n',ModStat.w0);
fprintf(fid,'Taxes_GDP_in_p             %5.6f \n',ModStat.taxes_total_toY_p);
fprintf(fid,'Taxes_all                  %5.15f \n',ModStat.taxes_total);
fprintf(fid,'Taxes_w                    %5.6f \n',ModStat.taxes_w);
fprintf(fid,'Taxes_e                    %5.6f \n',ModStat.taxes_e);
fprintf(fid,'Welfare_all                %5.6f \n',ModStat.Wtot);
fprintf(fid,'Welfare_entre              %5.6f \n',ModStat.We);
fprintf(fid,'Welfare_work               %5.6f \n',ModStat.Ww);
fprintf(fid,'frac_bor_constr            %5.6f \n',ModStat.frac_bor);
fprintf(fid,'frac_bor_bussiness         %5.6f \n',ModStat.frac_bor1);

fprintf(fid,'\n');
fprintf(fid,'Targets \n');
% Listing Targets (REVISION)
fprintf(fid,'Hours_worked               %5.6f \n',ModStat.model_targets(1));   % hours worked as a % of time endowment
fprintf(fid,'K/Y                        %5.6f \n',ModStat.model_targets(2));             
fprintf(fid,'Share_entre_inc_decl       %5.6f \n',ModStat.model_targets(3));
fprintf(fid,'Share_entre_inc_true       %5.6f \n',ModStat.share_inc_entre_true_p);
fprintf(fid,'Share_hiring_se            %5.6f \n',ModStat.model_targets(4)); %share of self-employed who hire external labor
fprintf(fid,'Exit_rate_(E_to_W)         %5.6f \n',ModStat.model_targets(5));
fprintf(fid,'Share_entre_inc_q1         %5.6f \n',ModStat.model_targets(6)); % share of entre by quintiles of total income  
fprintf(fid,'Share_entre_inc_q2         %5.6f \n',ModStat.model_targets(7)); % share of entre by quintiles of total income 
fprintf(fid,'Share_entre_inc_q3         %5.6f \n',ModStat.model_targets(8)); % share of entre by quintiles of total income 
fprintf(fid,'Share_entre_inc_q4         %5.6f \n',ModStat.model_targets(9)); % share of entre by quintiles of total income 
fprintf(fid,'Share_entre_inc_q5         %5.6f \n',ModStat.model_targets(10)); % share of entre by quintiles of total income 
fprintf(fid,'Fraction_entrepreneurs     %5.6f \n',ModStat.model_targets(11));
fprintf(fid,'Assets_entrepreneurs       %5.6f \n',ModStat.model_targets(12));
fprintf(fid,'Wealth_median_ratio_E_W    %5.6f \n',ModStat.wealth_med_ratio_E_W);
fprintf(fid,'Tax_evasion_q1             %5.6f \n',ModStat.model_targets(13));
fprintf(fid,'Tax_evasion_q2             %5.6f \n',ModStat.model_targets(14));
fprintf(fid,'Tax_evasion_q3             %5.6f \n',ModStat.model_targets(15));
fprintf(fid,'Tax_evasion_q4             %5.6f \n',ModStat.model_targets(16));
fprintf(fid,'Tax_evasion_q5             %5.6f \n',ModStat.model_targets(17));
fprintf(fid,'Income_Gap_all             %5.6f \n',ModStat.model_targets(18));
fprintf(fid,'Income_Gap_entre           %5.6f \n',ModStat.inc_gap_entre_p); % additional for a manual run
fprintf(fid,'Tax_Gap_entre              %5.6f \n',ModStat.tax_gap_entre_p); % additional for a manual run
fprintf(fid,'Tax_Gap_all                %5.6f \n',ModStat.tax_gap_all_p);   % additional for a manual run
fprintf(fid,'Taxes_GDP_in_p             %5.6f \n',ModStat.model_targets(19));
fprintf(fid,'Firm_size_hiring_b2        %5.6f \n',ModStat.model_targets(20)); % %percent of self-employed who hire 1-4 empl. conditional on hiring
fprintf(fid,'Firm_size_hiring_b3        %5.6f \n',ModStat.model_targets(21)); % 5-9
fprintf(fid,'Firm_size_hiring_b4        %5.6f \n',ModStat.model_targets(22)); % 10-19 
fprintf(fid,'Firm_size_hiring_b5        %5.6f \n',ModStat.model_targets(23)); % 20+
fprintf(fid,'r                          %5.15f \n',ModStat.r0);
fprintf(fid,'ED_(capital)               %5.15f \n',ModStat.ED);
fprintf(fid,'Mean_leverage_SE           %5.6f \n',ModStat.leverage_mean);

% Table with wealth and income concentration
fprintf(fid,'\n');
fprintf(fid,'Wealth_Concentration_All \n');
fprintf(fid,'Gini                       %5.6f \n',ModStat.model_wealth_all(1));
fprintf(fid,'Mean/Median                %5.6f \n',ModStat.model_wealth_all(2));
fprintf(fid,'Bottom_40                  %5.6f \n',ModStat.model_wealth_all(3));
fprintf(fid,'Top_20                     %5.6f \n',ModStat.model_wealth_all(4));
fprintf(fid,'Top_10                     %5.6f \n',ModStat.model_wealth_all(5));
fprintf(fid,'Top_1                      %5.6f \n',ModStat.model_wealth_all(6));

fprintf(fid,'\n');
fprintf(fid,'Wealth_Concentration_Entre \n');
fprintf(fid,'Gini                       %5.6f \n',ModStat.model_wealth_entre(1));
fprintf(fid,'Mean/Median                %5.6f \n',ModStat.model_wealth_entre(2));
fprintf(fid,'Bottom_40                  %5.6f \n',ModStat.model_wealth_entre(3));
fprintf(fid,'Top_20                     %5.6f \n',ModStat.model_wealth_entre(4));
fprintf(fid,'Top_10                     %5.6f \n',ModStat.model_wealth_entre(5));
fprintf(fid,'Top_1                      %5.6f \n',ModStat.model_wealth_entre(6));

fprintf(fid,'\n');
fprintf(fid,'Wealth_Concentration_Work \n');
fprintf(fid,'Gini                       %5.6f \n',ModStat.model_wealth_work(1));
fprintf(fid,'Mean/Median                %5.6f \n',ModStat.model_wealth_work(2));
fprintf(fid,'Bottom_40                  %5.6f \n',ModStat.model_wealth_work(3));
fprintf(fid,'Top_20                     %5.6f \n',ModStat.model_wealth_work(4));
fprintf(fid,'Top_10                     %5.6f \n',ModStat.model_wealth_work(5));
fprintf(fid,'Top_1                      %5.6f \n',ModStat.model_wealth_work(6));

fprintf(fid,'\n');
fprintf(fid,'Income_Concentration_All \n');
fprintf(fid,'Gini_true                  %5.6f \n',ModStat.gini_inc_true_p);
fprintf(fid,'Gini_decl                  %5.6f \n',ModStat.gini_inc_decl_p);
fprintf(fid,'Mean/Median                %5.6f \n',ModStat.model_income_all(2));
fprintf(fid,'Bottom_40                  %5.6f \n',ModStat.model_income_all(3));
fprintf(fid,'Top_20                     %5.6f \n',ModStat.model_income_all(4));
fprintf(fid,'Top_10                     %5.6f \n',ModStat.model_income_all(5));
fprintf(fid,'Top_1                      %5.6f \n',ModStat.model_income_all(6));

fprintf(fid,'\n');
fprintf(fid,'Income_Concentration_entre \n');
fprintf(fid,'Gini_true                  %5.6f \n',ModStat.gini_incse_p);
fprintf(fid,'Mean/Median                %5.6f \n',ModStat.incse_mean_to_median);
fprintf(fid,'Bottom_40                  %5.6f \n',ModStat.incse_perc_p(1));
fprintf(fid,'Top_20                     %5.6f \n',ModStat.incse_perc_p(2));
fprintf(fid,'Top_10                     %5.6f \n',ModStat.incse_perc_p(3));
fprintf(fid,'Top_1                      %5.6f \n',ModStat.incse_perc_p(4));

fprintf(fid,'\n');
fprintf(fid,'Income_Concentration_entre_decl \n');
fprintf(fid,'Gini_decl                  %5.3f \n',ModStat.gini_incse_decl_p);
fprintf(fid,'Mean/Median                %5.3f \n',ModStat.incse_mean_to_median_decl);
fprintf(fid,'Bottom_40                  %5.3f \n',ModStat.incse_perc_decl_p(1));
fprintf(fid,'Top_20                     %5.3f \n',ModStat.incse_perc_decl_p(2));
fprintf(fid,'Top_10                     %5.3f \n',ModStat.incse_perc_decl_p(3));
fprintf(fid,'Top_1                      %5.3f \n',ModStat.incse_perc_decl_p(4));

fprintf(fid,'\n');
fprintf(fid,'=============================================================');

fclose(fid);

end %END FUNCTION "txt_export"


    




