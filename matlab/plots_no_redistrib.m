%% This script is called by make_figures_paper_fine_v3.m
%{
It generates Figure (5)
Top panel: partial equilibrium
Bottom panel: general equilibrium
No fiscal redistribution
Notes: Updated by Alesandro Di Nola on September 15, 2020. %% Folder RED
%}

y_axis_low  = 99;
y_axis_high = 101;

switch method
    case 'PE'
        %% Figure 5(a): "Misreporting Rate, PE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),inc_gap_vec_norm(1:indmax),'-o','linewidth',ldw)
        xline(1.75,'--r','linewidth',ldw);
        grid on
        xlabel('Fine s','FontSize',fontw)
        %ylabel('Income Gap (%)','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([88 106])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1; title('The Impact of the Fine, Misrep. rate');end
        print([SaveDir 'fig5a'],fileExt)
        
        %% Figure 5(b): "Tax Revenues, PE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),taxes_all_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),taxes_e_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),taxes_w_vec_norm(1:indmax),'-^','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Total','Self-employed','Workers','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([94 113])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1;title('The Impact of the Fine, tax revenues');end
        hold off
        print([SaveDir 'fig5b'],fileExt)
        
        %% Figure 5(c): "Share of Self-employed, PE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),share_entre_p_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),cond_mean_theta_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Share','Average Productivity','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([98.5 y_axis_high])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1; title('The Impact of the Fine, share and prod. of SE'); end
        hold off
        print([SaveDir 'fig5c'],fileExt)
        
        %% Figure 5(d): "Production, PE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),capital_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),output_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),labor_vec_norm(1:indmax),'-^','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Total Capital','Total Output','Total Labor','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([96 102])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1; title('The Impact of the Fine, aggregates'); end
        hold off
        print([SaveDir 'fig5d'],fileExt)
        
        
        
    case 'GE'  % GENERAL EQUILIBRIUM, NO REDISTRIBUTION
        
        %% Figure 5(a): "Misreporting and Prices, GE" 
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),inc_gap_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),r_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),w_vec_norm(1:indmax),'-^','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Income Gap','r','w','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            %ylim([10.5 12.7])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1;title('The Impact of the Fine, Misrep. rate');end
        hold off
        print([SaveDir 'fig5a_ge'],fileExt)
                
        %% Figure 5(b): ""Tax Revenues, GE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),taxes_all_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),taxes_e_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),taxes_w_vec_norm(1:indmax),'-^','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Total','Self-employed','Workers','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([94 113])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1;title('The Impact of the Fine, tax revenues');end
        hold off
        print([SaveDir 'fig5b_ge'],fileExt)
        
        %% Figure 5(c): "Share of Self-Employed, GE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),share_entre_p_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),cond_mean_theta_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Share','Average Productivity','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([98.5 y_axis_high])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1;title('The Impact of the Fine, share and prod. of SE');end
        hold off
        print([SaveDir 'fig5c_ge'],fileExt)
        
        %% Figure 5(d): "Production, GE"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),capital_vec_norm(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),output_vec_norm(1:indmax),'-s','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),labor_vec_norm(1:indmax),'-^','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Total Capital','Total Output','Total Labor','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Normalized to benchmark','FontSize',fontw)
        if axis_tight==0
            ylim([96 102])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1;title('The Impact of the Fine, aggregates');end
        hold off
        print([SaveDir 'fig5d_ge'],fileExt)
        
        
        
    otherwise
        error('Selected method does not exist!!')
        
end
