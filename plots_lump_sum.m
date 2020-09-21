%% This script is called by make_figures_paper_fine_v3.m
%{
It generates Figures (6a) and (6b)
panel (a): All Households, Lump-Sum
panel (b): By Occupation, Lump-sum
Notes: Updated by Alesandro Di Nola on September 15, 2020. %% Folder RED
%}

switch method
    case 'PE'
        
        % Not relevant
        
    case 'GE'
        
        %% Figure 6(a): "All Households, Lump-Sum"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),cev_total_p(1:indmax),'-o','linewidth',ldw)
        xline(1.75,'--r','linewidth',ldw);
        grid on
        %legend('Total','location','best')
        xlabel('Fine s','FontSize',fontw)
        ylabel('Consumption equivalent units (%)','FontSize',fontw)
        if axis_tight==0
            %ylim([99 101])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1; title('The Impact of the Fine, Welfare'); end
        hold off
        print([SaveDir 'fig6a'],fileExt)
        
        %% Figure 6(b): "By Occupation, Lump-sum"
        fig_ct=fig_ct+1;
        figure(fig_ct)
        plot(s_vec(1:indmax),cev_total_se_p(1:indmax),'-o','linewidth',ldw)
        hold on
        plot(s_vec(1:indmax),cev_total_w_p(1:indmax),'-s','linewidth',ldw)
        hold on
        xline(1.75,'--r','linewidth',ldw);
        grid on
        legend('Self-employed','Workers','location','best','FontSize',fontw)
        xlabel('Fine s','FontSize',fontw)
        ylabel('Consumption equivalent units (%)','FontSize',fontw)
        if axis_tight==0
            %ylim([y_axis_low y_axis_high])
            xlim([s_vec(1) s_vec(end)])
        else
            axis tight
        end
        if ts==1;title('The Impact of the Fine, Welfare');end
        hold off
        print([SaveDir 'fig6b'],fileExt)
        
    otherwise
        error('Selected method does not exist!!')
        
end
