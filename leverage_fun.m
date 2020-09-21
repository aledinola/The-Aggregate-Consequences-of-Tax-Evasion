function [frac_bor1,leverage_percentiles,leverage_mean] = leverage_fun(policycapdet,agrid_dist,lambda,mu_se,nagrid_dist,negrid,ntgrid,dispp,agrid,policycap)

% DESCRIPTION:
%This function computes the fraction of SE who are borrowing constrained
%i.e. a<k=lambda*a. These people are borrowing the maximum allowed,lambda*a-a
%What is leverage for those who borrow i.e. k>a
% mean{max[(k-a)/a,0]}

% INPUTS:
% policycapdet: policy function for k on finer grid


frac_bor1=0; %cannot borrow for business



for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            if round(policycapdet(i,j,t),1) >= round(lambda*agrid_dist(i),1)
                frac_bor1 = frac_bor1 + mu_se(i,j,t);
            end
        end
    end
end

leverage = zeros(nagrid_dist,negrid,ntgrid);
for t=1:ntgrid
    for j=1:negrid
        for i=1:nagrid_dist
            temp = (policycapdet(i,j,t)-agrid_dist(i))/agrid_dist(i);
            if temp>0
                leverage(i,j,t) = temp;
            end
        end
    end
end

leverage_vec = leverage(:);
mu_se_vec    = mu_se(:);

q = [0.01,0.02,0.10,0.25,0.5,0.75,0.90,0.95,0.99];
leverage_percentiles = quantili(leverage_vec,mu_se_vec,q);
leverage_mean        = sum(leverage_vec.*mu_se_vec);


if dispp ==1
    fprintf('Share of borcon SE: %f \n',frac_bor1)
    fprintf('Mean of leverage SE: %f \n', leverage_mean)
    for i=1:length(q)
        fprintf('Percentiles of leverage SE: %f \t %f \n',q(i),leverage_percentiles(i) );
    end
end
%leverage_out = struct();

% %Debug
% if 0
% figure(1)
% plot(agrid_dist,agrid_dist,agrid_dist,lambda*agrid_dist,agrid_dist,policycapdet(:,4,10))
% legend('k=a','k=lambda*a','k(a,eps,theta)')
% figure(2)
% plot(agrid_dist,leverage(:,4,10))
% figure(3)
% plot(agrid,agrid,agrid,lambda*agrid,agrid,policycap(:,4,10))
% legend('k=a','k=lambda*a','k(a,eps,theta)')
% 
% figure(4)
% [leverage_vec_sort,ind_sort] = sort(leverage_vec);
% mu_se_vec_sort = mu_se_vec(ind_sort);
% nbins = 100;
% 
% [histw, vinterval] = histwc(leverage_vec, mu_se_vec, nbins);
% 
% histogram(leverage_vec,nbins,'Normalization','probability');
% title('Leverage Distribution in the Model')
% xlabel('Leverage, as max[(k-a)/a,0]')
% 
% end %END debug

end %END FUNCTION

