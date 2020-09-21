function [exit_ew,exit_we] = fun_exit_rate(Parameters,Grids,policy,flags,w0,r0,mu)
%-------------------------------------------------------------------------%
% DESCRIPTION:
% This function calls the external function trans_operator
%-------------------------------------------------------------------------%
% INPUTS: 
% Parameters: structure with model parameters
% Grids:      structure with grids
% policy:     structure with policy functions
% flags:      structure with flags
% w0,r0:      wage and interest rate
%-------------------------------------------------------------------------%
% OUTPUTS: 
% exit_ew: exit rate from E to W
% exit_we: exit rate from W to E (not targeted)
%-------------------------------------------------------------------------%

nagrid_dist   = Parameters.nagrid_dist;
negrid        = Parameters.negrid;
ntgrid        = Parameters.ntgrid;
occpoldet     = policy.occpoldet;

% First, we compute the transition operator  T(a,e,theta,a',e',theta')
% This requires a lot of RAM, hence it is done only if lowMemory = 0
if Parameters.mu_method      ==  1
    T = trans_operator(Parameters,Grids,policy,flags,w0,r0); % with lottery method for assets
elseif Parameters.mu_method  ==  2
    trans_operator1 % without lotteries
end
% It gives T, dimension(nagrid_dist,negrid,ntgrid,nagrid_dist,negrid,ntgrid)

%----------------- FROM ENTRE to WORK ------------------------------------%
% Now compute share of entre who become workers, using the transition T and
% the occupational choice decision rule
exit_ew_mat = zeros(nagrid_dist,negrid,ntgrid);

for ip=1:nagrid_dist % future a
    for jp=1:negrid % future shock eps
        for tp=1:ntgrid % future shock theta
            exit_ew_mat = exit_ew_mat + occpoldet.*T(:,:,:,ip,jp,tp)*(1-occpoldet(ip,jp,tp));
        end % future theta
    end % future eps
end % future a

% Normalize by the fraction of entre (to transform the number of exiting E
% into the share of exiting E out of total E)
share_entre = sum(sum(sum(occpoldet.*mu)));
exit_ew     = sum(exit_ew_mat(:).*mu(:)) / share_entre;

%----------------- FROM WORK to ENTRE ------------------------------------%

if nargout>1

exit_we_mat = zeros(nagrid_dist,negrid,ntgrid);

for ip=1:nagrid_dist % future a
    for jp=1:negrid % future shock eps
        for tp=1:ntgrid % future shock theta
            exit_we_mat = exit_we_mat + (1-occpoldet).*T(:,:,:,ip,jp,tp)*occpoldet(ip,jp,tp);
        end % future theta
    end % future eps
end % future a

% Normalize by the fraction of workers
exit_we = sum(exit_we_mat(:).*mu(:)) / (1-share_entre);

% Check that the exit rates computed above give rise to the perc. of
% workers and entrepreneurs in the stationary distrib.
entryexit_mat = [1-exit_we    exit_we;
                 exit_ew      1-exit_ew];

tempp   =   entryexit_mat^1000; % get invariant distribution of markov chain above
% disp('Check if numbers are the same 
% disp(tempp(1,:))        
% disp([1-share_entre, share_entre])

end

end %END function exit rate

