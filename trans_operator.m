function [T] = trans_operator(Parameters,Grids,policy,flags,w0,r0)

%{
DESCRIPTION: 
This script computes transition operator T(a,e,th,a',e',th')
with lottery for asset holdings (see Rios-Rull 1999).

INPUTS:
Parameters,Grids,policy,flags

OUTPUTS:
T, 6-D array: (a,e,th,a',e',th')

NOTES:
** Updated by Alessandro Di Nola on September 16, 2020.
%}

occpoldet     = policy.occpoldet;
policycapdet  = policy.policycapdet;
lepoldet      = policy.lepoldet;
policyndet    = policy.policyndet;
policyphidet  = policy.policyphidet;
apolse1det    = policy.apolse1det;
apolse0det    = policy.apolse0det;
apolwdet      = policy.apolwdet;

pn_1          = Parameters.pn_1;
pn_2          = Parameters.pn_2;
pn_3          = Parameters.pn_3;
delta         = Parameters.delta;
gamma         = Parameters.gamma;
vi            = Parameters.vi;
pflag         = Parameters.pflag;
nagrid_dist   = Parameters.nagrid_dist;
negrid        = Parameters.negrid;
ntgrid        = Parameters.ntgrid;

agrid_dist    = Grids.agrid_dist;
theta         = Grids.theta;
Peps          = Grids.Peps;
Ptheta        = Grids.Ptheta;
%-------------------------------------------------------------------------%

T = zeros(nagrid_dist,negrid,ntgrid,nagrid_dist,negrid,ntgrid);



for t=1:ntgrid % current theta
    for j=1:negrid % current eps
        for i=1:nagrid_dist % current a (assets)
            %-------------------------------------------------------------------------%
            %                          SELF - EMPLOYED there is a mistake here (work is
            %                          ok)
            %-------------------------------------------------------------------------%
            if occpoldet(i,j,t)==1 % self-employed
                % put p(k) mass on apolse1, 1-p(k) on apolse0
                kopt   = policycapdet(i,j,t);
                nopt   = policyndet(i,j,t);
                phiopt = policyphidet(i,j,t);
                leopt  = lepoldet(i,j,t);
                profitopt = theta(t)*prodfun(kopt,leopt,nopt,gamma,vi)-w0*nopt-delta*kopt-r0*kopt;
                profitopt_pos = max(profitopt,0); % ????
                p_k           = prob_audit(phiopt,profitopt_pos,theta(t), kopt, nopt, pflag, pn_1,pn_2, pn_3, leopt, gamma, vi);
                %                 [junk, kopt_ongrid] = min(abs(kopt-kgrid));
                %                 p_k = pk(kopt_ongrid);
                if p_k<0 || p_k>1
                    error('Wrong')
                end
                %% mass p_k of agents are detected ==> use apolse1det
                aopt1 = apolse1det(i,j,t);
                if aopt1<=agrid_dist(1) % all prob mass on FIRST grid point
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,1,jp,tp) =  T(i,j,t,1,jp,tp)+ p_k*Peps(j,jp)*Ptheta(t,tp);
                            
                        end
                    end
                elseif aopt1>=agrid_dist(end) % all prob mass on LAST grid point
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,end,jp,tp) = T(i,j,t,end,jp,tp)+ p_k*Peps(j,jp)*Ptheta(t,tp);
                            
                        end
                    end
                else
                    % determine the closest gridpoint below aopt1
                    jj_d = sum(agrid_dist<=aopt1);
                    % agrid(jjd) <= aopt1 < agrid(jjd+1)
                    weight_jj_d = (agrid_dist(jj_d+1)-aopt1)/(agrid_dist(jj_d+1)-agrid_dist(jj_d));
                    if ( weight_jj_d<0 ||  weight_jj_d>1); error('Wrong'); end
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,jj_d,jp,tp) = T(i,j,t,jj_d,jp,tp)+ weight_jj_d*p_k*Peps(j,jp)*Ptheta(t,tp);
                            T(i,j,t,jj_d+1,jp,tp) = T(i,j,t,jj_d+1,jp,tp)+ (1-weight_jj_d)*p_k*Peps(j,jp)*Ptheta(t,tp);
                            
                        end
                    end
                end % end if on aopt1
                %% mass (1-pk) are NOT detected ==> use apolse0det
                aopt0 = apolse0det(i,j,t);
                if aopt0<=agrid_dist(1)  % all prob mass on FIRST grid point
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,1,jp,tp) = T(i,j,t,1,jp,tp)+ (1-p_k)*Peps(j,jp)*Ptheta(t,tp);
                            
                        end
                    end
                    
                elseif aopt0>=agrid_dist(end) % all prob mass on LAST grid point
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,end,jp,tp) = T(i,j,t,end,jp,tp)+ (1-p_k)*Peps(j,jp)*Ptheta(t,tp);
                            
                            
                        end
                    end
                else
                    % determine the closest gridpoint below aopt0
                    jj_nd = sum(agrid_dist<=aopt0);
                    weight_jj_nd = (agrid_dist(jj_nd+1)-aopt0)/(agrid_dist(jj_nd+1)-agrid_dist(jj_nd));
                    if ( weight_jj_nd<0 ||  weight_jj_nd>1); error('Wrong'); end
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,jj_nd,jp,tp) = T(i,j,t,jj_nd,jp,tp)+ weight_jj_nd*(1-p_k)*Peps(j,jp)*Ptheta(t,tp);
                            T(i,j,t,jj_nd+1,jp,tp) = T(i,j,t,jj_nd+1,jp,tp)+ (1-weight_jj_nd)*(1-p_k)*Peps(j,jp)*Ptheta(t,tp);
                            
                            
                        end
                    end
                end % end if on aopt0
                %-------------------------------------------------------------------------%
                %                          WORKER
                %-------------------------------------------------------------------------%
            elseif occpoldet(i,j,t)==0 % worker
                aopt = apolwdet(i,j,t);
                if aopt<=agrid_dist(1) % all prob mass on FIRST grid point
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,1,jp,tp) = Peps(j,jp)*Ptheta(t,tp);
                        end
                    end
                elseif aopt>=agrid_dist(end) % all prob mass on LAST grid point
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,end,jp,tp) = Peps(j,jp)*Ptheta(t,tp);
                        end
                    end
                else
                    % determine the closest gridpoint below aopt
                    jj = sum(agrid_dist<=aopt);
                    weight_jj = (agrid_dist(jj+1)-aopt)/(agrid_dist(jj+1)-agrid_dist(jj));
                    % Now distribute prob. mass on points a(jj) and a(jj+1)
                    % on the finer grid
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,jj,jp,tp) = weight_jj*Peps(j,jp)*Ptheta(t,tp);
                            T(i,j,t,jj+1,jp,tp) = (1-weight_jj)*Peps(j,jp)*Ptheta(t,tp);
                        end
                    end
                end % end if on aopt
                %-------------------------------------------------------------------------%
            end % end if occupational choice
            
        end % end theta
    end % end epsilon
end % end assets a

% Check

if flags.debug_mode==1
    
    sumrowsT = zeros(nagrid_dist,negrid,ntgrid);
    for t=1:ntgrid % current theta
        for j=1:negrid % current eps
            for i=1:nagrid_dist % current a (assets)
                temp = squeeze(T(i,j,t,:,:,:));
                sumrowsT(i,j,t) = sum(temp(:));
                
            end
        end
    end
end



end %END function trans_operator
