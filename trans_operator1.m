% This script computes transition operator T(a,e,th,a',e',th')
% WITHOUT lottery for asset holdings

T = zeros(nagrid_dist,negrid,ntgrid,nagrid_dist,negrid,ntgrid);
apolwdet = zeros(nagrid_dist,negrid,ntgrid);
apolse1det = zeros(nagrid_dist,negrid,ntgrid);
apolse0det = zeros(nagrid_dist,negrid,ntgrid);

% Linear interp for policy functions for assets
for ie=1:negrid
    for it=1:ntgrid
        apolwdet(:,ie,it)=interp1(agrid',apolw(:,ie,it),agrid_dist');
        apolse1det(:,ie,it)=interp1(agrid',apolse1(:,ie,it),agrid_dist');
        apolse0det(:,ie,it)=interp1(agrid',apolse0(:,ie,it),agrid_dist');
    end
end


for i=1:nagrid_dist % current asset holdings
    for j=1:negrid % current shock eps
        for t=1:ntgrid % current shock theta
            % Compute pk
            kopt   = policycapdet(i,j,t);
            nopt   = policyndet(i,j,t);
            phiopt = policyphidet(i,j,t);
            profitopt = theta(t)*prodfun(kopt,le,nopt,gamma,vi)-w0*nopt-delta*kopt-r0*kopt;
            profitopt_pos = max(profitopt,0); % ????
            p_k           = prob_audit(phiopt,profitopt_pos,pn_1,pn_2);
%             [junk, kopt_ongrid] = min(abs(kopt-kgrid));
%             p_k = pk(kopt_ongrid);
            %% mass p_k of agents are detected ==> use apolse1det
            aopt1 = apolse1det(i,j,t);
            % determine the closest gridpoint
            [junkk,aopt1_ongrid] = min(abs(agrid_dist-aopt1));
            %% (1-pk) are not detected ==> use apolse0det
            aopt0 = apolse0det(i,j,t);
            % determine the closest gridpoint
            [junkk,aopt0_ongrid] = min(abs(agrid_dist-aopt0));
            %% some agents are worker ==> use apolwdet
            aoptw = apolwdet(i,j,t);
            % determine the closest gridpoint
            [junkk,aoptw_ongrid] = min(abs(agrid_dist-aoptw));
            % Now ready to distribute the mass
            %% workers
            if occpoldet(i,j,t)==0
                for jp=1:negrid
                    for tp=1:ntgrid
                        T(i,j,t,aoptw_ongrid,jp,tp) = Peps(j,jp)*Ptheta(t,tp);
                    end
                end
            %% Self-employed, detected and NOT detected    
            elseif occpoldet(i,j,t)==1
            
                if aopt0_ongrid==aopt1_ongrid
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,aopt1_ongrid,jp,tp) = Peps(j,jp)*Ptheta(t,tp);
                        end
                    end
                else % choice of D not equal to NOT D
                    for jp=1:negrid
                        for tp=1:ntgrid
                            T(i,j,t,aopt1_ongrid,jp,tp) = p_k*Peps(j,jp)*Ptheta(t,tp);
                            T(i,j,t,aopt0_ongrid,jp,tp) = (1-p_k)*Peps(j,jp)*Ptheta(t,tp);
                        end
                    end
                end % endif on "aopt0_ongrid==aopt1_ongrid"
                
            end
            
        end % end theta
    end % end epsilon
end % end assets a



% Check

if debug_mode==1
    sumrowsT = zeros(nagrid_dist,negrid,ntgrid);
    for i=1:nagrid_dist % current a (assets)
        for j=1:negrid % current eps
            for t=1:ntgrid % current theta
                temp = squeeze(T(i,j,t,:,:,:));
                sumrowsT(i,j,t) = sum(temp(:));
                
            end
        end
    end
end
















