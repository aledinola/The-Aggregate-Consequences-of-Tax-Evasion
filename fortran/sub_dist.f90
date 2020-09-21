subroutine sub_dist
	use mod_param
	use mod_baselib, only: prodfun, prob_audit,bracket 
	implicit none
    
    ! Declare local variables for distribution 
    real(dp) :: koptd, noptd, phioptd, profitoptd, profitopt_posd, leoptd
    
    !=================================================
    ! Invariant distribution
    !=================================================
    
    ! Compute distribution on original grid vs finer grid
    ! 2 methods: finer,finer2

    ! mu_method 1 = 'finer'  ! Distribution is computed on finer grid, size: nagrid_dist with lotteries a la Rios Rull
    ! mu_method 2 = 'finer2'  ! as in 'finer' but WITHOUT lotteries

    ! mu_method is an integer
 
    write(*,*) "============================================================="
    write(*,*) "STATIONARY DISTRIBUTION ITERATION..."
    write(*,*) "============================================================="
    
    select case (mu_method)
    
    case (1) ! method finer, our best so far
    
        mu = 1.0d0 / real(nagrid_dist*negrid*ntgrid,dp)
        
        ! mu dimension: (nagrid_dist,negrid,ntgrid)
        diff_d = 10.0d0
        iter_mu = 0
        
        do while (diff_d > tol_dist .and. iter_mu<=maxit_dist)
        
            iter_mu = iter_mu + 1
            
            mu_next = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
            
            do i=1,nagrid_dist ! current asset holdings
                do j=1,negrid ! current shock eps
                    do t=1,ntgrid ! current shock theta
                        
                        ! which policy for assets do I have to use?
                        if ( occpoldet(i,j,t) == 1 ) then ! self-employed
                            
                            ! put p_k mass on apolse1, 1-p_k on apolse0
                            koptd   = policycapdet(i,j,t)
                            noptd   = policyndet(i,j,t)
                            phioptd = policyphidet(i,j,t)
                            leoptd  = lepoldet(i,j,t)
                            !kopt_ongrid = minloc(abs(kopt-kgrid))
                            profitoptd = theta(t)*prodfun(koptd,leoptd,noptd)-w0*noptd-delta*koptd-r0*koptd ! scalar
                            profitopt_posd = max(profitoptd,0.0d0)
                            p_k = prob_audit(phioptd,profitopt_posd,theta(t),koptd,noptd,leoptd,pflag)
                            
                            !! mass p_k of agents are detected ==> use apolse1det
                            aopt1 = apolse1det(i,j,t)
                            
                            ! 3 cases to consider: aopt1 <= min(grid), aopt1 >= max(grid), aopt1 in between
                            if (aopt1<=agrid_dist(1) ) then ! all prob mass on FIRST grid point
                                
                                do jp=1,negrid
                                    do tp=1,ntgrid
                                
                                        mu_next(1,jp,tp) = &
                                        mu_next(1,jp,tp) + p_k*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                
                                    enddo
                                enddo
                        
                            elseif ( aopt1>=agrid_dist(nagrid_dist) ) then ! all prob mass on LAST grid point
                        
                                do jp=1,negrid
                                    do tp=1,ntgrid
                        
                                        mu_next(nagrid_dist,jp,tp) = &
                                        mu_next(nagrid_dist,jp,tp) + p_k*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                        
                                    enddo
                                enddo
                        
                            else
                        
                                ! determine the closest gridpoint below aopt1
                                call bracket(agrid_dist,aopt1,left,right)
                                jj_d = left
                                
                                ! agrid(jjd) <= aopt1 < agrid(jjd+1) 
                                weight_jj_d = &
                                (agrid_dist(jj_d+1)-aopt1) / (agrid_dist(jj_d+1) - agrid_dist(jj_d))
                            
                                  do jp=1,negrid
                                  do tp=1,ntgrid
                                
                                    mu_next(jj_d,jp,tp) = &
                                    mu_next(jj_d,jp,tp) + &
                                    weight_jj_d*p_k*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                    
                                    mu_next(jj_d+1,jp,tp) = &
                                    mu_next(jj_d+1,jp,tp) + &
                                    (1.0d0-weight_jj_d)*p_k*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                  
                                    enddo
                                enddo
     
                            endif ! endif on aopt1
                            
                            !! NOT detected ==> use apolse0det
                            aopt0 = apolse0det(i,j,t)
                    
                            ! 3 cases to consider: aopt0 <= min(grid), aopt0 >= max(grid), aopt0 in between
                            if (aopt0 <= agrid_dist(1)) then ! all prob mass on FIRST grid point
                            
                            do jp=1,negrid
                                do tp=1,ntgrid
                   
                                    mu_next(1,jp,tp) = &
                                    mu_next(1,jp,tp) + &
                                    (1.0d0-p_k)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                    
                                enddo
                            enddo

                            elseif (aopt0 >= agrid_dist(nagrid_dist) ) then ! all prob mass on LAST grid point
                        
                                do jp=1,negrid
                                    do tp=1,ntgrid
                                
                                        mu_next(nagrid_dist,jp,tp) = &
                                        mu_next(nagrid_dist,jp,tp) + &
                                        (1.0d0 - p_k)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                    
                                    enddo
                                enddo
                    
                            else
                            
                            ! determine the closest gridpoint below aopt0
                            call bracket(agrid_dist,aopt0,left,right)
                            jj_nd = left
                            weight_jj_nd = &
                            (agrid_dist(jj_nd+1)-aopt0)/(agrid_dist(jj_nd+1)-agrid_dist(jj_nd))
                            
                            do jp=1,negrid
                                do tp=1,ntgrid
                            
                                    mu_next(jj_nd,jp,tp) = &
                                    mu_next(jj_nd,jp,tp) + &
                                    weight_jj_nd*(1.0d0-p_k)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                    
                                    mu_next(jj_nd+1,jp,tp) = &
                                    mu_next(jj_nd+1,jp,tp) + &
                                    (1.0d0-weight_jj_nd)*(1.0d0-p_k)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                            
                                enddo
                            enddo
                    
                            endif ! endif on aopt0

                        elseif ( occpoldet(i,j,t)==0 ) then ! worker
                        
                            aopt = apolwdet(i,j,t)
                        
                            ! 3 cases to consider: aopt1 <= min(grid), aopt1 >= max(grid), aopt1 in between
                            if ( aopt <= agrid_dist(1) ) then ! all prob mass on FIRST grid point
                    
                                do jp=1,negrid
                                    do tp=1,ntgrid
                            
                                        mu_next(1,jp,tp) = &
                                        mu_next(1,jp,tp) + Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                        
                                    enddo
                                enddo
                            
                            elseif ( aopt>=agrid_dist(nagrid_dist) ) then ! all prob mass on LAST grid point
                    
                                do jp=1,negrid
                                
                                    do tp=1,ntgrid
                            
                                        mu_next(nagrid_dist,jp,tp) = &
                                        mu_next(nagrid_dist,jp,tp) + Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                            
                                    enddo
                                enddo

                            else
        
                                ! determine the closest gridpoint below aopt
                                call bracket(agrid_dist,aopt,left,right)
                                jj = left
                                weight_jj = (agrid_dist(jj+1)-aopt)/(agrid_dist(jj+1)-agrid_dist(jj))

                                ! Now distribute prob. mass on points a(jj) and a(jj+1)
                                ! on the finer grid
                                do jp=1,negrid
                                    do tp=1,ntgrid

                                        mu_next(jj,jp,tp) = &
                                        mu_next(jj,jp,tp) + &
                                        weight_jj*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                        
                                        mu_next(jj+1,jp,tp) = &
                                        mu_next(jj+1,jp,tp) + &
                                        (1.0d0-weight_jj)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                    
                                    enddo
                                enddo
            
                            endif ! end if on aopt
        
                        end if ! end if occupational choice

                    enddo ! end theta
                enddo ! end epsilon
            enddo ! end assets a
    
            
            !Check that also mu_next sum to 1
            check = sum(mu_next)
            mu_next = mu_next / check
            
            !Error for updating
            diff_d = maxval( ( abs( mu_next-mu )  ) )
            
            !Updating
            mu = mu_next
            
            if (display_mu==1) then
                write(*,*) "-----------------------------------"
                write(*,*) "Iteration: ", iter_mu
                write(*,"(' Check: ',f17.15)") check
                write(*,"(' Crit: ',f10.7)") diff_d
                write(*,*) "-----------------------------------"
            endif
            
 
    enddo !while loop

    case (2) ! without lotteries
        
            mu = 1.0d0/(nagrid_dist*negrid*ntgrid)
            
            diff_d = 10.0d0
            iter_mu = 0
        
            do while ( diff_d>tol_dist .and. iter_mu<=maxit_dist)
            
                iter_mu = iter_mu + 1
                mu_next = 0.0d0
            
                do i=1,nagrid_dist ! current asset holdings
                    do j=1,negrid ! current shock z
                        do t=1,ntgrid ! current shock theta
                
                            ! put p_k mass on apolse1, 1-p_k on apolse0
                            koptd   = policycapdet(i,j,t)
                            noptd   = policyndet(i,j,t)
                            phioptd = policyphidet(i,j,t)
                            leoptd  = lepoldet(i,j,t)
                            !kopt_ongrid = minloc(abs(kopt-kgrid))
                            profitoptd = theta(t)*prodfun(koptd,leoptd,noptd)-w0*noptd-delta*koptd-r0*koptd ! scalar
                            profitopt_posd = max(profitoptd,0.0d0)
                            p_k = prob_audit(phioptd,profitopt_posd,theta(t),koptd,noptd,leoptd,pflag)
                            
                            !! mass p_k of agents are detected ==> use apolse1det
                            aopt1 = apolse1det(i,j,t)
                            
                            ! determine the closest gridpoint
                            aopt1_ongrid = minloc(abs(agrid_dist-aopt1))
                            
                            !! (1-pk) are not detected ==> use apolse0det
                            aopt0 = apolse0det(i,j,t)
                            
                            ! determine the closest gridpoint
                            aopt0_ongrid = minloc(abs(agrid_dist-aopt0))
                            
                            !! some agents are worker ==> use apolwdet
                            aoptw = apolwdet(i,j,t)
                            
                            ! determine the closest gridpoint
                            aoptw_ongrid = minloc(abs(agrid_dist-aoptw))
                            
                            ! Now ready to distribute the mass
                    
                                do jp=1,negrid
                                    do tp=1,ntgrid
            
                                        ! workers
                                        mu_next(aoptw_ongrid(1),jp,tp) = &
                                        mu_next(aoptw_ongrid(1),jp,tp) + &
                                        (1-occpoldet(i,j,t))*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                        
                                        ! SE, detected
                                        mu_next(aopt1_ongrid(1),jp,tp) = &
                                        mu_next(aopt1_ongrid(1),jp,tp) + &
                                        p_k*occpoldet(i,j,t)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)
                                        
                                        ! SE, not detected
                                        mu_next(aopt0_ongrid(1),jp,tp) = &
                                        mu_next(aopt0_ongrid(1),jp,tp) + &
                                        (1-p_k)*occpoldet(i,j,t)*Peps(j,jp)*Ptheta(t,tp)*mu(i,j,t)

                                        enddo ! jp
                                    enddo ! tp
            
                        enddo ! end theta
                    enddo ! end epsilon
                enddo ! end assets a
            
                
                ! Check that also mu_next sum to 1
                check = sum(mu_next)
                mu_next = mu_next/check
                diff_d = maxval( ( abs( mu_next-mu )  ) )
                
                
                !Print only for the last iterations
                !if (diff_d < 2.0d0*tol_dist) then
                write(*,*) "-----------------------------------"
                write(*,*) "Iteration: ", iter_mu
                write(*,"(' Check: ',f17.15)") check
                write(*,"(' Crit: ',f10.7)") diff_d
                write(*,*) "-----------------------------------"
                !endif
                
           
                mu = mu_next
                
            enddo ! end while

    case default

        write(*,*) 'ERROR, mu_method must be either 1 or 2'
   
    end select ! end select case
    
    if (diff_d <= tol_dist) then
        write(*,*) 'Distrib. converged succesfully'
        exitflag_mu = 0
    else
        write(*,*) 'Distrib. did not converge succesfully'
        exitflag_mu = 1
    endif 
    
    if ( iter_mu < maxit_dist ) then 
        write(*,*) "============================================================="
        write(*,*) "STATIONARY DISTRIBUTION FOUND"
        write(*,*) "============================================================="
    else
        write(*,*) "============================================================="
        write(*,*) "MAXIMUM NUMBER OF ITERATIONS REACHED"
        write(*,*) "============================================================="
        
    endif
    
end subroutine sub_dist
