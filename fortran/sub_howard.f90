subroutine sub_howard

use mod_param
use mod_baselib, only: utilfun,linint,tax_work,tax_entre, prodfun, cost_evasion, prob_audit
implicit none

! Declare local variables for HOWARD 
    integer  :: acc, optw, loptw, aopt1ind, aopt0ind, kopth, phiopth, aopt1ind_bis, aopt0ind_bis, leopth
    real(dp) :: EVw, utilwopt, V0_interp_wopt, EV1, EV0, ywopt, cwopt, V0_interp_1opt, V0_interp_0opt, Vd_opt, Vn_opt
    real(dp) :: c1opt, c0opt, util1opt, util0opt, Vse1opt, incomeopt, incomeopt_pos, profitopt, profitopt_pos, prob_audit_opt
    real(dp) :: optwgrid, diff_howard, evaded_taxes_opt, evaded_taxes_opt_pos
    real(dp) :: Vhoward1(nagrid,negrid,ntgrid), Vhoward0(nagrid,negrid,ntgrid) 
    real(dp) :: Vwopt(nagrid,negrid,ntgrid), Vseopt(nagrid,negrid,ntgrid)
    
! Do Howard acceleration

    ! Start Howard improvement step from V0 
    ! and keep optimal policy fixed from VF step
    Vhoward0 = V
    Vhoward1 = 0.0d0 ! to be filled in during Howard step
    
    
    
    do acc = 1,accmax
        
        ! Start each Howard iter with clean arrays
        ! Optimal value of a worker
        Vwopt  = 0.0d0  !dimension: (nagrid, negrid, ntgrid)
        ! Optimal value of a self-employed
        Vseopt = 0.0d0  !dimension: (nagrid, negrid, ntgrid)
        Vse1opt = 0.0d0 !scalar
        Vd_opt = 0.0d0 !scalar
        Vn_opt = 0.0d0 !scalar
        
        do o = 1,ntgrid  ! current state - theta
            do j = 1,negrid  ! current state - epsilon
                do i=1,nagrid !current state - a
                    
                    ! Use old policy functions from the previous VFI step
                    ! These are all SHARED variables in the openMP implementation
                    ! WORKER, index for optimal a' and optimal labor supply
                    optw  = decisw(i,j,o)  !dim: scalar
                    loptw = decisl(i,j,o) !dim: scalar
                    ! Optimal indexes for capital, labor hiring, tax evasion and labor supply entre
                    kopth   = Icap(i,j,o)
                    nopt    = Inpol(i,j,o)
                    phiopth = Iphi(i,j,o)
                    leopth  = lepolind(i,j,o)
                    ! SELF-EMPLOYED, index for optimal a', detected and not detected
                    aopt1ind = apolse1ind(i,j,o)
                    aopt0ind = apolse0ind(i,j,o)
                    
                    !Workers - consumption, utility, value functions            
                    ywopt = w0*eps_grid(j)*lgrid(loptw)+ r0*agrid(i) !dim: scalar 
                    cwopt  = ywopt + agrid(i) - tax_work(ywopt)- agrid_fine(optw) !consumption of worker, dimension: scalar
                    !Utility - worker
                    utilwopt = utilfun(cwopt,lgrid(loptw))  ! scalar
                    
                    ! compute continuation value, for given eps and theta
                    EVw = 0.0d0 ! dimension: scalar
                    do m=1,negrid  ! future state - epsilon
                        do n=1,ntgrid  ! future state - theta
                            !write(*,*) size(agrid_fine(optw))
                            optwgrid = agrid_fine(optw)
                            V0_interp_wopt = linint(agrid,Vhoward0(:, m, n), optwgrid) ! scalar
                            !V0_interp_wopt = 1.0d0
                            EVw = EVw + V0_interp_wopt*Peps(j,m)*Ptheta(o,n) !go through all combinations
                        enddo  !future state - theta
                    enddo  !future state - epsilon
                        
                    !Value function iteration - worker
                    Vwopt(i,j,o) = utilwopt + beta*EVw !dimension: scalar
                        
                    ! Compute Vse1opt, given optimal choices for k and f
                    k = kopth
                    n = nopt
                    f = phiopth
                    l = leopth
                    
                    if (kgrid(k) <= lambda*agrid(i)) then
                    
                    !Entrepreneurs - consumption and utility
                    profitopt = theta(o)*prodfun(kgrid(k),legrid(l),ngrid(n))-w0*ngrid(n)-delta*kgrid(k)-r0*kgrid(k) ! scalar, pi greek
                    incomeopt = profitopt+r0*agrid(i)-cost_evasion(phigrid(f)) ! scalar, y^E
                    
                    incomeopt_pos = max(incomeopt,0.0d0) ! IMPORTANT!!!!!!
                    profitopt_pos = max(profitopt,0.0d0)
                            
                    ! Optimal continuation value, detected and not detected 
                    EV1 = 0.0d0 !dimension: scalar
                    EV0 = 0.0d0 !dimension: scalar
                    do m=1,negrid  !future state - epsilon
                        do n=1,ntgrid  !future state - theta
                            V0_interp_1opt = linint(agrid,Vhoward0(:, m, n),agrid_fine(aopt1ind)) ! scalar
                            EV1 = EV1 + V0_interp_1opt*Peps(j,m)*Ptheta(o,n) !go through all combinations
                            V0_interp_0opt = linint(agrid,Vhoward0(:, m, n),agrid_fine(aopt0ind)) ! scalar
                            EV0 = EV0 + V0_interp_0opt*Peps(j,m)*Ptheta(o,n) !go through all combinations
                        enddo  !future state - theta
                    enddo  !future state - epsilon
                    
                    evaded_taxes_opt = tax_entre(profitopt_pos + r0*agrid(i) ) - tax_entre(profitopt_pos*(1-phigrid(f)) + r0*agrid(i) )
                    evaded_taxes_opt_pos = max(evaded_taxes_opt,0.0d0)
                    ! c1 for being caught, dim: scalar
                    c1opt = incomeopt + agrid(i) - tax_entre(profitopt_pos*(1-phigrid(f)) + r0*agrid(i)) - agrid_fine(aopt1ind) - s*(evaded_taxes_opt_pos) 
                    ! c0 not caught, dim: scalar
                    c0opt = incomeopt + agrid(i) - tax_entre(profitopt_pos*(1-phigrid(f)) + r0*agrid(i)) - agrid_fine(aopt0ind)
                    
                    !Compute the optimal utility 
                    ! dim: scalar
                    if (c1opt>0.0d0) then
                        util1opt = utilfun(c1opt,legrid(l))
                    else
                        util1opt = large_negative
                    endif
                    if (c0opt>0.0d0) then
                        util0opt = utilfun(c0opt,legrid(l))
                    else
                        util0opt = large_negative
                    endif
                     
                    !Optimal saving choice for each phi and k
                    Vd_opt = util1opt + beta*EV1 !detected, dimension: scalar
                    Vn_opt = util0opt + beta*EV0 ! NOT detected: scalar
                    
                    !Create the expected value with respect to p(k)
                    prob_audit_opt = prob_audit(phigrid(f),profitopt_pos,theta(o),kgrid(k),ngrid(n),legrid(l),pflag)
                    Vse1opt = prob_audit_opt*Vd_opt + (1.0d0-prob_audit_opt)*Vn_opt
                    
                    ! Optimal value of being ENTRE
                    Vseopt(i,j,o) = Vse1opt
                    
                    
                    else ! collateral constraint is violated
                        
                       Vseopt(i,j,o) = large_negative
                        
                    endif ! END IF collateral
                                        
                    !Write down the current value function
                    Vhoward1(i,j,o) = max(Vwopt(i,j,o), Vseopt(i,j,o))
                    
                enddo  !current state - a
            enddo ! current state - epsilon
        enddo ! current state - theta
        
        !Error for VFI
        diff_howard = maxval((abs(Vhoward1-Vhoward0)))
        if (display_howard==1) then
            write(*,*) "-----------------------------------"
            write(*,'(a,f0.7)') " Diff Howard = ", diff_howard
            write(*,*) "-----------------------------------"
        endif
        
        ! Update 
        Vhoward0 = Vhoward1
        
    enddo ! end accelerations 1/accmax
    
    V = Vhoward1
    write(*,*) "-----------------------------------"
    write(*,*) "Howard acceleration completed"
    write(*,*) "-----------------------------------"

         
end subroutine sub_howard
   
