subroutine sub_vfi()
	use mod_param
	use mod_baselib, only: utilfun,linintv,tax_work,tax_entre, prodfun, cost_evasion, prob_audit, interp1
    use omp_lib
	implicit none
    
    ! Declare local variables for sub_vfi
    real(dp) :: kval,nval,leval,lwval,phival
    

    !=================================================
    ! Value function iteration
    !=================================================
    
    !call omp_set_num_threads(7)
    
    !$omp parallel if (par_fortran==1)
        write(*,*) 'Parallel hello to you!'
    !$omp end parallel   
   
    write(*,*) "============================================================="
    write(*,*) "VALUE FUNCTION ITERATION..."
    write(*,*) "============================================================="
    
    if (vfi_tricks == 1) then
        
        write(*,*) "-----------------------------------"
        write(*,*) "using speedup tricks on k and n"
        write(*,*) "-----------------------------------"
    
    endif
    
    if (concavity == 1 ) then
        
        write(*,*) "using speedup tricks on ap (concavity)"
        write(*,*) "-----------------------------------"
    
    endif 
    


    diff = 10.0d0
    itervf = 0
    ! Initialize policy functions for Howard acceleration
    !Time
    t1 = omp_get_wtime()
	!call cpu_time(t1)
    
    
    !Initial guess: V0(a,eps,theta)
    do while (diff > tolV .and. itervf <= maxitV)
        
        ! Initialize occupational choice polfun
        occpol = 0 ! integer values
        
        ! Value function of a worker
        !These initializations are probably useless
        Vw  = 0.0d0  !value function of a worker, dimension: (nagrid, negrid, ntgrid)
        decisw = 0 !index on the decision rule (saving) for worker, dimension(nagrid, negrid, ntgrid)
        decisl = 0 !index on the decision rule (labor supply) for worker, dimension(nagrid, negrid, ntgrid)
        decisw_temp = 0 !dimension(nlgrid)
        
        !$omp parallel if (par_fortran==1) default(none) private(j,o,EV,m,n,V0_interp,min_aw,min_a0,min_a1,i,il,lwval,yw,utilw,iap,cw,&
        !$ vtempW,RHSw,decisw_temp,tempw,Vse1,k,kval,ixn,nval,ile,leval,f,phival,profit,&
        !$ income, profit_pos,income_pos,evaded_taxes,evaded_taxes_pos,Vd,c1,vtemp1,&
        !$ Vn,c0,vtemp0,prob_audit_temp,Vse1_temp,min_k,min_n,max_k,max_n,decis1,decis0) &
        !$ shared(negrid,ntgrid,nagrid,nlgrid,eps_grid,&
        !$ agrid,agrid_fine,lgrid,beta,r0,w0,concavity,nagrid_fine,Peps,Ptheta,Vw,decisl,decisw,V0,&
        !$ kgrid,lambda,ngrid,nlegrid,legrid,phigrid,delta,Icap,Inpol,Iphi,lepolind,Vse,V,nkap,nngrid,nphi,&
        !$ theta,s,pflag,vfi_tricks,tricks_k_down,tricks_k_up,tricks_n_down,tricks_n_up,&
        !$ apolse1ind,apolse0ind,display_iter)
        !$omp do collapse(2)
        do j = 1,negrid  ! current state - epsilon
            do o = 1,ntgrid  ! current state - theta
                
                if (display_iter==2) then
                    write(*,'(1X,A,I2,4X,A,I2)') 'j = ', j, 'o = ',o 
                endif
                
                ! compute continuation value, for given eps and theta
                EV = 0.0d0 !dimension: (nagrid_fine)

                do m=1,negrid  !future state - epsilon
                    do n=1,ntgrid  !future state - theta
                        call linintv(agrid,V0(:, m, n),agrid_fine,V0_interp)
                        !V0_interp = interp1(agrid,V0(:, m, n),agrid_fine) 
                        EV = EV + V0_interp*Peps(j,m)*Ptheta(o,n) !go through all combinations

                    enddo  !future state - theta
                enddo  !future state - epsilon
                
                !Minimum and max assets choices (for workers and entre)
                min_aw = 1
                min_a0 = 1
                min_a1 = 1
                
                !!Minimum and max capital/labor choices for entrepreneurs 
                min_k = 1
                min_n = 1
                max_k = nkap
                max_n = nngrid
            
                !Introduce the state variable on assets
                do i=1,nagrid !current state - a
                    !----------------------------------------------------------------
                    !Start workers subproblem
                    decisw_temp = 0
                    
                    do il=1,nlgrid !labor choice (monotonicity TBA)
                        
                        lwval = lgrid(il)
                        yw  = w0*eps_grid(j)*lwval+ r0*agrid(i) !income of the worker (scalar) 
                        
                        !Max over next-period assets, given labor = lwval
                        utilw = large_negative ! utilw: scalar
                        do iap = min_aw,nagrid_fine   !future state a'
                            !consumption of worker, dimension: scalar
                            cw  = yw + agrid(i) - tax_work(yw)- agrid_fine(iap) 
                            if ( cw>0.0d0 ) then
                                vtempW = utilfun(cw,lwval) + beta*EV(iap)
                                if (vtempW > utilw) then 
                                    utilw    = vtempW
                                    RHSw(il) = utilw
                                    decisw_temp(il) = iap !best savings choice for each lwval
                                else
                                    if ( concavity==1 ) exit
                                endif !end max concavity a'
                            else
                                exit 
                            endif !end neg consumpt. check
                        enddo !end iap
                        
                    enddo !end il labor supply
                    
                    tempw     = maxloc(RHSw) !RHSw has dim [nlgrid]
                    Vw(i,j,o) = RHSw(tempw(1))
                    decisl(i,j,o) = tempw(1) !policy for labor supply
                    decisw(i,j,o) = decisw_temp(decisl(i,j,o)) !policy for a', given optimal l*
                    !End workers subproblem
                    !----------------------------------------------------------------
                    
                    !!----------------------------------------------------------------
                    !!Start entre subproblem
                    !!Assigning values to the value function VSE1(k,ixn,f,ile) and decision over assets (k,ixn,f,ile)
                    Vse1     = large_negative !scalar
                    decis1   = 1              !scalar
                    decis0   = 1              !scalar 
                    !
                    !!                                    
                    !!!Entrepreneurs - consumption and utility
                    do k=min_k,max_k   !1,nkap ! business capital - k
                        
                        kval = kgrid(k)
                        
                        ! Check collateral constraint k<=lambda*a(i)
                        if (kval <= lambda*agrid(i)) then ! collateral constr OK
                            
                            do ixn =min_n,max_n !1,nngrid !labor hiring
                                
                                nval = ngrid(ixn)
                                
                                do ile = 1,nlegrid !labor supply
                                    
                                    leval = legrid(ile)
                                    
                                    do f=1,nphi ! fraction hidden 
                                        
                                    phival = phigrid(f)
                                    
                                    ! income depends on (a,eps,theta,k,n,le) but not a'
                                    profit = theta(o)*prodfun(kval,leval,nval)-w0*nval-delta*kval-r0*kval ! scalar, pi greek
                                    income = profit+r0*agrid(i)-cost_evasion(phival) ! scalar, y^E
                                    
                                    profit_pos = max(profit,0.0d0) ! IMPORTANT!!!!!!
                                    income_pos = max(income,0.0d0) 
                                    
                                    evaded_taxes = &
                                    tax_entre(profit_pos + r0*agrid(i) ) - tax_entre(profit_pos*(1-phival) + r0*agrid(i) )
                                
                                    evaded_taxes_pos = max(evaded_taxes,0.0d0)
                                    
                                    Vd = large_negative
                                    ! Saving choice if detected
                                    do iap = min_a1,nagrid_fine
                                        ! c1 for being caught, dim: scalar
                                        c1 = income + agrid(i) - tax_entre(profit_pos*(1-phival) + r0*agrid(i) ) &
                                            - s*(evaded_taxes_pos) - agrid_fine(iap)
                                        if ( c1>0.0d0 ) then
                                            vtemp1 = utilfun(c1,leval) + beta*EV(iap) !detected, dimension: scalar
                                            !!! Check against current best
                                            if ( vtemp1 > Vd ) then
                                                Vd     = vtemp1
                                                decis1 = iap
                                            else
                                                if ( concavity==1 ) exit
                                            endif !end max concavity
                                        else
                                            ! if c1 is neg, exit loop over a'
                                            exit
                                        endif !end neg consumpt. check
                                    enddo !end iap
                                
                                    Vn = large_negative
                                    ! Saving choice if NOT detected
                                    do iap = min_a0,nagrid_fine
                                        ! c0 not caught, dim: scalar
                                        c0 = income + agrid(i) - tax_entre(profit_pos*(1-phival)+r0*agrid(i)) &
                                            - agrid_fine(iap)
                                        if ( c0>0.0d0 ) then
                                            vtemp0 = utilfun(c0,leval) + beta*EV(iap) !detected, dimension: scalar
                                            !!! Check against current best
                                            if ( vtemp0 > Vn ) then
                                                Vn     = vtemp0
                                                decis0 = iap
                                            else
                                                if ( concavity==1 ) exit
                                            endif !end max concavity
                                        else
                                            ! if c0 is neg, exit loop over a'
                                            exit
                                        endif !end neg consumpt. check
                                    enddo !end iap
                                
                                
                                    !Optimal saving choice for each k,n,phi
                                    ! decis1(k,ixn,f,ile) - detected
                                    ! decis0(k,ixn,f,ile) - not detected
                                    !Create the expected value with respect to p(k)
                            
                                    prob_audit_temp = prob_audit(phival,profit_pos,theta(o),kval,nval,leval,pflag)
                                    
                                    if ((prob_audit_temp<0.0d0) .or. (prob_audit_temp>1.0d0)) then
                                        write(*,*) 'pk is not in [0,1] '
                                        pause
                                    endif
                                    
                                    Vse1_temp       = prob_audit_temp*Vd + (1.0d0-prob_audit_temp)*Vn
                                    
                                    if (Vse1_temp > Vse1) then
                                        Vse1              = Vse1_temp
                                        Icap(i,j,o)       = k
                                        Inpol(i,j,o)      = ixn
                                        Iphi(i,j,o)       = f
                                        lepolind(i,j,o)   = ile
                                        apolse1ind(i,j,o) = decis1
                                        apolse0ind(i,j,o) = decis0
                                        
                                    endif 
                        
                                    enddo !END fraction hidden - phi
                                
                                enddo !END labor supply entre - ile
                            
                            enddo !END labor hiring - ixn
                        
                        else ! collateral constraint is violated
                            !Exit the do loop for k
                            exit
                        endif ! END IF collateral
                        
                    enddo  ! business capital - k
                    
                    !Store the optimal value of Vse1 in V^E
                    Vse(i,j,o) = Vse1
                    
                    !----------------------------------------------------------------
                    !END entre subproblem
                    !
                    
                    !Write down current occupational choice
                    ! either fixed or optimal
                    V(i,j,o) = max(Vw(i,j,o),Vse(i,j,o)) 
                    !if ( Vse(i,j,o) > Vw(i,j,o) ) then
                    !    occpol(i,j,o) = 1
                    !endif
                    
                    !!Update the minimum choices for assets for next iteration
                    !!Assumes policy function a'(a) is monotone 
                    !!VERY SAFE
                    min_aw = max(decisw(i,j,o),1)
                    min_a1 = max(apolse1ind(i,j,o),1)
                    min_a0 = max(apolse0ind(i,j,o),1)
                    !
                    if (vfi_tricks == 1) then
                        !----------------------------------------------------------------!
                        !Update the minimum and max choices for capital/labor for next iteration
                        !Assumes policy function k(a),n(a) are continuous 
                        !NOT SAFE
                        min_k = max(Icap(i,j,o)-int(tricks_k_down*nkap),1)
                        max_k = min(Icap(i,j,o)+int(tricks_k_up*nkap),nkap)
                        min_n = max(Inpol(i,j,o)-int(tricks_n_down*nngrid),1)
                        max_n = min(Inpol(i,j,o)+int(tricks_n_up*nngrid),nngrid)
                        !----------------------------------------------------------------!
                    endif
                                        
                enddo  !current state - a
    
            enddo ! current state - theta
    enddo ! current state - epsilon
        !$omp enddo
        !$omp end parallel

        
        
 !--------------------------------------------------  HOWARD ------------------------------------------------------------------------!       
    if (howard_flag == 1) then    
        write(*,*) "-----------------------------------"
        write(*,*) "Howard acceleration"
        write(*,*) "-----------------------------------"
        call sub_howard  
    else
        write(*,*) "NO Howard"
    endif ! if howard option
 !-------------------------------------------------- END HOWARD ---------------------------------------------------------------------!       

        
        !Number of iterations           
        itervf = itervf + 1
        !Error for VFI
        diff = maxval( ( abs( V-V0 )  ) )
        
        if (display_iter==1 .OR. display_iter==2) then
            write(*,*) "-----------------------------------"
            write(*,'(a,I4)') " iteration = ", itervf
            write(*,'(a,f0.7)') " norm_total = ", diff
            write(*,*) "-----------------------------------"
        endif
        
        !Update the guess for the value function             
        V0 = V
        
    enddo !while loop ends here
    
    if (display_iter==2) then
        write(*,*) "============================================================="
        write(*,*) "VALUE FUNCTION FOUND"
        write(*,*) "============================================================="
    endif 
    
    !call cpu_time(t2)
    t2 = omp_get_wtime()
    write(*,*) "============================================================="
    write(*,*) 'VFI runs for',real(t2-t1),'seconds.'
    write(*,*) "============================================================="
    
    if (diff <= tolV) then
        write(*,*) 'VFI converged succesfully'
        exitflag_vfi = 0
    else
        write(*,*) 'VFI did not converge succesfully'
        exitflag_vfi = 1
    endif 
    
    !=================================================
    ! Compute policy functions
    !=================================================
    
    !! Policy function for occupational choice
    if (fix_occpol == 0) then
        occpol = 0 ! dimension: (nagrid,negrid,ntgrid)
        where (Vse > Vw) occpol = 1 ! occpol=1 if agent is ENTRE
    endif 

    ! pol for asset a(t+1), workers
    ! decisw is defined on coarse grid and maps into fine grid
    apolw = 0.0d0
    lpolw = 0.0d0
    do i=1,nagrid
        do j=1,negrid
            do t=1,ntgrid
            
                apolw(i,j,t) = agrid_fine(decisw(i,j,t))
                lpolw(i,j,t) = lgrid(decisl(i,j,t))
                
            enddo
        enddo
    enddo
    
    !pol for consumption "cpolw", workers. Dim: (nagrid,negrid,ntgrid)
    
    cpolw = 0.0d0
    do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                yw           = w0*eps_grid(j)*lpolw(i,j,t)+r0*agrid(i)
                cpolw(i,j,t) = yw+agrid(i)-tax_work(yw)-apolw(i,j,t) 
            enddo
        enddo
    enddo
    
    
        
    ! Policy functions for a(t+1), self-employed, detected or not
    ! apolse1ind = 0 ! index, dimension: (nagrid,negrid,ntgrid) 
    ! apolse0ind = 0 ! index, dimension: (nagrid,negrid,ntgrid)
    ! Now compute the value on the appropriate grid

     do t=1,ntgrid
         do j=1,negrid
             do i=1,nagrid
                 apolse1(i,j,t) = agrid_fine(apolse1ind(i,j,t))
                 apolse0(i,j,t) = agrid_fine(apolse0ind(i,j,t))
             enddo
        enddo
     enddo

    ! Policy functions for business capital, labor hiring and tax evasion
    ! Icap ! index, dimension: (nagrid,negrid,ntgrid)
    ! Inpol ! index, dimension: (nagrid,negrid,ntgrid)
    ! Iphi ! index, dimension: (nagrid,negrid,ntgrid)
    ! lepolind ! index, dimension: (nagrid,negrid,ntgrid)
    ! Now compute the values on the grids kgrid and phigrid
    
    policycap = 0.0d0
    policyn   = 0.0d0
    policyphi = 0.0d0
    lepol     = 0.0d0
    
    do t = 1,ntgrid       !current state theta
        do j = 1,negrid   !current state epsilon
            do i=1,nagrid !current state assets
                policycap(i,j,t) = kgrid(Icap(i,j,t))
                policyphi(i,j,t) = phigrid(Iphi(i,j,t))
                policyn(i,j,t)   = ngrid(Inpol(i,j,t))
                lepol(i,j,t)     = legrid(lepolind(i,j,t))
            enddo
        enddo
    enddo
    
    !pol for consumption, self-employed, conditional on detection.
    ! cpolse1, dimension: (nagrid,negrid,ntgrid)
    ! cpolse0, dimension: (nagrid,negrid,ntgrid)
    
    cpolse1 = 0.0d0
    cpolse0 = 0.0d0
    do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                kval   = policycap(i,j,t)
                leval  = lepol(i,j,t)
                nval   = policyn(i,j,t)
                !phival = phigrid(f)   !!!!!!!!!
                phival = policyphi(i,j,t)
                profit = theta(t)*prodfun(kval,leval,nval)-w0*nval-delta*kval-r0*kval 
                income = profit+r0*agrid(i)-cost_evasion(phival) 
                profit_pos       = max(profit,0.0d0) 
                evaded_taxes     = tax_entre(profit_pos+r0*agrid(i))-tax_entre(profit_pos*(1-phival)+r0*agrid(i))
                evaded_taxes_pos = max(evaded_taxes,0.0d0)
                ! 0-Not detected
                cpolse0(i,j,t)   = income+agrid(i)-tax_entre(profit_pos*(1-phival)+r0*agrid(i))-apolse0(i,j,t)
                ! 1-Detected
                cpolse1(i,j,t)   = income+agrid(i)-tax_entre(profit_pos*(1-phival)+r0*agrid(i))-s*(evaded_taxes_pos)-apolse1(i,j,t)
            enddo
        enddo
    enddo
    

    
    !=================================================
    ! Refine policy functions 
    ! Policy functions are  computed on the finer grid 
    ! for the distribution iteration
    ! using interpolation
    ! ================================================

    apolwdet = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    lpolwdet = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    apolse1det = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    apolse0det = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    policycapdet = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    policyndet   = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    lepoldet     = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    policyphidet = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    Vsedet = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    Vwdet = 0.0d0 ! dimension: (nagrid_dist,negrid,ntgrid)
    occpoldet = 0 ! integer, dimension: (nagrid_dist,negrid,ntgrid)
    cpolwdet   = 0.0d0 !real, dimension: (nagrid_dist,negrid,ntgrid)
    cpolse0det = 0.0d0 !real, dimension: (nagrid_dist,negrid,ntgrid)
    cpolse1det = 0.0d0 !real, dimension: (nagrid_dist,negrid,ntgrid)
    
    ! Linear interp for policy functions for assets and consumption
    do it=1,ntgrid
        do ie=1,negrid
            apolwdet(:,ie,it) = interp1(agrid,apolw(:,ie,it),agrid_dist)
            lpolwdet(:,ie,it) = interp1(agrid,lpolw(:,ie,it),agrid_dist)
            apolse1det(:,ie,it) = interp1(agrid,apolse1(:,ie,it),agrid_dist)
            apolse0det(:,ie,it) = interp1(agrid,apolse0(:,ie,it),agrid_dist)
            !-------!
            cpolwdet(:,ie,it)   = interp1(agrid,cpolw(:,ie,it),agrid_dist)
            cpolse0det(:,ie,it) = interp1(agrid,cpolse0(:,ie,it),agrid_dist)
            cpolse1det(:,ie,it) = interp1(agrid,cpolse1(:,ie,it),agrid_dist)
        enddo
    enddo

    ! Linear interp for policy function of capital, evasion, labor hiring and labor supply entre
    do it=1,ntgrid
        do ie=1,negrid
            policycapdet(:,ie,it) = interp1(agrid,policycap(:,ie,it),agrid_dist)
            policyndet(:,ie,it)   = interp1(agrid,policyn(:,ie,it),agrid_dist)
            policyphidet(:,ie,it) = interp1(agrid,policyphi(:,ie,it),agrid_dist)
            lepoldet(:,ie,it)     = interp1(agrid,lepol(:,ie,it),agrid_dist)
        enddo
    enddo

    ! For value functions, NONlinear interp schemes don't work well
    do ie=1,negrid
        do it=1,ntgrid
    
            Vsedet(:,ie,it) = interp1(agrid,Vse(:,ie,it),agrid_dist)
            Vwdet(:,ie,it)  = interp1(agrid,Vw(:,ie,it),agrid_dist)
    
        enddo
    enddo

    if (fix_occpol == 0) then ! standard case (no decompositions)
        do ie=1,negrid
            do it=1,ntgrid
                do ia=1,nagrid_dist

                    if ( Vsedet(ia,ie,it) >= Vwdet(ia,ie,it) ) then

                        occpoldet(ia,ie,it) = 1

                    endif

                enddo
            enddo
        enddo
    else  ! decompositions, use fixed occpol
        do ie=1,negrid
            do it=1,ntgrid

                ! do interpolation
                occpoldet(:,ie,it) = interp1(agrid,dble(occpol(:,ie,it)),agrid_dist)
            enddo
        enddo

    endif
    

    
end subroutine sub_vfi
    