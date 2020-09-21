    program main
    use mod_param
    use mod_baselib
    use omp_lib
    
    !Authors: Di Nola, Kocharkov, Scholl, Tkhir. 
    !Updated on September 17, 2020.
    
	implicit none
	
	!Time
	!call cpu_time(t1)
    t_start = omp_get_wtime()
			
    !=================================================
    ! Reading in parameters and objects
    !=================================================
    
    !Reading in the dimensions:
	open(15, file='dimensions.txt', status='old')
	read(15,*) dimensions
	close(15)

    real_params_dim = dimensions(1)
	int_params_dim  = dimensions(2)
	
	!Allocate
	allocate(real_params(real_params_dim))
	allocate(int_params(int_params_dim))

    !Reading in the real parameters for the problems:
    open(15, file='real_params.txt', status='old')
	read(15,*) real_params
	close(15)
    
    !Reading in the integer parameters for the problems:
    open(15, file='int_params.txt', status='old')
	read(15,*) int_params
	close(15)
	
	!Reading in basic parameters
    beta         = real_params(1)
    delta        = real_params(2)
    tolV         = real_params(3)
    r0           = real_params(4)
    w0           = real_params(5)
    lambda_work  = real_params(6) 
    tau_work     = real_params(7)
    tau_s        = real_params(8)
    lambda_entre = real_params(9)
    tau_entre    = real_params(10)
    s            = real_params(11)
    lambda       = real_params(12)
    sigma        = real_params(13)
    tol_dist     = real_params(14)
    vi           = real_params(15)
    tricks_a_down     = real_params(16)
    tricks_a_up       = real_params(17)
    tricks_k_down     = real_params(18)
    tricks_k_up       = real_params(19)
    tricks_eval_k     = real_params(20)
    uncmean_eps       = real_params(21)
    b_work  = real_params(22)
    p_work  = real_params(23)
    s_work  = real_params(24)
    b_entre = real_params(25)
    p_entre = real_params(26)
    s_entre = real_params(27)
    taxrate_work  = real_params(28) 
    taxrate_entre = real_params(29) 
    sigma2        = real_params(30)
    psi           = real_params(31)
    gamma         = real_params(32)
    pn_1          = real_params(33) 
    pn_2          = real_params(34) 
    tricks_n_up   = real_params(35)
    tricks_n_down = real_params(36)
    pn_3          = real_params(37)
    cc0           = real_params(38)
    cc1           = real_params(39)
    cc2           = real_params(40)
    
    nagrid              = int_params(1)
    nagrid_fine         = int_params(2)
    nagrid_dist         = int_params(3)
    negrid              = int_params(4)
    ntgrid              = int_params(5)
    maxitV              = int_params(6)
    evade_only_profit   = int_params(7)
    nphi                = int_params(8)
    nkap                = int_params(9)
    mu_method           = int_params(10)
    vfi_tricks          = int_params(11)
    howard_flag         = int_params(12)
    accmax              = int_params(13)
    taxfunc             = int_params(14)
    maxit_dist          = int_params(15)
    fix_occpol          = int_params(16)
    fix_kpol            = int_params(17)
    fix_apol_nd         = int_params(18)
    nlgrid              = int_params(19) 
    nngrid              = int_params(20)
    concavity           = int_params(21)
    pflag               = int_params(22)
    nlegrid             = int_params(23)
    display_howard      = int_params(24)
    display_mu          = int_params(25)
    display_iter        = int_params(26)
    par_fortran         = int_params(27)
    

	!Allocate
	allocate(V0_0(nagrid*negrid*ntgrid))
	allocate(V0(nagrid,negrid,ntgrid))
    ! for decompositions
    allocate(occpol_fixed_0(nagrid*negrid*ntgrid))
	allocate(occpol_fixed(nagrid,negrid,ntgrid))
    allocate(policycap_fixed_0(nagrid*negrid*ntgrid))
	allocate(policycap_fixed(nagrid,negrid,ntgrid))
    allocate(Icap_fixed_0(nagrid*negrid*ntgrid))
	allocate(Icap_fixed(nagrid,negrid,ntgrid))
    
    
    !----
	allocate(agrid(nagrid))
    allocate(lgrid(nlgrid))
	allocate(agrid_fine(nagrid_fine))
	allocate(agrid_dist(nagrid_dist))	
	allocate(eps_grid(negrid))
	allocate(theta(ntgrid))
	allocate(phigrid(nphi))
	allocate(kgrid(nkap))
    allocate(ngrid(nngrid))
    allocate(legrid(nlegrid))
	allocate(Peps0(negrid*negrid))
	allocate(Peps(negrid,negrid))
	allocate(Ptheta0(ntgrid*ntgrid))
	allocate(Ptheta(ntgrid,ntgrid))
	
    !Reading routine
	call sub_read
	
	!Deallocate
	deallocate(real_params)
	deallocate(int_params)
	deallocate(V0_0)
    deallocate(occpol_fixed_0)
    deallocate(policycap_fixed_0)
	deallocate(Peps0)
	deallocate(Ptheta0)
					
    !=================================================
    ! Value function iteration
    !=================================================
    
    !Allocate
    allocate(decisw(nagrid,negrid,ntgrid))
    allocate(decisl(nagrid,negrid,ntgrid))
    allocate(occpol(nagrid,negrid,ntgrid))
    allocate(occpoldet(nagrid_dist,negrid,ntgrid))
    allocate(Iphi(nagrid,negrid,ntgrid))
    allocate(Icap(nagrid,negrid,ntgrid))
    allocate(Inpol(nagrid,negrid,ntgrid))
    
    allocate(apolse1ind(nagrid,negrid,ntgrid))
    allocate(apolse0ind(nagrid,negrid,ntgrid))
    allocate(Vw(nagrid,negrid,ntgrid))
    allocate(V(nagrid,negrid,ntgrid))
    allocate(Vse(nagrid,negrid,ntgrid))
    allocate(EV(nagrid_fine))
    allocate(V0_interp(nagrid_fine))
    
    allocate(apolwdet(nagrid_dist,negrid,ntgrid))
    allocate(lpolwdet(nagrid_dist,negrid,ntgrid))
    allocate(apolse1det(nagrid_dist,negrid,ntgrid))
    allocate(apolse0det(nagrid_dist,negrid,ntgrid))
    allocate(policycapdet(nagrid_dist,negrid,ntgrid))
    allocate(policyndet(nagrid_dist,negrid,ntgrid))
    allocate(lepoldet(nagrid_dist,negrid,ntgrid))
    allocate(policyphidet(nagrid_dist,negrid,ntgrid))
    allocate(Vsedet(nagrid_dist,negrid,ntgrid))
    allocate(Vwdet(nagrid_dist,negrid,ntgrid))
    
    allocate(cpolwdet(nagrid_dist,negrid,ntgrid))
    allocate(cpolse0det(nagrid_dist,negrid,ntgrid))
    allocate(cpolse1det(nagrid_dist,negrid,ntgrid))
    
    allocate(apolw(nagrid,negrid,ntgrid))
    allocate(lpolw(nagrid,negrid,ntgrid))
    allocate(apolse1(nagrid,negrid,ntgrid))
    allocate(apolse0(nagrid,negrid,ntgrid))
    allocate(policycap(nagrid,negrid,ntgrid))
    allocate(policyn(nagrid,negrid,ntgrid))
    allocate(lepolind(nagrid,negrid,ntgrid))
    allocate(lepol(nagrid,negrid,ntgrid))
    allocate(policyphi(nagrid,negrid,ntgrid))
    allocate(cpolw(nagrid,negrid,ntgrid))
    allocate(cpolse1(nagrid,negrid,ntgrid))
    allocate(cpolse0(nagrid,negrid,ntgrid))
    allocate(RHSw(nlgrid))
    allocate(decisw_temp(nlgrid))
    
    !call sub_vfi
    call sub_vfi()
    
    
    !Deallocate
    if (allocated(V0))         deallocate(V0)
    if (allocated(decisw))     deallocate(decisw)    
    if (allocated(Iphi))       deallocate(Iphi)
    if (allocated(apolse1ind)) deallocate(apolse1ind)
    if (allocated(apolse0ind)) deallocate(apolse0ind)
    if (allocated(util1))      deallocate(util1)
    if (allocated(util0))      deallocate(util0)
    if (allocated(EV))         deallocate(EV)
    if (allocated(V0_interp))  deallocate(V0_interp)
    if (allocated(Vsedet))     deallocate(Vsedet)
    if (allocated(Vwdet))      deallocate(Vwdet)
    if (allocated(k_anal))     deallocate(k_anal)
    if (allocated(RHSw))       deallocate(RHSw)
    if (allocated(decisw_temp)) deallocate(decisw_temp)
    if (allocated(RHS1))       deallocate(RHS1)
    if (allocated(RHS0))       deallocate(RHS0)
    
    
    
    !!=================================================
    !! Stationary distribution iteration
    !!=================================================
    !
    !!Allocate
    allocate(mu(nagrid_dist,negrid,ntgrid))
    allocate(mu_next(nagrid_dist,negrid,ntgrid))
        
    !Distribution iteration routine
    call sub_dist()
    
    !Deallocate
    deallocate(mu_next)
    deallocate(kgrid)
    deallocate(agrid_dist)
    deallocate(Peps)
    deallocate(Ptheta)
    
    
    !=================================================
    ! Writing out objects
    !=================================================
    
    !Allocate arrays for output
    allocate(mu_vec(nagrid_dist*negrid*ntgrid) )
    allocate(occpoldet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(policycapdet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(policyphidet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(Icap_vec(nagrid_dist*negrid*ntgrid) )
    allocate(occpol_vec(nagrid*negrid*ntgrid) )
    allocate(policycap_vec(nagrid*negrid*ntgrid) )
    allocate(policyn_vec(nagrid*negrid*ntgrid) )
    allocate(lepol_vec(nagrid*negrid*ntgrid) )
    allocate(policyndet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(lepoldet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(policyphi_vec(nagrid*negrid*ntgrid) )
    allocate(apolw_vec(nagrid*negrid*ntgrid) )
    allocate(lpolw_vec(nagrid*negrid*ntgrid) )
    allocate(lpolwdet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(apolwdet_vec(nagrid_dist*negrid*ntgrid) )
    allocate(apolse1_vec(nagrid*negrid*ntgrid) )
    allocate(apolse0_vec(nagrid*negrid*ntgrid) )
    allocate(apolse1det_vec(nagrid_dist*negrid*ntgrid) )
    allocate(apolse0det_vec(nagrid_dist*negrid*ntgrid) )
    allocate(Vw_vec(nagrid*negrid*ntgrid) )
    allocate(Vse_vec(nagrid*negrid*ntgrid) )
    allocate(V_vec(nagrid*negrid*ntgrid) )
    
    !!Writing routine
    call sub_write
    
    ! Deallocate
    deallocate(mu)
    deallocate(occpoldet)
    deallocate(policycapdet)
    deallocate(policyndet)
    deallocate(policyphidet)
    deallocate(apolw)
    deallocate(lpolw)
    deallocate(apolse1)
    deallocate(apolse0)
    deallocate(occpol)
    deallocate(policycap)
    deallocate(policyn)
    deallocate(lepol_vec)
    deallocate(lepoldet_vec)
    deallocate(policyphi)
    deallocate(Vw)
    deallocate(Vse)
    deallocate(V)
    deallocate(mu_vec)
    deallocate(occpoldet_vec)
    deallocate(policycapdet_vec)
    deallocate(policyphidet_vec)
    deallocate(occpol_vec)
    deallocate(policycap_vec)
    deallocate(policyn_vec)
    deallocate(policyndet_vec)
    deallocate(policyphi_vec)
    deallocate(apolw_vec)
    deallocate(lpolw_vec)
    deallocate(lpolwdet)
    deallocate(lpolwdet_vec)
    deallocate(apolwdet_vec)
    deallocate(apolse1_vec)
    deallocate(apolse0_vec)
    deallocate(apolse1det)
    deallocate(apolse0det)
    deallocate(apolse1det_vec)
    deallocate(apolse0det_vec)
    deallocate(Vw_vec)
    deallocate(Vse_vec)
    deallocate(V_vec)
    deallocate(agrid)
    deallocate(lgrid)
    deallocate(agrid_fine)
    deallocate(eps_grid)
    deallocate(theta)
    deallocate(phigrid)
    deallocate(cpolw)
    deallocate(cpolse1)
    deallocate(cpolse0)
    deallocate(cpolwdet)
    deallocate(cpolse1det)
    deallocate(cpolse0det)
    deallocate(lepoldet)
    deallocate(apolwdet)
    deallocate(Icap)
    deallocate(Inpol)
    deallocate(lepol)
    deallocate(lepolind)
    deallocate(decisl)

    !Time
    !call cpu_time(t2)
    t_end = omp_get_wtime()
    
    write(*,*) "============================================================="
    !write(*,*) 'FORTRAN executable runs for',real(t2-t1),'seconds.'
    write(*,*) 'FORTRAN executable runs for',real(t_end-t_start),'seconds.'
    write(*,*) "============================================================="
    
	!pause
        	
end program main
