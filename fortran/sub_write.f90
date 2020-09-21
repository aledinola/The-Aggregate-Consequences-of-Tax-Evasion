subroutine sub_write
	use mod_param
	implicit none
	
    ! Write the outputs as txt files (to be imported into matlab) 
    
    ! Vectorize all arrays
    mu_vec = reshape(mu, [nagrid_dist*negrid*ntgrid] )
    
    occpoldet_vec    = reshape(occpoldet, [nagrid_dist*negrid*ntgrid] )
    policycapdet_vec = reshape(policycapdet, [nagrid_dist*negrid*ntgrid] )
    policyphidet_vec = reshape(policyphidet, [nagrid_dist*negrid*ntgrid] )
    
    apolw_vec   = reshape(apolw, [nagrid*negrid*ntgrid] )
    lpolw_vec   = reshape(lpolw, [nagrid*negrid*ntgrid] )
    apolse1_vec = reshape(apolse1, [nagrid*negrid*ntgrid] )
    apolse0_vec = reshape(apolse0, [nagrid*negrid*ntgrid] )
    
    apolse1det_vec = reshape(apolse1det, [nagrid_dist*negrid*ntgrid] )
    apolse0det_vec = reshape(apolse0det, [nagrid_dist*negrid*ntgrid] )
    policyndet_vec = reshape(policyndet, [nagrid_dist*negrid*ntgrid] )
    lepoldet_vec   = reshape(lepoldet, [nagrid_dist*negrid*ntgrid] )
    lpolwdet_vec   = reshape(lpolwdet, [nagrid_dist*negrid*ntgrid] )
    apolwdet_vec   = reshape(apolwdet, [nagrid_dist*negrid*ntgrid] )
    
    occpol_vec    = reshape(occpol, [nagrid*negrid*ntgrid] )
    policycap_vec = reshape(policycap, [nagrid*negrid*ntgrid] )
    policyn_vec   = reshape(policyn, [nagrid*negrid*ntgrid] )
    policyphi_vec = reshape(policyphi, [nagrid*negrid*ntgrid] )
    lepol_vec     = reshape(lepol, [nagrid*negrid*ntgrid] )
    Icap_vec      = reshape(Icap, [nagrid*negrid*ntgrid] )
    Vw_vec        = reshape(Vw, [nagrid*negrid*ntgrid] )
    Vse_vec       = reshape(Vse, [nagrid*negrid*ntgrid] )
    V_vec         = reshape(V, [nagrid*negrid*ntgrid] )
    
    !Write
	open(15, file='mu.txt', status='replace')
    do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) mu_vec(i)
    enddo
	close(15)
    
    open(15, file='occpoldet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) occpoldet_vec(i)
    enddo
	close(15)
    
    !---------------- Capital policy functions -----------------------!
    !Integer
    open(15, file='Icap.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) Icap_vec(i)
    enddo
	close(15)
    !Real
    open(15, file='policycap.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) policycap_vec(i)
    enddo
	close(15)
    !Real, finer grid
    open(15, file='policycapdet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) policycapdet_vec(i)
    enddo
	close(15)
    
    !---------------- Labor hiring policy functions -----------------------!
    !Integer
    open(15, file='Inpol.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                write(15,*) Inpol(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    !Real
    open(15, file='policyn.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) policyn_vec(i)
    enddo
	close(15)
    !Real, finer grid
    open(15, file='policyndet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) policyndet_vec(i)
    enddo
	close(15)
    
    !---------------- Labor supply SE -----------------------!
    !Integer
    open(15, file='lepolind.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                write(15,*) lepolind(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    !Real
    open(15, file='lepol.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) lepol_vec(i)
    enddo
	close(15)
    !Real, finer grid
    open(15, file='lepoldet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) lepoldet_vec(i)
    enddo
	close(15)
    
    !---------------- Labor supply Workers -----------------------!
    !Integer
    open(15, file='decisl.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                write(15,*) decisl(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    !Real
    open(15, file='lpolw.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) lpolw_vec(i)
    enddo
	close(15)
    !Finer grid
    open(15, file='lpolwdet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) lpolwdet_vec(i)
    enddo
	close(15)
    
    open(15, file='policyphidet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) policyphidet_vec(i)
    enddo
	close(15)
    
    open(15, file='apolwdet.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) apolwdet_vec(i)
    enddo
	close(15)
    
    open(15, file='occpol.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) occpol_vec(i)
    enddo
	close(15)
    
    open(15, file='policyphi.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) policyphi_vec(i)
    enddo
	close(15)
    
    open(15, file='apolw.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) apolw_vec(i)
    enddo
	close(15)
    
    open(15, file='apolse1.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) apolse1_vec(i)
    enddo
	close(15)
    
    open(15, file='apolse0.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) apolse0_vec(i)
    enddo
	close(15)
    
    open(15, file='apolse1det.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) apolse1det_vec(i)
    enddo
	close(15)
    
    open(15, file='apolse0det.txt', status='replace')
	do i=1,nagrid_dist*negrid*ntgrid
	write(15,*) apolse0det_vec(i)
    enddo
	close(15)
    
    open(15, file='Vw.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) Vw_vec(i)
    enddo
	close(15)
    
    open(15, file='Vse.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) Vse_vec(i)
    enddo
	close(15)
    
    open(15, file='V.txt', status='replace')
	do i=1,nagrid*negrid*ntgrid
	write(15,*) V_vec(i)
    enddo
	close(15)
    
    open(15, file='cpolw.txt', status='replace')
    !Vectorize cpolw in column order (as using (:) in Matlab) 
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                write(15,*) cpolw(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    
    open(15, file='cpolse0.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                write(15,*) cpolse0(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    
    open(15, file='cpolse1.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid
                write(15,*) cpolse1(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    
    open(15, file='cpolwdet.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid_dist
                write(15,*) cpolwdet(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    
    open(15, file='cpolse1det.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid_dist
                write(15,*) cpolse1det(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    
    open(15, file='cpolse0det.txt', status='replace')
	do t=1,ntgrid
        do j=1,negrid
            do i=1,nagrid_dist
                write(15,*) cpolse0det(i,j,t)
            enddo
        enddo
    enddo
    close(15)
    
    open(15, file='exitflag_vfi.txt', status='replace')
	write(15,*) exitflag_vfi
	close(15)
    
    open(15, file='exitflag_mu.txt', status='replace')
	write(15,*) exitflag_mu
	close(15)
              
end subroutine sub_write
