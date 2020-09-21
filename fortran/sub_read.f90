subroutine sub_read
	use mod_param
	implicit none
	
   !Guess for the value function:
    open(15, file='V0.txt', status='old')
	read(15,*) V0_0
	close(15)
	V0 = reshape(V0_0,[nagrid, negrid, ntgrid])
    
   !Used for decompositions:
    if (fix_occpol == 1) then
        open(15, file='occpol_fixed_note.txt', status='old')
	    read(15,*) occpol_fixed_0
	    close(15)
	    occpol_fixed = reshape(occpol_fixed_0,[nagrid, negrid, ntgrid])
    endif
    
    if (fix_kpol == 1) then
        
        open(15, file='policycap_fixed_note.txt', status='old')
	    read(15,*) policycap_fixed_0
	    close(15)
	    policycap_fixed = reshape(policycap_fixed_0,[nagrid, negrid, ntgrid])
        
        open(15, file='Icap_fixed_note.txt', status='old')
	    read(15,*) Icap_fixed_0
	    close(15)
	    Icap_fixed = reshape(Icap_fixed_0,[nagrid, negrid, ntgrid])
    
    endif
    
    if (fix_apol_nd == 1) then
        
        open(15, file='apolse0_fixed_note.txt', status='old')
	    read(15,*) apolse0_fixed_note_0
	    close(15)
	    apolse0_fixed_note = reshape(apolse0_fixed_note_0,[nagrid, negrid, ntgrid])
        
    endif
    
    
    !------------
    
    
   	!Reading in the 3 grids for assets:
    open(15, file='agrid.txt', status='old')
	read(15,*) agrid
	close(15)

    open(15, file='agrid_fine.txt', status='old')
	read(15,*) agrid_fine
	close(15)
	
    open(15, file='agrid_dist.txt', status='old')
	read(15,*) agrid_dist
	close(15)

    !Reading in the grid for epsilon:
    open(15, file='eps_grid.txt', status='old')
	read(15,*) eps_grid
	close(15)

    !Reading in the grid for theta:
    open(15, file='theta.txt', status='old')
	read(15,*) theta
	close(15)

    !Reading in the grid for phi (fraction of income taxes evaded):
    open(15, file='phigrid.txt', status='old')
	read(15,*) phigrid
	close(15)

    !Reading in the grid for kapital (business capital):
    open(15, file='kgrid.txt', status='old')
	read(15,*) kgrid
	close(15)
	
	!Reading in transition matrix for epsilon:
	open(15, file='Peps.txt', status='old')
	read(15,*) Peps0
	close(15)
    Peps = reshape(Peps0,[negrid, negrid])
	
	!Reading in the transition matrix for theta:
	open(15, file='Ptheta.txt', status='old')
	read(15,*) Ptheta0
	close(15)
    Ptheta = reshape(Ptheta0,[ntgrid, ntgrid])
	
    !Reading in the grid for labor supply
	open(15, file='lgrid.txt', status='old')
	read(15,*) lgrid
	close(15)
    
    !Reading in the grid for labor hirings:
	open(15, file='ngrid.txt', status='old')
	read(15,*) ngrid
	close(15)

    !Reading in the grid for labor supply entre:
	open(15, file='legrid.txt', status='old')
	read(15,*) legrid
	close(15)
    

end subroutine sub_read
