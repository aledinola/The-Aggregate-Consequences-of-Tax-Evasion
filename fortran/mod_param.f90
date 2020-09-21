module mod_param
	implicit none
    !This version: v15
    
	!Declaring variables:
	!--------------------

	!Declaring parameters:	
    integer, parameter :: dp = kind(0.0d0)
    real(dp), parameter :: large_negative = -1d100 !-10.0d0**25
    
    !Decalring basic dimensions for input:
    integer :: dimensions(2), real_params_dim, int_params_dim
    
    !Declaring the grids:
    integer :: nagrid, nagrid_fine, nagrid_dist, negrid, ntgrid, nphi, nkap, nlegrid
    
    !Integers:
	integer :: maxitV, maxit_dist, evade_only_profit, mu_method, temp(4), tempw(1), temp_phi(2), max_eval_k
	integer :: i, j, ii, o, k, f, l, itervf, m, ixn, ie, it, ia, t, ix, iy, iz, il, ixl, iap, n, ile
    integer :: iter_mu, jp, tp, jj_d, jj_nd, jj, left, right, vfi_tricks, howard_flag, accmax
    integer :: min_aw, min_a0, min_a1, max_aw, max_a0, max_a1, min_k, min_n, max_k, max_n
    integer :: kopt_ongrid(1), aoptw_ongrid(1), aopt1_ongrid(1), aopt0_ongrid(1), display_howard, display_mu, display_iter
    integer :: fix_occpol, fix_kpol, fix_apol_nd, nlgrid, nngrid, concavity, exitflag_vfi, exitflag_mu, pflag
    integer :: par_fortran
 
	!Real variables:
    real(dp) :: tolV, r0, w0, s, lambda, t1, t2, t_start, t_end, Vd, Vn, tricks_eval_k, Vse1, Vse1_temp
    real(dp) :: sigma, beta, tol_dist, vi, delta, tricks_a_down, tricks_a_up, tricks_k_down, tricks_k_up, tricks_n_up, tricks_n_down
    real(dp) :: yw, income, profit, income_pos, profit_pos, evaded_taxes, evaded_taxes_pos, diff, kopt, nopt, phiopt
    real(dp) :: c0, c1, vtemp1, vtemp0, cw, utilw, vtempW, prob_audit_temp
    real(dp) :: diff_d, p_k, aopt1,aopt0, aoptw, aopt, check, weight_jj_d, weight_jj_nd, weight_jj
    real(dp) :: sigma2, psi, gamma, pn_1, pn_2, pn_3, cc0, cc1, cc2, decis1, decis0

    !Decalring allocatable integer variables:
    integer, allocatable :: int_params(:)
    integer, allocatable :: decisw(:,:,:), decisl(:,:,:), occpol(:,:,:), occpoldet(:,:,:)
    integer, allocatable :: Iphi(:,:,:), Icap(:,:,:), Inpol(:,:,:)
    integer, allocatable :: apolse1ind(:,:,:), apolse0ind(:,:,:)
    integer, allocatable :: occpoldet_vec(:), occpol_vec(:), decisw_temp(:)
    integer, allocatable :: occpol_fixed_0(:), occpol_fixed(:,:,:)
    integer, allocatable :: Icap_fixed_0(:), Icap_fixed(:,:,:), Icap_vec(:)
            
    !Declaring allocatable real variables:
    real(dp), allocatable :: real_params(:), V0_0(:), V0(:,:,:), V(:,:,:)
    real(dp), allocatable :: agrid(:), agrid_fine(:), agrid_dist(:)
    real(dp), allocatable :: eps_grid(:), theta(:), phigrid(:), kgrid(:), ngrid(:)
    real(dp), allocatable :: Peps0(:), Peps(:,:), Ptheta0(:), Ptheta(:,:) !, pk(:)
    real(dp), allocatable :: lgrid(:), hgrid(:), legrid(:), cpolw(:,:,:), cpolse1(:,:,:), cpolse0(:,:,:)
    real(dp), allocatable :: cpolwdet(:,:,:), cpolse1det(:,:,:), cpolse0det(:,:,:)
    
    real(dp), allocatable :: Vw(:,:,:), Vse(:,:,:)
    real(dp), allocatable :: util1(:), util0(:)
    real(dp), allocatable :: EV(:), V0_interp(:), apolwdet(:,:,:) 
    real(dp), allocatable :: apolse1det(:,:,:), apolse0det(:,:,:), lpolwdet(:,:,:)
    real(dp), allocatable :: policyphidet(:,:,:), policycapdet(:,:,:), policyndet(:,:,:), lepoldet(:,:,:)
    real(dp), allocatable :: Vsedet(:,:,:), Vwdet(:,:,:)
    real(dp), allocatable :: apolw(:,:,:), lpolw(:,:,:), apolse1(:,:,:), apolse0(:,:,:)
	real(dp), allocatable :: policycap(:,:,:), policyn(:,:,:), policyphi(:,:,:), k_anal(:,:)
    real(dp), allocatable :: RHSw(:), RHS1(:), RHS0(:), lepolind(:,:,:), lepol(:,:,:)
    real(dp), allocatable :: policycap_fixed_0(:), policycap_fixed(:,:,:) 
    real(dp), allocatable :: apolse0_fixed_note_0(:,:,:), apolse0_fixed_note(:,:,:)
    
    
    real(dp), allocatable :: mu(:,:,:), mu_next(:,:,:), mu_vec(:)

    real(dp), allocatable :: policycapdet_vec(:), policyphidet_vec(:), lpolwdet_vec(:), apolwdet_vec(:)
    real(dp), allocatable :: policycap_vec(:), policyn_vec(:), policyphi_vec(:), policyndet_vec(:), lepoldet_vec(:)
    real(dp), allocatable :: apolw_vec(:), lpolw_vec(:), apolse1_vec(:), apolse0_vec(:) 
    real(dp), allocatable :: apolse1det_vec(:), apolse0det_vec(:), lepol_vec(:)
    real(dp), allocatable :: Vse_vec(:), Vw_vec(:), V_vec(:)
       
end module mod_param
