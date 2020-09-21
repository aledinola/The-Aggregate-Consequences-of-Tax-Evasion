module mod_baselib
    
    use mod_param
    
    implicit none
    
    !Declaring parameters:	
    integer, parameter :: dp1 = kind(0.0d0)
    integer   :: taxfunc
    real(dp1) :: lambda_entre, tau_entre, lambda_work, tau_work, tau_s, uncmean_eps
    real(dp1) :: b_work, p_work, s_work, b_entre, p_entre, s_entre, taxrate_work, taxrate_entre
         
    contains
    
FUNCTION linint(x,y,xi)
! linear interpolation of function y on grid x at interpolation point xi
	IMPLICIT NONE
	REAL(dp1), DIMENSION(:), INTENT(IN) :: x,y
	REAL(dp1), INTENT(IN) :: xi
	REAL(dp1) :: linint
	REAL(dp1) :: a,b,d
	INTEGER :: n,i
	n=size(x)
	IF (size(y)/=n) THEN
		PRINT *, 'linint: x and y must be of the same size'
		STOP 'program terminated by linint'
	END IF
	i=max(min(locate(x,xi),n-1),1)
	d=x(i+1)-x(i)
	IF (d == 0.0) STOP 'bad x input in splint'
	a=(x(i+1)-xi)/d
	b=(xi-x(i))/d
	linint=a*y(i)+b*y(i+1)
END FUNCTION linint
!******************************************************************************!

SUBROUTINE linintv(x,y,xi,yi)
! linear interpolation of function y on grid x at interpolation vector xi
	IMPLICIT NONE
	REAL(dp1), DIMENSION(:), INTENT(IN)  :: x,y,xi
	REAL(dp1), DIMENSION(:), INTENT(OUT) :: yi
	REAL(dp1) :: a,b,d
	INTEGER :: m,n,i,j
	n=size(x)
	IF (size(y)/=n) THEN
		PRINT *, 'linintv: x and y must be of the same size'
		STOP 'program terminated by linintv'
	END IF
	m=size(xi)
	IF (size(yi)/=m) THEN
		PRINT *, 'linintv: xi and yi must be of the same size'
		STOP 'program terminated by linintv'
	END IF
	DO j=1,m
		i=max(min(locate(x,xi(j)),n-1),1)
		d=x(i+1)-x(i)
		IF (d == 0.0) THEN
			STOP 'bad x input in linintv'
		END IF
		a=(x(i+1)-xi(j))/d
		b=(xi(j)-x(i))/d
		yi(j)=a*y(i)+b*y(i+1)
	END DO
END SUBROUTINE linintv
!******************************************************************************!
   
PURE FUNCTION locate(xx,x)
	IMPLICIT NONE
	REAL(dp1), DIMENSION(:), INTENT(IN) :: xx
	REAL(dp1), INTENT(IN) :: x
	INTEGER :: locate
	INTEGER :: n,il,im,iu
	n=size(xx)
	il=0
	iu=n+1
	do
		if (iu-il <= 1) exit
		im=(iu+il)/2
		if (x >= xx(im)) then
			il=im
		else
			iu=im
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=il
	end if
END FUNCTION locate
!******************************************************************************!


    !******************************************************************************!
    pure function interp1(x,y,xi) result(yi)
        ! original code by  John Burkardt
        real(dp1), intent(in), dimension(:) :: x, y, xi
        real(dp1), dimension(:), allocatable :: yi
        integer :: j, nx, l, r
    
        nx = size(xi)
        allocate(yi(nx))
        doj:do j = 1,nx
           call bracket (x, xi(j), l, r)
           yi(j)= y(l)+(y(r)-y(l))*(xi(j)-x(l))/(x(r)-x(l))
        end do doj
    end function interp1
    
    !******************************************************************************!
    pure function interp1h(x,y,xi) result(yi)
        ! original code by  John Burkardt
        real(dp1), intent(in), dimension(:) :: x, y
        real(dp1), intent(in) :: xi  ! it is a scalar!!
        real(dp1) :: yi
        integer :: j, nx, l, r
    
        call bracket (x, xi, l, r)
        yi = y(l)+(y(r)-y(l))*(xi-x(l))/(x(r)-x(l))
       
    end function interp1h
    
    !******************************************************************************!
    pure subroutine bracket(x, xval, l, r)
        ! original code by  John Burkardt
        real(dp1), intent(in), dimension(:) :: x
        real(dp1), intent(in) :: xval
        integer, intent(out) :: l, r
        integer :: i, n
    
        n=size(x)
        do i = 2, n - 1
           if ( xval < x(i) ) then
              l = i - 1
              r = i
              return
           end if
        end do
        l = n - 1
        r = n
    end subroutine bracket
    
    !******************************************************************************!
    pure function utilfun(c,l) result(F)
        ! NOTE: maybe write also a separate utilfun for entre
        implicit none

        !Declare inputs:
        real(dp1), intent(in) :: c, l ! consumption and labor
        !Declare outputs:
        real(dp1) :: F 
        
        F = c**(1.0d0-sigma)/(1.0d0-sigma) - psi*l**(1.0d0+sigma2)/(1.0d0+sigma2)
        
       
    end function utilfun
    
    !******************************************************************************!
    pure function prodfun(k,le,n) result(F)
        ! NOTE: prodfun is defined without theta shock
        implicit none

        !Declare inputs:
        real(dp1), intent(in) :: k, le, n ! capital, labor supply, labor hirings
        !Declare outputs:
        real(dp1) :: F 
        
        F = (k**gamma*(le+n)**(1-gamma))**vi
        
       
    end function prodfun
    !******************************************************************************!
    
    pure function cost_evasion(phi) result(F)
        ! NOTE: prodfun is defined without theta shock
        implicit none

        !Declare inputs:
        real(dp1), intent(in) :: phi ! capital, labor supply, labor hirings
        !Declare outputs:
        real(dp1) :: F
        
        if (phi>0d0) then
            F = cc0+cc1*phi+cc2*phi**2
        else
            F = 0d0
        endif
        
       
    end function cost_evasion
    
    !******************************************************************************!
    pure function prob_audit(phi,profit,theta,k,n,le,flag) result(pk)
        ! INPUTS
        ! phi: fraction of misreported business income
        ! profit: business income
        implicit none

        !Declare inputs:
        real(dp1), intent(in) :: phi, profit, theta, k, n, le
        integer, intent(in) :: flag
        !Declare outputs:
        real(dp1) :: pk 
        ! Declare locals:
        real(dp1) :: temp ! 
        
        SELECT CASE (flag)
            
        CASE (1) !misreported income
            temp = phi*profit
            
        CASE (2) !f(theta,k,n)
            temp = theta*prodfun(k,le,n)
            
        CASE (3) !f(k,n)
            temp = prodfun(k,le,n)
        
        CASE (4) !k
            temp = k
            
        CASE (5) !n hired labor
            temp = n
            
        END SELECT
        
        pk   = pn_3 + (1.0d0-pn_3)/(1.0d0 + pn_1*exp(-pn_2*temp ))
        
       
    end function prob_audit
    
    !******************************************************************************!
    
    pure function tax_entre(income)
        implicit none

        !Declare inputs:
        real(dp1), intent(in) :: income ! income in levels !!!!!!
        !Declare outputs:
        real(dp1) :: tax_entre ! total taxes paid 
        ! Declare locals:
        real(dp1) :: y ! normalized income !!!!!!
        
        SELECT CASE (taxfunc)
            
        CASE (1) ! HSV
            ! Impose normalization
            y = income /(w0*uncmean_eps)
            tax_entre = y - lambda_entre*y**(1.0d0-tau_entre) + tau_s
        
        CASE (2) ! GS 
            tax_entre = income*(b_entre - b_entre*(s_entre*income**p_entre + 1)**(-1.0d0/p_entre)) + tau_s
            
        CASE (3) ! flat
            tax_entre = taxrate_entre*income + tau_s
            
        END SELECT
        
        return
    end function tax_entre 
    !******************************************************************************!
    pure function tax_work(income)
        implicit none

        !Declare inputs:
        real(dp1), intent(in) :: income ! income in levels 
        !Declare outputs:
        real(dp1) :: tax_work ! total taxes paid 
        ! Declare locals:
        real(dp1) :: y 
        
        SELECT CASE (taxfunc)
        
        CASE (1) ! HSV
            ! Impose normalization
            y = income /(w0*uncmean_eps)
            tax_work = y - lambda_work*y**(1.0d0-tau_work) + tau_s
            
        CASE (2) ! GS 
            tax_work = income*(b_work - b_work*(s_work*income**p_work + 1)**(-1.0d0/p_work)) + tau_s
            
        CASE (3) ! flat
            tax_work = taxrate_work*income + tau_s
            
        END SELECT
        
        return
    end function tax_work     
!******************************************************************************!

end module mod_baselib
