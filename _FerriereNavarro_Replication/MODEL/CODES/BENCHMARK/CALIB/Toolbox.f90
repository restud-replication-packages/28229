module Toolbox
    
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!
!!!!!     Search points in a vector     !!!!!
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vsearchprm(y,v,vprm) returns an integer with the two closest points in v to y
!
! Input:
!   y    = Ny x 1 vector
!   v    = Nv x 1 vector
!   vprm = real - curvature parameter in v
!
! Output:
!   ind = Ny x 1 vector
!
! Routine assumes v is constructed as:
!			v = vmin + (vmax-vmin)*(x**1/vprm)
! with vmin = v(1), vmax = v(Nv), where x is an Nv by 1 vector uniformly distributed between 0 and 1.
! vprm is the curvature: vprm = 0 for L-shaped, vprm = 1 for linear
!
!
! Note:
! if y(i) <= v(1) : ind(i) = 1
! if y(i) >= v(p) : ind(i) = p-1
! else            : v(ind(i)) <= y(i) < v(ind(i)+1)
!
! Please contact us if you find any errors:
! Axelle Ferriere, Paris School of Economics: axelle.ferriere@psemail.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: January 2021
! This version:  January 2021
function vsearchprm_FET(y,v,vprm) result(ind)

	implicit none
	real(8), intent(in)           :: y(:), v(:), vprm
	integer, dimension(size(y,1)) :: ind

	real(8), dimension(size(y,1)) :: yhat, nhat
	real(8) :: Dx, vmin, vmax
	integer :: Ny, Nv

	Ny = size(y,1); Nv = size(v,1);
	Dx = 1.0D0/dble(Nv-1)

	vmin = v(1); vmax = v(Nv)

	yhat = ((y-vmin)/(vmax-vmin))**vprm
	nhat = 1.0D0 + (yhat/Dx)

	ind = floor(nhat)

	where (y >= vmax) ind = Nv-1
	where (y <= vmin) ind = 1

end function vsearchprm_FET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!
!!!!!        Equally spaced grid        !!!!!
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function LINSPACE_FET(xlb,xub,Nx)  result(xout)

	implicit none

	real(8), intent(in) :: xlb, xub
	integer, intent(in) :: Nx

	real(8)             :: xout(Nx)

	real(8) :: dx
	integer :: ix

	dx = (xub-xlb)/(dble(Nx)-1.0D0)

	xout = [(xlb + dx*dble(ix-1), ix = 1, Nx)]


end function LINSPACE_FET


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!
!!!!!        Equally spaced grid        !!!!!
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Tauchen(rhox,sgx,Nx,mx,xvec,Px)

    implicit none
    real(8), intent(in)  :: rhox, sgx, mx
    integer, intent(in)  :: Nx
    real(8), intent(out) :: xvec(Nx), Px(Nx,Nx)

    real(8) :: dx, zz, xu, xl
    integer :: ix, ix_np

    real(8), external :: normcdf

    xvec(Nx) = mx*sqrt((sgx**2.0D0)/(1.0D0-(rhox**2.0D0)))
    xvec(1)  = - xvec(Nx)
    dx       = (xvec(Nx)-xvec(1))/dble(Nx-1)

    do ix = 2, Nx
        xvec(ix) = xvec(ix-1) + dx
    enddo

    ! Computing probabilities
    ! j = ix
    ! k = ix_np
    do ix = 1, Nx
	    zz        = (xvec(1) + (dx/2.0D0) - rhox*xvec(ix))/sgx
	    Px(ix,1)  = normcdf(zz)

	    do ix_np = 2, Nx-1
		    xu           = (xvec(ix_np) + (dx/2.0D0) - rhox*xvec(ix))/sgx
            xl           = (xvec(ix_np) - (dx/2.0D0) - rhox*xvec(ix))/sgx
            Px(ix,ix_np) = normcdf(xu) - normcdf(xl)
        enddo

        zz        = (xvec(Nx) - (dx/2.0D0) - rhox*xvec(ix))/sgx
	    Px(ix,Nx) = 1.0D0 - normcdf(zz)

		Px(ix,:) = Px(ix,:)/sum(Px(ix,:))
    enddo

    xvec = exp(xvec)

end subroutine Tauchen

! Notation
! Px(ix,ix_np) = move from ix to ix_np
    
end module Toolbox

