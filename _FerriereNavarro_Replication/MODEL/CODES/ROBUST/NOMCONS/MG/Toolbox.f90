module Toolbox
    
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                       !!!!!
!!!!!       Matrix Inverse Computation      !!!!!
!!!!!                                       !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! inv(A) computes the inverse of matrix A
!
! Input:
!   A = n x n matrix
!
! Output:
!   Ainv = n x n matrix, inverse of A
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: July 2016
! NEEDS UPDATE

function inv_FET(A) result(Ainv)

  real(8), dimension(:,:), intent(in)     :: A  
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for MKL
  integer, dimension(size(A,1)) :: ipiv  ! pivot indices
  integer :: n, info    

  ! Store A in Ainv to prevent it from being overwritten by MKL
  Ainv = A
  n    = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)
  ! call SGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)
  ! call SGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  
end function inv_FET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                       !!!!!   
!!!!!      Linear Least Square Routine      !!!!!
!!!!!                                       !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! coef_sol(B,A) solves for C in the system of equations A*C = B
! Input: 
!   A : an m x n matrix, with m >= n
!   B : an m x 1 vector, the rhs of the system
!
! Output:
!   C : an n x 1 vector
! The function uses the LAPACK routine DGELS().  See details here: https://software.intel.com/en-us/node/469160#EC9BE639-8638-4AF2-A4AC-74C9E0334883
! 
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: July 2016
! This version : Feb 2017
! NEEDS TESTING

subroutine coef_sol_FET(C,B,A)

    real(8), intent(in)  :: B(:), A(:,:)   
    real(8), intent(out) :: C(size(A,2))
    
    real(8) :: Aaux(size(A,1),size(A,2)) , Baux(size(B,1))    
    real(8) :: work( min(size(A,1),size(A,2))*10 + max(size(A,1),size(A,2))*65) 
    integer :: M, N, lwork, info, lwmax
    
    M = size(A,1)
    N = size(A,2)        
    
    lwmax = min(M,N)*10 + max(M,N)*65
    
    if (M < N) then
        stop 'ERROR: more unknowns than equations, least squares solution could not be computed!'
    end if    
    
    Baux = B; Aaux = A
    
    ! Query the optimal workspace    
    lwork = -1
    call DGELS('N', M, N, 1, Aaux, M, Baux, M, work, lwork, info)                
    lwork = min(lwmax, int(work(1)))
    ! print *, 'int(work(1)) = ', int(work(1))
    ! print *, 'lwmax = ', lwmax
    ! pause
    ! if ( lwork == lwmax) print *, 'increase lwmax in coef_sol'
    
    ! Solve for X in AX = B            
    call DGELS('N', M, N, 1, Aaux, M, Baux, M, work, lwork, info)    
    
    
     if (info > 0 ) then
         print *, 'ERROR: The algorithm computing SVD failed to converge;'
         print *, 'the least squares solution could not be computed!'
         stop
     else
         C = Baux         
     end if 

end subroutine coef_sol_FET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                               !!!!!
!!!!!       Kronecker Product       !!!!!
!!!!!                               !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! kron(A,B) computes kronecker product between A and B
!
! Input:
!   A = mA x nA matrix (real)
!   B = mB x nB matrix (real)
!
! Output:
!   K = (mA*mB) x (nA*nB) matrix (real)
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: July 2016
! This version : Feb 2017

subroutine kron_FET(K,A,B) 

    real(8), intent(in)  :: A(:,:), B(:,:)
    real(8), intent(out) :: K(size(A,1)*size(B,1),size(A,2)*size(B,2)) 
    
    integer :: mA, nA, mB, nB
    integer :: rA, cA, rB, cB, rK_l, rK_r, cK_a, cK_b
    
    mA = size(A,1)
    nA = size(A,2)
    mB = size(B,1)
    nB = size(B,2)
    
    K = 0.0D0
    
    do ra = 1, mA
    do ca = 1, nA
        
        rk_l = 1 + (ra-1)*mB
        rk_r = rk_l + mB - 1
        
        ck_a = 1 + (ca-1)*nB
        ck_b = ck_a + nB - 1
        
        K(rk_l:rk_r,ck_a:ck_b) = A(ra,ca)*B;
        
    end do
    end do
    
end subroutine kron_FET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                              !!!!!
!!!!!     Vectorized Kronecker     !!!!!
!!!!!                              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! kron_vec(A,B) computes the kronecker product of A and B row by row
! 
! Input:
!   A = rA x cA matrix (real)
!   B = rB x cB matrix (real)
!   need rA = rB, otherwise an error occurs
!
! Output:
!   K = rA x (cA*cB) matrix (real)
!   K(i,:) = kron(A(i,:),B(i,:))
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: July 2016
! This version : Feb 2017

subroutine kron_vec_FET(K,A,B)

    implicit none

    real(8), intent(in)  :: A(:,:), B(:,:)    
    real(8), intent(out) :: K(size(A,1),size(A,2)*size(B,2))     
    
    integer :: rA, cA, rB, cB
    integer :: icA, icB, ic, cK_l, cK_r    
    
    
    rA = size(A,1)
    cA = size(A,2)
    rB = size(B,1)
    cB = size(B,2)
        
    if (rA /= rB) then
        stop 'ERROR: size(A,1) \= size(B,1), cannot perform vectorized kronecker!'
    end if
    
    do icA = 1, cA                
        do icB = 1, cB                        
            ic = (icA - 1)*cB + icB       
            call vdmul( rA, A(:,icA), B(:,icB), K(:,ic) )
            ! K(:,ic) = A(:,icA)*B(:,icB)            
        end do        
    end do

end subroutine kron_vec_FET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                              !!!!!
!!!!!     element x element matrix multiplication  !!!!!
!!!!!                                              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mat_dot(K,A,B) computes the element x element multiplication og A and B and stores it in K
! 
! Input:
!   A = rA x cA matrix (real)
!   B = rB x cB matrix (real)
!   need rA = rB and cA = cB, otherwise an error occurs
!
! Output:
!   K = rA x cA (= rB x cB)  matrix (real)
!   K(i,j) = A(i,j)*B(i,j)
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: Feb 2017
! This version : Feb 2017


subroutine mat_dot_FET(K,A,B)

    implicit none
    real(8), intent(in)  :: A(:,:), B(:,:)
    real(8), intent(out) :: K(size(A,1),size(A,2))
    
    real(8) :: Avec(size(A,1)*size(A,2),1), Bvec(size(B,1)*size(B,2),1), Kvec(size(A,1)*size(A,2),1)
    
    integer :: rA, cA, rB, cB
    integer :: icA, icB, ic, cK_l, cK_r    
    
    
    rA = size(A,1)
    cA = size(A,2)
    rB = size(B,1)
    cB = size(B,2)
        
    if ( (rA /= rB) .or. (cA /= cB) ) then
        stop 'ERROR: size(A) \= size(B), cannot perform element x element multiplication !'
    end if
    
    Avec = reshape(A ,  (/rA*cA,1/) )
    Bvec = reshape(B ,  (/rB*cB,1/) )
    
    call vdmul( rA*cA , Avec , Bvec, Kvec)
    
    K = reshape(Kvec, (/rA, cA/))
    
    ! do ic = 1, cA
    !     
    !     call vdmul( rA, A(:,ic), B(:,ic), K(:,ic) )
    !     
    ! end do


end subroutine mat_dot_FET


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                              !!!!!
!!!!!              matrix multiplication           !!!!!
!!!!!                                              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! matmul_FET(A,B,C) computes the multiplication of A times B and stores it in C
! 
! Input:
!   A = rA x cA matrix (real)
!   B = rB x cB matrix (real)
!   need cA = rB otherwise an error occurs
!
! Output:
!   C = rA x cB   matrix (real)
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: Oct 2018
! This version : Oct 2018
! 

subroutine matmul_FET(A,B,C)

    implicit none
    real(8), intent(in)  :: A(:,:), B(:,:)
    real(8), intent(out) :: C(size(A,1),size(B,2))
    
    integer :: rA, cA, rB, cB
    ! integer :: icA, icB, ic, cK_l, cK_r    
    
    
    rA = size(A,1)
    cA = size(A,2)
    rB = size(B,1)
    cB = size(B,2)
    
    if (cA .NE. rB ) then
        stop 'ERROR: columns(A) \= rows(B), cannot multiply AxB!'        
    end if
    
    call dgemm('N', 'N', rA, cB, cA, 1.0D0, A, rA, B, rB, 0.0D0, C, rA)

end subroutine matmul_FET


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                              !!!!!
!!!!!       matrix x vector multiplication         !!!!!
!!!!!                                              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! matvec_FET(A,V,Y) computes the multiplication of matrix A times vector V and stores it in vector Y
! 
! Input:
!   A = rA x cA matrix (real)
!   V = rV x 1  vector (real)
!   need cA = rV otherwise an error occurs
!
! Output:
!   Y = rA x 1   vector (real)
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
!
! First version: Oct 2018
! This version : Oct 2018
! 

subroutine matvec_FET(A,V,Y)

    implicit none
    real(8), intent(in)  :: A(:,:), V(:)
    real(8), intent(out) :: Y(size(A,1))
    
    integer :: rA, cA, rV
    
    rA = size(A,1)
    cA = size(A,2)
    rV = size(V)
    
    if (cA .NE. rV ) then
        stop 'ERROR: columns(A) \= rows(V), cannot multiply AxV!'        
    end if
    
    call dgemv('N', rA, cA, 1.0D0, A, rA, V, 1, 0.0D0, Y, 1)
    
    ! call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    

end subroutine matvec_FET




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!    
!!!!!     Search points in a vector     !!!!!   
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! vsearch(X,v) returns an integer with the two closest points in v to X
!
! Input: 
!   X = Nx x 1 vector
!   v = p  x 1 vector
!
! Output:
!   ind = Nx x 1 vector
!
! Note:
! if X(i) <= v(1) : ind(i) = 1
! if X(i) >= v(p) : ind(i) = p-1
! else            : v(ind(i)) <= X(i) < v(ind(i)+1)
!
! Please contact us if you find any errors:
! Axelle Ferriere, European University Institute: axelle.ferriere@eui.eu
! Gaston Navarro, Federal Reserve Board: gaston.m.navarro@frb.gov
! 
! First version: July 2016
! NEEDS UPDATE

function vsearch_FET(X,v) result(ind)

    real(8), intent(in)           :: X(:), v(:)
    integer, dimension(size(X,1)) :: ind
    
    integer, dimension(size(X,1)) :: i_lb, i_ub, i_mid, flag
    real(8), dimension(size(X,1)) :: v_mid 
    integer :: Nx, p, ssss
    
    p = size(v,1)
    
    flag = 1
    ssss=0d0
    where( X >= v(p) )
        
   
        ind  = p  -1
        flag = 0 
    end where
    where(X <= v(1))
        ind  = 1
        flag = 0
    end where    
    
    i_lb = 1
    i_ub = p
    i_mid = floor( dble(i_lb + i_ub)/2.0d0 )    
    
    v_mid = v(i_mid)    
   ! ssss=0
    do while(maxval(flag)>0) 
    ! ssss=   ssss+1
    ssss=ssss+1
    Sdrama=0d0
    if (ssss >100) stop 'Something is wrong in vsearch, maybe NAN'   
    ! print*, ssss
        where(flag == 1 .and. X < v_mid) 
            i_ub  = i_mid
            i_mid = ceiling( dble(i_lb+i_ub)/2.0d0 ) 
            v_mid = v(i_mid)
        end where
        where(flag == 1 .and. X >= v_mid) 
            i_lb  = i_mid
            i_mid = floor( dble(i_lb+i_ub)/2.0d0)
            v_mid = v(i_mid)
        end where
        where(i_ub == i_lb + 1) 
            flag = 0
            ind = i_lb
        end where

    end do
    
end function vsearch_FET    

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                   !!!!!
!!!!!         compute percentiles       !!!!!
!!!!!                                   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function percofx(x,omga,Nx,pvec,Np) result(xout)

    implicit none
    real(8), intent(in)  :: x(:), omga(:), pvec(:)
    integer, intent(in)  :: Np, Nx

    real(8)              :: xout(Np)

    real(8) :: omgax(Nx), indx(Nx), xmax(1), xmin(1), xlb, xub, xm, pp, Fx, errp, tolp, tolx
    integer :: ip, itp, maxitp

    omgax = omga/sum(omga(:))

    xmin = minval(x); xmax = maxval(x)
    xlb  = xmin(1)  ; xub  = xmax(1)

    tolp = 1e-5; maxitp = 1000; tolx = 1e-8
    do ip = 1, Np
        print *, 'perc ',ip,'/',Np

        pp = pvec(ip)
        do itp = 1, maxitp
            xm = 0.5D0*(xlb+xub)

            indx = 0.0D0
            where(x <= xm)
                indx = 1.0D0;
            end where
            Fx = sum(indx*omgax)

            errp = Fx - pp
            ! print *, 'itp = ', itp,', errp = ', errp
            if ( (abs(errp)<tolp) .or. (xub-xlb<tolx) ) then
                xout(ip) = xm
                
                !---bounds for next iteration
                xub = xmax(1)
                exit
            else
                if (Fx < pp) then   ! low mass, need to increase xm
                    xlb = xm
                else
                    xub = xm
                endif
            endif
            
            if (itp == maxitp) then
                print *, 'cannot find pp = ', pp
                print *, 'xlb = ', xlb, ', xub = ', xub
                print *, 'Fx = ', Fx
                xout(ip) = xm
                STOP
            endif
        enddo

    enddo

end function percofx

    
end module Toolbox

