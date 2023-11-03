module ModuleTRANSITION

	use GLOBALS
	use ModuleQERRORS_TR
	use ModuleSAVE
	use Toolbox
	
contains



    subroutine Compute_sequential_decomp

    implicit none

    real(8) :: fkp1(Nprm), xk(Nprm), errf(1), tolf

    integer, parameter :: maxit = 5500
    integer :: it, itsave, itax, tt, is, ib, ih
    real(8) :: wwx, dbs, err_t(maxit), wwy, cnew, cold
    real(8) :: ymn, yl(NS), rba(NS), yk(NS), ymnDATA, cvMtoD, dy, mpc(NS,Nbeta)

    ibar = -1.0D0 + (1.0D0+rk)*Pibar

    ! Steady-state

    xk(0*T_TR+1:1*T_TR) = rk
    xk(1*T_TR+1:2*T_TR) = Pibar
    xk(2*T_TR+1:3*T_TR) = wgeH
    xk(3*T_TR+1:4*T_TR) = lbd
    xk(4*T_TR+1:5*T_TR) = 0d0   ! dividends



    if (index_decomp==1) then
        open(1,  file = '../CG/OUTPUT/wgeH_TR.txt',    status = 'unknown')
        open(2,  file = '../CG/OUTPUT/Pi_TR.txt',    status = 'unknown')
        read(1, *) wgeH_TR
        read(2, *) Pi_TR
        close(1); close(2);
        open(1,  file = '../CG/OUTPUT/rk_TR.txt',    status = 'unknown')
        read(1, *) rk_TR
        close(1)
        open(1,  file = '../CG/OUTPUT/div_TR.txt',    status = 'unknown')
        read(1, *) div_TR
        close(1)
        open(1, file = '../CG/OUTPUT/lbd_TR.txt', status = 'unknown')
        open(2, file = '../CG/OUTPUT/gma_TR.txt', status = 'unknown')
        read(1, *) lbd_TR
        read(2, *) gma_TR
        close(1); close(2);

        open(1, file = '../CG/OUTPUT/G_TR.txt', status = 'unknown')
        read(1, *) G_TR
        close(1)

        xk(0*T_TR+1:1*T_TR) = rk_TR
        xk(1*T_TR+1:2*T_TR) = Pi_TR
        xk(2*T_TR+1:3*T_TR) = wgeH_TR
        xk(3*T_TR+1:4*T_TR) = lbd_TR
        xk(4*T_TR+1:5*T_TR) = div_TR! dividends

        fkp1 = ED_TR(xk)

        call SaveAggregates_TR_decomp
        call SaveVFandPolicies_TR

    elseif (index_decomp == 2) then


        open(1,  file = '../CG/OUTPUT/wgeH_TR.txt',    status = 'unknown')
        open(2,  file = '../CG/OUTPUT/Pi_TR.txt',    status = 'unknown')
        read(1, *) wgeH_TR
        read(2, *) Pi_TR
        close(1); close(2);
        open(1,  file = '../CG/OUTPUT/rk_TR.txt',    status = 'unknown')
        read(1, *) rk_TR
        close(1)
        open(1,  file = '../CG/OUTPUT/div_TR.txt',    status = 'unknown')
        read(1, *) div_TR
        close(1)

        G_TR = G
        gma_TR = gma
        lbd_TR = lbd

        xk(0*T_TR+1:1*T_TR) = rk_TR
        xk(1*T_TR+1:2*T_TR) = Pi_TR
        xk(2*T_TR+1:3*T_TR) = wgeH_TR
        xk(3*T_TR+1:4*T_TR) = lbd  ! xk(2*T_TR+1:3*T_TR) = xlbd
        xk(4*T_TR+1:5*T_TR) = div_TR   ! dividends

        fkp1 = ED_TR(xk)

        call SaveAggregates_TR_decomp
    

    elseif (index_decomp == 3) then


        open(1, file = '../CG/OUTPUT/lbd_TR.txt', status = 'unknown')
        open(2, file = '../CG/OUTPUT/gma_TR.txt', status = 'unknown')
        read(1, *) lbd_TR
        read(2, *) gma_TR
        close(1); close(2);


        xk(0*T_TR+1:1*T_TR) = rk
        xk(1*T_TR+1:2*T_TR) = Pibar
        xk(2*T_TR+1:3*T_TR) = wgeH
        xk(3*T_TR+1:4*T_TR) = lbd_TR  ! xk(2*T_TR+1:3*T_TR) = xlbd
        xk(4*T_TR+1:5*T_TR) = 0d0

        fkp1 = ED_TR(xk)

        call SaveAggregates_TR_decomp

    endif 

    end subroutine Compute_sequential_decomp
	

end module ModuleTRANSITION
