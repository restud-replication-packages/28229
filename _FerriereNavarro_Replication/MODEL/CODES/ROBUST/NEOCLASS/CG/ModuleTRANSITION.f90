module ModuleTRANSITION

	use GLOBALS
	use ModuleQERRORS_TR
	use ModuleSAVE
	use Toolbox
	
contains


    subroutine Compute_sequential

    implicit none

    real(8) :: fkp1(Nprm), xk(Nprm), errf(1), tolf

    integer, parameter :: maxit = 9500
    integer :: it, itsave, itax, tt, index_load
    real(8) :: wwx, dbs, err_t(maxit), wwy

    ibar = -1.0D0 + (1.0D0+rk)*Pibar

    index_load = 1 ! load previous guess

    !---initial point
    xk(0*T_TR+1:1*T_TR) = rk
    xk(1*T_TR+1:2*T_TR) = lbd  ! xk(2*T_TR+1:3*T_TR) = xlbd
    xk(2*T_TR+1:3*T_TR) = 0d0   ! dividends
    xk(3*T_TR+1:4*T_TR) = 1.0d0   ! qk

    if (index_load == 1) then
        open(1, file = 'OUTPUT/xkLOAD.txt',     status = 'unknown')
        read(1, *) xk
        close(1)
    endif

    !---itearation
    wwx = 0.99D0; tolf = 1e-5; itsave = 5; wwy = 0.9995D0

    do it = 1, maxit

        !---evaluate guess
        fkp1 = ED_TR(xk)

        errf   = maxval(abs(fkp1))
        print *, '****************************'
        print *, '***  it, err = ', it, errf
        print *, '****************************'

        if (errf(1) < 4e-5) then
            wwx = 0.99D0
!            wwy = 0.995D0
        endif

        if (errf(1) < 5e-5) then
            wwx = 0.995D0
 !           wwy = 0.999D0
        endif

        if (errf(1) < tolf ) then
            print *, 'equilibrium found'
            call SaveAggregates_TR
            exit
        elseif (errf(1) < 7e-5 .and. it > 12000) then
            print *, 'exiting with errf(1) = ', errf(1)
            call SaveAggregates_TR
            exit
        else
            !--- update xk

            xk(0*T_TR+1:1*T_TR) = wwx*rk_TR  + (1d0-wwx)*rkb_TR
            xk(1*T_TR+1:2*T_TR) = wwx*lbd_TR + (1d0-wwx)*lbdb_TR
            xk(2*T_TR+1:3*T_TR) = wwx*div_TR + (1d0-wwx)*divb_TR
            xk(3*T_TR+1:4*T_TR) = wwy*qk_TR  + (1d0-wwy)*qkb_TR



        endif

        err_t(it) = errf(1)
        if (itsave <= 0    ) then

            print *, '*** saving results ***'

            call SaveAggregates_TR

            !---save errors
            open(1,  file = 'OUTPUT/err_t.txt',     status = 'unknown')
            do itax = 1, it
                write(1, '(*(f18.6))') ( err_t(itax) )
            enddo
            close(1)

            !!! save current guess
            open(2,  file = 'OUTPUT/xk.txt',     status = 'unknown')
            do tt = 1, Nprm
                write(2, '(*(f18.6))') ( xk(tt) )
            enddo
            close(2)
            open(1, file = 'OUTPUT/xkLOAD.txt',     status = 'unknown')
            write(1, *) xk
            close(1)
            itsave = 5
        endif
        itsave = itsave-1

    enddo ! end loop it


    end subroutine Compute_sequential
	

end module ModuleTRANSITION
