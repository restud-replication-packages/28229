module ModuleTRANSITION

	use GLOBALS
	use ModuleQERRORS_TR
	use ModuleSAVE
	use Toolbox
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!! 				Quasi-Newton Step				!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine Compute_QNEWTON
	
	implicit none
	real(8) :: xk(Nprm), xkp1(Nprm), xik(Nprm), yk(Nprm), fk(Nprm), fik(Nprm), fkp1(Nprm), Bjinvk(Nprm,Nprm), Bjinvkp1(Nprm,Nprm), Btf(Nprm), dx(Nprm), eax(1), ek, eik, eikm, ekp1, errf(1), tolf
	
	real(8)          :: JACx(Nprm,JxC)
	real(8)          :: dmult(T_TR)
	! real(8)          :: rkMAT(T_TR,JxC), wgeMAT(T_TR,JxC), qkMAT(T_TR,JxC), MCMAT(T_TR,JxC), YMAT(T_TR,JxC), YsMAT(T_TR,JxC), LMAT(T_TR,JxC), KdMAT(T_TR,JxC), BCMAT(T_TR,JxC), AMAT(T_TR,JxC), DMAT(T_TR,JxC)
	! real(8)          :: rkMATs(T_TR,Nprm), wgeMATs(T_TR,Nprm), qkMATs(T_TR,Nprm), MCMATs(T_TR,Nprm), YMATs(T_TR,Nprm), YsMATs(T_TR,Nprm), LMATs(T_TR,Nprm), KdMATs(T_TR,Nprm), BCMATs(T_TR,Nprm), AMATs(T_TR,Nprm), DMATs(T_TR,Nprm)
	integer          :: ic, ccin, ccou, loadxk
	character(len=5) :: ics
	
	integer, parameter :: maxit = 1500
	integer :: it, itbs, maxitbs, itsave, JACit, itax, tt, ivar, iprm, ipx
	real(8) :: wwx, dbs, err_t(maxit)
	
	!---steady-state objects
	lbdlb = 0.001D0; lbdub = 0.999D0
	
	ibar = -1.0D0 + (1.0D0+rk)*Pibar
	
	Dk   = dlta*Kagg
	
    !---initial point
    xk(0*T_TR+1:1*T_TR) = rk
    xk(1*T_TR+1:2*T_TR) = Pibar
    xk(2*T_TR+1:3*T_TR) = wgeH
    xk(3*T_TR+1:4*T_TR) = lbd  ! xk(2*T_TR+1:3*T_TR) = xlbd
    xk(4*T_TR+1:5*T_TR) = 0d0   ! dividends
	
	
	!---compute Jacobian
	loadxk = 1
	!---read case
	if ( COMMAND_ARGUMENT_COUNT() /= 1) then
		print *, 'ERROR: A command line argument specifying the parameter is required'
		STOP
	endif
	print *, 'about to read'
	call get_command_argument(1, iWxs)
	print *, 'done reading'
	read(iWxs,*) iWx
	print *, 'done writing iWx'
	
	if (iWx > 0) then
		call Compute_JACOBIAN_xc 	! call Compute_JACOBIAN_T
		! call Compute_JACOBIAN_full
		STOP
	else
		if (loadxk == 1) then
		
			open(1, file = 'OUTPUT/BjinvkLOAD.txt',    status = 'unknown')
			read(1, *) Bjinvk
			close(1)
			
			open(1, file = 'OUTPUT/xkLOAD.txt',    status = 'unknown')
			read(1, *) xk
			close(1)
			
			ED_ss = ED_TR(xk); fk = ED_ss
			print *, ' maxval(abs(fk)) = ', maxval(abs(fk))
		
		else
			ED_ss = ED_TR(xk); fk = ED_ss
			print *, ' maxval(abs(fk)) = ', maxval(abs(fk))
			
			ccin = 1
			do ic = 1, NC
				if (ic < 10) then
					write (ics,'(I1)') ic
				elseif (ic < 100) then
					write (ics,'(I2)') ic
				elseif (ic < 1000) then
					write (ics,'(I3)') ic
				else
					print *, 'need ics - at ic = ', ic
					STOP
				endif
				print *,'********************'
				print *, 'ics = ', ics
				print *, 'ic  = ', ic
				print *,'********************'
				open(1, file = 'OUTPUT/JACmatxcLOAD_'//trim(ics)//'.txt',    status = 'unknown')
				read(1, *) JACx
				close(1)
				
				ccou = ccin + JxC - 1
			
				JACmat(:,ccin:ccou) = JACx
			
				ccin = ccou+1
			enddo
			
			Bjinvk  = inv_FET(JACmat)
			
		
			open(1, file = 'OUTPUT/JACmat_kx.txt',    status = 'unknown')
			do iprm = 1, Nprm
				! print *, 'saving at iprm = ', iprm
				write(1, '(*(f18.6))') (JACmat(iprm,ipx), ipx = 1, Nprm)
			enddo
			close(1)
			
			open(1, file = 'OUTPUT/JACmat_kx.txt',    status = 'unknown')
			do iprm = 1, Nprm
				! print *, 'saving at iprm = ', iprm
				write(1, '(*(f18.6))') (JACmat(iprm,ipx), ipx = 1, Nprm)
			enddo
			close(1)

			open(1, file = 'OUTPUT/Bjinvk_kx.txt',    status = 'unknown')
			do iprm = 1, Nprm
				! print *, 'saving at iprm = ', iprm
				write(1, '(*(f18.6))') (Bjinvk(iprm,ipx), ipx = 1, Nprm)
			enddo
			close(1)

			open(1, file = 'OUTPUT/ED_ss.txt',    status = 'unknown')
			do iprm = 1, Nprm
				! print *, 'saving at iprm = ', iprm
				write(1, '(*(f18.6))') (fk(iprm) )
			enddo
			close(1)
		endif
		
		
	endif
	! STOP
	
	
	! ED_ss = ED_TR(xk)
	
	fk      = ED_ss
	eax     = maxval(abs(fk));
	ek      = eax(1)
		
	!---itearation
	wwx = 0.99D0; tolf = 4e-5
	itsave = 5; maxitbs = 5; dbs = 10.0D0; JACit = 1
	
	mult_lp = 0.0D0
	do it = 1, maxit
	
		!---new guess
		call matvec_FET(Bjinvk,fk,Btf) ! matvec_FET(A,V,Y)
		xkp1 = xk - Btf
		xkp1 = wwx*xk + (1.0D0-wwx)*xkp1
		
		dx   = xkp1 - xk
		
		eikm = 99900000000.0D0
		
		!---backstep check
		do itbs = 1, maxitbs
			xik = xk + dx
			fik = ED_TR(xik)
			
			eax = maxval(abs(fik))
			eik = eax(1)
			
			if ( eik < ek ) then 			! this was a good step, move forward
				! print *, ''
				exit
			elseif ( eikm < eik ) then		! last step was better, undo last shrinking step and exit
				dx  = dbs*dx
				exit
			endif
			
			! if you reached here, (eik > ek) but (eik < eikm) -> keep shrinking
			eikm = eik
			dx   = dx/dbs
			
		enddo ! end loop itbs
		
		xik = xk + dx
		
		xkp1 = xik; fkp1 = fik; ekp1 = eik;
		
		errf   = maxval(abs(fkp1))
		print *, '****************************'
		print *, '***  it, itbs = ', it, itbs
		print *, '*** errf = ', errf
		print *, '****************************'
				
		dmult   = mult - mult_lp
		mult_lp = mult
		print *, '--- err mult(1:40): ', maxval(abs(dmult(1:40))), ' in : ', maxloc(abs(dmult(1:40)))
		
		if (errf(1) < 1e-3) then
			wwx = 0.995D0
			maxitbs = 5
			JACit = 1
		endif
		
		if (errf(1) < 1e-4) then
!			JACit = 0
		endif
		
		if (errf(1) < tolf ) then
			print *, 'equilibrium found'
			call SaveAggregates_TR
			exit
		elseif (errf(1) < 7e-5 .and. it > 1200) then
			print *, 'exiting with errf(1) = ', errf(1)
			call SaveAggregates_TR

            exit
		else
			!--- update x and B
			if (JACit==1) then
				Bjinvkp1 = Update_JACOBIAN(xk,xkp1,fk,fkp1,Bjinvk)
			else
				Bjinvkp1 = Bjinvk
			endif

			!--- rename stuff
			xk     = xkp1
			fk     = fkp1
			Bjinvk = Bjinvkp1
			ek     = ekp1
			
		endif
		
		err_t(it) = errf(1)
		if (itsave <= 0	) then
		
			print *, '*** saving results ***'
		
			call SaveAggregates_TR
			
			!---save errors
			open(1,  file = 'OUTPUT/err_t.txt',     status = 'unknown')
			do itax = 1, it
				write(1, '(*(f18.6))') ( err_t(itax) )
			enddo
			close(1)
			
			!!! save current guess
			open(1,  file = 'OUTPUT/fk.txt',     status = 'unknown')
			open(2,  file = 'OUTPUT/xk.txt',     status = 'unknown')
			do tt = 1, Nprm
			    write(1, '(*(f18.6))') ( fk(tt) )
				write(2, '(*(f18.6))') ( xk(tt) )
			enddo
			close(1); close(2)
			open(1, file = 'OUTPUT/xkLOAD.txt',     status = 'unknown')
			open(2, file = 'OUTPUT/BjinvkLOAD.txt', status = 'unknown')
			write(1, *) xk
			write(2, *) Bjinvk
			close(1); close(2)
			itsave = 5
		endif
		itsave = itsave-1
	
	enddo ! end loop it
	
	
	end subroutine Compute_QNEWTON
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Transition computation, variables: {Y_t,Pi_t,lbd_t,Dk_t}
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine Compute_JACOBIAN_xc
	
	implicit none
	real(8) :: xk(Nprm), xx_t(Nprm), ED_t(Nprm), Jdstep
	real(8) :: JACx(Nprm,JxC)
	integer :: ic, cc, inn, tt, colm, iv1, iv2, iprm, ipx
	
	real(8) :: rkMAT(T_TR,JxC), wgeMAT(T_TR,JxC), qkMAT(T_TR,JxC), MCMAT(T_TR,JxC), YMAT(T_TR,JxC), YsMAT(T_TR,JxC), LMAT(T_TR,JxC), KdMAT(T_TR,JxC), BCMAT(T_TR,JxC), AMAT(T_TR,JxC), DMAT(T_TR,JxC)
	
	!---check if we have enough computers
	if (Nprm == NC*JxC) then
		print *, 'Nprm - NC*JxC = ', Nprm - NC*JxC
		print *, 'moving to Jacobian'
	else
		print *, 'Nprm - NC*JxC = ', Nprm - NC*JxC
		print *, 'something is wrong -- stopping'
		STOP
	endif
	
	!---steady-state as initial point
    xk(0*T_TR+1:1*T_TR) = rk
    xk(1*T_TR+1:2*T_TR) = Pibar
    xk(2*T_TR+1:3*T_TR) = wgeH
    xk(3*T_TR+1:4*T_TR) = lbd  ! xk(2*T_TR+1:3*T_TR) = xlbd
    xk(4*T_TR+1:5*T_TR) = 0d0   ! dividends
	
	ED_ss = ED_TR(xk)
	print *, 'maxval(abs(ED_ss(:))) = ', maxval(abs(ED_ss(:)))
		
	ic = iWx
	
	print *, '*************************************'
	print *, '*** JAC on computer ic  = ', ic, '/', NC, '***'
	print *, '*************************************'
	
	Jdstep = 1e-5; colm = 1
	do cc = 1, JxC
		print *, '*** JAC cc = ', cc,'/', JxC ,'***'
		inn = cc + JxC*(ic-1)   ! column for Jacobian

        xx_t      = xk

        if ((inn>=4*T_TR+1) .and. (inn<=5*T_TR)) then
            xx_t(inn) = xk(inn) +Jdstep
            ED_t      = ED_TR(xx_t)

            JACx(:,colm) = (ED_t-ED_ss)/Jdstep
        else
            xx_t(inn) = xk(inn)*(1.0D0+Jdstep)
            ED_t      = ED_TR(xx_t)

            JACx(:,colm) = (ED_t-ED_ss)/(Jdstep*xk(inn))
        endif
		colm = colm+1
	enddo
	
	open(1, file = 'OUTPUT/JACmatxcLOAD_'//trim(iWxs)//'.txt',    status = 'unknown')
    write(1, *) JACx
    close(1)
	
	
	end subroutine Compute_JACOBIAN_xc
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Transition computation, variables: {Y_t,Pi_t,lbd_t,Dk_t}
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine Compute_JACOBIAN_full
	
	implicit none
	real(8) :: xk(Nprm), xx_t(Nprm), ED_t(Nprm), Jdstep
	real(8) :: JACxF(Nprm,Nprm)
	integer :: inn, tt, iprm, ipx
	
	
	!---steady-state as initial point
    xk(0*T_TR+1:1*T_TR) = rk
    xk(1*T_TR+1:2*T_TR) = Pibar
    xk(2*T_TR+1:3*T_TR) = wgeH
    xk(3*T_TR+1:4*T_TR) = lbd  ! xk(2*T_TR+1:3*T_TR) = xlbd
    xk(4*T_TR+1:5*T_TR) = 0d0   ! dividends
	
	ED_ss = ED_TR(xk)
	print *, 'maxval(abs(ED_ss(:))) = ', maxval(abs(ED_ss(:)))
	
	Jdstep = 1e-5;
	do inn = 1, Nprm
		print *, '*** JAC full inn = ', inn,'/', Nprm ,'***'
		
		xx_t      = xk
		xx_t(inn) = xk(inn)*(1.0D0+Jdstep)
		ED_t      = ED_TR(xx_t)

		JACxF(:,inn) = (ED_t-ED_ss)/(Jdstep*xk(inn))
		
	enddo
	
	open(1, file = 'OUTPUT/JACmat_kx.txt',    status = 'unknown')
	do iprm = 1, Nprm
	! print *, 'saving at iprm = ', iprm
		write(1, '(*(f18.6))') (JACxF(iprm,ipx), ipx = 1, Nprm)
	enddo
	close(1)
			
	open(1, file = 'OUTPUT/JACmatxcLOAD_FULL.txt',    status = 'unknown')
    write(1, *) JACxF
    close(1)
	
	end subroutine Compute_JACOBIAN_full
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!---Update Jacobian
	function Update_JACOBIAN(xk,xkp1,fk,fkp1,Bk)  result(Bkp1)
	
	implicit none
	real(8), intent(in) :: xk(Kuvar), xkp1(Kuvar), fk(Kuvar), fkp1(Kuvar), Bk(Kuvar,Kuvar)
	real(8)             :: Bkp1(Kuvar,Kuvar)
	
	real(8) :: dfk(Kuvar), ukvec(Kuvar), dkvec(Kuvar), Btf(Kuvar), Bnum1(Kuvar,1), Bnum2(1,Kuvar), Bnum(Kuvar,Kuvar), dkvecT(1,Kuvar), denB
	
	!!! update x and B
	dfk   = fkp1 - fk
	call matvec_FET(Bk,dfk,ukvec)
	dkvec = xkp1 - xk

	!!! numerator computations
	Bnum1(:,1)  = dkvec - ukvec  ! Bnum1 = dk-uk -> a Nx1 "matrix"

	dkvecT(1,:) = dkvec
	call matmul_FET(dkvecT,Bk,Bnum2) 	   ! this is Bnum2 = dk^T * Bk  -> a 1xN "matrix"   ~~ matmul_FET(A,B,C)
	call matmul_FET(Bnum1,Bnum2,Bnum)      ! this is Bnum = (dk-uk)*(dk^T*Bk) -> a NxN matrix

	!!! denominator computations
	denB  = sum(dkvec*ukvec)

	!!! update Binv
	Bkp1 = Bk + (Bnum/denB)
	
	end function Update_JACOBIAN


end module ModuleTRANSITION
