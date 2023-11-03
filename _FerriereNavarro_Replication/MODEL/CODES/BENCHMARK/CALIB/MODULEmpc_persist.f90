module MODULEmpc_persist

	use GLOBALS
	use FUNCTIONS_lpe ! Use the same files as for the LPE to compute a mini-transition
	use ModuleCOMPUTATIONS_lpe

contains

	!-------------------------------------------	
	!---lpe computation
	!-------------------------------------------
	subroutine COMPUTEmpc_persist
	
	implicit none

    integer, parameter :: maxit = 5500
    integer :: it, itsave, itax, tt, is, ib, ih, ix
    real(8) :: wwx, dbs, err_t(maxit), wwy, cnew, cold, rhoG
    real(8) :: ymn, yl(NS), rba(NS), yk(NS), ymnDATA, cvMtoD, dy, mpc(NS,Nbeta)

    print *, 'Compute MPC_persistent'

	!---prices and taxes
	wgeH_TR = wgeH
	r_TR    = rk
	qk_TR   = qk
	tauk_TR = tauk
	lbd_TR  = lbd
	gma_TR  = gma	
	TF_TR   = TF

    !---compute ymean
    ymn = 0.0D0
    do ib = 1, Nbeta
    do ih = 1, Nh
    do is = 1, NS
        yl(is)  = wgeH*hvec(ih)*S(is,2)
        rba(is) = rk
        yk(is)  = rba(is)*S(is,1)
        
        ymn = ymn + (yl(is) + yk(is))*hpol(is,ib,ih)*mu(is,ib)
    enddo
    enddo
    enddo

    ymnDATA = 67000.0D0
    cvMtoD  = ymnDATA/ymn

    !---compute transfer dC and new policies
    dy = 500.0D0/cvMtoD            ! mpc out of dC dollars

    TF_TR(1) = TF + dy
    rhoG = 0.90D0
    do tt = 2, T_TR
        TF_TR(tt) = (1.0D0-rhoG)*TF + rhoG*TF_TR(tt-1)
    enddo

    DTAX_TR = 0d0

    !---compute transition
    call Model_VF_TR
    call ComputeMEASURE_TR
    call COMPUTE_HH_TR

    tt=1
    do is = 1, NS
    do ib = 1, Nbeta
        cnew = (hpol_TR(is,ib,tt,1)*cpol_TR(is,ib,tt,1) + hpol_TR(is,ib,tt,2)*cpol_TR(is,ib,tt,2))
        cold = (hpol(is,ib,1)*cpol(is,ib,1) + hpol(is,ib,2)*cpol(is,ib,2))
        mpc(is,ib) =  (cnew-cold)/dy
    enddo
    enddo

    print *, 'MPC out of a persistent transfer = ', (Cagg_TR(1:10)-Cagg)/(dy)

    open(1,  file = 'OUTPUT/mpc_persist.txt',   status = 'unknown')
    do is = 1, NS
        write(1, '(*(f18.6))') ( mpc(is,ib) , ib = 1,Nbeta )
    enddo
    close(1);
	
	end subroutine COMPUTEmpc_persist
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       				Compute VF       					  !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Model_VF_TR
    
    implicit none
    integer :: tt, ib, ih, is, iloc(NS), dcase4(NS)
    real(8) :: xl(NS), xr(NS), yl(NS), yk(NS), rba(NS), TAX(NS)   		
    integer :: clck_counts_beg, clck_counts_end, clck_rate
    
    ! initialize last period --> value and policies
    Vwrk_TR(:,:,T_TR)     = Vwrk
	Vnow_TR(:,:,T_TR)     = Vnow	
	do ih = 1, Nh
		hpol_TR(:,:,T_TR,ih)   = hpol(:,:,ih)
		apol_TR(:,:,T_TR,ih)   = apol(:,:,ih)
		cpol_TR(:,:,T_TR,ih)   = cpol(:,:,ih)
	enddo
    
    
    
    call system_clock ( clck_counts_beg, clck_rate )
    do tt = T_TR-1,1,-1
        ! if ( tt == 1 .or. tt == T_TR -1) then
        !    print *, 'tt = ', tt
        ! endif
        
        ! compute Ve		
        call ComputeVe_TR(tt+1)
        
        ! compute optimal policy -- scalar
        do ib = 1, Nbeta
            
            !$OMP PARALLEL DO DEFAULT(SHARED)
			do is = 1, NS				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!---case 1: no work -> hh = 0				
				yl(is)  = wgeH_TR(tt)*0.0D0*S(is,2)
				rba(is) = RbFUN_TR(S(is,1),tt)	
				yk(is)  = rba(is)*S(is,1)
				TAX(is) = TAXFUN_TR(yl(is),yk(is),tt)
				! TAX(is) = TAX(is) + DTAX_TR(tt)*abs(TAX(is)) 
				xr(is)	= yl(is) + yk(is) + S(is,1) + TF_TR(tt) - TAX(is)
				xl(is)  = amin				

				apol_TR(is,ib,tt,1) = SolveBella_TR(xl(is),xr(is),0.0D0,is,ib,tt)
				Vnow_TR(is,ib,tt)   = VFa_eval_TR(apol_TR(is,ib,tt,1),0.0D0,is,ib,tt,1)								
				
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!---case 2: work -> hh = hbar				
				yl(is)  = wgeH_TR(tt)*hbar*S(is,2)
				rba(is) = RbFUN_TR(S(is,1),tt)	
				yk(is)  = rba(is)*S(is,1)
				TAX(is) = TAXFUN_TR(yl(is),yk(is),tt)
				! TAX(is) = TAX(is) + DTAX_TR(tt)*abs(TAX(is)) 
				xr(is)	= yl(is) + yk(is) + S(is,1) + TF_TR(tt) - TAX(is)				
				xl(is)  = amin

				apol_TR(is,ib,tt,2) = SolveBella_TR(xl(is),xr(is),hbar,is,ib,tt)				
				Vwrk_TR(is,ib,tt)   = VFa_eval_TR(apol_TR(is,ib,tt,2),hbar,is,ib,tt,1)				
				
			enddo ! end loop is	            
            !$OMP END PARALLEL DO 			
			
        enddo  ! end loop ib                
       
        
    enddo  ! end loop tt
	
	! compute Ve tt=1 -> to compute hpol
	tt=1
    call ComputeVe_TR(tt)

    call system_clock ( clck_counts_end, clck_rate )
    ! print *, 'VF_TR seconds = ' , (clck_counts_end - clck_counts_beg) / real (clck_rate)
    
    end subroutine Model_VF_TR

end module MODULEmpc_persist
