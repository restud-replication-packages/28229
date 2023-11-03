module MODULElpe

	use GLOBALS
	use FUNCTIONS_lpe
	use ModuleCOMPUTATIONS_lpe

contains

	!-------------------------------------------	
	!---lpe computation
	!-------------------------------------------
	subroutine COMPUTElpe
	
	implicit none
	real(8) :: lpeH_agg, lpeL_agg
	real(8) :: DTshock, rhoDT	
	integer :: tt

    print *, 'Compute LPE_tau'
	
	!---tax shock
	DTshock = 0.01D0; rhoDT = 0.9D0
	DTAX_TR(1) = DTshock 
	do tt = 2, T_TR
		DTAX_TR(tt) = rhoDT*DTAX_TR(tt-1)
	enddo
	
	!---prices and taxes
	wgeH_TR = wgeH
	r_TR    = rk
	qk_TR   = qk
	tauk_TR = tauk
	lbd_TR  = lbd
	gma_TR  = gma	
	TF_TR   = TF
	
	
	!---compute transition
	call Model_VF_TR
	call ComputeMEASURE_TR
	call COMPUTE_HH_TR
	
	call COMPUTE_lpex
	
	
	!---compute lpes
	lpeH_agg = (100.0D0*(Hagg_TR(1)-Hagg)/Hagg)/(100.0d0*DTAX_TR(1))
	lpeL_agg = (100.0D0*(Lagg_TR(1)-Lagg)/Lagg)/(100.0d0*DTAX_TR(1))
	
!	print *, 'Hagg, Hagg_TR(1)', Hagg, Hagg_TR(1)
	print *, 'lpeH_agg = ', lpeH_agg
!	print *, 'Lagg, Lagg_TR(1)', Lagg, Lagg_TR(1)
	print *, 'lpeL_agg = ', lpeL_agg

    print *, 'Reproduce Erosa et al.'

    !---tax shock
    DTAX_TR = 0.0D0

    !---prices and taxes
    wgeH_TR = wgeH
    wgeH_TR(1) = wgeH*1.01D0

    r_TR    = rk
    qk_TR   = qk
    tauk_TR = tauk
    lbd_TR  = lbd
    gma_TR  = gma
    TF_TR   = TF


    !---compute transition
    call Model_VF_TR
    call ComputeMEASURE_TR
    call COMPUTE_HH_TR


    !---compute lpes in response to 1% increase in pre-tax wage
    lpeH_agg = (100.0D0*(Hagg_TR(1)-Hagg)/Hagg)
    lpeL_agg = (100.0D0*(Lagg_TR(1)-Lagg)/Lagg)

    !    print *, 'Hagg, Hagg_TR(1)', Hagg, Hagg_TR(1)
    print *, 'lpeH_agg = ', lpeH_agg
    !    print *, 'Lagg, Lagg_TR(1)', Lagg, Lagg_TR(1)
    ! print *, 'lpeL_agg = ', lpeL_agg


	
	end subroutine COMPUTElpe
	
	
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

end module MODULElpe
