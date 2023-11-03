module ModuleQERRORS_TR

	use GLOBALS
	use ModuleCOMPUTATIONS_TR
	use FUNCTIONS_TR
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!         Compute Excess Demand - given Yagg_TR and Pi_TR         !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ED_TR(xxin)  result(ED_out) ! function ED_TR(Yin_TR,Piin_TR,xlbdin_TR)  result(ED_out)
	
	implicit none
	real(8), dimension(Nprm), intent(in) :: xxin	
	real(8), dimension(Nprm)             :: ED_out	
	
	real(8), dimension(T_TR)             :: albdin_TR
	real(8), dimension(T_TR)             :: ED_L, ED_A
	
	integer :: maxitlbd_t, itlbd_t, errloc(1), errlocTF(1), dlbd, tt
	real(8) :: errlbd_t(1), errTF_t(1), tollbd_t, wwlbd, wwTF, deltax
	
	maxitlbd_t = 1000; tollbd_t = 1e-5; wwlbd = 0.95D0; wwTF = 0.95D0  
	
	!---assign x to variables
	rk_TR         = xxin(0*T_TR+1:1*T_TR)
    Pi_TR        = xxin(1*T_TR+1:2*T_TR)
	wgeH_TR       = xxin(2*T_TR+1:3*T_TR)
	lbd_TR        = xxin(3*T_TR+1:4*T_TR)   ! xlbd_TR       = xxin(2*T_TR+1:3*T_TR)
    div_TR        = xxin(4*T_TR+1:5*T_TR)

    !---compute div and Tdx_TR
    do tt = 1, T_TR
        !---compute Tdx_TR
        deltax       = div_TR(tt)/Ex
        Tdx_TR(:,tt) = deltax*S(:,2)
    enddo


    ! Compute i_tr given Pi_t
    i_TR(1) = ibar
    do tt = 2, T_TR
        i_TR(tt) = ((Pi_TR(tt-1)/Pibar)**phiPI)*(1d0+ibar)-1d0
    enddo
    r_TR = (1d0+i_TR)/Pi_TR -1d0

    ! Compute qk_TR

    ! Compute qtk using rk and rtk
    qk_TR(T_TR) = qk
    do tt = T_TR-1,1,-1
        qk_TR(tt) = (qk_TR(tt+1)+rk_TR(tt+1))/(1d0+r_TR(tt+1))
    enddo


	!---compute HH problem
	call Model_VF_TR	
	
	!---compute measure and government's budget	
	call ComputeMEASURE_TR
	call Update_HHandTAXES_TR
    
	
	ED_out(0*T_TR+1:1*T_TR) = 1.0D0*(rk_TR - rkb_TR) ! Lagg_TR - Ld_TR
    ED_out(1*T_TR+1:2*T_TR) = qk_TR - qkb_TR
	ED_out(2*T_TR+1:3*T_TR) = 1.0d0*(wgeH_TR(1:T_TR) - wgeHb_TR) ! rk_TR - rb_TR
	ED_out(3*T_TR+1:4*T_TR) = BC_TR*1.0D0
    ED_out(4*T_TR+1:5*T_TR) = div_TR - divb_TR

	
	end function ED_TR
	
	
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
    Vwrk_TR(:,:,T_TR+1)     = Vwrk
	Vnow_TR(:,:,T_TR+1)     = Vnow	
	do ih = 1, Nh
		hpol_TR(:,:,T_TR+1,ih)   = hpol(:,:,ih)
		apol_TR(:,:,T_TR+1,ih)   = apol(:,:,ih)
		cpol_TR(:,:,T_TR+1,ih)   = cpol(:,:,ih)
	enddo
    
    
    call system_clock ( clck_counts_beg, clck_rate )
    do tt = T_TR,1,-1
        ! if ( tt == 1 .or. tt == T_TR -1) then
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
				xr(is)	= min(amax,yl(is) + yk(is) + S(is,1) + TF_TR(tt) + Tdx_TR(is,tt) - TAX(is))
				xl(is)  = amin				

				apol_TR(is,ib,tt,1) = SolveBella_TR(xl(is),xr(is),0.0D0,is,ib,tt)
				Vnow_TR(is,ib,tt)   = VFa_eval_TR(apol_TR(is,ib,tt,1),0.0D0,is,ib,tt,1)								
				
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!---case 2: work -> hh = hbar				
				yl(is)  = wgeH_TR(tt)*hbar*S(is,2)
				rba(is) = RbFUN_TR(S(is,1),tt)	
				yk(is)  = rba(is)*S(is,1)
				TAX(is) = TAXFUN_TR(yl(is),yk(is),tt)
				xr(is)	= min(amax,yl(is) + yk(is) + S(is,1) + TF_TR(tt) + Tdx_TR(is,tt) - TAX(is))
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

end module ModuleQERRORS_TR
