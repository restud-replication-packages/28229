module ModuleCOMPUTATIONS_TR

	use GLOBALS
	use FUNCTIONS_TR
	use Toolbox
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!  		Given Yagg_TR, Pi_TR, Dk_TR -> compute aggregates  		!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       					Compute MEASURE       				!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ComputeMEASURE_TR
    
    implicit none
    real(8) :: mm(NS,Nbeta), a_np(NS), da(NS)
    integer :: aind_np(NS)
    integer :: tt, ib, ih, is, ix_np, is_np
    
    
    
    mu_TR(:,:,1) = mu
    do tt = 1, T_TR-1
        
        mm = 0.0D0
        
        do ib = 1, Nbeta
			do ih = 1, Nh
				a_np    = min( apol_TR(:,ib,tt,ih) , amax)				
				aind_np = vsearchprm_FET(a_np,avec,avecpar)						
				da      = (a_np - avec(aind_np))/(avec(aind_np+1) - avec(aind_np))
            
				do is = 1, NS
					do ix_np = 1,Nx
						is_np = aind_np(is) + Na*( ix_np - 1)
						mm(is_np,ib)   = mm(is_np,ib)   + (1.0D0-da(is))*Px(xind(is),ix_np)*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
						mm(is_np+1,ib) = mm(is_np+1,ib) +    da(is)     *Px(xind(is),ix_np)*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
					enddo
				enddo
			enddo   
!            mm(:,ib)         = Pbeta(ib)*mm(:,ib)/sum(mm(:,ib))     
        enddo

        mu_TR(:,:,tt+1) = matmul(mm,Pbmat) !mm(:,ib)

     
    enddo
    
    end subroutine ComputeMEASURE_TR
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!  					Update HH side and Taxes		  			!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine Update_HHandTAXES_TR
	
	implicit none
	real(8) :: yl(NS), yk(NS), rba, taxl, taxlss, GTOT, deltax, I0, Iss, aty_tr, aty_ss, Ixx(T_TR)
	integer :: tt, ib, ih, is
	
	!---update HH side	
	Cagg_TR = 0.0D0
    Aagg_TR = 0.0D0    
    Lagg_TR = 0.0D0	
	! Anp     = 0.0D0
	
	Aagg_TR(1) = Aagg
	
	do tt = 1, T_TR
		do ib = 1, Nbeta
		do ih = 1, Nh
			Cagg_TR(tt) = Cagg_TR(tt) + sum(cpol_TR(:,ib,tt,ih)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt))
			Lagg_TR(tt) = Lagg_TR(tt) + sum(hvec(ih)*S(:,2)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt))

			Aagg_TR(tt+1) = Aagg_TR(tt+1) + sum(apol_TR(:,ib,tt,ih)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt)) 			
		enddo
		enddo
	enddo

	
	
	!---update governments budget
	RevK_TR     = 0.0D0    
    RevL_TR     = 0.0D0    
    Rev_TR      = 0.0D0     
	DB_TR       = 0.0D0	
	
	! INTEGRAL_TR = 0.0D0	 ! INTEGRAL_TR(tt) = INTEGRAL_TR(tt) + (yl**(1.0D0-gma_TR(tt)))*hprob_TR(1,is,ib,tt)*mu_TR(is,ib,tt)                
    DB_TR(1) = DB
    index_tax_TR = 0d0
	
	do tt = 1, T_TR
        I0=0d0
        Iss = 0d0
		do ib = 1, Nbeta
		do ih = 1, Nh
			do is = 1, NS
				yl(is)  = wgeH_TR(tt)*hvec(ih)*S(is,2)
				rba     = RbFUN_TR(S(is,1),tt)	
				yk(is)  = rba*S(is,1)
				
				RevK_TR(tt) = RevK_TR(tt) + tauk_TR(tt)*yk(is)*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
			
				
				if (yl(is)>0.0D0) then

                    if (tt>50) then
                        aty_tr = lbd_TR(tt)*(yl(is)**(1.0D0-gma_TR(tt)))

                        taxl            = yl(is) - aty_tr
                        RevL_TR(tt)     = RevL_TR(tt) + taxl*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
                        I0  = I0 + ((yl(is)**(1.0D0-gma_TR(tt))) )*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)

                    else
                        aty_tr = lbd_TR(tt)*(yl(is)**(1.0D0-gma_TR(tt)))
                        aty_ss = lbd*(yl(is)**(1.0D0-gma))

                        taxl            = yl(is) - min(aty_tr,aty_ss)
                        RevL_TR(tt)     = RevL_TR(tt) + taxl*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
                        if (aty_tr<=aty_ss) then ! normal case, I actually pay the TR case
                            I0  = I0 + ((yl(is)**(1.0D0-gma_TR(tt))) )*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
                            if (tt<=4) then
                            index_tax_TR(is,ib) = index_tax_TR(is,ib) + 1d0
                            endif
                        else
                            Iss = Iss + ((yl(is)**(1.0D0-gma)) )*hpol_TR(is,ib,tt,ih)*mu_TR(is,ib,tt)
                        endif

                    endif

				endif
			enddo
		enddo
		enddo
		
		!---deficit before labor taxes
		DEFnoL_TR(tt) = G_TR(tt) + (1.0D0+r_TR(tt))*DB_TR(tt) + TF_TR(tt) - RevK_TR(tt)
		
		DB_TR(tt+1)   = DB + (1.0D0-DB_adj)*(DEFnoL_TR(tt)-DEFnoL)

        !---compute deficit
        Rev_TR(tt) = RevK_TR(tt) + RevL_TR(tt)
        GTOT       = G_TR(tt)    + (1.0D0+r_TR(tt))*DB_TR(tt) + TF_TR(tt)
        BC_TR(tt)  = Rev_TR(tt)  + DB_TR(tt+1) - GTOT

        !lbdb_TR(tt) = (wgeH_TR(tt)*Lagg_TR(tt) - GTOT + DB_TR(tt+1) + RevK_TR(tt))/I0

    enddo


    !---asset market clearing given qk
    do tt = 1, T_TR+1
        if (tt == 1) then
            Kagg_TR(tt) = (Aagg_TR(tt) - DB_TR(tt))/qk
        else
            Kagg_TR(tt) = (Aagg_TR(tt) - DB_TR(tt))/qk_TR(tt-1)
        endif
        Kagg_TR(tt) = max(Kagg_TR(tt),1e-4)
    enddo

    Yagg_TR = zbar*(Kagg_TR(1:T_TR)**(1d0-alpha))*(Lagg_TR**alpha)



    !---marginal cost MC_t from Phillips
    do tt = 1, T_TR-1
        MC_TR(tt) = (Pi_TR(tt)-Pibar)*Pi_TR(tt) + ((epsPi-1.0D0)/thetaPi) - (1.0D0/(1.0D0+r_TR(tt+1)))*(Pi_TR(tt+1)-Pibar)*Pi_TR(tt+1)*(Yagg_TR(tt+1)/Yagg_TR(tt))
    enddo
    MC_TR(T_TR) = (Pi_TR(T_TR)-Pibar)*Pi_TR(T_TR) + ((epsPi-1.0D0)/thetaPi)
    MC_TR = MC_TR*(thetaPi/epsPi)


    ! Compute the new capital and investment and implied qkhat
    do tt = 1, T_TR
        Dk_TR(tt) = Kagg_TR(tt+1) - (1.0D0-dlta)*Kagg_TR(tt)
    enddo
    Iagg_TR = Dk_TR + (phiK/2.0D0)*( ((Dk_TR/Kagg_TR(1:T_TR)) - dlta)**2.0D0 ) * Kagg_TR(1:T_TR)
    qkb_TR  = 1.0D0 + phiK*( (Dk_TR/Kagg_TR(1:T_TR)) - dlta )

    ! Update rkb and wage
    rkb_TR = (1d0-alpha)*( MC_TR*zbar )*((Lagg_TR/Kagg_TR(1:T_TR))**alpha)-dlta*qk_TR

    !---update wge_TR
    wge_TR = alpha*( MC_TR*zbar )*( (Kagg_TR(1:T_TR)/Lagg_TR )**(1.0D0-alpha) )

    !---Compute Piw_TR
!    Piw_TR(2:T_TR) = Pi_TR(2:T_TR)*wge_TR(2:T_TR)/wge_TR(1:T_TR-1)
!    Piw_TR(1) = Pi_TR(1)*wge_TR(1)/wge
!
!    ! Compute wgeHb_TR from the Wage Philipps curve
!    do tt = 1, T_TR-1
!        wgeHb_TR(tt) = (Piw_TR(tt)-Piwbar)*Piw_TR(tt) + ((epswPi-1.0D0)/thetawPi)*wge_TR(tt) - (1.0D0/(1.0D0+r_TR(tt+1)))*(Piw_TR(tt+1)-Piwbar)*Piw_TR(tt+1)*(Lagg_TR(tt+1)/Lagg_TR(tt))
!    enddo
!    wgeHb_TR(T_TR) = (Piw_TR(T_TR)-Piwbar)*Piw_TR(T_TR) + ((epswPi-1.0D0)/thetawPi)*wge_TR(T_TR)
!    wgeHb_TR = wgeHb_TR*(thetawPi/epswPi)

    wgeHb_TR = ((epswPi-1d0)/epswPi)*wge_TR


    ! And finally dividends
    !---compute div and Tdx_TR
    do tt = 1, T_TR
        !---divY_TR
        theta_TR(tt) = (thetaPi/2.0D0)*((Pi_TR(tt)-Pibar)**2.0D0)*Yagg_TR(tt)
        divY_TR(tt)  = (1.0D0-MC_TR(tt))*Yagg_TR(tt) - 1.0D0*theta_TR(tt) - FC

        ! divN_TR
        thetaw_TR(tt) = 0d0 !(thetawPi/2.0d0)*((Piw_TR(tt)-Piwbar)**2.0D0)*Lagg_TR(tt)
        divw_TR(tt)  = (wge_TR(tt)-wgeH_TR(tt))*Lagg_TR(tt) - 1.0D0*thetaw_TR(tt) - FCW

        !---divK_TR
        divK_TR(tt) = qk_TR(tt)*Dk_TR(tt) - Iagg_TR(tt)
    !        divK_TR(tt)  = qk_TR(tt)*PHI_TR(tt)*Kagg_TR(tt) - Iagg_TR(tt) - FCk

        !---divBANK_TR
        divBank_TR(tt) = Aagg_TR(tt+1) + (qk_TR(tt)+rk_TR(tt))*Kagg_TR(tt) + (1.0D0+r_TR(tt))*DB_TR(tt) - qk_TR(tt)*Kagg_TR(tt+1) -  DB_TR(tt+1) - (1.0D0+r_TR(tt))*Aagg_TR(tt)

        !---div_TR
        divb_TR(tt)   = divY_TR(tt) + divK_TR(tt) + divW_TR(tt) + divBank_TR(tt)
        
    enddo


	end subroutine Update_HHandTAXES_TR
        
end module ModuleCOMPUTATIONS_TR
