module ModuleCOMPUTATIONS_lpe

	use GLOBALS
	use FUNCTIONS_lpe
	use Toolbox
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       					Compute MEASURE       				!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ComputeMEASURE_TR
    
    implicit none
    real(8) :: mm(NS,Nbeta), a_np(NS), da(NS), muaux(NS,Nbeta)
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
        enddo	! end loop beta

        mu_TR(:,:,tt+1) =  matmul(mm,Pbmat)

    enddo
    
    end subroutine ComputeMEASURE_TR
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!  					Update HH side    				  			!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine COMPUTE_HH_TR
	
	implicit none
	real(8) :: yl(NS), yk(NS), rba, taxl, taxlss, GTOT, deltax
	integer :: tt, ib, ih, is
	
	!---update HH side	
	Cagg_TR = 0.0D0
    Aagg_TR = 0.0D0    
    Lagg_TR = 0.0D0	
	Hagg_TR = 0.0D0
	! Anp     = 0.0D0
	
	Aagg_TR(1) = Aagg
	
	do tt = 1, T_TR
		do ib = 1, Nbeta
		do ih = 1, Nh
			Cagg_TR(tt) = Cagg_TR(tt) + sum(cpol_TR(:,ib,tt,ih)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt))
			Lagg_TR(tt) = Lagg_TR(tt) + sum(hvec(ih)*S(:,2)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt))
			Hagg_TR(tt) = Hagg_TR(tt) + sum(hvec(ih)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt))

			Aagg_TR(tt+1) = Aagg_TR(tt+1) + sum(apol_TR(:,ib,tt,ih)*hpol_TR(:,ib,tt,ih)*mu_TR(:,ib,tt)) 			
		enddo
		enddo
	enddo
	
	
	end subroutine COMPUTE_HH_TR
	
	!------------------------------------------------------------------------------
	!------------------------------------------------------------------------------
	!-----lpe x (s,beta)
	subroutine COMPUTE_lpex
	
	implicit none
	real(8) :: lpe(NS,Nbeta), yls(NS,Nbeta)
	integer :: tt, ih, ib, is
	
	
	tt=1; ih=2;
	
	do ib = 1, Nbeta
	do is = 1, NS
		lpe(is,ib) = ((hpol_TR(is,ib,tt,ih) - hpol(is,ib,ih))/hpol(is,ib,ih))/(100d0*DTAX_TR(1))
		yls(is,ib) = wgeH*hvec(ih)*hpol(is,ib,ih)*S(is,2)
	enddo
	enddo
	
	open(1,  file = 'OUTPUT/lpe.txt',   status = 'unknown')    	
	open(2,  file = 'OUTPUT/yls.txt',   status = 'unknown')    	
    do is = 1, NS
        write(1, '(*(f18.6))') ( lpe(is,ib) , ib = 1,Nbeta )        
		write(2, '(*(f18.6))') ( yls(is,ib) , ib = 1,Nbeta )        
	enddo		
	close(1); close(2); 

	
	
	end subroutine COMPUTE_lpex
	
	
end module ModuleCOMPUTATIONS_lpe
