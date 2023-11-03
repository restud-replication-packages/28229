module ModuleSAVE
    
    use GLOBALS
    
contains    
    
    !--------------------------------------------------
    !---save grids and parameters
    !--------------------------------------------------
    subroutine SaveGRIDS
    
    implicit none
    real(8) :: paramN(3), prmsv(3)
    integer :: is, ia, ix, ix_np, ii
    
    open(1,  file = 'OUTPUT/avec.txt',    status = 'unknown')    
    open(2,  file = 'OUTPUT/xvec.txt',    status = 'unknown')            
    open(3,  file = 'OUTPUT/params.txt',   status = 'unknown')    
    do ia = 1, Na;     write(1, *) avec(ia); enddo    
    do ix = 1, Nx;     write(2, *) xvec(ix); enddo        
    write(3, '(*(f18.6))') [dble(Na), dble(Nx), dble(Nbeta)]
    close(1); close(2); close(3); close(4)
	
	open(1,  file = 'OUTPUT/Svec.txt',    status = 'unknown')
    do is = 1, NS
		write(1, *) (S(is,ii), ii = 1, 2) 	
	enddo
	close(1)
	
	
	open(1,  file = 'OUTPUT/Px.txt',   status = 'unknown')
	open(2,  file = 'OUTPUT/pix.txt',  status = 'unknown')
	do ix = 1, Nx
		write(1, '(*(f18.6))') ( Px(ix,ix_np), ix_np = 1, Nx )
		write(2, '(*(f18.6))') ( pix(ix) )		
	enddo
	close(1); close(2)
	
	
	prmsv = [alpha, dlta, rhow]
	open(1,  file = 'OUTPUT/prmsv.txt',   status = 'unknown')
	write(1, '(*(f18.6))') (prmsv(ii), ii = 1, 3)
	close(1)
    end subroutine SaveGRIDS
	
	!--------------------------------------------------
    !---save policies
    !--------------------------------------------------
    subroutine SaveVFandPolicies
	
	implicit none
	integer :: is, ib, ih, ip
	
	open(1, file = 'OUTPUT/V.txt',    status = 'unknown')
	open(2, file = 'OUTPUT/Ve.txt',    status = 'unknown')	
	open(3, file = 'OUTPUT/mu.txt',    status = 'unknown')	
	do is = 1, NS
		write(1, '(*(f18.6))') ( V(is,ib),  ib = 1, Nbeta)
		write(2, '(*(f18.6))') ( Ve(is,ib), ib = 1, Nbeta)		
		write(3, '(*(f18.6))') ( mu(is,ib), ib = 1, Nbeta)		
	enddo	
	close(1); close(2); close(3)	

	
	open(1,  file = 'OUTPUT/apol.txt',   status = 'unknown')
    open(2,  file = 'OUTPUT/cpol.txt',   status = 'unknown')
    open(3,  file = 'OUTPUT/hpol.txt',   status = 'unknown')
	do ib = 1, Nbeta
    do is = 1, NS
        write(1, '(*(f18.6))') ( apol(is,ib,ih) , ih = 1,2 )
        write(2, '(*(f18.6))') ( cpol(is,ib,ih) , ih = 1,2 )
        write(3, '(*(f18.6))') ( hpol(is,ib,ih) , ih = 1,2 )		
	enddo	
	enddo
	close(1); close(2); close(3); 
	

	
	! for the initial guess
    open(1,  file = 'OUTPUT/VF.txt',      status = 'unknown')
    write(1, *) V
    close(1)  
	end subroutine SaveVFandPolicies
	
	!--------------------------------------------------
    !---save aggregates
    !--------------------------------------------------
    subroutine SaveAggregates 
	
	implicit none
	real(8) :: CALIB(2), guess(2), bbounds(2), prices(3), AGG(9)
	integer :: ii
	
	!--guess: rk and lbd
	guess = [rk, lbd]
	open(1,  file = 'OUTPUT/guess.txt',      status = 'unknown')
    write(1, *) (guess(ii), ii = 1, 2)
    close(1)
	
	!---prices	
	prices = [wge, rk, qk]
	open(1,  file = 'OUTPUT/prices.txt',      status = 'unknown')
    write(1, '(*(f18.6))') (prices(ii), ii = 1, 3)
    close(1)
	
	! bta bounds
    bbounds = [blb, bub]
    open(1,  file = 'OUTPUT/bbounds.txt',      status = 'unknown')
    write(1, *) ( bbounds(ii), ii = 1, 2)
    close(1)
	
	!---calibrated parameters
    CALIB = [Yagg, B0]
    open(1,  file = 'OUTPUT/CALIB_YandB_save.txt',      status = 'unknown')
    write(1, *) (CALIB(ii), ii = 1, 2)
    close(1)
	
	open(1, file = 'OUTPUT/VeLOAD.txt',    status = 'unknown')
	open(2, file = 'OUTPUT/muLOAD.txt',    status = 'unknown')
    write(1, *) Ve
	write(2, *) mu
    close(1); close(2)
	
	AGG = [lbd, G, Cagg, Yagg, Lagg, EMP, Kagg, DB, FC]
	open(1,  file = 'OUTPUT/AGGREGATES.txt',      status = 'unknown')
    write(1, '(*(f18.6))') (AGG(ii), ii = 1, 9)
    close(1)
	

	
	end subroutine SaveAggregates
	
	!--------------------------------------------------
    !---Save Transition Aggregates     
    !--------------------------------------------------
    subroutine SaveAggregates_TR
    
    implicit none
    
    integer :: is, tt, ttax
    
    open(1,  file = 'OUTPUT/G_TR.txt',     status = 'unknown')
    open(2,  file = 'OUTPUT/lbd_TR.txt',   status = 'unknown')
    open(3,  file = 'OUTPUT/gma_TR.txt',   status = 'unknown')
    open(4,  file = 'OUTPUT/wge_TR.txt',   status = 'unknown')    
    open(5,  file = 'OUTPUT/rk_TR.txt',    status = 'unknown')	
	open(6,  file = 'OUTPUT/qk_TR.txt',    status = 'unknown')
    open(7,  file = 'OUTPUT/tauk_TR.txt',  status = 'unknown')
    open(8,  file = 'OUTPUT/Yagg_TR.txt',  status = 'unknown')    
    open(9, file = 'OUTPUT/Cagg_TR.txt',   status = 'unknown')
    open(10, file = 'OUTPUT/Lagg_TR.txt',  status = 'unknown')
	open(11, file = 'OUTPUT/Iagg_TR.txt',  status = 'unknown')	
	open(12, file = 'OUTPUT/Dk_TR.txt',    status = 'unknown')	
    open(13, file = 'OUTPUT/div_TR.txt',   status = 'unknown')
    open(14, file = 'OUTPUT/wgeH_TR.txt',  status = 'unknown')
    do tt = 1, T_TR
         write(1 , '(*(f18.6))') ( G_TR(tt)   )
         write(2 , '(*(f18.6))') ( lbd_TR(tt) )
         write(3 , '(*(f18.6))') ( gma_TR(tt) )
         write(4 , '(*(f18.6))') ( wge_TR(tt) )
         write(5 , '(*(f18.6))') ( rk_TR(tt) )         
		 write(6 , '(*(f18.6))') ( qk_TR(tt) )
         write(7 , '(*(f18.6))') ( tauk_TR(tt) )
         write(8 , '(*(f18.6))') ( Yagg_TR(tt) )         
         write(9 , '(*(f18.6))') ( Cagg_TR(tt) )
         write(10,'(*(f18.6))')  ( Lagg_TR(tt) )
		 write(11,'(*(f18.6))')  ( Iagg_TR(tt) )		 
		 write(12,'(*(f18.6))')  ( Dk_TR(tt) )                  
         write(13,'(*(f18.6))')  ( div_TR(tt) )                  
         write(14,'(*(f18.6))')  ( wgeH_TR(tt) )
    enddo
    close(1); close(2); close(3) ; close(4) ; close(5) ; close(6) ; close(7); 
	close(8); close(9); close(10); close(11); close(12); close(13); close(14)
	
	open(1,  file = 'OUTPUT/Pi_TR.txt',   status = 'unknown')    
    open(2,  file = 'OUTPUT/i_TR.txt',    status = 'unknown')    	
    do tt = 1, T_TR
         write(1, '(*(f18.6))') ( Pi_TR(tt)   )         
         write(2, '(*(f18.6))') ( i_TR(tt) )                
    end do
    close(1); close(2);
    
	open(1,  file = 'OUTPUT/DB_TR.txt',   status = 'unknown')
	open(2,  file = 'OUTPUT/Kd_TR.txt',   status = 'unknown')
	open(3,  file = 'OUTPUT/DB_TR.txt',   status = 'unknown')	
	open(4,  file = 'OUTPUT/Aagg_TR.txt', status = 'unknown')		
	do tt = 1, T_TR+1
         write(1, '(*(f18.6))') ( DB_TR(tt)   )
		 write(2, '(*(f18.6))') ( Kagg_TR(tt)   )
		 write(3, '(*(f18.6))') ( DB_TR(tt)   )		 
		 write(4, '(*(f18.6))') ( Aagg_TR(tt)   )		 
	enddo	 
	close(1); close(2); close(3); close(4);
    
	open(1,  file = 'OUTPUT/Rev_TR.txt',    status = 'unknown')
    open(2,  file = 'OUTPUT/RevL_TR.txt',   status = 'unknown')
    open(3,  file = 'OUTPUT/RevK_TR.txt',   status = 'unknown')    
	open(4,  file = 'OUTPUT/TF_TR.txt',     status = 'unknown')    	
	open(5,  file = 'OUTPUT/BC_TR.txt',     status = 'unknown')	
    do tt = 1, T_TR
         write(1, '(*(f18.6))') ( Rev_TR(tt)  )
         write(2, '(*(f18.6))') ( RevL_TR(tt) )
         write(3, '(*(f18.6))') ( RevK_TR(tt) )         
		 write(4, '(*(f18.6))') ( TF_TR(tt)   )         		 
		 write(5, '(*(f18.6))') ( BC_TR(tt)   )
    enddo
    close(1); close(2); close(3); close(4);  close(5); 
		
	
	open(1,  file = 'OUTPUT/Tdx_TR_LOAD.txt',      status = 'unknown')
	write(1, *) Tdx_TR
	close(1)	
    
    end subroutine SaveAggregates_TR

    !--------------------------------------------------
    !---save policies for transitin
    !--------------------------------------------------
    subroutine SaveVFandPolicies_TR

    implicit none
    integer :: tt, is, ib, ih

    open(1, file = 'OUTPUT/Vwrk_TR.txt',    status = 'unknown')
    open(2, file = 'OUTPUT/Vnow_TR.txt',    status = 'unknown')
    open(3, file = 'OUTPUT/mu_TR.txt',    status = 'unknown')
    do tt = 1, T_TR
    do ib = 1, Nbeta
    do is = 1, NS
        write(1, '(*(f18.6))') ( Vwrk_TR(is,ib,tt) )
        write(2, '(*(f18.6))') ( Vwrk_TR(is,ib,tt) )
        write(3, '(*(f18.6))') ( mu_TR(is,ib,tt)   )
    enddo
    enddo
    enddo
    close(1); close(2); close(3)

    open(1,  file = 'OUTPUT/apol_TR.txt',   status = 'unknown')
    open(2,  file = 'OUTPUT/cpol_TR.txt',   status = 'unknown')
    open(3,  file = 'OUTPUT/hpol_TR.txt',   status = 'unknown')
    do ih = 1, Nh
    do tt = 1, T_TR
    do ib = 1, Nbeta
    do is = 1, NS
        write(1, '(*(f18.6))') ( apol_TR(is,ib,tt,ih)  )
        write(2, '(*(f18.6))') ( cpol_TR(is,ib,tt,ih)  )
        write(3, '(*(f18.6))') ( hpol_TR(is,ib,tt,ih)  )
    enddo
    enddo
    enddo
    enddo
    close(1); close(2); close(3);

    end subroutine SaveVFandPolicies_TR

    
end module ModuleSAVE
