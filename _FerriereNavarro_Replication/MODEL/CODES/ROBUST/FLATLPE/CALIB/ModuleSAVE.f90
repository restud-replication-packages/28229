4module ModuleSAVE
    
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
    open(4,  file = 'OUTPUT/hpolw.txt',   status = 'unknown')
	do is = 1, NS
		write(1, '(*(f18.6))') ( V(is,ib),  ib = 1, Nbeta)
		write(2, '(*(f18.6))') ( Ve(is,ib), ib = 1, Nbeta)		
		write(3, '(*(f18.6))') ( mu(is,ib), ib = 1, Nbeta)		
        write(4, '(*(f18.6))') ( hpol(is,ib,2), ib = 1, Nbeta)
	enddo
	close(1); close(2); close(3); close(4)

	
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
	
	
	! open(1, file = 'OUTPUT/cpol.txt',   status = 'unknown')
	! open(2, file = 'OUTPUT/apol.txt',   status = 'unknown')
	! open(3, file = 'OUTPUT/kpol.txt',   status = 'unknown')
	! open(4, file = 'OUTPUT/Vwp.txt',    status = 'unknown')
	! open(5, file = 'OUTPUT/pprob.txt',  status = 'unknown')
	! do ih = 1, Nh	
	! do is = 1, NS
	! 	write(1, '(*(f18.6))') ( cpol(is,ih,ip) )
	! 	write(2, '(*(f18.6))') ( bpol(is,ih,ip) )
	! 	write(3, '(*(f18.6))') ( kpol(is,ih,ip) )
	! 	write(4, '(*(f18.6))') ( Vwp(is,ih,ip)  )
	 !	write(5, '(*(f18.6))') ( pprob(is,ih,ip)  )		
	! enddo
	! enddo
	! enddo
	! close(1); close(2); close(3); close(4); close(5)
	
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
	real(8) :: CALIB(2), guess(2), bbounds(2), prices(3), taxes(3)
	integer :: ii
	
	!--guess: rk and lbd
	guess = [rk, lbd]
	open(1,  file = 'OUTPUT/guess.txt',      status = 'unknown')
    write(1, *) (guess(ii), ii = 1, 2)
    close(1)
	
	!---prices	
	prices = [wge, rk, qk]
	open(1,  file = 'OUTPUT/prices.txt',      status = 'unknown')
    write(1, '(*(f18.6))') (prices(ii), ii = 1, 4)
    close(1)
    
    taxes = [lbd,gma,TF]
    open(1, file='OUTPUT/tax.txt', status = 'unknown')
    write(1, '(*(f18.6))') (taxes(ii), ii = 1, 3)
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

	
	end subroutine SaveAggregates
    
end module ModuleSAVE
