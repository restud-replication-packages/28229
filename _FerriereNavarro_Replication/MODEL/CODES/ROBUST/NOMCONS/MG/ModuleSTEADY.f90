module ModuleSTEADY

	use GLOBALS
    use FUNCTIONS
    use ModuleCOMPUTATIONS
!	use ModuleMPC
    use ModuleSAVE
    use Toolbox
!    use ModuleSIMULATION_PANEL
    
contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!          GE Loop           !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Model_STEADYSTATE
    
    implicit none
	real(8) :: tolGE, errGE, tolBC, wwb, wwlbd, guess(2)
    integer :: itb, maxitb, itlbd, maxitlbd, itsave, adjparm, ib, ia
    real(8) :: prices(3), bbounds(2), Aggvec(4)
    
    adjparm = 0      ! =1 to adjust paramteres across iterations
    
    tolGE = 1e-5; maxitb   = 1500; wwb  = 0.75D0
    tolBC = 1e-7; maxitlbd = 25;  wwlbd = 0.50D0
    itsave = -1
	
	!---beta vec
	bub = 0.9995D0*1.0D0/(1.0D0+rk*(1.0D0-tauk)) ! = 0.9918
    blb = 0.90D0   

	
	open(1,  file = 'OUTPUT/bbounds.txt',      status = 'unknown')
    read(1, *) bbounds
    close(1)
    blb = bbounds(1); bub = bbounds(2);
	
	do ib = 1, Nbeta
		bvec(ib) = btaH - dble(Nbeta-ib)*Dbta
	enddo
	
	!---initial prices	
	MC	= (epsPi-1.0D0)/epsPi
	qk  = 1.0D0
	
	lbd = 0.683D0
	
	! guess = [rk, lbd]
	open(1,  file = 'OUTPUT/guess.txt',      status = 'unknown')
    read(1, *) guess
    close(1)
	lbd = guess(2)
	
	wge = alpha*( (MC*zbar)**(1.0D0/alpha) )*( ((1.0D0-alpha)/(rk+(qk*dlta)))**((1.0D0-alpha)/alpha) )

    wgeH = ((epswPi-1d0)/epswPi)*wge
	
	Vnow = 0.0D0; Vwrk = 0.0D0
	
	!---mu
    ! mu = 1.0D0
    ! do ib = 1, Nbeta
    !     mu(:,ib) = Pbeta(ib)*mu(:,ib)/sum(mu(:,ib))
    ! enddo 	
	open(1, file = 'OUTPUT/muLOAD.txt',    status = 'unknown')    
	read(1, *) mu
    close(1)
	
	
	!***********************************************************************************
	!---LOOP
	!***********************************************************************************
	do itb = 1, maxitb
	
		btaH = 0.5D0*(blb + bub)
		! btaH = 0.993700D0
		Dbta = 0.0575D0/dble(Nbeta-1);
		
		do ib = 1, Nbeta
			bvec(ib) = btaH - dble(Nbeta-ib)*Dbta
		enddo
	
		print *, '*****--------------------------*****'
		print *, '***** itb = ', itb, ', btaH = ', btaH,' *****'
		print *, '*****--------------------------*****'	
		
	
		do itlbd = 1, maxitlbd
            
			call Model_VF
			call ComputeMeasure
			call ComputeAggregates
            
            print *, '** itlbd = ', itlbd, '  **'
            print *, '** lbd = ', lbd, ', lbdnew = ', lbdnew            
            print *, '** BC = ', BC, '  **'
            
            if (abs(BC) < tolBC ) then
                print *, 'lambda found '
                exit
            else
                lbd = wwlbd*lbd + (1.0D0-wwlbd)*lbdnew
            endif
			!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! adjust parameters !!!
            print *, '*********************************'
            print *, 'EMP   = ', EMP			
            print *, 'DB/Yagg = ', DB/Yagg
            print *, 'G/Yagg = ', G/Yagg
            print *, 'TF/Yagg = ', TF/Yagg
            print *, '*********************************'
            
            if (adjparm == 1) then
                G  = GYcal*Yagg  ! G  = 0.15D0*Yagg
                TF = TFcal*Yagg  ! TF = 0.05D0*Yagg
				DB = DBcal*Yagg  ! DB = 3.60D0*Yagg ! 2.40D0*Yagg
                
                if (EMP > EMPcal + EMPtol)     then ! if (EMP > 0.82D0 ) then                    
					B0 = 1.05D0*B0
                elseif (EMP < EMPcal - EMPtol) then  ! elseif (EMP < 0.78D0)  then                    
					B0 = 0.95D0*B0
                endif                
                
            endif            
            !!!!!!!!!!!!!!!!!!!!!!!!!
		enddo !end itlbd	
		
		
		call SaveVFandPolicies
		call SaveAggregates 
		
		errGE = Aagg - (DB+Kagg) ! r - rnew 
        print *, '*********************************'
        print *, 'itb = ', itb, ', errGE = ', errGE
        print *, 'Yagg = ', Yagg
        print *, 'EMP  = ', EMP
		print *, 'Kagg/Yagg = ', Kagg/Yagg
		print *, 'KYxx = ', KYxx
        print *, 'DB/Yagg = ', DB/Yagg
        print *, 'G/Yagg = ', G/Yagg
        print *, 'TF/Yagg = ', TF/Yagg
        print *, '*********************************'
		
		!---update btaH bounds
		if (abs(errGE) < tolGE) then
			print *, 'bta found' ! print *, 'bta and rk found'
			call SaveVFandPolicies
			call SaveAggregates 
            exit
		else
			if (errGE > 0.0D0) then  !too much saving, decrease btaH
				bub = wwb*bub + (1.0D0-wwb)*btaH				
			else
				blb = wwb*blb + (1.0D0-wwb)*btaH
			endif	
			
		endif	
	enddo !end itb
	

	
	
	end subroutine Model_STEADYSTATE
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Value Function       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Model_VF
    
    implicit none 
    
    real(8) :: Vnownew(NS,Nbeta), Vwrknew(NS,Nbeta), xl(NS), xr(NS), yl(NS), yk(NS), rba(NS), TAX(NS)    
    real(8) :: errVvec(NS,Nbeta), errVnow(1), errVwrk(1), errV(1), wwV, tolV
    integer :: ib, is, ih, iloc(NS), itV, maxitV, itH, maxitH, showV
    integer :: clck_counts_beg, clck_counts_end, clck_rate
    
    maxitV = 2000; wwV = 0.25D0; tolV = 1e-9; maxitH = 40
    showV  = 10
	
	do itv = 1, maxitV
        ! call system_clock ( clck_counts_beg, clck_rate )
		
		!---compute Ve        
        call ComputeVe
        
		!---update V
		do ib = 1, Nbeta
			!$OMP PARALLEL DO 
            do is = 1, NS    
			
				!---case 1: no work -> h = 0
				yl(is)  = wgeH*0.0D0*S(is,2)
				rba(is) = RbFUN(S(is,1))	
				yk(is)  = rba(is)*S(is,1)
				TAX(is) = TAXFUN(yl(is),yk(is))
				xr(is)	= min(amax,yl(is) + yk(is) + qk*S(is,1) + TF - TAX(is))
				xl(is)  = amin
				
				apol(is,ib,1)  = SolveBella(xl(is),xr(is),0.0D0,is,ib)
				Vnownew(is,ib) = VFa_eval(apol(is,ib,1),0.0D0,is,ib,1)										
			
				!---case 2: work -> h = hbar
				yl(is)  = wgeH*hbar*S(is,2)
				rba(is) = RbFUN(S(is,1))	
				yk(is)  = rba(is)*S(is,1)
				TAX(is) = TAXFUN(yl(is),yk(is))
				xr(is)	= min(amax,yl(is) + yk(is) + qk*S(is,1) + TF - TAX(is))
				xl(is)  = amin
				
				apol(is,ib,2)  = SolveBella(xl(is),xr(is),hbar,is,ib)
				Vwrknew(is,ib) = VFa_eval(apol(is,ib,2),hbar,is,ib,1)
				
			enddo !end loop is
			!$OMP END PARALLEL DO 
		enddo  !end loop ib
		
		!---update V
        ! errVvec = abs(V-Vnew)
		errVwrk  = maxval(abs(Vwrk-Vwrknew))
		errVnow  = maxval(abs(Vnow-Vnownew))
        
		errV     = max(errVnow, errVwrk)
        
        ! print *, '***itV = ', itV 
		! print *, '***errVwrk, errVnow = ', errVwrk, errVnow
        if (showV <= 0) then 
            print *, 'itV = ', itV
			print *, 'errVwrk, errVnow = ', errVwrk, errVnow
            showV = 50 + 1
        end if
        showV = showV -1 
        
        
        if (errV(1) < tolV ) then
            print *, 'Value functions converged at itV = ', itV     
            exit
        else
			Vnow = wwV*Vnow + (1.0D0-wwV)*Vnownew
            Vwrk = wwV*Vwrk + (1.0D0-wwV)*Vwrknew
        endif

        do ib = 1, Nbeta
        do ih = 1, Nh
            aind_np_ho(:,ib,ih) = vsearchprm_FET(apol(:,ib,ih),avec,avecpar)
        enddo
        enddo
		
		do itH = 1, maxitH			
			call ComputeVe
            do ib = 1, Nbeta
				Vnownew(:,ib) = VFvec_itH(0.0D0,ib,aind_np_ho(:,ib,1))
                Vwrknew(:,ib) = VFvec_itH(hbar,ib,aind_np_ho(:,ib,2))
            enddo

            errVwrk  = maxval(abs(Vwrk-Vwrknew))
            errVnow  = maxval(abs(Vnow-Vnownew))
            errV     = max(errVnow, errVwrk)

            if (errV(1) < tolV ) then
                ! print *, 'Value functions converged at itV = ', itV
                exit
            else
                Vnow = Vnownew
                Vwrk = Vwrknew
            endif
			
        enddo
		
	enddo  !end loop itv	
	
	call SaveVFandPolicies
	
	end subroutine Model_VF
	
end module ModuleSTEADY
