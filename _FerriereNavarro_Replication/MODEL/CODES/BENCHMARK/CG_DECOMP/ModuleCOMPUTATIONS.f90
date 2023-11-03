module ModuleCOMPUTATIONS
    
    use GLOBALS
	use FUNCTIONS
    use Toolbox
    
contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Stationary Measure       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ComputeMeasure
    
    implicit none
    
    real(8) :: munew(NS,Nbeta), a_np(NS), da(NS), da_bh(NS,Nbeta,Nh), muaux(NS,Nbeta)
    integer :: aind_np(NS), aind_np_bh(NS,Nbeta,Nh), ib, ih, is, ix_np, is_np
    real(8) :: tolmu, wwmu, errmu(1)
    integer :: itmu, maxitmu, showmu
    
    maxitmu = 7500; tolmu = 1e-13; wwmu = 0.0D0
    showmu  = 0    
	
	! print *, 'mu = ', mu	
	! print *, 'hprob = ', hprob
	
	do ib = 1, Nbeta
	do ih = 1, Nh
		a_np    = min( apol(:,ib,ih) , amax)
		aind_np = vsearchprm_FET(a_np,avec,avecpar)		
		da      = (a_np - avec(aind_np))/(avec(aind_np+1) - avec(aind_np))
		
		aind_np_bh(:,ib,ih) = aind_np
		da_bh(:,ib,ih)      = da	
	enddo
	enddo
    
    do itmu = 1, maxitmu
        
        munew = 0.0D0
		
		do ib = 1, Nbeta
			do ih = 1, Nh
				a_np    = min( apol(:,ib,ih) , amax)
				aind_np = aind_np_bh(:,ib,ih)
				da      = da_bh(:,ib,ih)
				
				do is = 1, NS
					do ix_np = 1, Nx
						is_np = aind_np(is) + Na*( ix_np - 1)
						munew(is_np,ib)   = munew(is_np,ib)   + (1.0D0-da(is))*Px(xind(is),ix_np)*hpol(is,ib,ih)*mu(is,ib)
						munew(is_np+1,ib) = munew(is_np+1,ib) +    da(is)     *Px(xind(is),ix_np)*hpol(is,ib,ih)*mu(is,ib)					
					enddo
				
				enddo
			enddo
			!munew(:,ib) = Pbeta(ib)*munew(:,ib)/sum(munew(:,ib))
		enddo
        muaux = munew
        
        munew = matmul(muaux,Pbmat)

        errmu = maxval(abs(mu - munew))
        
        if (errmu(1)<tolmu) then
            print *, 'measure converged itmu = ', itmu
            ! print *, 'sum(mu) = ', sum(mu)            
            exit
        else            
            mu = wwmu*mu + (1.0D0-wwmu)*munew
            mu = mu/sum(mu)
        end if        
        
        if (showmu <=0 ) then
            print *, 'itmu ', itmu, 'errmu = ', errmu
            showmu = 500
        end if
        showmu = showmu-1
        

        
    enddo !end loop itmu
    
    end subroutine ComputeMeasure
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Compute Aggregates       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    subroutine ComputeAggregates
    
    implicit none
    
    real(8) :: GTOT, yl(NS), yk(NS), rba, hh(NS), hg0(NS), EMPv2
    integer :: is, ib, ih
	
	Cagg = 0.0D0
	Aagg = 0.0D0
	Hagg = 0.0D0
	Lagg = 0.0D0
	
	do ib = 1, Nbeta
	do ih = 1, Nh
		Cagg = Cagg + sum(cpol(:,ib,ih)*hpol(:,ib,ih)*mu(:,ib))
		Aagg = Aagg + sum(apol(:,ib,ih)*hpol(:,ib,ih)*mu(:,ib))
		Hagg = Hagg + sum(hpol(:,ib,ih)*hpol(:,ib,ih)*mu(:,ib))
		Lagg = Lagg + sum(hvec(ih)*S(:,2)*hpol(:,ib,ih)*mu(:,ib))		
	enddo
	enddo
	
	Kagg = (zbar*MC*(1.0D0-alpha))/(rk+(qk*dlta))
	Kagg = (Kagg**(1.0D0/alpha)) * Lagg
	Kdmn = ((1.0D0-alpha)/alpha) * (wge/(rk+(qk*dlta))) * Lagg
	Kagg = max(Kagg,1e-4)
	
	Yagg = zbar*(Kagg**(1.0D0-alpha)) * (Lagg**alpha)

	KYxx = zbar*(MC*(1.0D0-alpha))/(rk+(qk*dlta))
	
	!!! employment + average hours
	EMP   = 0.0D0 ! 	Hmean = 0.0D0
	do ib = 1, Nbeta
		EMP = EMP + sum(hpol(:,ib,2)*mu(:,ib))		
	enddo 
		
    
    ! FCw = (wge - wgeH*(1.0D0-swH))*Lagg
	! FC = (zbar - wgeH)*Lagg
	FC = Yagg*(1.0D0-MC)
    FCW = (wge-wgeH)*Lagg
	
	!---update lbd
	RevK = 0.0D0;  RevL = 0.0D0; INTEGRAL = 0.0D0;
	do ib = 1, Nbeta
	do ih = 1, Nh
		do is = 1, NS
			yl(is)  = wgeH*hvec(ih)*S(is,2)
			rba     = RbFUN(S(is,1))	
			yk(is)  = rba*S(is,1)
			
			RevK = RevK + tauk*rba*S(is,1)*hpol(is,ib,ih)*mu(is,ib)
			
			if (yl(is) > 0.0D0 ) then
				RevL     = RevL     + ( yl(is) - lbd*(yl(is)**(1.0D0-gma)) ) *hpol(is,ib,ih)*mu(is,ib)
				INTEGRAL = INTEGRAL + ( yl(is)**(1.0D0-gma) )*hpol(is,ib,ih)*mu(is,ib)
			endif
		enddo
	enddo
	enddo
	
	Rev  = RevK + RevL
    GTOT = G + rk*DB + TF 
    BC   = Rev - GTOT 
    
    lbdnew = RevK + wgeH*Lagg - GTOT
    lbdnew = lbdnew/INTEGRAL    
    
	DEFnoL = G + (1.0D0+rk)*DB + TF - RevK

    
    end subroutine ComputeAggregates
	

	
    
    
    
end module ModuleCOMPUTATIONS
