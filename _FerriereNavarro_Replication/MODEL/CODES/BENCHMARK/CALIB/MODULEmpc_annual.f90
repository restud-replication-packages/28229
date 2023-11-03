module MODULEmpc_annual

	use GLOBALS
	use FUNCTIONS	
	use Toolbox
	
contains

	!-------------------------------------------	
	!---mpc annual computation
	!-------------------------------------------
	subroutine COMPUTEmpc_annual
	
	implicit none
	integer, parameter :: Tsim = 4
	
	real(8) :: ymn, ymnDATA, cvMtoD, yl(NS), yk(NS), rba(NS), dy, ady(NS), a_np, da
	real(8) :: cpolt1(NS,Nbeta,Nh), apolt1(NS,Nbeta,Nh), hpolt1(NS,Nbeta,Nh), Vnowdy(NS,Nbeta), Vwrkdy(NS,Nbeta), Dnow(NS), Dwrk(NS), corr(NS), cslope(NS,Nbeta,Nh)
	real(8) :: GMAt(NS,Nbeta), GMAtp1(NS,Nbeta), GGx(NS,Nbeta), cseh(NS,Nbeta), cseht1(NS,Nbeta), csimax
	real(8) :: csim(NS,Nbeta), csimdy(NS,Nbeta), mpc_annual(NS,Nbeta), caggsim(Tsim), caggsimdy(Tsim)
	integer :: ib, ih, is, tt, iz, iz0, ia0, ix0, ib0, is0, izt, iat, ixt, ibt, ist, ib_np, ix_np, is_np, iz_np, aind_dy(NS), Sind_dy(NS), aidy, Sidy, aindv(1), aind_np, itshow, dois		
	
	!---compute ymean
	ymn = 0.0D0
	do ib = 1, Nbeta
	do ih = 1, Nh
		do is = 1, NS
			yl(is)  = wgeH*hvec(ih)*S(is,2)
			rba(is) = RbFUN(S(is,1))	
			yk(is)  = rba(is)*S(is,1)
			
			ymn = ymn + (yl(is) + yk(is))*hpol(is,ib,ih)*mu(is,ib)
		enddo
	enddo
	enddo
	
	ymnDATA = 67000.0D0
	cvMtoD  = ymnDATA/ymn
	
	!---compute transfer dC and new policies
	dy = 500.0D0/cvMtoD			! mpc out of dC dollars
	
	!---compute c and a at a+dy <--> policy at t=1 -- interpolate policies	
	ady     = S(:,1) + dy
	aind_dy = vsearchprm_FET(ady,avec,avecpar)
	Sind_dy = aind_dy + Na*(xind-1)   
	
	do ib = 1, Nbeta
	do ih = 1, Nh
		do is = 1, NS
			aidy = aind_dy(is)
			Sidy = Sind_dy(is)
			
			cslope(is,ib,ih) = (ady(is) - avec(aidy)) * ( cpol(Sidy+1,ib,ih) - cpol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))			
		
			cpolt1(is,ib,ih) = cpol(Sidy,ib,ih) + (ady(is) - avec(aidy)) * ( cpol(Sidy+1,ib,ih) - cpol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))			
			apolt1(is,ib,ih) = apol(Sidy,ib,ih) + (ady(is) - avec(aidy)) * ( apol(Sidy+1,ib,ih) - apol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))			
		enddo
	enddo
	enddo
	
	!---compute values and hpol at a+dy <--> interpolate values
	do ib = 1, Nbeta
		do is = 1, NS
			aidy = aind_dy(is)
			Sidy = Sind_dy(is)
		
			!---value of not working
			Vnowdy(is,ib) = Vnow(Sidy,ib) + (ady(is) - avec(aidy)) * ( Vnow(Sidy+1,ib) - Vnow(Sidy,ib) )/(avec(aidy+1) - avec(aidy))			
			
			!---value of working
			Vwrkdy(is,ib) = Vwrk(Sidy,ib) + (ady(is) - avec(aidy)) * ( Vwrk(Sidy+1,ib) - Vwrk(Sidy,ib) )/(avec(aidy+1) - avec(aidy))						
		enddo
		
		Dnow = Vnowdy(:,ib);  Dwrk = Vwrkdy(:,ib);  
		corr = max( (Dnow/rhow) , (Dwrk/rhow) )
		
		hpolt1(:,ib,1) = exp((Dnow/rhow)-corr)/(exp((Dnow/rhow)-corr)+exp((Dwrk/rhow)-corr))  ! prob of h=0
		hpolt1(:,ib,2) = 1.0D0 - hpolt1(:,ib,1) 					
	enddo	
	
	! do ib = 1, Nbeta
	! do is = 1, NS
	! 	print *, 'hpolt1 = ', hpolt1(is,ib,1), hpolt1(is,ib,2)
	! enddo
	! enddo
	! STOP
	
	
	!----------------------------------------------------------
	!---compute consumption for t=1,..,4 with no shock	
	cseh = 0.0D0
	do ib = 1, Nbeta
	do is = 1, NS
		! iz = is + NS*(ib-1)
		cseh(is,ib) = ( hpol(is,ib,1)*cpol(is,ib,1) + hpol(is,ib,2)*cpol(is,ib,2) )
	enddo
	enddo			
	
	print *, '***********************'	
	print *, '*** no shock'
	
	csim = 0.0D0; caggsim = 0.0D0;	itshow = 0;	dois   = 1;
	if (dois == 1) then
	do ib0 = 1, Nbeta	
	do is0 = 1, NS
	
		GMAt          = 0.0D0;	
		GMAt(is0,ib0) = 1.0D0; 
		
		! tt = 1; 		
		! csimax        = cseh(is0,ib0)*GMAt(is0,ib0)
		! csim(is0,ib0) = csim(is0,ib0) + csimax
		! caggsim(tt) = caggsim(tt) + csimax*mu(is0,ib0)
		
		do tt = 1, Tsim
			GMAtp1 = 0.0D0
			
			do ibt = 1, Nbeta
			do ist = 1, NS			
				ixt = xind(ist)
				
				csim(is0,ib0) = csim(is0,ib0) + (cseh(ist,ibt)*GMAt(ist,ibt))
				caggsim(tt)   = caggsim(tt)   + (cseh(ist,ibt)*GMAt(ist,ibt))*mu(is0,ib0)				
				
				do ih = 1, Nh
					a_np    = min( apol(ist,ibt,ih) , amax)
					aindv   = vsearchprm_FET([a_np],avec,avecpar)		
					aind_np = aindv(1)
					da      = (a_np - avec(aind_np))/(avec(aind_np+1) - avec(aind_np))
					
					do ib_np = 1, Nbeta
					do ix_np = 1, Nx
						is_np = aind_np + Na*(ix_np - 1)
			
						GMAtp1(is_np  ,ib_np) = GMAtp1(is_np  ,ib_np) + (1.0D0-da)*Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)
						GMAtp1(is_np+1,ib_np) = GMAtp1(is_np+1,ib_np) +    da     *Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)		
					enddo
					enddo
				
				enddo !end loop ih 
			enddo !end loop ist
			enddo !end loop ibt
			
			! print *, 'sum( GMAtp1 ) = ', sum( GMAtp1 )
			GMAtp1 = GMAtp1/sum(GMAtp1)			
			GMAt   = GMAtp1
		enddo !ned loop tt
		
		if (itshow<=0) then
			print *, 'is0 = ', is0,'/',NS,', ib0 = ', ib0,'/',Nbeta
			itshow = 1500
		endif
		itshow = itshow-1
	
	enddo !end loop is0	
	enddo !end loop ib0	
	endif
	
	
	
!	print *, 'Cagg = ', Cagg
!	print *, 'caggsim = ', caggsim
	! print *, 'csim(1:10,:) = ', csim(1:10,:)
	!----------------------------------------------------------
	!----------------------------------------------------------
	
	! STOP
	!----------------------------------------------------------
	!---compute consumption for t=1,..,4 with shock
	
	cseht1 = 0.0D0
	do ib = 1, Nbeta
	do is = 1, NS
		! iz = is + NS*(ib-1)
		cseht1(is,ib) = ( hpolt1(is,ib,1)*cpolt1(is,ib,1) + hpolt1(is,ib,2)*cpolt1(is,ib,2) )
	enddo
	enddo
	
	print *, '***********************'	
	print *, '*** with shock'
	
	csimdy = 0.0D0; caggsimdy = 0.0D0;	itshow = 0;	dois   = 1;
	if (dois == 1) then
	do ib0 = 1, Nbeta	
	do is0 = 1, NS
	
		GMAt          = 0.0D0;	
		GMAt(is0,ib0) = 1.0D0; 		
		
		do tt = 1, Tsim
			GMAtp1 = 0.0D0
			
			do ibt = 1, Nbeta
			do ist = 1, NS			
				ixt = xind(ist)
				
				if (tt == 1) then
					csimdy(is0,ib0) = csimdy(is0,ib0) + (cseht1(ist,ibt)*GMAt(ist,ibt))
					caggsimdy(tt)   = caggsimdy(tt)   + (cseht1(ist,ibt)*GMAt(ist,ibt))*mu(is0,ib0)				
				else
					csimdy(is0,ib0) = csimdy(is0,ib0) + (cseh(ist,ibt)*GMAt(ist,ibt))
					caggsimdy(tt)   = caggsimdy(tt)   + (cseh(ist,ibt)*GMAt(ist,ibt))*mu(is0,ib0)				
				endif
				
				
				do ih = 1, Nh
					if (tt==1) then
						a_np  = max( min( apolt1(ist,ibt,ih) , amax), amin)
					else
						a_np  = min( apol(ist,ibt,ih) , amax)
					endif
					
					aindv   = vsearchprm_FET([a_np],avec,avecpar)		
					aind_np = aindv(1)
					da      = (a_np - avec(aind_np))/(avec(aind_np+1) - avec(aind_np))
					
					do ib_np = 1, Nbeta
					do ix_np = 1, Nx
						is_np = aind_np + Na*(ix_np - 1)
						
						if (tt==1) then
							GMAtp1(is_np  ,ib_np) = GMAtp1(is_np  ,ib_np) + (1.0D0-da)*Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpolt1(ist,ibt,ih)*GMAt(ist,ibt)
							GMAtp1(is_np+1,ib_np) = GMAtp1(is_np+1,ib_np) +    da     *Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpolt1(ist,ibt,ih)*GMAt(ist,ibt)		
						else
							GMAtp1(is_np  ,ib_np) = GMAtp1(is_np  ,ib_np) + (1.0D0-da)*Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)
							GMAtp1(is_np+1,ib_np) = GMAtp1(is_np+1,ib_np) +    da     *Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)		
						endif
			
						
					enddo
					enddo
				
				enddo !end loop ih 
			enddo !end loop ist
			enddo !end loop ibt
			
			! print *, 'sum( GMAtp1 ) = ', sum( GMAtp1 )
			GMAtp1 = GMAtp1/sum(GMAtp1)		
			GMAt   = GMAtp1	
			
		enddo !ned loop tt
		
		if (itshow<=0) then
			print *, 'is0 = ', is0,'/',NS,', ib0 = ', ib0,'/',Nbeta
			itshow = 1500
		endif
		itshow = itshow-1
	
	enddo !end loop is0	
	enddo !end loop ib0	
	endif
	
!	print *, 'Cagg = ', Cagg
!	print *, 'caggsim   = ', caggsim
!	print *, 'caggsimdy = ', caggsimdy
	
	
	! print *, 'csim(1:10,:) = ', csim(1:10,:)
	!----------------------------------------------------------
	!----------------------------------------------------------
	
	!----------------------------------------------------------
	!---ya nos vamos
	
	mpc_annual = (csimdy - csim)/dy
	
	print *, 'mpc annual = ', sum(mpc_annual*mu)
	
	open(1,  file = 'OUTPUT/mpc_annual.txt',   status = 'unknown')    	
	open(2,  file = 'OUTPUT/csim.txt'      ,   status = 'unknown')    	
	open(3,  file = 'OUTPUT/csimdy.txt'    ,   status = 'unknown')    	
    do is = 1, NS
        write(1, '(*(f18.6))') ( mpc_annual(is,ib) , ib = 1,Nbeta )        
		write(2, '(*(f18.6))') ( csim(is,ib)       , ib = 1,Nbeta )        
		write(3, '(*(f18.6))') ( csimdy(is,ib)     , ib = 1,Nbeta )        
	enddo		
	close(1); close(2); close(3);
	!----------------------------------------------------------
	!----------------------------------------------------------
	
	
	end subroutine COMPUTEmpc_annual

end module MODULEmpc_annual
