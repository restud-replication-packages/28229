module ModuleMPC

	use GLOBALS
	use FUNCTIONS
	use Toolbox

contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!           compute mpc            !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine COMPUTEmpc
	
	implicit none
	integer :: ib, ih, is, aind_dy(NS), Sind_dy(NS), aidy, Sidy
	real(8) :: ymn, ymnDATA, cvMtoD, yl(NS), yk(NS), rba(NS), dy, ady(NS), Caggx, Caggdy, mpc(NS,Nbeta), mpc_mn, mpc_dCagg 
	real(8) :: cdy(NS,Nbeta,Nh), Vnowdy(NS,Nbeta), Vwrkdy(NS,Nbeta), Dnow(NS), Dwrk(NS), corr(NS), hpoldy(NS,Nbeta,Nh), cpolE(NS,Nbeta), cpolEdy(NS,Nbeta)
	
	!---compute ymean
	ymn = 0.0D0
	do ib = 1, Nbeta
	do ih = 1, Nh
		do is = 1, NS
			yl(is)  = wge*hvec(ih)*S(is,2)
			rba(is) = RbFUN(S(is,1))	
			yk(is)  = rba(is)*S(is,1)
			
			ymn = ymn + (yl(is) + yk(is))*hpol(is,ib,ih)*mu(is,ib)
		enddo
	enddo
	enddo
	
	ymnDATA = 62000.0D0
	cvMtoD  = ymnDATA/ymn
	
	!---compute transfer dC and new policies
	dy = 500.0D0/cvMtoD			! mpc out of dC dollars
	
	!---compute consumption at a+dy <--> interpolate policies	
	ady     = S(:,1) + dy
	aind_dy = vsearchprm_FET(ady,avec,avecpar)
	Sind_dy = aind_dy + Na*(xind-1)   
	
	do ib = 1, Nbeta
	do ih = 1, Nh
		do is = 1, NS
			aidy = aind_dy(is)
			Sidy = Sind_dy(is)
		
			cdy(is,ib,ih) = cpol(Sidy,ib,ih) + (ady(is) - avec(aidy)) * ( cpol(Sidy+1,ib,ih) - cpol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))			
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
		
		hpoldy(:,ib,1) = exp((Dnow/rhow)-corr)/(exp((Dnow/rhow)-corr)+exp((Dwrk/rhow)-corr))  ! prob of h=0
		hpoldy(:,ib,2) = 1.0D0 - hpoldy(:,ib,1) 	
		
		
	enddo
	
	!---average cons -- integrate over h_eps
	cpolE = 0.0D0; cpolEdy = 0.0D0;  
	Caggx = 0.0D0; Caggdy  = 0.0D0; mpc_mn = 0.0D0 
	do ib = 1, Nbeta
		do ih = 1, Nh		
			cpolE(:,ib)   = cpolE(:,ib)   + hpol(:,ib,ih)  *cpol(:,ib,ih)		
			
			cpolEdy(:,ib) = cpolEdy(:,ib) + hpoldy(:,ib,ih)*cdy(:,ib,ih)		
		enddo
		
		Caggx  = Caggx  + sum(cpolE(:,ib)  *mu(:,ib))
		Caggdy = Caggdy + sum(cpolEdy(:,ib)*mu(:,ib))
		
		mpc(:,ib) = (cpolEdy(:,ib)-cpolE(:,ib))/dy	
		
		mpc_mn = mpc_mn + sum(mpc(:,ib)*mu(:,ib))
	enddo
		
	mpc_dCagg = (Caggdy-Caggx)/dy
	
	print *, 'Cagg = ', Cagg, ', Caggx = ', Caggx
	print *, 'dy = ', dy
	print *, 'Caggdy = ', Caggdy
	print *, 'Caggdy-Caggx = ', Caggdy-Caggx	
	print *, 'mpc_mn = ', mpc_mn
	print *, 'mpc_dCagg = ', mpc_dCagg
	
	print *, 'annual mpc_mn = ', 1.0D0 - ((1.0D0-mpc_mn)**4.0D0)
	
	
	
	open(1,  file = 'OUTPUT/mpc.txt',   status = 'unknown')    	
    do is = 1, NS
        write(1, '(*(f18.6))') ( mpc(is,ib) , ib = 1,Nbeta )        
	enddo		
	close(1); 
	
	end subroutine COMPUTEmpc
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!         VF eval a scalar         !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function VFa_eval_MPC(a_in,h_in,is_in,ib_in,dy,svdum)  result(V_out)
	
	implicit none
	real(8), intent(in)    :: a_in, h_in, dy
	integer, intent(in)    :: is_in, ib_in, svdum
	real(8)                :: V_out
	
	real(8) :: aa, rba, yl, yk, TAX, INC, anp, hh, cc, uu, hg0, V_np
	integer :: aindv_np(1), aind_np, Sind_np
	
	aa   = S(is_in,1)
	anp  = a_in
	hh   = h_in
	
	hg0 = 0.0D0
	if (hh > 0.0D0) then
		hg0 = 1.0D0
	endif
	
	rba  = RbFUN(aa)	
	yl   = wge*hh*S(is_in,2)
    yk   = rba*S(is_in,1)
    TAX  = TAXFUN(yl,yk)
	
	INC  = yl + qk*yk + qk*S(is_in,1) + TF + Dy - TAX        
	
	cc   = INC - qk*anp
	cc   = max(cc,tiny)
	
	! uu = log(cc) - B*( (hh**(1.0D0+varphi)) / (1.0D0+varphi) ) - B0*hg0
	uu = log(cc) - B0*hg0
	
	aindv_np = vsearchprm_FET([anp],avec,avecpar)
	aind_np  = aindv_np(1)
    Sind_np  = aind_np + Na*(xind(is_in)-1)    
	
	 
    V_np = Ve(Sind_np,ib_in) + (anp - avec(aind_np)) * ( Ve(Sind_np+1,ib_in) - Ve(Sind_np,ib_in) )/(avec(aind_np+1) - avec(aind_np))
    
	V_out = uu + bvec(ib_in)*V_np  
	
	if (svdum == 1) then
		if (hh == 0.0D0) then
			cpol(is_in,ib_in,1) = cc     !	hpol(1,is_in,ib_in) = hh
		elseif (hh > 0.0D0) then
			cpol(is_in,ib_in,2) = cc     !  hpol(2,is_in,ib_in) = hh
		else
			print *, 'wrong hh = ', hh
		endif
		
	endif	
	
	end function VFa_eval_MPC	

end module ModuleMPC