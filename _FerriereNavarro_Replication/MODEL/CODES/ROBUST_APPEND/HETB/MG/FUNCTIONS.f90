module FUNCTIONS
    
    use GLOBALS
    use Toolbox
	! use Toolbox_CE
contains
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Compute Ve       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ComputeVe
	
    integer :: ib
    integer :: ix_np, is_np(NS)
	real(8) :: Dwrk(NS), Dnow(NS), corr(NS), Pbt(Nbeta,Nbeta)

    do ib = 1, Nbeta
	
	!---integrate Gumbel shock
	Dnow = Vnow(:,ib);  Dwrk = Vwrk(:,ib);
	corr = max( (Dnow/rhow) , (Dwrk/rhow) )
	
	Vtl(:,ib) = rhow*log( exp((Dnow/rhow)-corr) + exp((Dwrk/rhow)-corr) ) + rhow*corr
	
	
	!---prob of working
	hpol(:,ib,1) = exp((Dnow/rhow)-corr)/(exp((Dnow/rhow)-corr)+exp((Dwrk/rhow)-corr))  ! prob of h=0
	hpol(:,ib,2) = 1.0D0 - hpol(:,ib,1)
	
	!---integragte x
    Ve(:,ib) = 0.0D0
    do ix_np = 1, Nx
        is_np = aind + Na*( ix_np - 1)
        
        Ve(:,ib) = Ve(:,ib) + Px(xind,ix_np)*Vtl(is_np,ib)
    enddo

    enddo

    Pbt = TRANSPOSE(Pbmat)
    Ve  = MATMUL(Ve,Pbt)

    end subroutine ComputeVe
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!         VF eval a scalar         !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function VFa_eval(a_in,h_in,is_in,ib_in,svdum)  result(V_out)
	
	implicit none
	real(8), intent(in)    :: a_in, h_in
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
	yl   = wgeH*hh*S(is_in,2)
    yk   = rba*S(is_in,1)
    TAX  = TAXFUN(yl,yk)
	
	INC  = yl + yk + S(is_in,1) + TF - TAX
	
	cc   = INC - anp
	cc   = max(cc,tiny)
	
	! uu = log(cc) - B*( (hh**(1.0D0+varphi)) / (1.0D0+varphi) ) - B0*hg0
	! uu = log(cc) - B0*hg0
	uu = log(cc) - uBvec(ib_in)*hg0
	
	
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
	
	end function VFa_eval
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!         VF eval for policies     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function VFvec_itH(hh,ib_in,aind_np)  result(V_out)

    implicit none

    integer, intent(in) :: ib_in, aind_np(NS)
    real(8), intent(in) :: hh
    real(8)             :: V_out(NS), cin(NS), ain(NS)
    integer             :: Sind_np(NS)

    real(8), dimension(NS) ::  uh, uu, V_np, hg0

    if (hh == 0.0D0) then
        cin = cpol(:,ib_in,1)
        ain = apol(:,ib_in,1)
        hg0 = 0.0D0
    elseif (hh > 0.0D0) then
        cin = cpol(:,ib_in,2)
        ain = apol(:,ib_in,2)
        hg0 = 1.0D0
    else
        print *, 'wrong hh = ', hh
    endif

    ! uu = log(cin) - B0*hg0
	uu = log(cin) - uBvec(ib_in)*hg0
	

    !    aind_np = vsearchprm_FET(ain,avec,avecpar)
    Sind_np = aind_np + Na*(xind-1)

    V_np = Ve(Sind_np,ib_in) + (ain - avec(aind_np)) * ( Ve(Sind_np+1,ib_in) - Ve(Sind_np,ib_in) )/(avec(aind_np+1) - avec(aind_np))

    V_out = uu + bvec(ib_in)*V_np

    end function VFvec_itH
	
	
	!----------------------------------------
    !---Taxes
    !----------------------------------------
	function TAXFUN(yl,yk) result(T_out)
    
    implicit none    
    real(8), intent(in) :: yl, yk    
    real(8) 			:: T_out
    
    T_out = tauk*yk
    if (yl > 0.0D0) then
        T_out = T_out + yl - lbd*(yl**(1.0D0-gma))
    endif
    
    end function TAXFUN	
	
	!----------------------------------------
    !---Rb function
    !----------------------------------------
	function RbFUN(bin) result(rb_out)
	
	implicit none
	real(8), intent(in) :: bin
	real(8)             :: rb_out
	
!	if (bin>=0.0D0) then
		rb_out = rk
!	else
!		rb_out = rk + rbcost
!	endif
	end function RbFUN
	
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Golden Search       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function SolveBella(xl1,xr1,h_in,is_in,ib_in) result(x_out)    

    implicit none

    real(8), intent(in) :: xl1, xr1, h_in
    integer, intent(in) :: is_in, ib_in
    real(8)             :: x_out
    real(8)             :: xl, xr
    
    real(8)            :: alpha1, alpha2
    real(8)            :: d, x1, x2, f1, f2, x1new, x2new, f1new, f2new
    real(8), parameter :: tol = 1e-8
    
    
    alpha1 = 0.5D0*(3.0D0 - sqrt(5.0D0))
    alpha2 = 0.5D0*(sqrt(5.0D0) - 1.0D0)
    
    xl = xl1
    xr = xr1
    
    if(xr < xl) xr = xl   
    
    d  = xr - xl
    
    x1 = xl + alpha1*d  
    x2 = xl + alpha2*d  
    
    f1 = VFa_eval(x1,h_in,is_in,ib_in,0)
    f2 = VFa_eval(x2,h_in,is_in,ib_in,0)
        
    d  = alpha1*alpha2*d
    
    f1new = f1
    f2new = f2
    x1new = x1
    x2new = x2
    
    do while ( d > tol) 
        
        ! print *, 'maxval(d) = ', maxval(d)
        
        f1 = f1new
        f2 = f2new
        x1 = x1new
        x2 = x2new
        
        d = d*alpha2
        
        if (f2 < f1 ) then
            
            x2new = x1
            x1new = x1 - d
            
            f2new = f1
            f1new = VFa_eval(x1-d,h_in,is_in,ib_in,0)            
            
        else
        
            x2new = x2 + d
            x1new = x2
            
            f2new = VFa_eval(x2+d,h_in,is_in,ib_in,0)            
            f1new = f2
        
        endif
        
    end do
    
    if ( f2 < f1) then
        x_out = x1
    else
        x_out = x2
    endif
    
    end function SolveBella

    
end module FUNCTIONS
