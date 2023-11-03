module FUNCTIONS_lpe
    
    use GLOBALS
    use Toolbox

contains
    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Compute Ve       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ComputeVe_TR(tt_in)
	
    integer :: tt_in
    integer :: ix_np, is_np(NS), ib
	real(8) :: Dwrk(NS), Dnow(NS), corr(NS), Pbt(Nbeta,Nbeta)
	
    do ib = 1, Nbeta

        !---integrate Gumbel shock
        Dnow = Vnow_TR(:,ib,tt_in);  Dwrk = Vwrk_TR(:,ib,tt_in);
        corr = max( (Dnow/rhow) , (Dwrk/rhow) )
        
        
        Vtl_TR(:,ib,tt_in) = rhow*log( exp((Dnow/rhow)-corr) + exp((Dwrk/rhow)-corr) ) + rhow*corr
        
        
        !---prob of working
        hpol_TR(:,ib,tt_in,1) = exp((Dnow/rhow)-corr)/(exp((Dnow/rhow)-corr)+exp((Dwrk/rhow)-corr))  ! prob of h=0
        hpol_TR(:,ib,tt_in,2) = 1.0D0 - hpol_TR(:,ib,tt_in,1)
        
        !---integragte x
        Ve_TR(:,ib,tt_in) = 0.0D0
        do ix_np = 1, Nx
            is_np = aind + Na*( ix_np - 1)
            
            Ve_TR(:,ib,tt_in) = Ve_TR(:,ib,tt_in) + Px(xind,ix_np)*Vtl_TR(is_np,ib,tt_in)
        enddo

    enddo

    Pbt = TRANSPOSE(Pbmat)
    Ve_TR(:,:,tt_in)  = MATMUL(Ve_TR(:,:,tt_in),Pbt)
	
    end subroutine ComputeVe_TR
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!         VF eval a scalar         !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function VFa_eval_TR(a_in,h_in,is_in,ib_in,tt_in,svdum)  result(V_out)
	
	implicit none
	real(8), intent(in)    :: a_in, h_in
	integer, intent(in)    :: is_in, ib_in, tt_in, svdum
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
	
	yl   = wgeH_TR(tt_in)*hh*S(is_in,2)
	rba  = RbFUN_TR(aa,tt_in)		
    yk   = rba*S(is_in,1)
    TAX  = TAXFUN_TR(yl,yk,tt_in)
	! TAX  = TAX + DTAX_TR(tt_in)*abs(TAX)
	
	INC  = yl + yk + S(is_in,1) + TF_TR(tt_in) - TAX        
	
	cc   = INC - anp
	cc   = max(cc,tiny)
	
	! uu = log(cc) - B*( (hh**(1.0D0+varphi)) / (1.0D0+varphi) ) - B0*hg0
	! uu = log(cc) - B0*hg0
	uu = log(cc) - uBvec(ib_in)*hg0
	
	
	aindv_np = vsearchprm_FET([anp],avec,avecpar)
	aind_np  = aindv_np(1)
    Sind_np  = aind_np + Na*(xind(is_in)-1)    
	
	 
    V_np = Ve_TR(Sind_np,ib_in,tt_in+1) + (anp - avec(aind_np)) * ( Ve_TR(Sind_np+1,ib_in,tt_in+1) - Ve_TR(Sind_np,ib_in,tt_in+1) )/(avec(aind_np+1) - avec(aind_np))
    
	V_out = uu + bvec(ib_in)*V_np  
	
	if (svdum == 1) then
		if (hh == 0.0D0) then
			cpol_TR(is_in,ib_in,tt_in,1) = cc     !	hpol(1,is_in,ib_in) = hh
		elseif (hh > 0.0D0) then
			cpol_TR(is_in,ib_in,tt_in,2) = cc     !  hpol(2,is_in,ib_in) = hh
		else
			print *, 'wrong hh = ', hh
		endif
		
	endif	
	
	end function VFa_eval_TR
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Taxes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function TAXFUN_TR(yl,yk,tt_in) result(T_out)
    
    implicit none    
    real(8), intent(in) :: yl, yk    
	integer, intent(in) :: tt_in
    real(8) 			:: T_out, taxl
    
    T_out = tauk_TR(tt_in)*yk
    if (yl > 0.0D0) then
		taxl  = yl - lbd_TR(tt_in)*(yl**(1.0D0-gma_TR(tt_in)))
		taxl  = taxl + DTAX_TR(tt_in)*abs(taxl)
		T_out = T_out + taxl 
        ! T_out = T_out + yl - lbd_TR(tt_in)*(yl**(1.0D0-gma_TR(tt_in)))
    endif
    
    end function TAXFUN_TR	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Rb function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	function RbFUN_TR(bin,tt_in) result(rb_out)
	
	implicit none
	real(8), intent(in) :: bin
	integer, intent(in) :: tt_in
	real(8)             :: rb_out
	
	rb_out = r_TR(tt_in)

	end function RbFUN_TR
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Golden Search       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function SolveBella_TR(xl1,xr1,h_in,is_in,ib_in,tt_in) result(x_out)    

    implicit none

    real(8), intent(in) :: xl1, xr1, h_in
    integer, intent(in) :: is_in, ib_in, tt_in
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
    
    f1 = VFa_eval_TR(x1,h_in,is_in,ib_in,tt_in,0)
    f2 = VFa_eval_TR(x2,h_in,is_in,ib_in,tt_in,0)
        
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
            f1new = VFa_eval_TR(x1-d,h_in,is_in,ib_in,tt_in,0)            
            
        else
        
            x2new = x2 + d
            x1new = x2
            
            f2new = VFa_eval_TR(x2+d,h_in,is_in,ib_in,tt_in,0)            
            f1new = f2
        
        endif
        
    end do
    
    if ( f2 < f1) then
        x_out = x1
    else
        x_out = x2
    endif
    
    end function SolveBella_TR

    
end module FUNCTIONS_lpe
