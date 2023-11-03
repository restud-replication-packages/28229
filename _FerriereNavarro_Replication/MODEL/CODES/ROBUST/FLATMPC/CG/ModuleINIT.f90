module ModuleINIT
    
    use GLOBALS
    use ModuleSAVE
    use Toolbox    
contains

    !------------------------------
    !---Grids
    !------------------------------
    subroutine BuildGrids
    	
    implicit none    	
    integer :: is, ia, ix, itx, iav(1), iax, ib
	real(8) :: pixnew(Nx), errxv(Nx), errx(1), gxa(Na), avx(Na), dax, errxb(Nbeta), Pbetanew(Nbeta)
	
    
    !---productivity
    call Tauchen(rhox,sgx,Nx,mx,xvec,Px)
	
	!---ergodic x
	pix = 1.0D0/dble(Nx); errx = 8.0D0
	do itx = 1, 500
		pixnew = matmul(pix,Px)
		errxv  = pix - pixnew 
		errx   = maxval(abs(errxv))
		if (errx(1) < 1e-8) then
            Ex = sum(pix*xvec)

			exit
		else
			pix = pixnew
		endif	
    enddo
	
	
    !---assets grid	
    gxa  = LINSPACE_FET(0.0D0,1.0D0,Na)
	!-usual grid

    avec = amin + (amax-amin)*(gxa**(1.0D0/avecpar))
	
    
    !---Svec
    is=1
    do ix = 1, Nx
    do ia = 1, Na       
        S(is,:) = [avec(ia), xvec(ix)]
        
        xind(is) = ix
        aind(is) = ia            
        
        is=is+1 
    enddo
    enddo      
    ! print *, 'kvec = ', kvec
    
    !---h vec
    hvec = [0.0D0, hbar]
    
    !---phivec    
    ! phivec = [0.0D0, phip]

    ! Beta-bec

    Pbmat = 0d0
    Pbmat(1,2) = 1d0/200d0
    Pbmat(2,1) = 0.5d0/200d0
    Pbmat(2,3) = Pbmat(2,1)
    Pbmat(3,2) = Pbmat(1,2)

    do ib = 1, Nbeta
        Pbmat(ib,ib) = 1d0-1d0/200d0
    enddo

    !---ergodic nbeta
    Pbeta = 1.0D0/dble(Nbeta); errx = 8.0D0
    do itx = 1, 1500
        Pbetanew = matmul(Pbeta,Pbmat)
        errxb  = Pbeta - Pbetanew
        errx   = maxval(abs(errxb))
        if (errx(1) < 1e-8) then
!            print *, 'Nbeta = ', Pbeta
            exit
        else
        Pbeta =     Pbetanew
        endif
    enddo
    
    call SaveGRIDS
    end subroutine BuildGrids
    
    
end module ModuleINIT
