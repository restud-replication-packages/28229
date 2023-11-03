program MAIN

	use GLOBALS
    use ModuleINIT
	use ModuleSTEADY
	use ModuleTRANSITION

	implicit none    
    real(8) :: CALIB(2)
    integer :: tt
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Parameters      !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !---tax parameters
    gma  = 0.10D0    
    tauk = 0.35D0
    
    !---calibrtaed parameters
    B    = 2.054D0 ! 50.93D0
    DB   = 3.46D0
    G    = 0.22D0
    TF   = (1.0D0/3.0D0)*G
    B0   = 0.5D0
    ! read parameters if you have    
    ! open(1,  file = 'OUTPUT/CALIBparm.txt',      status = 'unknown')
    ! read(1, *) CALIB 					
    ! close(1)        
	! G  = CALIB(1); DB = CALIB(2); TF = CALIB(3); B  = CALIB(4); B0 = CALIB(5); 	  ! G  = CALIB(1); DB = CALIB(2); TF = CALIB(3); B  = CALIB(4)		
	
    open(1,  file = 'OUTPUT/CALIB_YandB_save.txt',      status = 'unknown')
    read(1, *) CALIB
    close(1)
	Yagg = CALIB(1); B0 = CALIB(2);
	G  = GYcal*Yagg  
	TF = TFcal*Yagg  
	DB = DBcal*Yagg  
    
    call BuildGrids    
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Steady-State       !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    call Model_STEADYSTATE
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!       Transition      !!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    
    !---select policies
    gma_TR   = gma    
    tauk_TR  = tauk  
    TF_TR    = TF
    DB_adj	 = 0.5D0  ! =1 for no deficit financing
	
    
    !---G shock
    Gshock  = 1.0D0*0.01D0; rhoG = 0.90D0
    G_TR    = 0.0D0; 
    G_TR(1) = G*(1.0D0+Gshock)
    do tt = 2, T_TR
        G_TR(tt) = (1.0D0-rhoG)*G + rhoG*G_TR(tt-1)
    enddo  

	!---gma path
	! gma_TR(1) = 0.11D0; rhogma = 0.80D0
	! do tt = 2, T_TR
	! 	gma_TR(tt) = (1.0D0-rhogma)*gma + rhogma*gma_TR(tt-1)
	! end do
    
    ! call Compute_JACOB 
	! call Compute_LINEARSOL
	call Compute_QNEWTON
!    call Compute_sequential

end program MAIN
