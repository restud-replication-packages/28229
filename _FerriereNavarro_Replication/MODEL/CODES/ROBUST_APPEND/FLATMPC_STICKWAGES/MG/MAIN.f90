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
    
    !---calibated parameters
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
    tauk_TR  = tauk
    TF_TR    = TF
    DB_adj	 = 0.5D0  ! =1 for no deficit financing
	
    !gma_TR   = gma
    gma_TR(1) = 0.11d0
    rhoG = 0.9D0
    do tt = 2, T_TR
        gma_TR(tt) = (1.0D0-rhoG)*gma + rhoG*gma_TR(tt-1)
    enddo


    !---G shock
    Gshock  = 1.0D0*0.01D0; rhoG = 0.90D0
    G_TR    = 0.0D0; 
    G_TR(1) = G*(1.0D0+Gshock)
    do tt = 2, T_TR
        G_TR(tt) = (1.0D0-rhoG)*G + rhoG*G_TR(tt-1)
    enddo  

	call Compute_QNEWTON

end program MAIN
