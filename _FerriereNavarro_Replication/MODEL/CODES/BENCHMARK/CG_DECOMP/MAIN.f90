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
    
    !---calibrated parameters
    open(1,  file = '../CG/OUTPUT/CALIB_YandB_save.txt',      status = 'unknown')
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
    Gshock  = 1.0D0*0.01D0;
    rhoG = 0.90D0
    G_TR    = 0.0D0; 
    G_TR(1) = G*(1.0D0+Gshock)
    do tt = 2, T_TR
        G_TR(tt) = (1.0D0-rhoG)*G + rhoG*G_TR(tt-1)
    enddo  

    index_decomp = 1
    ! index = 1 to feed all sequences and save _TR policies at the hh level;
    ! index = 2 to feed all but taxes
    ! index = 3 to feed only taxes

    call Compute_sequential_decomp

end program MAIN
