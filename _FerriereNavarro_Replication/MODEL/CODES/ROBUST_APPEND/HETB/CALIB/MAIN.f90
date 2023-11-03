program MAIN

	use GLOBALS
    use ModuleINIT
	use ModuleSTEADY

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

end program MAIN
