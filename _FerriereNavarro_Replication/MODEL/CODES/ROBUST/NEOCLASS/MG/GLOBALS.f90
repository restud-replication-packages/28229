module GLOBALS
    
    !-------------------------------------------------------
    !---Parameters
    !-------------------------------------------------------
    !---household
    real(8), parameter :: sgc    = 2.0D0                        ! utility consumption curvature 
    ! real(8)            :: bta    = 0.985D0                      ! discount factor -> calibration
    real(8), parameter :: hbar   = 0.33D0                       ! time endowment
    real(8), parameter :: cbar   = 0.00D0                       ! subsitence consumption
    
    real(8), parameter :: rhox = 0.939D0, sgx = 0.287D0       ! persistence and std for productivity        	
    
    !---Technology    
    real(8), parameter :: alpha = 0.66D0                        ! labor exopnent on prod function 
    real(8), parameter :: dlta  = 0.0235D0                       ! depreciation rate
	
	!---betas
	real(8) :: btaH, Dbta
	
	real(8), parameter :: rhow = 0.066D0      ! for hours
	
	!---taxes
    real(8) :: gma, tauk
	
	!---technology for NK block
	real(8), parameter :: epsPi   = 7.0D0                       ! elasticity of substitution across nk
    real(8), parameter :: thetaPi = 0d0 !200.0D0                     ! cost of adjusting prices -> only for transition
    real(8)            :: swH                                   ! subsidy     
    real(8)            :: FC, FCW                               ! fixed production cost -> should add it to calibration
	real(8)            :: MC									! marginal cost in setady-state
	! real(8), parameter :: zbar = epsPi/(epsPi-1.0D0)    		! igp productivity
	real(8), parameter :: zbar = 1.0D0     		! igp productivity
	
	!---capital good producer
	
	!---monetary block
    real(8), parameter :: Pibar = 1.0D0			! Inflation target
	real(8), parameter :: phiPI = 1.5D0
	real(8)            :: ibar

    ! Sticky wages
    real(8), parameter :: Piwbar = 1.0d0, phiwPi = 1.5d0, epswPi = 7.0d0, thetawPi = 0d0 !200.0D0
	
    !---parameters to calibrate(ish)
    real(8) :: G, DB, B, B0, TF    
	
	!---calibration targets
	real(8), parameter :: GYcal  = 0.10D0
	real(8), parameter :: TFcal  = 0.081D0
	real(8), parameter :: DBcal  = 1.0D0
	real(8), parameter :: EMPcal = 0.75D0,      EMPtol = 0.01D0
	real(8), parameter :: Hcal   = 1.0D0/3.0D0, Htol   = 0.02D0
    real(8)            :: Yaggcal 
    real(8), parameter :: Ytol = 0.01D0 
    
    !---prices
    real(8)            :: wge, qk, lbd, wgeH
    real(8), parameter :: rk = -1.0D0 + (1.035D0**(1.0D0/4.0D0)) ! rk = 0.005D0
	
	
	!---aggregate quantities
	real(8) :: Cagg, Aagg, Kagg, Kdmn, KYxx, Hagg, Lagg, Yagg, EMP, ylmean
	real(8) :: BC, Rev, RevK, RevL, INTEGRAL, DEFnoL, Dk	
	
	!---iteration	
	real(8) :: rknew, lbdnew, Ynew, what	
	real(8) :: blb, bub, rklb, rkub, lbdlb, lbdub, Yagglb, Yaggub, B0lb, B0ub
	real(8) :: albd, xlbd
	
	!---additional stuff
    real(8), parameter :: tiny  = 1E-13                 ! lower bound for consumption    
	real(8), parameter :: Vlow  = -950000 

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Grids
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer, parameter :: Na = 351, Nx = 15, Nbeta = 3
    real(8)            :: avec(Na), xvec(Nx), bvec(Nbeta)	
    integer, parameter :: Nh = 2
    real(8)            :: hvec(Nh)
	real(8), parameter :: avecpar = 0.3D0
	
    real(8), parameter :: mx = 3
    real(8)            :: x_lb, x_ub
    real(8)            :: Px(Nx,Nx), pix(Nx), Ex ! , pixcdf(Nx), phidx(Nx), mupx
	
	real(8)            :: Pbeta(Nbeta), Pbmat(Nbeta,Nbeta)

    ! real(8)            :: bmin = -2.0D0,     bmax = 40.0D0
	real(8)            :: amin = 0.0D0,    amax = 700.0D0
    
    integer, parameter :: NS = Na*Nx
    real(8)            :: S(NS,2)  ! vectorized states
    integer            :: aind(NS), xind(NS)
    ! real(8) :: Sbkx(Nb,Nk,Nx)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Value functions and policies
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), dimension(NS,Nbeta)    :: V, Vwrk, Vnow, Vtl, Ve, mu ! , Vn, Va, mu  	
	real(8), dimension(NS,Nbeta,Nh) :: hpol, apol, cpol
	real(8), dimension(NS)          :: Tdx
    
    integer, dimension(Ns,Nbeta,Nh) :: aind_np_ho


    real(8) :: Vwnew(NS)
    
    integer :: savedum	

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!         				 Transition        					 !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, parameter :: T_TR = 252
    
	!---Shock and gma_TR
	real(8) :: Gshock, rhoG
    real(8) :: Tshock, rhoT
	real(8) :: rhogma
	real(8) :: DB_adj 
	
	!---fiscal
	real(8), dimension(T_TR)   :: tauk_TR, gma_TR, TF_TR, G_TR
	real(8), dimension(T_TR)   :: RevL_TR, RevK_TR, Rev_TR
	real(8), dimension(T_TR+1) :: DB_TR 
	
	!---monetary
	real(8), dimension(T_TR) :: Pi_TR, i_TR
	
	!---Aggregates	
	real(8), dimension(T_TR)   :: wge_TR, rk_TR, r_TR, qk_TR, lbd_TR, lbdb_TR, xlbd_TR, albd_TR, div_TR, divb_TR, divY_TR, divK_TR, divBANK_TR, theta_TR, rkb_TR, wgeb_TR, qkb_TR
    real(8), dimension(T_TR)   :: wgeH_TR, wgeHb_TR, Piw_TR, thetaw_TR, divw_TR
	real(8), dimension(T_TR)   :: Cagg_TR, Lagg_TR, Yagg_TR, Ld_TR, Ysup_TR, Dk_TR, Iagg_TR, MC_TR
	real(8), dimension(T_TR+1) :: Aagg_TR, Kagg_TR
	
	!---eqilibrium computation
	real(8), dimension(T_TR) :: BC_TR, DEFnoL_TR
	real(8)                  :: Tdx_TR(NS,T_TR)
		
	
	!---Newton set-up
	integer, parameter :: Nvar = 3               ! number of variables to solve for: {r,w,lbd_t, div,q}
	integer, parameter :: Nprm = Nvar*T_TR       ! Pi_t  from t=1 to t=T_TR
	integer, parameter :: Kuvar = Nprm
	real(8)            :: JACmat(Nprm,Nprm), ED_ss(Nprm)
	real(8)            :: yylb_TR(Nvar),yyub_TR(Nvar)
	
	
	!---value function and policies
	real(8), dimension(NS,Nbeta,T_TR+1)    :: Vwrk_TR, Vnow_TR, Vtl_TR, Ve_TR
	real(8), dimension(NS,Nbeta,T_TR+1,Nh) :: hpol_TR, apol_TR, cpol_TR
    real(8), dimension(NS,Nbeta,T_TR)      :: mu_TR
    real(8), dimension(NS,Nbeta) :: index_tax_TR
	
	!---to read inputs and save
	character(len=5)     :: iWxs
	integer              :: iWx 
	character(len = 250) :: wvS
	
	integer, parameter :: NC  = 12  ! NC  = # of computers you're calling
	integer, parameter :: JxC = 63  ! JxC = # of columns ran by each computer --> need NC*JxC = Nprm !!  
		


    !---capital good producer
    real(8), parameter :: qk0  = 1.00D0   ! steady-state price of capital
    real(8), parameter :: phiK = 0d0 !15.30D0

	
end module GLOBALS
