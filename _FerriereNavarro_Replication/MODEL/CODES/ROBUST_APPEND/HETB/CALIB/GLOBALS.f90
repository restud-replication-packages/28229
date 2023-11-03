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
	
    !---preference shocks
	real(8), parameter :: rhow = 0.0605D0      ! for hours

	!---taxes
    real(8) :: gma, tauk
	
	!---technology for NK block
	real(8), parameter :: epsPi   = 7.0D0                       ! elasticity of substitution across nk
    real(8), parameter :: thetaPi = 200.0D0                     ! cost of adjusting prices -> only for transition
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
    real(8), parameter :: Piwbar = 1.0d0, phiwPi = 1.5d0, epswPi = 7.0d0, thetawPi = 200.0D0
	
    !---parameters to calibrate
    real(8) :: G, DB, B0, TF
	
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
    real(8)            :: avec(Na), xvec(Nx), bvec(Nbeta), uBvec(Nbeta)	
    integer, parameter :: Nh = 2
    real(8)            :: hvec(Nh)
	real(8), parameter :: avecpar = 0.30D0
	
    real(8), parameter :: mx = 3
    real(8)            :: x_lb, x_ub
    real(8)            :: Px(Nx,Nx), pix(Nx), Ex ! , pixcdf(Nx), phidx(Nx), mupx
	
	real(8)            :: Pbeta(Nbeta), Pbmat(Nbeta,Nbeta)

	real(8)            :: amin = 0.0D0,    amax = 700.0D0
    
    integer, parameter :: NS = Na*Nx
    real(8)            :: S(NS,2)  ! vectorized states
    integer            :: aind(NS), xind(NS)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---Value functions and policies
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), dimension(NS,Nbeta)    :: V, Vwrk, Vnow, Vtl, Ve, mu ! , Vn, Va, mu  	
	real(8), dimension(NS,Nbeta,Nh) :: hpol, apol, cpol
	real(8), dimension(NS)          :: Tdx

    integer, dimension(Ns,Nbeta,Nh) :: aind_np_ho

    
    real(8) :: Vwnew(NS)
    
    integer :: savedum	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !---calibration objects
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(8) :: mborr  !  pmuadj(Np),	
	
	
    !---for lpe computation
    integer, parameter :: T_TR = 20

    real(8), dimension(T_TR) :: DTAX_TR

    real(8), dimension(T_TR)   :: wge_TR, wgeH_TR, rk_TR, r_TR, qk_TR, tauk_TR, lbd_TR, gma_TR, TF_TR
    real(8), dimension(T_TR)   :: Cagg_TR, Lagg_TR, Hagg_TR
    real(8), dimension(T_TR+1) :: Aagg_TR

    !---value function and policies
    real(8), dimension(NS,Nbeta,T_TR+1)    :: Vwrk_TR, Vnow_TR, Vtl_TR, Ve_TR
    real(8), dimension(NS,Nbeta,T_TR+1,Nh) :: hpol_TR, apol_TR, cpol_TR
    real(8), dimension(NS,Nbeta,T_TR)      :: mu_TR

    real(8), dimension(NS,T_TR)      :: Tdx_TR = 0d0

    !---capital good producer
    real(8), parameter :: qk0  = 1.00D0   ! steady-state price of capital
    real(8), parameter :: phiK = 15.30D0
			
end module GLOBALS
