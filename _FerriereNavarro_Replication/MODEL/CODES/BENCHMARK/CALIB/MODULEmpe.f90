module MODULEmpe

	use GLOBALS
	use FUNCTIONS	
	use Toolbox
	
contains

	!-------------------------------------------	
	!---mpc annual computation
	!-------------------------------------------
	subroutine COMPUTEmpe
	
	implicit none
	integer, parameter :: Tsim = 20

    integer :: index_mpe
	
	real(8) :: ymn, ymnDATA, cvMtoD, yl(NS), yk(NS), rba(NS), dy, ady(NS), a_np, da,  yw(NS), ynw(NS)
	real(8) :: cpolt1(NS,Nbeta,Nh), apolt1(NS,Nbeta,Nh), hpolt1(NS,Nbeta,Nh), Vnowdy(NS,Nbeta), Vwrkdy(NS,Nbeta), Dnow(NS), Dwrk(NS), corr(NS), cslope(NS,Nbeta,Nh)
	real(8) :: GMAt(NS,Nbeta), GMAtp1(NS,Nbeta), GGx(NS,Nbeta), yseh(NS,Nbeta), yseht1(NS,Nbeta)
    real(8) :: ysim(NS,Nbeta), ysimdy(NS,Nbeta), mpe_annual(NS,Nbeta), yaggsim(Tsim), yaggsimdy(Tsim)
    integer :: ib, ih, is, tt, iz, iz0, ia0, ix0, ib0, is0, izt, iat, ixt, ibt, ist, ib_np, ix_np, is_np, iz_np, aind_dy(NS), Sind_dy(NS), aidy, Sidy, aindv(1), aind_np, itshow, dois


    do index_mpe = 1, 3
	
        !---compute ymean
        ymn = 0.0D0
        do ib = 1, Nbeta
        do ih = 1, Nh
            do is = 1, NS
                yl(is)  = wgeH*hvec(ih)*S(is,2)
                rba(is) = RbFUN(S(is,1))
                yk(is)  = rba(is)*S(is,1)
                
                ymn = ymn + (yl(is) + yk(is))*hpol(is,ib,ih)*mu(is,ib)
            enddo
        enddo
        enddo

        yw = lbd*(wgeH*hbar*S(:,2))**(1d0-gma) + TF
        ynw = 0d0  + TF

        
        ymnDATA = 67000.0D0
        cvMtoD  = ymnDATA/ymn
        
        !---compute transfer dC and new policies

        if (index_mpe==1) then
            print *, 'MPE SMALL SIZE'
            dy = 165000.0D0/cvMtoD            ! MPE out of dC dollars
        elseif (index_mpe==2) then
            print *, 'MPE MEDIUM SIZE'
            dy = 650000.0D0/cvMtoD            ! MPE out of dC dollars
        elseif (index_mpe==3) then
            print *, 'MPE LARGE SIZE'
            dy = 2000000.0D0/cvMtoD            ! MPE out of dC dollars
        endif

        !---compute c and a at a+dy <--> policy at t=1 -- interpolate policies
        ady     = S(:,1) + dy
        aind_dy = vsearchprm_FET(ady,avec,avecpar)
        Sind_dy = aind_dy + Na*(xind-1)
        
        do ib = 1, Nbeta
        do ih = 1, Nh
            do is = 1, NS
                aidy = aind_dy(is)
                Sidy = Sind_dy(is)
                
                cslope(is,ib,ih) = (ady(is) - avec(aidy)) * ( cpol(Sidy+1,ib,ih) - cpol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))
            
                cpolt1(is,ib,ih) = cpol(Sidy,ib,ih) + (ady(is) - avec(aidy)) * ( cpol(Sidy+1,ib,ih) - cpol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))
                apolt1(is,ib,ih) = apol(Sidy,ib,ih) + (ady(is) - avec(aidy)) * ( apol(Sidy+1,ib,ih) - apol(Sidy,ib,ih) )/(avec(aidy+1) - avec(aidy))
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
            
            hpolt1(:,ib,1) = exp((Dnow/rhow)-corr)/(exp((Dnow/rhow)-corr)+exp((Dwrk/rhow)-corr))  ! prob of h=0
            hpolt1(:,ib,2) = 1.0D0 - hpolt1(:,ib,1)
        enddo
        

        
        
        !----------------------------------------------------------
        !---compute income for t=1,..,20 with no shock
        do ib = 1, Nbeta
        do is = 1, NS
            ! iz = is + NS*(ib-1)
            yseh(is,ib) = ( hpol(is,ib,1)*ynw(is) + hpol(is,ib,2)*yw(is) )
        enddo
        enddo
        
        print *, '***********************'
        print *, '*** no shock'
        
        itshow = 0;	dois   = 1;
        ysim = 0.0D0; yaggsim = 0.0D0
        if (dois == 1) then
        do ib0 = 1, Nbeta
        do is0 = 1, NS
        
            GMAt          = 0.0D0;
            GMAt(is0,ib0) = 1.0D0;
            
            
            do tt = 1, Tsim
                GMAtp1 = 0.0D0
                
                do ibt = 1, Nbeta
                do ist = 1, NS
                    ixt = xind(ist)
                                    
                    ysim(is0,ib0) = ysim(is0,ib0) + (yseh(ist,ibt)*GMAt(ist,ibt))
                    yaggsim(tt)   = yaggsim(tt)   + (yseh(ist,ibt)*GMAt(ist,ibt))*mu(is0,ib0)

                    do ih = 1, Nh
                        a_np    = min( apol(ist,ibt,ih) , amax)
                        aindv   = vsearchprm_FET([a_np],avec,avecpar)
                        aind_np = aindv(1)
                        da      = (a_np - avec(aind_np))/(avec(aind_np+1) - avec(aind_np))
                        
                        do ib_np = 1, Nbeta
                        do ix_np = 1, Nx
                            is_np = aind_np + Na*(ix_np - 1)
                
                            GMAtp1(is_np  ,ib_np) = GMAtp1(is_np  ,ib_np) + (1.0D0-da)*Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)
                            GMAtp1(is_np+1,ib_np) = GMAtp1(is_np+1,ib_np) +    da     *Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)
                        enddo
                        enddo
                    
                    enddo !end loop ih
                enddo !end loop ist
                enddo !end loop ibt
                
                ! print *, 'sum( GMAtp1 ) = ', sum( GMAtp1 )
                GMAtp1 = GMAtp1/sum(GMAtp1)
                GMAt   = GMAtp1
            enddo !ned loop tt
            
            if (itshow<=0) then
                print *, 'is0 = ', is0,'/',NS,', ib0 = ', ib0,'/',Nbeta
                itshow = 1500
            endif
            itshow = itshow-1
        
        enddo !end loop is0
        enddo !end loop ib0
        endif
        
        !----------------------------------------------------------
        !----------------------------------------------------------
        
        ! STOP
        !----------------------------------------------------------
        !---compute consumption for t=1,..,4 with shock
        
        yseht1 = 0.0D0
        do ib = 1, Nbeta
        do is = 1, NS
            ! iz = is + NS*(ib-1)
            yseht1(is,ib) = ( hpolt1(is,ib,1)*ynw(is) + hpolt1(is,ib,2)*yw(is) )
        enddo
        enddo
        
        print *, '***********************'
        print *, '*** with shock'
        
        itshow = 0;	dois   = 1;
        ysimdy = 0.0D0; yaggsimdy = 0.0D0;
        if (dois == 1) then
        do ib0 = 1, Nbeta
        do is0 = 1, NS
        
            GMAt          = 0.0D0;
            GMAt(is0,ib0) = 1.0D0;
            
            do tt = 1, Tsim
                GMAtp1 = 0.0D0
                
                do ibt = 1, Nbeta
                do ist = 1, NS
                    ixt = xind(ist)
                    
                    if (tt == 1) then
                        
                        ysimdy(is0,ib0) = ysimdy(is0,ib0) + (yseht1(ist,ibt)*GMAt(ist,ibt))
                        yaggsimdy(tt)   = yaggsimdy(tt)   + (yseht1(ist,ibt)*GMAt(ist,ibt))*mu(is0,ib0)

                    else
                        
                        ysimdy(is0,ib0) = ysimdy(is0,ib0) + (yseh(ist,ibt)*GMAt(ist,ibt))
                        yaggsimdy(tt)   = yaggsimdy(tt)   + (yseh(ist,ibt)*GMAt(ist,ibt))*mu(is0,ib0)

                    endif
                    
                    
                    do ih = 1, Nh
                        if (tt==1) then
                            a_np  = max( min( apolt1(ist,ibt,ih) , amax), amin)
                        else
                            a_np  = min( apol(ist,ibt,ih) , amax)
                        endif
                        
                        aindv   = vsearchprm_FET([a_np],avec,avecpar)
                        aind_np = aindv(1)
                        da      = (a_np - avec(aind_np))/(avec(aind_np+1) - avec(aind_np))
                        
                        do ib_np = 1, Nbeta
                        do ix_np = 1, Nx
                            is_np = aind_np + Na*(ix_np - 1)
                            
                            if (tt==1) then
                                GMAtp1(is_np  ,ib_np) = GMAtp1(is_np  ,ib_np) + (1.0D0-da)*Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpolt1(ist,ibt,ih)*GMAt(ist,ibt)
                                GMAtp1(is_np+1,ib_np) = GMAtp1(is_np+1,ib_np) +    da     *Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpolt1(ist,ibt,ih)*GMAt(ist,ibt)
                            else
                                GMAtp1(is_np  ,ib_np) = GMAtp1(is_np  ,ib_np) + (1.0D0-da)*Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)
                                GMAtp1(is_np+1,ib_np) = GMAtp1(is_np+1,ib_np) +    da     *Px(ixt,ix_np)*Pbmat(ibt,ib_np)*hpol(ist,ibt,ih)*GMAt(ist,ibt)
                            endif
                
                            
                        enddo
                        enddo
                    
                    enddo !end loop ih
                enddo !end loop ist
                enddo !end loop ibt
                
                ! print *, 'sum( GMAtp1 ) = ', sum( GMAtp1 )
                GMAtp1 = GMAtp1/sum(GMAtp1)
                GMAt   = GMAtp1
                
            enddo !ned loop tt
            
            if (itshow<=0) then
                print *, 'is0 = ', is0,'/',NS,', ib0 = ', ib0,'/',Nbeta
                itshow = 1500
            endif
            itshow = itshow-1
        
        enddo !end loop is0
        enddo !end loop ib0
        endif
        
!        print *, 'yaggsim   = ', yaggsim
!        print *, 'yaggsimdy = ', yaggsimdy
!        
        ! print *, 'csim(1:10,:) = ', csim(1:10,:)
        !----------------------------------------------------------
        !----------------------------------------------------------
        
        !----------------------------------------------------------
        !---ya nos vamos
        
        mpe_annual = (ysimdy - ysim)/dy
        
        if (index_mpe == 1) then

            print *, 'mpe small size = ', sum(mpe_annual*mu)/dble(Tsim/4) ! Average per year over 5 years

            open(2,  file = 'OUTPUT/mpe_smallsize.txt',   status = 'unknown')
            open(5,  file = 'OUTPUT/ysim_smallsize.txt'      ,   status = 'unknown')
            do is = 1, NS
                write(2, '(*(f18.6))') ( mpe_annual(is,ib) , ib = 1,Nbeta )
                write(5, '(*(f18.6))') ( ysim(is,ib)       , ib = 1,Nbeta )
            enddo
            close(2); close(5)
            !----------------------------------------------------------
            !----------------------------------------------------------

        elseif (index_mpe == 2) then

            print *, 'mpe medium size = ', sum(mpe_annual*mu)/dble(Tsim/4) ! Average per year over 5 years

            open(2,  file = 'OUTPUT/mpe_mediumsize.txt',   status = 'unknown')
            open(5,  file = 'OUTPUT/ysim_mediumsize.txt'      ,   status = 'unknown')
            do is = 1, NS
                write(2, '(*(f18.6))') ( mpe_annual(is,ib) , ib = 1,Nbeta )
                write(5, '(*(f18.6))') ( ysim(is,ib)       , ib = 1,Nbeta )
            enddo
            close(2); close(5)

        elseif (index_mpe==3) then


            print *, 'mpe large size = ', sum(mpe_annual*mu)/dble(Tsim/4) ! Average per year over 5 years

            open(2,  file = 'OUTPUT/mpe_largesize.txt',   status = 'unknown')
            open(5,  file = 'OUTPUT/ysim_largesize.txt'      ,   status = 'unknown')
            do is = 1, NS
                write(2, '(*(f18.6))') ( mpe_annual(is,ib) , ib = 1,Nbeta )
                write(5, '(*(f18.6))') ( ysim(is,ib)       , ib = 1,Nbeta )
            enddo
            close(2); close(5)

        else
            
            print *, 'something is going wrong in MPE model'

        endif

    enddo
	
	end subroutine COMPUTEmpe

end module MODULEmpe
