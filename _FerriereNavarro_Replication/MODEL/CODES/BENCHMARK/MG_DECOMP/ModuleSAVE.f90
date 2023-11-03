4module ModuleSAVE
    
    use GLOBALS
    
contains    
    
    !--------------------------------------------------
    !---save grids and parameters
    !--------------------------------------------------
    subroutine SaveGRIDS
    
    implicit none

    end subroutine SaveGRIDS
	
	!--------------------------------------------------
    !---save policies
    !--------------------------------------------------
    subroutine SaveVFandPolicies
	
	implicit none
	
	end subroutine SaveVFandPolicies
	
	!--------------------------------------------------
    !---save aggregates
    !--------------------------------------------------
    subroutine SaveAggregates 
	
	implicit none
	
	end subroutine SaveAggregates
	
	!--------------------------------------------------
    !---Save Transition Aggregates     
    !--------------------------------------------------
    subroutine SaveAggregates_TR
    
    implicit none

    end subroutine SaveAggregates_TR


    subroutine SaveAggregates_TR_decomp

    implicit none

    integer :: is, tt, ttax
    character(len=5) :: ics

    write (ics,'(I1)') index_decomp

    open(1, file = 'OUTPUT_'//trim(ics)//'/G_TR.txt',     status = 'unknown')
    open(2, file = 'OUTPUT_'//trim(ics)//'/Yagg_TR.txt',  status = 'unknown')
    open(3, file = 'OUTPUT_'//trim(ics)//'/Lagg_TR.txt',  status = 'unknown')
    open(4, file = 'OUTPUT_'//trim(ics)//'/Kagg_TR.txt',  status = 'unknown')

    open(5, file = 'OUTPUT_'//trim(ics)//'/Avgetax.txt',  status = 'unknown')
    open(6, file = 'OUTPUT_'//trim(ics)//'/Cagg_TR.txt',   status = 'unknown')
    open(7, file = 'OUTPUT_'//trim(ics)//'/Iagg_TR.txt',   status = 'unknown')
    open(8, file = 'OUTPUT_'//trim(ics)//'/theta_TR.txt',  status = 'unknown')
    open(9, file = 'OUTPUT_'//trim(ics)//'/thetaw_TR.txt', status = 'unknown')

    open(10, file = 'OUTPUT_'//trim(ics)//'/lbd_TR.txt',   status = 'unknown')
    open(11, file = 'OUTPUT_'//trim(ics)//'/gma_TR.txt',   status = 'unknown')
    open(12, file = 'OUTPUT_'//trim(ics)//'/r_TR.txt',     status = 'unknown')
    open(13, file = 'OUTPUT_'//trim(ics)//'/rk_TR.txt',    status = 'unknown')
    open(14, file = 'OUTPUT_'//trim(ics)//'/qk_TR.txt',    status = 'unknown')
    open(15, file = 'OUTPUT_'//trim(ics)//'/wgeH_TR.txt',  status = 'unknown')
    open(16, file = 'OUTPUT_'//trim(ics)//'/wge_TR.txt',   status = 'unknown')
    open(17, file = 'OUTPUT_'//trim(ics)//'/div_TR.txt',   status = 'unknown')

    open(20, file = 'OUTPUT_'//trim(ics)//'/div_nobank_TR.txt',   status = 'unknown')
    open(21, file = 'OUTPUT_'//trim(ics)//'/Pi_TR.txt',   status = 'unknown')

    do tt = 1, T_TR
         write(1 , '(*(f18.6))') ( G_TR(tt)   )
         write(2 , '(*(f18.6))') ( Yagg_TR(tt) )
         write(3 , '(*(f18.6))') ( Lagg_TR(tt) )
         write(4 , '(*(f18.6))') ( Kagg_TR(tt) )

         write(5 , '(*(f18.6))') ( Avgetax(tt) )
         write(6 , '(*(f18.6))') ( Cagg_TR(tt) )
         write(7 , '(*(f18.6))') ( Iagg_TR(tt) )
         write(8 , '(*(f18.6))') ( theta_TR(tt) )
         write(9 , '(*(f18.6))') ( thetaw_TR(tt) )

         write(10,'(*(f18.6))')  ( lbd_TR(tt) )
         write(11,'(*(f18.6))')  ( gma_TR(tt) )
         write(12,'(*(f18.6))')  ( r_TR(tt) )
         write(13,'(*(f18.6))')  ( rk_TR(tt) )
         write(14,'(*(f18.6))')  ( qk_TR(tt) )
         write(15,'(*(f18.6))')  ( wgeH_TR(tt) )
         write(16,'(*(f18.6))')  ( wge_TR(tt) )
         write(17,'(*(f18.6))')  ( div_TR(tt) )

    !     write(18, '(*(f18.6))') ( Avgetax_LPE(tt) )
    !     write(19, '(*(f18.6))') ( Avgetax_LPE_xinc(tt))

        write(20,'(*(f18.6))')  ( div_TR(tt) - divBank_TR(tt) )
        write(21,'(*(f18.6))')  ( Pi_TR(tt) )

    enddo
    close(1); close(2); close(3) ; close(4) ; close(5) ; close(6) ; close(7);
    close(8); close(9); close(10); close(11); close(12); close(13); close(14);
    close(15); close(16); close(17); ! close(18); close(19);
    close(20); close(21)


    open(1,  file = 'OUTPUT_'//trim(ics)//'/DB_TR.txt',   status = 'unknown')
    open(4,  file = 'OUTPUT_'//trim(ics)//'/Aagg_TR.txt', status = 'unknown')
    do tt = 1, T_TR+1
         write(1, '(*(f18.6))') ( DB_TR(tt)   )
         write(4, '(*(f18.6))') ( Aagg_TR(tt)   )
    enddo
    close(1); close(2); close(3); close(4);

    end subroutine SaveAggregates_TR_decomp


    !--------------------------------------------------
    !---save policies for transitin
    !--------------------------------------------------
    subroutine SaveVFandPolicies_TR

    implicit none
    integer :: tt, is, ib, ih

    character(len=5) :: ics

    write (ics,'(I1)') index_decomp


    open(1, file = 'OUTPUT_'//trim(ics)//'/Vwrk_TR.txt',    status = 'unknown')
    open(2, file = 'OUTPUT_'//trim(ics)//'/Vnow_TR.txt',    status = 'unknown')
    open(3, file = 'OUTPUT_'//trim(ics)//'/mu_TR.txt',    status = 'unknown')
    do tt = 1, T_TR
    do ib = 1, Nbeta
    do is = 1, NS
        write(1, '(*(f18.6))') ( Vwrk_TR(is,ib,tt) )
        write(2, '(*(f18.6))') ( Vwrk_TR(is,ib,tt) )
        write(3, '(*(f18.6))') ( mu_TR(is,ib,tt)   )
    enddo
    enddo
    enddo
    close(1); close(2); close(3)

    open(1,  file = 'OUTPUT_'//trim(ics)//'/apol_TR.txt',   status = 'unknown')
    open(2,  file = 'OUTPUT_'//trim(ics)//'/cpol_TR.txt',   status = 'unknown')
    open(3,  file = 'OUTPUT_'//trim(ics)//'/hpol_TR.txt',   status = 'unknown')
    do ih = 1, Nh
    do tt = 1, T_TR
    do ib = 1, Nbeta
    do is = 1, NS
        write(1, '(*(f18.6))') ( apol_TR(is,ib,tt,ih)  )
        write(2, '(*(f18.6))') ( cpol_TR(is,ib,tt,ih)  )
        write(3, '(*(f18.6))') ( hpol_TR(is,ib,tt,ih)  )
    enddo
    enddo
    enddo
    enddo
    close(1); close(2); close(3);

    end subroutine SaveVFandPolicies_TR

    end module ModuleSAVE
