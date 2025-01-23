
      TYPE phon_typ

        sequence 

        integer                 :: par        
        integer                 :: jj      
        integer                 :: tz 
        double precision   :: enf
!	integer*8 :: us
!	integer*8 :: ihom

      END TYPE phon_typ


      TYPE amp_typ

        sequence

        integer                 :: isp
        integer                 :: ila
        double precision        :: am

      END TYPE amp_typ

      TYPE level_typ

        sequence

        integer                 :: npr
        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j
        integer                 :: nr
        integer                 :: tz
        double precision        :: en        ! s.p.-energy

      END TYPE level_typ

      TYPE amp_ph_typ

        sequence

        integer                 :: par
        integer                 :: hol
        double precision        :: am

      END TYPE amp_ph_typ







