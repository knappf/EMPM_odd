      TYPE level_typ

        sequence
        
	integer                 :: npr
        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j        
        integer                 :: tz
        double precision        :: en        ! s.p.-energy
      
      END TYPE level_typ

      
      TYPE phon_typ

        sequence 

        integer                 :: par        
        integer                 :: j      
        double precision        :: enf
    	integer                 :: us
        integer                 :: tz

      END TYPE phon_typ

      TYPE phonbase_typ

        sequence 

        integer                 :: isp        
        integer                 :: ila      
        integer                 :: tz
      END TYPE phonbase_typ

      TYPE amp_typ

        sequence 

        integer                 :: par        
        integer                 :: hol      
        double precision        :: am

      END TYPE amp_typ

      TYPE rho_typ

        sequence 

        integer*2                 :: ilap
        integer*2                 :: j
        integer*2                 :: i1        
        integer*2                 :: i2      
        double precision        :: ro
!        real(kind=4)              :: ro

     END TYPE rho_typ

     TYPE ro_typ

       sequence

       integer                    :: ilap
       integer*2                    :: j
       integer*2                    :: i1
       integer*2                    :: i2
!       real                       :: rho
      double precision        :: ro
    END TYPE ro_typ

     TYPE ro2_typ
       sequence
       integer                   :: ilap
       integer (kind=1)           :: j
       integer (kind=2)           :: i1
       integer (kind=2)           :: i2
       real(kind=4)        :: ro

      END TYPE ro2_typ


     











     

      
