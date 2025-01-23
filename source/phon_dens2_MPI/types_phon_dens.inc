      TYPE level_typ

        sequence

        integer                 :: npr
        integer                 :: n
        integer                 :: l
        integer                 :: j         ! 2*j
        integer                 :: tz
        double precision        :: en        ! s.p.-energy
       

      END TYPE level_typ


      TYPE amp2_typ

        sequence
      
        integer                 :: is
        integer                 :: ig     
        double precision        :: am
!        real                    :: am

      END TYPE amp2_typ

      TYPE ro_typ

       sequence

       integer                      :: ib
       integer*2                    :: is
       integer*2                    :: i1
       integer*2                    :: i2
!       real                        :: rho
      double precision              :: rho

      END TYPE ro_typ

      TYPE roc_typ

       sequence

       integer                    :: ib
       integer                    :: is
       integer                    :: i1
       integer                    :: i2
       real                       :: rho
!      double precision        :: rho

      END TYPE roc_typ

      
      TYPE phon_typ

        sequence 

        integer                 :: par        
        integer                 :: jj      
        integer                 :: tz 
        double precision   :: enf
!	integer*8 :: us
!	integer*8 :: ihom

      END TYPE phon_typ

      TYPE rho_typ

      sequence

      integer*2                 :: ilap
      integer*2                 :: j
      integer*2                 :: i1
      integer*2                 :: i2
      double precision        :: ro
!      real                      :: ro

      END TYPE rho_typ

      TYPE ro2_typ
       sequence
       integer                   :: ib
       integer (kind=1)           :: is
       integer (kind=2)           :: i1
       integer (kind=2)           :: i2
       real(kind=4)        :: rho

      END TYPE ro2_typ

!      TYPE ro2_real_typ
       
!       sequence
!       integer                   :: ib
!       integer (kind=1)           :: is
!       integer (kind=2)           :: i1
!       integer (kind=2)           :: i2
!       real(kind=4)        :: rho
       
!      END TYPE ro2_real_typ






     

      
