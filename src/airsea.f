
        SUBROUTINE AIRSEA( Faw, Achem, ii, ROW, COL, HOUR, SALT, TW,
     &                     DEPTH, AREA, HIS, FRICE, VERBOSE )

        USE CPARAM
        USE CIN
        USE COUT

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'C_model_inc'

!#######################################################################################
! DESC: Air-sea exchange of elemental and dimethly mercury
!       Writes directly to Tc( :,HGDEM3D ) and ATMIN( ROW,COL,HOUR,MCHG0 )
!
! HISTORY: MERCY v2.0 - 13.10.2021
!
! AUTHOR:  johannes.bieser@hereon.de
!
! LICENSE: GNU Version 3
!
!########################################################################################

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        REAL                          , INTENT(OUT)   :: Faw     ! Flux air-water Hg0                 [ng/(m² h)]
        REAL, DIMENSION( NDREI,NCHEM ), INTENT(INOUT) :: Achem
        INTEGER                       , INTENT(IN)    :: ii
        INTEGER                       , INTENT(IN)    :: ROW !i [1,m]
        INTEGER                       , INTENT(IN)    :: COL !k [1,n]
        INTEGER                       , INTENT(IN)    :: HOUR ! hour index for ATMIN array			[1,24]

        REAL, DIMENSION( NDREI )      , INTENT(IN)    :: SALT   ! Salinity					[PSU]
        REAL, DIMENSION( NDREI )      , INTENT(IN)    :: TW     ! Water temperature				[°C]
        REAL, DIMENSION( NDREI )      , INTENT(IN)    :: DEPTH  ! Depth of upper most ocean grid cell 		[m]
        REAL, DIMENSION( NDREI )      , INTENT(IN)    :: AREA
        REAL, INTENT(IN)    :: HIS      ! height of ice shield				[m]
        REAL, INTENT(IN)    :: FRICE    ! fraction of ice cover				[0.0-1.0]

        INTEGER, INTENT(IN) :: VERBOSE  ! Triggers debug information			[0,1,2]

C...........   LOCAL VARIABLES and their descriptions:
        REAL, PARAMETER     :: CMIN   = 1.0E-08
        REAL, PARAMETER     :: H_DMHG = 0.31     !Caq = Cair / H_DMHg
        REAL, PARAMETER     :: Kdmhg  = 3.48E-06 !reaction rate [s-1]
        REAL, PARAMETER     :: R      = 8.314462 ![gas constant [J mol-1 K-1]

        INTEGER, PARAMETER :: ITS10 = ITS/10

        INTEGER NNC

!	C......... Atmospheric physical variables
        REAL    U10     ! 10m wind speed	 [m/s]
        REAL    HEIGHT  ! grid cell height (atm) [m]

!	C......... local variables used for air sea exchange
        REAL    S       ! Salinity capped at 35              [psu]
        REAL    Sc, Sc35, Sc00 ! Schmidt number              [dimensionless]
        REAL    DHg_f   ! Diffusivity of Hg in fresh water   [cm²/s]
        REAL    DHg_s   ! Diffusivity of Hg in salt water    [cm²/s]
        REAL    v       ! kinematic viscosity                [cm²/s]
        REAL    H       ! Henry's law constant 		     [dimensionless]
        REAL    k, k600 ! Transfer velocity                  [m/s]
        REAL    kw,ka   ! Mass transfer coefficient water    [m/h]
        REAL    dc      ! eqilibrium exchange rate 	     [ng/m³]
        REAL    Kol     ! Resistance parameter 		     [m/h]
        REAL    Flux    ! effective water to air flux        [ng/m²]
C...........   LOCAL PARAMETERS
        CHARACTER( 256 )    MESG                ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'AIRSEA' ! program name

C***********************************************************************
C   begin body of subroutine AIRSEA
        DO NNC = 1,NHGTRANS
          IF( Achem( ii,NNC ) .LT. 0. ) THEN
            print*, "ERROR 107: asx.f -> negative ",
     &              VNAMOUT( NNC ), Achem( ii,NNC )
            Achem( ii,NNC ) = CMIN
          END IF
        END DO

!	C. Set helper variables
        Faw = 0.
        U10   = ATMIN( ROW, COL, HOUR, MWSDP10 )                ! [m/s]

!       C. Height of lowest atmospheric grid cell
!	HEIGHT = ATMIN( ROW, COL, HOUR, MZF )
        HEIGHT = 3000.                                           ! [m]

!	C. Transform wind speed at z[m] to wind speed at 10[m]
!       U10 = ( 10.4 / ( LOG( z ) + 8.1 ) ) * Uz                ! [m/s]

!	C. Henry constant [dimensionless] (temperature dependent) (0°C H=0.15) -> (30°C H=0.36)
        H   = EXP( -2404.3 / ( TW(ii) + 273.15 ) + 6.915 )      ! [] (Andersson et al., 2008)

!       C. Schmidt number (Kuss et al., 2014; Kuss et al., 2018)
!          For Sc of O2 see Wanninkhof 2014 & Rovelli et al., 2016
!       DHg_s = 0.0011 * EXP( -11.06 / ( R * ( 273.15 + TW(ii) ) ) ) ! [cm2/s] Kuss et al., 2018
!       DHg_f = 0.0335 * EXP( -18.63 / ( R * ( 273.15 + TW(ii) ) ) ) ! [cm2/s]
        S = SALT(ii)                                            ! [psu]
        IF( S .GT. 35. ) S = 35.
        Sc35 = - 0.0398 * TW(ii)**3 + 3.3910 * TW(ii)**2
     &         - 118.02 * TW(ii)    + 1948.2                    ! [] Kuss et al., 2014
        Sc00 = - 0.0304 * TW(ii)**3 + 2.7457 * TW(ii)**2
     &         - 118.13 * TW(ii)    + 2226.2                    ! [] Kuss et al., 2014
        Sc   = ( Sc35 * S + Sc00 * ( 35. - S ) ) / 35.          ! []

!       C. Transfer velocity    ! See Benallal et al., 2017 for a review
        k600 = 0.222 * U10**2 + 0.333 * U10                     ! [cm/hr] Nightinggale et al., 2000; Kuss et al., 2009
!       k600 = 0.251 * U10**2                                   ! [cm/hr] Wanninkhof, 2014
        k    = k600 * ( Sc / 600. )**(-0.5)                     ! [cm/hr]
        k    = k / 3600. / 100.                                 ! [m/s]
        kol  = 1. / k                                           ! [s/m]

        dc    = ( ATMIN( ROW,COL,HOUR,MCHG0 ) / H )              ! [ng/m3]
     &        - Achem( ii,HGDEM3D ) * 1000

!       C. Mercury flux
        Faw   = dc / kol * ITS10                                ! [ng/m2]

        Flux  = Faw / ( DEPTH(ii) * 1000. )                     ! [ng/L]

        IF( HIS .GT. 0.1 ) THEN
            Faw = Faw * ( 1 - FRICE )                           ! [ng/m2]
        END IF

        IF( -1.0*Flux .GE. Achem( ii,HGDEM3D ) ) THEN
            print*, 'ERROR 111: asx overshoot ',
     &              Achem( ii,HGDEM3D ), Flux, ii,
     &              kol, dc, Faw, ATMIN( ROW,COL,HOUR,MCHG0 )
            print*, Faw, Flux, DEPTH( ii )
            Flux = -1.0 * ( Achem( ii,HGDEM3D ) - CMIN )
        END IF

        ATMIN( ROW, COL, HOUR, MCHG0 ) =                        ![ng/m³]
     &  ATMIN( ROW, COL, HOUR, MCHG0 ) - Faw / HEIGHT

        Achem( ii,HGDEM3D ) =                                   ![pg/m³] 
     &  Achem( ii,HGDEM3D ) + Flux

        HGOUT2D( ROW, COL, 1, HGSAX2D ) =                       ![kg]
     &  HGOUT2D( ROW, COL, 1, HGSAX2D ) + Faw * AREA( ii ) * 1.E-12

        ! DMHg exchange
        Achem( ii,DMHG3D ) = Achem( ii,DMHG3D ) * EXP( -ITS10 * Kdmhg )

C.........  End of subroutine:
        END SUBROUTINE AIRSEA
