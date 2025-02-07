        SUBROUTINE SEDIMENT(ii, ii_i, ii_k, ii_j, Achem, Abio, VERBOSE)

	USE CPARAM
	USE CDATA
	USE COUT
        USE CIN

	IMPLICIT NONE

C...........
!
! DESC:         Sedimentation, resuspention and burial
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........

C...........   INCLUDES:
        INCLUDE 'C_model_inc'

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        INTEGER                          , INTENT(IN)     :: ii         ! cell number           []
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_i       ! m = row               [1-M]
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_k       ! n = col               [1-N]
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_j       ! layer index           [1-ILO] (ilo<=19)
        REAL, DIMENSION( ndrei,nchem )   , INTENT(INOUT)  :: Achem      ! chemistry array       [ng/L]
        REAL, DIMENSION( ndrei,3:ninbio ), INTENT(IN)     :: Abio       ! biology array         [mgC/L]
        INTEGER                          , INTENT(IN)     :: VERBOSE    ! print debug information [0,1,2]

C...........   PARAMETERS and their descriptions:
        LOGICAL, PARAMETER  :: ONOFF = .FALSE.

C...... Hg2+ --> Hg0
        REAL  Vrd
        REAL, PARAMETER :: krd = 6.00E-7        ![s-1]        !Kuss et al., 2015

C...... mono-methylation
C...... km1 : Hg2+ --> MMHg
        REAL  Vm1
        REAL, PARAMETER :: km1 = 4.40E-07       ![s-1]        !Monperrus, 2007 max

C...... demethylation
C...... kdm : MMHg --> Hg2+
        REAL  Vdm
        REAL, PARAMETER :: kdm = 6.94E-07       ![s-1]        !Monperrus, 2007

C...... (R20) reductive de-methylation
C...... krm : MMHg --> Hg0
        REAL  Vrm
        REAL, PARAMETER :: krm = 2.22E-09       ![s-1]        !


        REAL, PARAMETER     :: Ksed = 5.7870E-05        ! Sedimentation rate (5 meter per day)  [m/s]  (was 3.5 d-1)
        REAL, PARAMETER     :: Kres = 2.8935E-04        ! Resuspention rate  (25 per day)	[s-1]
        REAL, PARAMETER     :: Kbur = 1.1575E-10        ! Burial rate	     (0.00001 per day)  [s-1]
        REAL, PARAMETER     :: UCRIT = 0.01             ! Threshold for resuspension		[m/s]
!       REAL, PARAMETER     :: DPS = 1./( 24. * 3600. ) ! Days per second 			[d/s]

C...........   LOCAL VARIABLES and their descriptions:
        REAL    USTAR   ! Critical velocity	[m/s]
        REAL    R_SED, R_RES, R_BUR
        REAL    F2_SED_HG, F2_SED_ME, F3_SED_HG, F3_SED_ME
        REAL    F2_RES_HG, F2_RES_ME, F3_RES_HG, F3_RES_ME
        REAL    F2_BUR_HG, F2_BUR_ME
        REAL    RED, MET, DEM, REM
        REAL    lz

!	C. Helper variables
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'SEDIMENT' ! program name

C***********************************************************************
C   begin body of subroutine SEDIMENT

!	Switch this module on/off
!      IF( ONOFF ) THEN

!       Might be a problem not sure yet
        lz = MAX( 0., Tdept( ii ) )
        IF( lz .LT. 1 ) THEN
            print*, 'WARNING: Bottom depth ', lz
            lz = 1.
        END IF

!       1. Critical velocity
        USTAR = SQRT( Tvn(ii_i(ii),ii_k(ii),ii_j(ii))**2
     &              + Tun(ii_i(ii),ii_k(ii),ii_j(ii))**2 )

!	2. Sedimentation
        IF( USTAR .LT. UCRIT ) THEN
          R_SED = MIN( Ksed * ITS, lz ) / lz                            ! [] relative settling [0:1]
        ELSE
          R_SED = 0.
        END IF
        
        F3_SED_HG = R_SED * Achem( ii,HGPOC3D )                         ! [ng/dm^3] 3d sedimentation flux
        F3_SED_ME = R_SED * Achem( ii,MEHGPOC3D )
        F2_SED_HG = 1000. * F3_SED_HG * lz                              ! [ng/m^2] 2d sedimentation flux
        F2_SED_ME = 1000. * F3_SED_ME * lz

!	3. Resuspension
        IF( USTAR .GE. UCRIT ) THEN
          R_RES = Kres * ITS                                            ! [] relative resuspention [0:1]
        ELSE
          R_RES = 0.
        END IF

        F2_RES_HG = R_RES * Ts( ii_i(ii),ii_k(ii),1,HGSED2D )           ! [ng/m^2] 2d resuspension flux
        F2_RES_ME = R_RES * Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D )         
        F3_RES_HG = F2_RES_HG / ( 1000. * lz )                          ! [ng/dm^3] 3d resuspension flux
        F3_RES_ME = F2_RES_ME / ( 1000. * lz )

!	4. Burial
        IF( USTAR .LT. UCRIT ) THEN
          R_BUR = Kbur * ITS                                            ! [] = relative burial
        ELSE
          R_BUR = 0.
        END IF

        F2_BUR_HG = R_BUR * Ts( ii_i(ii),ii_k(ii),1,HGSED2D )           ! [ng/m^2] 2d burial flux
        F2_BUR_ME = R_BUR * Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D ) 

!	5. Apply changes to data fields
        Achem( ii,HGPOC3D ) =
     &  Achem( ii,HGPOC3D ) - F3_SED_HG + F3_RES_HG

        Achem( ii,MEHGPOC3D ) =
     &  Achem( ii,MEHGPOC3D ) - F3_SED_ME + F3_RES_ME 

        Ts( ii_i(ii),ii_k(ii),1,HGSED2D ) =
     &  Ts( ii_i(ii),ii_k(ii),1,HGSED2D ) + F2_SED_HG
     &                                    - F2_RES_HG
     &                                    - F2_BUR_HG

        Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D ) = 
     &  Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D ) + F2_SED_ME
     &                                      - F2_RES_ME
     &                                      - F2_BUR_ME

        Ts( ii_i(ii),ii_k(ii),1,HGBUR2D ) =
     &  Ts( ii_i(ii),ii_k(ii),1,HGBUR2D ) + F2_BUR_HG

        Ts( ii_i(ii),ii_k(ii),1,MEHGBUR2D ) =
     &  Ts( ii_i(ii),ii_k(ii),1,MEHGBUR2D ) + F2_BUR_ME

        Ts( ii_i(ii),ii_k(ii),1,HGRES2D ) =
     &  Ts( ii_i(ii),ii_k(ii),1,HGRES2D ) + F3_RES_HG - F3_SED_HG

        Ts( ii_i(ii),ii_k(ii),1,MEHGRES2D ) =
     &  Ts( ii_i(ii),ii_k(ii),1,MEHGRES2D ) + F3_RES_ME - F3_SED_ME

!	6. Sediment chemistry
        Vrd = 1. - EXP( -ITS * krd )
        RED = Vrd * Ts( ii_i(ii),ii_k(ii),1,HGSED2D )
        Vm1 = 1. - EXP( -ITS * Km1 )
        MET = Vm1 * Ts( ii_i(ii),ii_k(ii),1,HGSED2D )
        Vdm = 1. - EXP( -ITS * kdm )
        DEM = Vdm * Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D )
        Vrm = 1. - EXP( -ITS * krm )
        REM = Vrm * Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D )

        Ts( ii_i(ii),ii_k(ii),1,HGSED2D )   =
     &  Ts( ii_i(ii),ii_k(ii),1,HGSED2D )   - MET + DEM - RED
        Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D ) =
     &  Ts( ii_i(ii),ii_k(ii),1,MEHGSED2D ) + MET - DEM - REM
        Achem( ii,HGDEM3D ) =
     &  Achem( ii,HGDEM3D ) + ( RED + REM ) / ( 1000. * lz )

!	8. Exchange of dissolved species
!	TODO: HGDOC interaction with sediment

        END SUBROUTINE SEDIMENT

