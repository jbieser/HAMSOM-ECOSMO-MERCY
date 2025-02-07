
        SUBROUTINE PARTPART( Achem, Abio, VERBOSE, ii )

        USE CPARAM
        USE CIN
        USE COUT

        IMPLICIT NONE
C...........
!
! DESC:         Three-way particle partitioning
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........

        INCLUDE 'C_model_inc'

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        INTEGER                          , INTENT(IN)     :: ii           ! cell number
        REAL, DIMENSION( ndrei,nchem )   , INTENT(INOUT)  :: Achem        ! chemistry fields
        REAL, DIMENSION( ndrei,3:ninbio ), INTENT(IN)     :: Abio         ! biology fields
        INTEGER                          , INTENT(IN)     :: VERBOSE      ! print debug information [0,1,2]

C.........  Local fields
        REAL, PARAMETER :: CMIN = 1.0E-6

!       REAL, PARAMETER :: PSR = 0.16   ![]      !Wozniak, 2010 (2.3%-42.0%)
        REAL, PARAMETER :: PSR = 0.10   ![]      !Sharif, 2014  (9.3%-10.8%)
!       REAL, PARAMETER :: PSR = 1.00   ![]      !No SPM case (POC/SPM ratio)

!       REAL, PARAMETER :: kd0 = 199526.23 ! log(kd) = 5.30 (4.2 - 6.9)  [l/kg] (Batrakova, 2014)
        REAL, PARAMETER :: kd0 = 251188.64 ! log(kd) = 5.40
!       REAL, PARAMETER :: kl0 = 316277.77 ! log(kl) = 5.50 [l/kg] (Allison and Allison, 2005)
        REAL, PARAMETER :: kl0 = 398107.17 ! log(kl) = 5.60 (5.4:  3.0 - 6.0) [l/kg] (Allison and Allison, 2005)
        
        REAL, PARAMETER :: kd0_mehg = 79432.83  ! log(kd) = 4.9 (4.2 - 6.2) [l/kg]  for MeHg (Allison and Allison, 2005)
        REAL, PARAMETER :: kl0_mehg = 100000.0  ! log(kl) = 5.0 (2.8 - 5.5) [l/kg]  for MeHg (Allison and Allison, 2005)

!       log(ks) sediment/water Hg = 4.9 MeHg = 3.6 (Allison and Allison, 2005)

        REAL, PARAMETER :: kfast = 1.92541E-4      ! [s-1] equilibrium half life time of 1 day
        REAL, PARAMETER :: Vfast = 1. - EXP( -ITS * kfast )
        REAL, PARAMETER :: kslow = 8.01980E-6      ! [s-1] equilibrium half life time of 1 hour
        REAL, PARAMETER :: Vslow = 1. - EXP( -ITS * kslow )
        REAL, PARAMETER :: Fred  = 0.6

        REAL               kd, kl                ! partitioning coefficients
        REAL               HgW, HgD, HgP         ! fractions in dissolved, colloidal, and particulate phase
        REAL               eqP, eqD, eqT
        REAL               adP1, adP2, adD1, adD2
        REAL               reP, reD
        REAL               HgT                   ! total mercury concentration
        REAL               MeHgT                 ! total methyl mercury concnetration
        REAL               TDISS, TPART          ! total particle and dissolved organic matter
        REAL               TOT
        REAL               diff

!       Partitioning coefficients
        TDISS = Abio( ii,DOM ) + Achem( ii,DTOM3D )             ! [mgC/m³]
        Kl = Kl0 * TDISS * 10.0 * 1.0E-9 * 1.72                       ! []

        TPART = Abio( ii,DET ) + Achem( ii,PTOM3D )       ! [mgC/m³]
        Kd = kd0 * TPART * 10.0 * 1.0E-9 * 1.72                       ! []

!       inorganic mercury speciation
        HgT = Achem( ii,HGDIS3D )
     &      + Achem( ii,HGDOC3D )
     &      + Achem( ii,HGPOC3D )

!       Achem( ii,HGDIS3D ) = HgT / ( 1 + kd + kl )
!       Achem( ii,HGDOC3D ) = Achem( ii,HGDIS3D ) * kl
!       Achem( ii,HGPOC3D ) = Achem( ii,HGPOC3D ) * kd 

        eqP = HgT - HgT / (kd + 1.)
        eqD = Achem( ii,HGDIS3D ) * kl
        eqT = eqP + eqD

        IF( eqT .GT. HgT ) THEN
            eqP = eqP * HgT / eqT
            eqD = eqD * HgT / eqT
        END IF

!       instant eq.
        Achem( ii,HGPOC3D ) = eqP
        Achem( ii,HGDOC3D ) = eqD
        Achem( ii,HGDIS3D ) = HgT - eqP - eqD

!       organic mercruy speciation
        Kd = kd0_mehg * TPART * 10.0 * 1.0E-9 * 1.72                  ! []
        Kl = Kl0_mehg * TDISS * 10.0 * 1.0E-9 * 1.72                  ! []

        MeHgT = Achem( ii,MEHGDIS3D )
     &        + Achem( ii,MEHGDOC3D )
     &        + Achem( ii,MEHGPOC3D )

        eqP = MeHgT - MeHgT / (kd + 1.)
        eqD = Achem( ii,MEHGDIS3D ) * kl
        eqT = eqP + eqD

        IF( eqT .GT. MeHgT ) THEN
            eqP = eqP * MeHgT / eqT
            eqD = eqD * MeHgT / eqT
        END IF

!       instant eq.
        Achem( ii,MEHGPOC3D ) = eqP
        Achem( ii,MEHGDOC3D ) = eqD
        Achem( ii,MEHGDIS3D ) = MeHgT - eqP - eqD

!       Check mass consistency - again
        IF( Achem( ii,HGDIS3D ) .LT. 0. ) THEN
            Achem( ii,HGDIS3D ) = 0.
        END IF
        IF( Achem( ii,HGDOC3D ) .LT. 0. ) THEN
            Achem( ii,HGDOC3D ) = 0.
        END IF
        IF( Achem( ii,HGPOC3D ) .LT. 0. ) THEN
            Achem( ii,HGPOC3D ) = 0.
        END IF
        IF( Achem( ii,MEHGDIS3D ) .LT. 0. ) THEN
            Achem( ii,MEHGDIS3D ) = 0.
        END IF
        IF( Achem( ii,MEHGDOC3D ) .LT. 0. ) THEN
            Achem( ii,MEHGDOC3D ) = 0.
        END IF
        IF( Achem( ii,MEHGPOC3D ) .LT. 0. ) THEN
            Achem( ii,MEHGPOC3D ) = 0.
        END IF

        TOT = Achem( ii,HGDIS3D )
     &      + Achem( ii,HGDOC3D )
     &      + Achem( ii,HGPOC3D )

        IF( TOT .NE. HgT ) THEN
            Achem( ii,HGDIS3D ) = Achem( ii,HGDIS3D ) * HgT / TOT
            Achem( ii,HGDOC3D ) = Achem( ii,HGDOC3D ) * HgT / TOT
            Achem( ii,HGPOC3D ) = Achem( ii,HGPOC3D ) * HgT / TOT
        END IF

          TOT = Achem( ii,MEHGDIS3D )
     &        + Achem( ii,MEHGDOC3D )
     &        + Achem( ii,MEHGPOC3D )

        IF( TOT .NE. MeHgT ) THEN
            Achem( ii,MEHGDIS3D ) = Achem( ii,MEHGDIS3D ) * MeHgT / TOT
            Achem( ii,MEHGDOC3D ) = Achem( ii,MEHGDOC3D ) * MeHgT / TOT
            Achem( ii,MEHGPOC3D ) = Achem( ii,MEHGPOC3D ) * MeHgT / TOT
        END IF

        END SUBROUTINE PARTPART
