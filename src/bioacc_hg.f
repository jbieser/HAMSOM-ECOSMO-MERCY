
        SUBROUTINE BIOACC_HG( Achem, Afish, Abio, Aflu, TEMP, SALT,
     &             VERBOSE, ii, ii_i, ii_k, ii_j, ii_j_max, ii_ldep )

        USE CPARAM
        USE CIN
        USE COUT
        USE CDATA

        IMPLICIT NONE

C...........
!
! DESC:! Bioaccumulation of inorganic Hg
!
! HISTORY: MERCY v2.0 - 13.10.2021
!
! AUTHOR:  johannes.bieser@hereon.de
!
! LICENSE: GNU Version 3
!
C...........

C...........   INCLUDES:

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INCLUDE 'C_model_inc'

C...........   PARAMETERS and their descriptions:
        REAL LOL, ROFL
        REAL, PARAMETER :: SEC_PER_DAY = 3600 * 24                      ! Seconds per day                 [s]
        REAL, PARAMETER :: LIM         = 0.99                           ! flux limiter                    []
        REAL, PARAMETER :: FRR         = 0.40                           ! Reminerilazation fraction       []

        REAL, PARAMETER :: CMIN        = 1.0E-10        ! minimum concentration                           [ng/m**3]
        REAL, PARAMETER :: FMAX        = 1.0            ! Maximum feeding fraction >=0.0  <=1.0           []

        REAL, PARAMETER :: pHgCl2   = 7.2E-06   ! Permeability of CH3HgCl2  [m s-1]
        REAL, PARAMETER :: pHgOHCl  = 6.0E-06   ! Permeability of CH3HgOHCl [m s-1]
        REAL, PARAMETER :: pHgOH2   = 1.1E-06   ! Permeability of CH3HgOH2  [m s-1]
        REAL, PARAMETER :: pCH3HgCl = 7.4E-06   ! Permeability of CH3HgCl   [m s-1]
        REAL, PARAMETER :: pCH3HgOH = .29E-06   ! Permeability of CH3HgOH   [m s-1]

        REAL, PARAMETER :: A_dia = 5.3375728E-03 ! Surface area [m**2 mgC-1]
        REAL, PARAMETER :: A_fla = 2.5228876E-03 ! Surface area [m**2 mgC-1]
        REAL, PARAMETER :: A_cya = 6.0300368E-03 ! Surface area [m**2 mgC-1]

        REAL, PARAMETER :: Ru_dia = 3.95E-08
        REAL, PARAMETER :: Ru_fla = 1.87E-08
        REAL, PARAMETER :: Ru_cya = 4.46E-08

        REAL, PARAMETER :: Ru_zos = 1.94E-10 !3.0E-10!1.94E-10 ! Surface area [m**3 mgC-1 s-1]
        REAL, PARAMETER :: Ru_zol = 1.94E-10 !3.0E-10!1.94E-10 ! Surface area [m**3 mgC-1 s-1]
        REAL, PARAMETER :: Ru_mac = 1.90E-10 ! Surface area [m**3 mgC-1 s-1]
        REAL, PARAMETER :: Ru_fsh = 3.88E-12 ! Surface area [m**3 mgC-1 s-1]

        REAL, PARAMETER :: Rr_dia = 6.58E-04 ! release rate [s-1] = 6.6 [day-1]
        REAL, PARAMETER :: Rr_fla = 3.11E-04 ! release rate [s-1] = 4.0 [day-1]
        REAL, PARAMETER :: Rr_cya = 7.43E-04 ! release rate [s-1] = 15.0 [day-1]

        REAL, PARAMETER :: Rre_zos =  6.94E-06  ! release rate [s-1] = 0.05   [day-1]
        REAL, PARAMETER :: Rre_zol =  6.94E-06  ! release rate [s-1] = 0.05   [day-1]
        REAL, PARAMETER :: Rre_mac =  6.94E-06  ! release rate [s-1] = 0.05   [day-1]
        REAL, PARAMETER :: Rre_fsh =  3.47E-07  ! release rate [s-1] = 0.0287 [day-1]

        REAL, PARAMETER :: Rri_zos =  5.78704E-07!*0.5  ! release rate [s-1] = 0.05   [day-1]
        REAL, PARAMETER :: Rri_zol =  5.78704E-07!*0.5  ! release rate [s-1] = 0.05   [day-1]
        REAL, PARAMETER :: Rri_mac =  5.78704E-07!*0.5  ! release rate [s-1] = 0.05   [day-1]
        REAL, PARAMETER :: Rri_fsh =  6.44676E-07!*0.5  ! release rate [s-1] = 0.0557 [day-1]

        REAL, PARAMETER :: E_det =  0.45    ! transfer efficiency [] ! [4,5,9]
        REAL, PARAMETER :: E_fla =  0.45    ! transfer efficiency [] ! [1]
        REAL, PARAMETER :: E_dia =  0.45    ! transfer efficiency [] ! [2]
        REAL, PARAMETER :: E_cya =  0.45    ! transfer efficiency [] ! [3]

        REAL, PARAMETER :: Ei_zos =  0.45   ! transfer efficiency [] ! [6]
        REAL, PARAMETER :: Ei_zol =  0.45   ! transfer efficiency [] ! [7]
        REAL, PARAMETER :: Ei_mac =  0.45   ! transfer efficiency [] ! [8]
        REAL, PARAMETER :: Ei_fsh =  0.158  ! transfer efficiency [] 

        REAL, PARAMETER :: Ee_zos =  0.45   ! transfer efficiency []
        REAL, PARAMETER :: Ee_zol =  0.45   ! transfer efficiency []
        REAL, PARAMETER :: Ee_mac =  0.45   ! transfer efficiency []
        REAL, PARAMETER :: Ee_fsh =  0.158  ! transfer efficiency []

        REAL :: FRAC_CL = 0.8               ! fraction of chlorinated mercury       []
        REAL :: FRAC_OH = 0.2               ! fraction of hydroxy mercury species   []
        REAL DEP

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        REAL, DIMENSION( ndrei,nchem )   , INTENT(INOUT)  :: Achem
        REAL, DIMENSION( ndrei,2 )       , INTENT(IN)     :: Afish
        REAL, DIMENSION( ndrei,3:ninbio ), INTENT(IN)     :: Abio
        REAL, DIMENSION( ndrei,nflu )    , INTENT(IN)     :: Aflu

        REAL, DIMENSION( ndrei )         , INTENT(IN)    :: TEMP     ! Water temperature        [°C]
        REAL, DIMENSION( ndrei )         , INTENT(IN)    :: SALT     ! Salinity                 [PSU]

        INTEGER, INTENT(IN) :: VERBOSE  ! print debug information				[0,1,2]

        INTEGER                          , INTENT(IN)     :: ii      ! cell number
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_i    ! m = row               [1-M]
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_k    ! n = col               [1-N]
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_j    ! layer index           [1-ILO] (ilo<=19)
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_j_max! number of layers in column
        INTEGER, DIMENSION( ndrei )      , INTENT(IN)     :: ii_ldep ! depth of lowest layer [m]

        INTEGER         NNC
        REAL            HGTOT, FNORM

        REAL vP
        REAL UP_FLA, UP_DIA, UP_CYA, UP_ZOS, UP_ZOL, UP_MAC, UP_FSH     !uptake rates
        REAL  R_FLA,  R_DIA,  R_CYA, RE_ZOS, RE_ZOL, RE_MAC, RE_FSH     !release rates
        REAL                         RI_ZOS, RI_ZOL, RI_MAC, RI_FSH     !release_rates
        REAL FS1, FS2, FS3, FS4, FS5                                    ! feeding rate of ZOS
        REAL FL1, FL2, FL3, FL4, FL5, FLe6, FLi6                        ! feeding rate of ZOL
        REAL                FI4, FI5, FIe6, FIi6, FIe7, FIi7, FIe8, FIi8! feeding rate of FSH
        REAL MI1, MI2,      MI4, MI5, MIe6, MIi6, MIe7, MIi7,       MI9 ! feeding rate of MAC
        REAL AS1, AS2, AS3, AS4, AS5                                    ! Feeding fraction accumulated
        REAL AL1, AL2, AL3, AL4, AL5, ALe6, ALi6            
        REAL                AI4, AI5, AIe6, AIi6, AIe7, AIi7, AIe8, AIi8
        REAL AM1, AM2,      AM4, AM5, AMe6, AMi6, AMe7, AMi7,       AM9
        REAL ES1, ES2, ES3, ES4, ES5                                    ! Feeding fraction excreted
        REAL EL1, EL2, EL3, EL4, EL5, ELe6, ELi6, ELe7, ELi7
        REAL                EI4, EI5, EIe6, EIi6, EIe7, EIi7, EIe8, EIi8
        REAL EM1, EM2,      EM4, EM5, EMe6, EMi6, EMe7, EMi7,       EM9
        REAL            MFLA, MDIA, MCYA, MZOSe, MZOLe, MFSHe, MMACe        ! Mortality rate
        REAL                              MZOSi, MZOLi, MFSHi, MMACi
        REAL            XFLA, XDIA, XCYA, XZOS, XZOL, XFSH, XMAC        ! mortality fraction to DET
        REAL            YFLA, YDIA, YCYA, YZOS, YZOL, YFSH, YMAC        ! mortality fraction to DOC
        REAL            FREM                                            ! Remineralization
        REAL            PHY_TOT, FLAFRAC, DIAFRAC, CYAFRAC
        REAL            ZOO_TOT, ZOSFRAC, ZOLFRAC

        REAL :: RDET = 0.
        REAL :: RDOC = 0.
        REAL SEDFRAC
        REAL HGMACe, HGMACi, HGSED                                      ! macrobenthos Hg2+ concentration in [ng/dm^3]

        REAL            Vm1, Vdm, MET
        REAL, PARAMETER :: km1 = 2.21E-07      ![s-1] !Monperrus, 2007 anoxic methylation
        REAL, PARAMETER :: km2 = 4.63E-10      ![s-1] !Lehnherr, oxic mehtylation
        REAL, PARAMETER :: kdm = 6.94E-07      ![s-1] !Monperrus, 2007

        REAL            Vbur, HGBUR
        REAL, PARAMETER :: Kbur = 1.1575E-10        ! Burial rate (0.00001 per day)  [s-1]

        INTEGER         I

C***********************************************************************
C   begin body of subroutine BIOACC
        DEP = MAX( 5,ii_ldep(ii) )
        if( ii_j(ii) .EQ. 1 ) DEP = 5.

!       total mercury
        HGTOT = 0.

        DO NNC = 1,NHGTRANS
          IF( Achem( ii,NNC ) .LT. 0. ) THEN
            print*, "ERROR 949: bioacc.f -> neative ",
     &              VNAMOUT( NNC ), Achem( ii,NNC )
            Achem( ii,NNC ) = CMIN
          END IF
          HGTOT = HGTOT + Achem( ii,NNC )
        END DO
       IF( Achem( ii,HGMACe3D ) .LT. 0. ) print*, 'MACe ',
     &          Achem( ii,HGMACe3D )
       IF( Achem( ii,HGMACe3D ) .LT. 0. ) print*, 'MACi ',
     &          Achem( ii,HGMACi3D )
!       Achem( ii,HGMACe3D ) = MAX( 0.,Achem( ii,HGMACe3D ) )
!       Achem( ii,HGMACi3D ) = MAX( 0.,Achem( ii,HGMACi3D ) )


        HGMACe = 0.
        HGMACi = 0.
        HGSED  = 0.

        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
            HGMACe = Achem( ii,HGMACe3D ) / 1000. / DEP ! [ng/dm^3} = [ng/m^2] / [dm^3/m^3] / [m]
            HGMACi = Achem( ii,HGMACi3D ) / 1000. / DEP
            HGSED  = Ts( ii_i(ii),ii_k(ii),1,HGSED2D ) / 1000. / DEP

        ! Set diagnostic 2d bottom output fields
            Ts( ii_i(ii),ii_k(ii),1,HGMACe2D ) = Achem( ii,HGMACe3D )
            Ts( ii_i(ii),ii_k(ii),1,HGMACi2D ) = Achem( ii,HGMACi3D )
            Ts( ii_i(ii),ii_k(ii),1,MAC2D   ) = Abio( ii,MAC )
            Ts( ii_i(ii),ii_k(ii),1,STOT2D  ) = Abio( ii,STOT )
        END IF

        HGTOT = HGTOT + Achem( ii,HGFISHi3D ) + Achem( ii,HGFISHe3D )
        HGTOT = HGTOT + HGMACe + HGMACi
        HGTOT = HGTOT + HGSED

        FRAC_CL = 0.8
        FRAC_OH = 0.2
        IF( SALT( ii ) .LT. 20. ) THEN
            FRAC_OH = 0.2 + 1. - EXP( -( 20. - SALT( ii ) ) / 14. )
            FRAC_CL = 1.0 - FRAC_OH
        END IF
        IF( SALT( ii ) .LE. 1E-6 ) THEN
            FRAC_OH = 0.96
            FRAC_CL = 0.04
        END IF
        IF( FRAC_OH .GT. 0.96 ) THEN
            FRAC_OH = 0.96
            FRAC_CL = 0.04
        ELSE IF( FRAC_OH .LT. 0.2 ) THEN
            FRAC_OH = 0.2
            FRAC_CL = 0.8
        END IF        

!       passive uptake
        vP = FRAC_CL * pHgCl2 + FRAC_OH * pHgOH2

!       UP_FLA = vP * A_fla * Abio( ii,FLA ) * Achem( ii,HGDIS3D ) * ITS
!       UP_DIA = vP * A_dia * Abio( ii,DIA ) * Achem( ii,HGDIS3D ) * ITS
!       UP_CYA = vP * A_cya * Abio( ii,CYA ) * Achem( ii,HGDIS3D ) * ITS
        UP_DIA = Ru_dia * Abio( ii,DIA ) * Achem( ii,HGDIS3D ) * ITS
        UP_FLA = Ru_fla * Abio( ii,FLA ) * Achem( ii,HGDIS3D ) * ITS
        UP_CYA = Ru_cya * Abio( ii,CYA ) * Achem( ii,HGDIS3D ) * ITS
        UP_ZOS = Ru_zos * Abio( ii,ZOS ) * Achem( ii,HGDIS3D ) * ITS
        UP_ZOL = Ru_zol * Abio( ii,ZOL ) * Achem( ii,HGDIS3D ) * ITS
        UP_MAC = Ru_mac * Abio( ii,MAC ) * Achem( ii,HGDIS3D ) * ITS
        UP_FSH = Ru_fsh * Abio( ii,FSH ) * Achem( ii,HGDIS3D ) * ITS

!       release
        R_FLA = Rr_fla * Achem( ii,HGFLA3D ) * ITS
        R_DIA = Rr_dia * Achem( ii,HGDIA3D ) * ITS
        R_CYA = Rr_cya * Achem( ii,HGCYA3D ) * ITS

!       external release
        RE_ZOS = Rre_zos * Achem( ii,HGZOSe3D ) * ITS
        RE_ZOL = Rre_zol * Achem( ii,HGZOLe3D ) * ITS
        RE_MAC = Rre_mac * HGMACe * ITS
        RE_FSH = Rre_fsh * Achem( ii,HGFISHe3D ) * ITS

!       internal release
        RI_ZOS = Rri_zos * Achem( ii,HGZOSi3D ) * ITS
        RI_ZOL = Rri_zol * Achem( ii,HGZOLi3D ) * ITS
        RI_MAC = Rri_mac * HGMACi * ITS
        RI_FSH = Rri_fsh * Achem( ii,HGFISHi3D ) * ITS

        RE_MAC = MAX( 0.,RE_MAC )
        RE_MAC = MIN( 0.5*HGMACi,RE_MAC )
        RI_MAC = MAX( 0.,RI_MAC )
        RI_MAC = MIN( 0.5*HGMACi,RI_MAC )

C........ Feeding
        IF( Abio( ii,FLA ) .LE. 1E-6 ) THEN
          FS1 = 0.
        ELSE
          FS1 = Achem( ii,HGFLA3D ) / Abio( ii,FLA )
     &        * Aflu( ii,ZOS_ON_FLA ) / SEC_PER_DAY * ITS
        END IF
        FS1 = MIN( LIM * Achem( ii,HGFLA3D ), FS1 )

        IF( Abio( ii,DIA ) .LE. 1E-6 ) THEN
          FS2 = 0.
        ELSE
          FS2 = Achem( ii,HGDIA3D ) / Abio( ii,DIA ) 
     &        * Aflu( ii,ZOS_ON_DIA ) / SEC_PER_DAY * ITS
        END IF
        FS2 = MIN( LIM * Achem( ii,HGDIA3D ), FS2 )

        IF( Abio( ii,CYA ) .LE. 1E-6 ) THEN
          FS3 = 0.
        ELSE
          FS3 = Achem( ii,HGCYA3D ) / Abio( ii,CYA )
     &        * Aflu( ii,ZOS_ON_CYA ) / SEC_PER_DAY * ITS
        END IF
        FS3 = MIN( LIM * Achem( ii,HGCYA3D ), FS3 )

        IF( Abio( ii,DET ) .LE. 1E-6 ) THEN
          FS4 = 0.
          FS5 = 0.
        ELSE
!         FS4 = Achem( ii,HGDET3D ) / Abio( ii,DET )
!    &        * Aflu( ii,ZOS_ON_DET ) / SEC_PER_DAY * ITS
          FS5 = Achem( ii,HGPOC3D ) / Abio( ii,DET )
     &        * Aflu( ii,ZOS_ON_DET ) / SEC_PER_DAY * ITS
        END IF
!       FS4 = MIN( LIM * Achem( ii,HGDET3D ), FS4 )
        FS4 = 0.
        FS5 = MIN( LIM * Achem( ii,HGPOC3D ), FS5 )

        IF( Abio( ii,FLA ) .LE. 1E-6 ) THEN
          FL1 = 0.
        ELSE
          FL1 = Achem( ii,HGFLA3D ) / Abio( ii,FLA )
     &        * Aflu( ii,ZOL_ON_FLA ) / SEC_PER_DAY * ITS
        END IF
        FL1 = MIN( LIM * Achem( ii,HGFLA3D ), FL1 )

        IF( Abio( ii,DIA ) .LE. 1E-6 ) THEN
          FL2 = 0.
        ELSE
          FL2 = Achem( ii,HGDIA3D ) / Abio( ii,DIA )
     &        * Aflu( ii,ZOL_ON_DIA ) / SEC_PER_DAY * ITS
        END IF
        FL2 = MIN( LIM * Achem( ii,HGDIA3D ), FL2 )

        IF( Abio( ii,CYA ) .LE. 1E-6 ) THEN
          FL3 = 0.
        ELSE
          FL3 = Achem( ii,HGCYA3D ) / Abio( ii,CYA )
     &        * Aflu( ii,ZOL_ON_CYA ) / SEC_PER_DAY * ITS
        END IF
        FL3 = MIN( LIM * Achem( ii,HGCYA3D ), FL3 )

        IF( Abio( ii,DET ) .LE. 1E-6 ) THEN
          FL4 = 0.
          FL5 = 0.
        ELSE
!         FL4 = Achem( ii,HGDET3D ) / Abio( ii,DET ) 
!    &        * Aflu( ii,ZOL_ON_DET ) / SEC_PER_DAY * ITS
          FL5 = Achem( ii,HGPOC3D ) / Abio( ii,DET )
     &        * Aflu( ii,ZOL_ON_DET ) / SEC_PER_DAY * ITS
        END IF
!       FL4 = MIN( LIM * Achem( ii,HGDET3D ), FL4 )
        FL4 = 0.
        FL5 = MIN( LIM * Achem( ii,HGPOC3D ), FL5 )

        IF( Abio( ii,ZOS ) .LE. 1E-6 ) THEN
          FLe6 = 0.
          FLi6 = 0.
        ELSE
          FLe6 = Achem( ii,HGZOSe3D ) / Abio( ii,ZOS )
     &        * Aflu( ii,ZOL_ON_ZOS ) / SEC_PER_DAY * ITS
          FLi6 = Achem( ii,HGZOSi3D ) / Abio( ii,ZOS )
     &        * Aflu( ii,ZOL_ON_ZOS ) / SEC_PER_DAY * ITS
        END IF
        FLe6 = MIN( LIM * Achem( ii,HGZOSe3D ), FLe6 )
        FLi6 = MIN( LIM * Achem( ii,HGZOSi3D ), FLi6 )

        IF( Abio( ii,DET ) .LE. 1E-6 ) THEN
          FI4 = 0.
          FI5 = 0.
        ELSE
!         FI4 = Achem( ii,HGDET3D ) / Abio( ii,DET )
!    &        * Aflu( ii,FSH_ON_DET ) / SEC_PER_DAY * ITS
          FI5 = Achem( ii,HGPOC3D ) / Abio( ii,DET )
     &        * Aflu( ii,FSH_ON_DET ) / SEC_PER_DAY * ITS
        END IF
!       FI4 = MIN( LIM * Achem( ii,HGDET3D ), FI4 )
        FI4 = 0.
        FI5 = MIN( LIM * Achem( ii,HGPOC3D ), FI5 )

        IF( Abio( ii,ZOS ) .LE. 1E-6 ) THEN
          FIe6 = 0.
          FIi6 = 0.
        ELSE
          FIe6 = Achem( ii,HGZOSe3D ) / Abio( ii,ZOS )
     &        * Aflu( ii,FSH_ON_ZOS ) / SEC_PER_DAY * ITS
          FIi6 = Achem( ii,HGZOSi3D ) / Abio( ii,ZOS )
     &        * Aflu( ii,FSH_ON_ZOS ) / SEC_PER_DAY * ITS
        END IF
        FIe6 = MIN( LIM * Achem( ii,HGZOSe3D ), FIe6 )
        FIi6 = MIN( LIM * Achem( ii,HGZOSi3D ), FIi6 )

        IF( Abio( ii,ZOL ) .LE. 1E-6 ) THEN
          FIe7 = 0.
          FIi7 = 0.
        ELSE
          FIe7 = Achem( ii,HGZOLe3D ) / Abio( ii,ZOL )
     &        * Aflu( ii,FSH_ON_ZOL ) / SEC_PER_DAY * ITS
          FIi7 = Achem( ii,HGZOLi3D ) / Abio( ii,ZOL )
     &        * Aflu( ii,FSH_ON_ZOL ) / SEC_PER_DAY * ITS
        END IF
        FIe7 = MIN( LIM * Achem( ii,HGZOLe3D ), FIe7 )
        FIi7 = MIN( LIM * Achem( ii,HGZOLi3D ), FIi7 )

!       Processes in benthic grid cells
        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN

          IF( Abio( ii,MAC ) .LE. 1E-6 ) THEN
            FIe8 = 0.
            FIi8 = 0.
          ELSE
            FIe8 = HGMACe / Abio( ii,MAC ) 
     &          * Aflu( ii,FSH_ON_MAC ) / SEC_PER_DAY * ITS
            FIi8 = HGMACi / Abio( ii,MAC )
     &          * Aflu( ii,FSH_ON_MAC ) / SEC_PER_DAY * ITS
          END IF
          FIe8 = MIN( LIM * HGMACe, FIe8 )
          FIi8 = MIN( LIM * HGMACi, FIi8 )

          PHY_TOT = Abio( ii,FLA ) + Abio( ii,DIA )
          IF( PHY_TOT .LE. 1E-6 ) THEN
            FLAFRAC = 0.5
            DIAFRAC = 0.5
          ELSE
            FLAFRAC = Abio( ii,FLA ) / PHY_TOT
            DIAFRAC = Abio( ii,DIA ) / PHY_TOT
          END IF

          IF( Abio( ii,FLA ) .LE. 1E-6 ) THEN
            MI1 = 0.
          ELSE
            MI1 = Achem( ii,HGFLA3D ) / Abio( ii,FLA ) / DEP     ! [ng/dm^3] / [mgC/m^3] / [m]
     &          * FLAFRAC * Aflu( ii,MAC_ON_PHY ) / SEC_PER_DAY * ITS  ! * [] * [mgC/m^2] / [s] * [s] = [ng/dm^3]

          END IF
          MI1 = MIN( LIM * Achem( ii,HGFLA3D ), MI1 )

          IF( Abio( ii,DIA ) .LE. 1E-6 ) THEN
            MI2 = 0.
          ELSE
            MI2 = Achem( ii,HGDIA3D ) / Abio( ii,DIA ) / DEP
     &          * DIAFRAC * Aflu( ii,MAC_ON_PHY ) / SEC_PER_DAY * ITS  
          END IF
          MI2 = MIN( LIM * Achem( ii,HGDIA3D ), MI2 )

          IF( Abio( ii,DET ) .LE. 1E-6 ) THEN
            MI4 = 0.
            MI5 = 0.
          ELSE
!           MI4 = Achem( ii,HGDET3D ) / Abio( ii,DET ) / DEP
!    &          * Aflu( ii,MAC_ON_DET ) / SEC_PER_DAY * ITS
            MI5 = Achem( ii,HGPOC3D ) / Abio( ii,DET ) / DEP
     &          * Aflu( ii,MAC_ON_DET ) / SEC_PER_DAY * ITS
          END IF
!         MI4 = MIN( LIM * Achem( ii,HGDET3D ), MI4 )
          MI4 = 0.
          MI5 = MIN( LIM * Achem( ii,HGPOC3D ), MI5 )

          ZOO_TOT = 0.2 * Abio( ii,ZOS ) + 0.3 * Abio( ii,ZOL )
          IF( ZOO_TOT .LE. 1E-6 ) THEN
            ZOSFRAC = 0.4
            ZOLFRAC = 0.6
          ELSE
            ZOSFRAC = 0.2 * Abio( ii,ZOS ) / ZOO_TOT
            ZOLFRAC = 0.3 * Abio( ii,ZOL ) / ZOO_TOT
          END IF

          IF( Abio( ii,ZOS ) .LE. 1E-6 ) THEN
            MIe6 = 0.
            MIi6 = 0.
          ELSE
            MIe6 = Achem( ii,HGZOSe3D ) / Abio( ii,ZOS ) / DEP
     &          * ZOSFRAC * Aflu( ii,MAC_ON_ZOO ) / SEC_PER_DAY * ITS
            MIi6 = Achem( ii,HGZOSi3D ) / Abio( ii,ZOS ) / DEP
     &          * ZOSFRAC * Aflu( ii,MAC_ON_ZOO ) / SEC_PER_DAY * ITS
          END IF
          MIe6 = MIN( LIM * Achem( ii,HGZOSe3D ), MIe6 )
          MIi6 = MIN( LIM * Achem( ii,HGZOSi3D ), MIi6 )

          IF( Abio( ii,ZOL ) .LE. 1E-6 ) THEN
            MIe7 = 0.
            MIi7 = 0.
          ELSE
            MIe7 = Achem( ii,HGZOLe3D ) / Abio( ii,ZOL ) / DEP
     &          * ZOLFRAC * Aflu( ii,MAC_ON_ZOO ) / SEC_PER_DAY * ITS
            MIi7 = Achem( ii,HGZOLi3D ) / Abio( ii,ZOL ) / DEP
     &          * ZOLFRAC * Aflu( ii,MAC_ON_ZOO ) / SEC_PER_DAY * ITS
          END IF
          MIe7 = MIN( LIM * Achem( ii,HGZOLe3D ), MIe7 )
          MIi7 = MIN( LIM * Achem( ii,HGZOLi3D ), MIi7 )

        ! Feeding on sediment
          IF( Abio( ii,STOT ) .LE. 1E-6 ) THEN
            MI9 = 0.
          ELSE

            SEDFRAC = Aflu( ii,MAC_ON_SED ) / Abio( ii,STOT )   ! [mgC/m^2] / [mgC/m^2] / [s] * [s] = []
     &              / SEC_PER_DAY * ITS                         ! dimensionless feeding fraction [0:1]

            IF( SEDFRAC .GE. 1.0 ) THEN
                print*, 'ERROR 230: MAC_ON_SEDFRAC >= 1.0',
     &                  Aflu( ii,MAC_ON_SED ), Abio( ii,STOT )
                SEDFRAC = 0.9999
            END IF

            MI9 = SEDFRAC * HGSED !Ts( ii_i(ii),ii_k(ii),1,HGSED2D )  ! [mgC m-2]
            MI9 = MIN( LIM * HGSED, MI9 )
          END IF
        ELSE
          FIe8 = 0.
          FIi8 = 0.
          MI1 = 0.
          MI2 = 0.
          MI4 = 0.
          MI5 = 0.
          MIe6 = 0.
          MIi6 = 0.
          MIe7 = 0.
          MIi7 = 0.
          MI9 = 0.
        END IF

!       Apply fraction for accumulation and excretion
        AS1 = E_FLA * FS1
        ES1 = ( 1 - E_FLA ) * FS1
        AS2 = E_DIA * FS2
        ES2 = ( 1 - E_DIA ) * FS2
        AS3 = E_CYA * FS3
        ES3 = ( 1 - E_CYA ) * FS3
        AS4 = E_DET * FS4
        ES4 = ( 1 - E_DET ) * FS4
        AS5 = E_DET * FS5
        ES5 = ( 1 - E_DET ) * FS5
        AL1 = E_FLA * FL1
        EL1 = ( 1 - E_FLA ) * FL1
        AL2 = E_DIA * FL2
        EL2 = ( 1 - E_DIA ) * FL2
        AL3 = E_CYA * FL3
        EL3 = ( 1 - E_CYA ) * FL3
        AL4 = E_DET * FL4
        EL4 = ( 1 - E_DET ) * FL4
        AL5 = E_DET * FL5
        EL5 = ( 1 - E_DET ) * FL5
        ALi6 = Ei_ZOS * FLi6
        ELi6 = ( 1 - Ei_ZOS ) * FLi6
        ALe6 = Ee_ZOS * FLe6
        ELe6 = ( 1 - Ee_ZOS ) * FLe6

        AI4 = E_DET * FI4
        EI4 = ( 1 - E_DET ) * FI4
        AI5 = E_DET * FI5
        EI5 = ( 1 - E_DET ) * FI5
        AIi6 = Ei_ZOS * FIi6
        EIi6 = ( 1 - Ei_ZOS ) * FIi6
        AIe6 = Ee_ZOS * FIe6
        EIe6 = ( 1 - Ee_ZOS ) * FIe6
        AIi7 = Ei_ZOL * FIi7
        EIi7 = ( 1 - Ei_ZOL ) * FIi7
        AIe7 = Ee_ZOL * FIe7
        EIe7 = ( 1 - Ee_ZOL ) * FIe7

!       Macrobenthos section
        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
          AIe8 = Ee_MAC * FIe8
          EIe8 = ( 1 - Ee_MAC ) * FIe8
          AIi8 = Ei_MAC * FIi8
          EIi8 = ( 1 - Ei_MAC ) * FIi8
          FIe8 = MIN( FIe8,HGMACe )
          FIi8 = MIN( FIi8,HGMACi )
          IF( FIi8 .GE. 1.0E-8 .AND. FIi8 .GE. HGMACi ) THEN
            print*, "ERROR 234: Mac_i over feeding",
     &              FIi8, DEP, HGMACi
          END IF
          IF( FIe8 .GE. 1.0E-8 .AND. FIe8 .GE. HGMACe ) THEN
            print*, "ERROR 234: Mac_e over feeding",
     &              FIe8, DEP, HGMACe
          END IF

          AM1 = E_FLA * MI1
          EM1 = ( 1 - E_FLA ) * MI1
          AM2 = E_DIA * MI2
          EM2 = ( 1 - E_DIA ) * MI2
          AM4 = E_DET * MI4 
          EM4 = ( 1 - E_DET ) * MI4
          AM5 = E_DET * MI5 
          EM5 = ( 1 - E_DET ) * MI5
          AMi6 = Ei_ZOS * MIi6 
          EMi6 = ( 1 - Ei_ZOS ) * MIi6
          AMe6 = Ee_ZOS * MIe6
          EMe6 = ( 1 - Ee_ZOS ) * MIe6
          AMi7 = Ei_ZOL * MIi7 
          EMi7 = ( 1 - Ei_ZOL ) * MIi7
          AMe7 = Ee_ZOL * MIe7 
          EMe7 = ( 1 - Ee_ZOL ) * MIe7
          AM9 = E_DET * MI9
          EM9 = ( 1 - E_DET ) * MI9
        ELSE
          AIe8 = 0.
          EIe8 = 0.
          AIi8 = 0.
          EIi8 = 0.
          AM1 = 0.
          EM1 = 0.
          AM2 = 0.
          EM2 = 0.
          AM4 = 0.
          EM4 = 0.
          AM5 = 0.
          EM5 = 0.
          AMe6 = 0.
          EMe6 = 0.
          AMi6 = 0.
          EMi6 = 0.
          AMe7 = 0.
          EMe7 = 0.
          AMi7 = 0.
          EMi7 = 0.
          AM9 = 0.
          EM9 = 0.
        END IF

!        IF( ii .EQ. 60559 ) THEN
        IF( 1 .EQ. 1 ) THEN
!         UP_MAC = 0.
!         RE_MAC = 0.
!         RI_MAC = 0.
          MMACe = 0.
          MMACi = 0.
          FIe8 = 0.
          AIe8 = 0.
          EIe8 = 0.
          FIi8 = 0.
          AIi8 = 0.
          EIi8 = 0.
          MI1 = 0.
          AM1 = 0.
          EM1 = 0.
          MI2 = 0.
          AM2 = 0.
          EM2 = 0.
          MI4 = 0.
          AM4 = 0.
          EM4 = 0.
          MI5 = 0.
          AM5 = 0.
          EM5 = 0.
          MIe6 = 0.
          AMe6 = 0.
          EMe6 = 0.
          MIi6 = 0.
          AMi6 = 0.
          EMi6 = 0.
          MIe7 = 0.
          AMe7 = 0.
          EMe7 = 0.
          MIi7 = 0.
          AMi7 = 0.
          EMi7 = 0.
          MI9 = 0.
          AM9 = 0.
          EM9 = 0.
        END IF

C........ Mortality
!      ng/L = ng/L * day-1 * day s-1 *  s
        IF( Abio( ii,FLA ) .LE. 1E-6 ) THEN
          MFLA = 0.
        ELSE
          MFLA = Achem( ii, HGFLA3D ) / Abio( ii,FLA )
     &         * Aflu( ii,MORT_FLA ) / SEC_PER_DAY * ITS
        END IF
        MFLA = MIN( 0.99 * Achem( ii,HGFLA3D ), MFLA )

        IF( Abio( ii,DIA ) .LE. 1E-6 ) THEN
          MDIA = 0.
        ELSE
          MDIA = Achem( ii, HGDIA3D ) / Abio( ii,DIA ) 
     &         * Aflu( ii,MORT_DIA ) / SEC_PER_DAY * ITS
        END IF
        MDIA = MIN( 0.99 * Achem( ii,HGDIA3D ), MDIA )

        IF( Abio( ii,CYA ) .LE. 1E-6 ) THEN
          MCYA = 0.
        ELSE
          MCYA = Achem( ii, HGCYA3D ) / Abio( ii,CYA )
     &         * Aflu( ii,MORT_CYA ) / SEC_PER_DAY * ITS
        END IF
        MCYA = MIN( 0.99 * Achem( ii,HGCYA3D ), MCYA )


        IF( Abio( ii,ZOS ) .LE. 1E-6 ) THEN
          MZOSe = 0.
          MZOSi = 0.
        ELSE
          MZOSe = Achem( ii, HGZOSe3D ) / Abio( ii,ZOS ) 
     &         * Aflu( ii,MORT_ZOS ) / SEC_PER_DAY * ITS
          MZOSi = Achem( ii, HGZOSi3D ) / Abio( ii,ZOS )
     &         * Aflu( ii,MORT_ZOS ) / SEC_PER_DAY * ITS
        END IF
        MZOSe = MIN( 0.99 * Achem( ii,HGZOSe3D ), MZOSe )
        MZOSi = MIN( 0.99 * Achem( ii,HGZOSi3D ), MZOSi )

        IF( Abio( ii,ZOL ) .LE. 1E-6 ) THEN
          MZOLe = 0.
          MZOLi = 0.
        ELSE
          MZOLe = Achem( ii, HGZOLe3D ) / Abio( ii,ZOL ) 
     &         * Aflu( ii,MORT_ZOL ) / SEC_PER_DAY * ITS
          MZOLi = Achem( ii, HGZOLi3D ) / Abio( ii,ZOL )
     &         * Aflu( ii,MORT_ZOL ) / SEC_PER_DAY * ITS
        END IF
        MZOLe = MIN( 0.99 * Achem( ii,HGZOLe3D ), MZOLe )
        MZOLi = MIN( 0.99 * Achem( ii,HGZOLi3D ), MZOLi )

        IF( Abio( ii,FSH ) .LE. 1E-6 ) THEN
          MFSHe = 0.
          MFSHi = 0.
        ELSE
          MFSHe = Achem( ii, HGFISHe3D )      / Abio( ii,FSH ) 
     &         * Aflu( ii,MORT_FSH ) / SEC_PER_DAY * ITS
          MFSHi = Achem( ii, HGFISHi3D )     / Abio( ii,FSH )
     &         * Aflu( ii,MORT_FSH ) / SEC_PER_DAY * ITS
        END IF
        MFSHe = MIN( 0.99 * Achem( ii,HGFISHe3D ), MFSHe )
        MFSHi = MIN( 0.99 * Achem( ii,HGFISHi3D ), MFSHi )

!       Fraction going to DOM (40%)
        XFLA = FRR * MFLA
        XDIA = FRR * MDIA
        XCYA = FRR * MCYA
        XZOS = FRR * ( MZOSe + MZOSi )
        XZOL = FRR * ( MZOLe + MZOLi )
        XFSH = FRR * ( MFSHe + MFSHi )

!       Fraction goint to Abio( ii,DET ) (60%)
        YFLA = ( 1. - FRR ) * MFLA
        YDIA = ( 1. - FRR ) * MDIA
        YCYA = ( 1. - FRR ) * MCYA
        YZOS = ( 1. - FRR ) * ( MZOSe + MZOSi )
        YZOL = ( 1. - FRR ) * ( MZOLe + MZOLi )
        YFSH = ( 1. - FRR ) * ( MFSHe + MFSHi )

!       Macrobenthos mortality [ng/m^2]
        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
          IF( Abio( ii,MAC ) .LE. 1E-6 ) THEN
            MMACe = 0.
            MMACi = 0.
          ELSE
            MMACe = HGMACe / Abio( ii,MAC )
     &           * Aflu( ii,MORT_MAC ) / SEC_PER_DAY * ITS
            MMACi = HGMACi / Abio( ii,MAC )
     &           * Aflu( ii,MORT_MAC ) / SEC_PER_DAY * ITS
          END IF
          MMACe = MIN( 0.99 * HGMACe, MMACe )
          MMACi = MIN( 0.99 * HGMACi, MMACi )
          XMAC = FRR * ( MMACe + MMACi )
          YMAC = ( 1. - FRR ) * ( MMACe + MMACi )
        ELSE
          MMACe = 0.
          MMACi = 0.
          XMAC = 0.
          YMAC = 0.
        END IF

!        IF( ii .EQ. 60559 ) THEN
        IF( 1 .EQ. 1 ) THEN
          UP_FSH = 0.
          MFSHe = 0.
          MFSHi = 0.
          AI4 = 0.
          AI5 = 0.
          AIe6 = 0.
          AIi6 = 0.
          AIe7 = 0.
          AIi7 = 0.
          AIe8 = 0.
          AIi8 = 0.
          MMACe = 0.
          MMACi = 0.
          XMAC = 0.
          YMAC = 0.
        END IF

        MI9 = 0
        XMAC = 0
        YMAC = 0
        IF( YMAC .LT. 0.) print*, 'NEGATIVE YMAC', YMAC,
     &      Aflu( ii,MORT_MAC ), Achem( ii,HGMACe3D ), Abio( ii,MAC )
        IF( XMAC .LT. 0.) print*, 'NEGATICE XMAC', XMAC

!       Apply bioaccumulation tendencies
        Achem( ii,HGDIS3D ) =
     &  Achem( ii,HGDIS3D ) + RDET   + RDOC 
     &                      - UP_FLA - UP_DIA - UP_CYA 
     &                      - UP_ZOS - UP_ZOL - UP_MAC - UP_FSH
     &                      + R_FLA  + R_DIA  + R_CYA 
     &                      + RE_ZOS + RE_ZOL + RE_MAC + RE_FSH
     &                      + RI_ZOS + RI_ZOL + RI_MAC + RI_FSH

        Achem( ii,HGDIS3D ) = ! used to go to HGDET
     &  Achem( ii,HGDIS3D ) - RDET - FS4 - FL4 - FI4 - MI4
     &               + YFLA + YDIA + YCYA + YZOS + YZOL + YFSH !+ YMAC

        Achem( ii,HGPOC3D ) =
     &  Achem( ii,HGPOC3D )        - FS5 - FL5 - FI5 - MI5

        Achem( ii,HGDIS3D ) = ! used to go to HGDOC
     &  Achem( ii,HGDIS3D ) - RDOC
     &                      + ES1 + ES2 + ES3 + ES4 + ES5
     &                      + EL1 + EL2 + EL3 + EL4 + EL5 + ELe6 + ELi6
     &                      + EI4 + EI5 + EIe6 + EIi6 + EIe7 + EIi7
     &                      + EIe8 + EIi8 
     &                      + EM1 + EM2 + EM4 + EM5 + EMe6 + EMi6 
     &                      + EMe7 + EMi7 + EM9
     &               + XFLA + XDIA + XCYA + XZOS + XZOL + XFSH + XMAC

        Achem( ii,HGFLA3D ) = 
     &  Achem( ii,HGFLA3D ) + UP_FLA - R_FLA - FS1 - FL1 - MI1 - MFLA

        Achem( ii,HGDIA3D ) = 
     &  Achem( ii,HGDIA3D ) + UP_DIA - R_DIA - FS2 - FL2 - MI2 - MDIA

        Achem( ii,HGCYA3D ) = 
     &  Achem( ii,HGCYA3D ) + UP_CYA - R_CYA - FS3 - FL3       - MCYA

        Achem( ii,HGZOSe3D ) = 
     &  Achem( ii,HGZOSe3D ) + UP_ZOS
     &                       - RE_ZOS - FLe6 - FIe6 - MIe6 - MZOSe

        Achem( ii,HGZOLe3D ) = 
     &  Achem( ii,HGZOLe3D ) + UP_ZOL
     &                       - RE_ZOL        - FIe7 - MIe7 - MZOLe

        Achem( ii,HGZOSi3D ) =
     &  Achem( ii,HGZOSi3D ) - RI_ZOS - FLi6 - FIi6 - MIi6 - MZOSi
     &                       + AS1 + AS2 + AS3 + AS4 + AS5

        Achem( ii,HGZOLi3D ) =
     &  Achem( ii,HGZOLi3D ) - RI_ZOL - FIi7 - MIi7 - MZOLi
     &                       + AL1 + AL2 + AL3 + AL4 + AL5 + ALe6 + ALi6

        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
          HGMACe = HGMACe + UP_MAC - RE_MAC - FIe8 - MMACe

          Achem( ii,HGMACe3D ) = HGMACe
  
          HGMACi = HGMACi          - RI_MAC - FIi8 - MMACi
     &         + AM1 + AM2 + AM4 + AM5 + AMe6 + AMi6 + AMe7 + AMi7 + AM9

          Achem( ii,HGMACi3D ) = HGMACi
        END IF

        Achem( ii,HGFISHe3D ) =
     &  Achem( ii,HGFISHe3D ) + UP_FSH - MFSHe 

        Achem( ii,HGFISHi3D ) =
     &  Achem( ii,HGFISHi3D ) +AI4 +AI5 +AIe6 +AIi6 +AIe7 +AIi7 +AIe8 +AIi8
     &                - MFSHi

        Ts( ii_i( ii ),ii_k( ii ),1,HGSED2D ) =
     &  Ts( ii_i( ii ),ii_k( ii ),1,HGSED2D ) - MI9  * 1000. * DEP
     &                                        + YMAC * 1000. * DEP

!move this to sediment.f
        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
            Achem( ii,SEDHG3D ) = Ts(ii_i( ii ),ii_k( ii ),1,HGSED2D )
        END IF

!...... Enforce minimum concentration and mass conservation
        FNORM = 0.
        DO NNC = 1,NHGTRANS
          IF( Achem( ii,NNC ) .LT. 0. ) THEN
            print*, "ERROR 950: bioacc.f -> negative ",
     &              VNAMOUT( NNC ), Achem( ii,NNC )
          END IF

          Achem( ii,NNC ) = MAX( CMIN,Achem( ii,NNC ) )
          FNORM = FNORM + Achem( ii,NNC )
        END DO

        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
            HGMACe = Achem( ii,HGMACe3D ) / 1000. / DEP ! [ng/dm^3} = [ng/m^2] / [dm^3/m^3] / [m]
            HGMACi = Achem( ii,HGMACi3D ) / 1000. / DEP
            HGSED  = Ts( ii_i(ii),ii_k(ii),1,HGSED2D ) / 1000. / DEP
        END IF

        FNORM = FNORM + Achem( ii, HGFISHe3D ) + Achem( ii, HGFISHi3D )
        FNORM = FNORM + HGMACe + HGMACi
        FNORM = FNORM + HGSED
        LOL = HGSED
! Set diagnostic 2d bottom output fields
        IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
            Ts( ii_i(ii),ii_k(ii),1,HGMACe2D ) = Achem( ii,HGMACe3D )
            Ts( ii_i(ii),ii_k(ii),1,HGMACi2D ) = Achem( ii,HGMACi3D )
            Ts( ii_i(ii),ii_k(ii),1,MAC2D   ) = Abio( ii,MAC2D )
            Ts( ii_i(ii),ii_k(ii),1,STOT2D  ) = Abio( ii,STOT2D )
        END IF

        ROFL = FNORM
        IF( FNORM .GT. 0 ) THEN
            FNORM = HGTOT / FNORM
        END IF

        IF( FNORM .GT. 1.005 .OR. FNORM .LT. 0.995 ) THEN
            print*, 'WARNING 123: bioacc_hg.f(809) FNORM=',
     &              FNORM,ii,ii_i(ii),ii_k(ii),ii_j(ii)
                print*, 'FNORM  = ', FNORM
                print*, 'BEFORE = ', ROFL
                print*, 'HGTOT  = ', HGTOT
                print*, 'DEPTH  = ', DEP
                print*, 'HGDIS  = ', Achem( ii,HGDIS3D )
                print*, ' UP_MAC = ', UP_MAC
                print*, 'HGMACe = ', HGMACe
                print*, ' RE_MAC = ', RE_MAC
                print*, 'HGMACi = ', HGMACi
                print*, ' RI_MAC = ', RI_MAC
                print*, '  FIi8 = ', FIi8
                print*, '  AM1  = ', AM1
                print*, '  AM2  = ', AM2
                print*, '  AM4  = ', AM4
                print*, '  AM5  = ', AM5
                print*, '  AMe6 = ', AMe6
                print*, '  AMi6 = ', AMi6
                print*, '  AMe7 = ', AMe7
                print*, '  AMi7 = ', AMi7
                print*, '  AM9  = ', AM9
                print*, ' -- '
                print*, 'FLA        = ', Abio( ii,FLA )
                print*, 'DIA        = ', Abio( ii,DIA )
                print*, 'MAC_ON_PHY = ', Aflu( ii,MAC_ON_PHY )
                print*, ' MI1 = ', MI1
                print*, ' MI2 = ', MI2
                print*, 'ZOS    = ', Abio( ii,ZOS )
                print*, 'ZOL    = ', Abio( ii,ZOL )
                print*, ' MIe6 = ', MII6
                print*, ' MIi6 = ', MIi6
                print*, ' MIe7 = ', MIe7
                print*, ' MIi7 = ', MIi7
                print*, 'MAC_ON_ZOO = ', Aflu( ii,MAC_ON_ZOO )
                print*, 'DET    = ', Abio( ii,ZOL )
                print*, ' MIe6 = ', MIe6
                print*, 'MAC_ON_DET = ', Aflu( ii,MAC_ON_DET )
                print*, 'SED    = ', Abio( ii,STOT )
                print*, ' MI4 = ', MI4
                print*, ' MI5 = ', MI5
                print*, 'SED before = ', LOL
                print*, 'SED after  = ', HGSED
                print*, 'MAC_ON_SED = ', Aflu( ii,MAC_ON_SED )
                print*, ' MI9 = ', MI9
                print*, 'MMACe = ', MMACe
                print*, 'MMACi = ', MMACi
                print*, 'FISH = ', Abio( ii,FSH )
                print*, 'FSH_ON_MAC = ', Aflu( ii,FSH_ON_MAC )
                print*, ' FIe8 = ', FIe8
                print*, ' FIi8 = ', FIi8       
                print*,'----------------------------------------------'
                print*, ''

            print*, 'RI_MAC ', RI_MAC, 'FIi8 ', FIi8, 'MMACi ', MMACi
            IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
!               print*, 'WARNING 123: bioacc.f(809) FNORM=',FNORM,ii
                print*, ' |-> benthic layer'
            END IF
        END IF

!       Check
        IF( FNORM .LT. 0. ) THEN
          DO NNC = 1,NHGTRANS
            print*, VNAMOUT(NNC), Achem( ii,NNC )
          END DO
            print*, 'MAC', Abio( ii,MAC )
            print*, 'HGMACe', Achem( ii,HGMACe3D )
            print*, 'HGMACi', Achem( ii,HGMACi3D )
            print*, 'FSH', Abio( ii,FSH )
            print*, 'HGFISHe', Achem( ii,HGFISHe3D )
            print*, 'HGFISHi', Achem( ii,HGFISHi3D )

            print*, 'HGSED', Ts( ii_i(ii),ii_k(ii),1,HGSED2D )
            print*, 'HGSED', Ts( ii_i(ii),ii_k(ii),1,HGSED2D )
        END IF

        DO NNC = 1, NHGTRANS
!          IF( ii_j( ii ) .NE. ii_j_max( ii ) ) THEN
             Achem( ii,NNC ) = Achem( ii,NNC ) * FNORM
!          END IF
        END DO
        Achem( ii,HGMACe3D ) = Achem( ii,HGMACe3D ) * FNORM
        Achem( ii,HGMACi3D ) = Achem( ii,HGMACi3D ) * FNORM

C.........  End of subroutine: 
        END SUBROUTINE BIOACC_HG
