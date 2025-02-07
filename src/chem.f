
        SUBROUTINE CHEM( Achem, Abio, TEMP, SALT, RADI, VERBOSE, ii )

        USE CPARAM
        USE CIN
        USE COUT

        IMPLICIT NONE

C...........
!
! DESC:         Chemical mechanism
!
! HISTORY:      13.10.2021    Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........


C...........   INCLUDES:

        INCLUDE 'C_model_inc'

C...........   EXTERNAL FUNCTIONS and their descriptions:

C...........   PARAMETERS and their descriptions:
        REAL, PARAMETER :: SEC_PER_DAY = 3600. * 24.
        REAL, PARAMETER :: CMIN = 1E-10         ! Minimum concentration
        REAL, PARAMETER :: RADI_TO_PAR = 0.5211 ! Jacovides et al., 2004
        REAL, PARAMETER :: E_TO_W = 4.57        ! McCree, 1972

        REAL  Vpam
!       REAL, PARAMETER :: kpam = 4.63E-7
        REAL, PARAMETER :: kpam = 1.1547E-7

C...... (R1) dark/chemical oxidation
C...... kod : Hg0 --> Hg2+
        REAL  Vod
        REAL, PARAMETER :: kod = 2.60E-06       ![s-1]        !Kuss et al., 2015

C...... (R2) dark/chemical reduction
C...... Hg2+ --> Hg0
        REAL  Vrd
        REAL  krd_T
!       REAL, PARAMETER :: krd = 0.0
        REAL, PARAMETER :: krd = 6.00E-7        ![s-1]        !Kuss et al., 2015

C...... (R3) biogenic reduction
C...... krb : Hg2+ --> Hg0
        REAL  Vbr
!       REAL, PARAMETER :: krb = 0.0
        REAL, PARAMETER :: krb = 8.06E-09       ![m3 mgC-1 s-1] !Kuss et al., 2015

C...... (R4) photolytic reduction
C...... kpr : Hg2+ + phot --> Hg0
        REAL  Vpr
!       REAL, PARAMETER :: kpr = 0.0
        REAL, PARAMETER :: kpr = 1.00E-08       ![s-1 W-1 m2] !Kuss et al., 2015

C...... (R5) photolytic oxidation
C...... kpo : Hg0 + phot --> Hg*
        REAL  Vpo
!       REAL, PARAMETER :: kpo = 0.0
        REAL, PARAMETER :: kpo = 0.24E-08       ![s-1 W-1 m2] !Kuss et al., 2015

C...... (R6) HgS formation
C...... ksf : Hg2+(aq) + H2S -> HgS(s) + 2H+
        REAL  Vhgf
!       REAL, PARAMETER :: ksf = 0.0
        REAL, PARAMETER :: ksf = 4.90E-4        ! [s-1]       !Slowey et al., 2010

C...... (R7) HgS-DOC dissolution
C...... HgS(p) + DOM(aq) -> HgS-DOM(aq)
        REAL Vhgd
!       REAL, PARAMETER :: ksd = 0.0
        REAL, PARAMETER :: ksd = 5.78E-6        ! [s-1]       !Ping Jiang et al., 2016

C...... (R8) HgS-DOC re-crystallisation
C...... HgS-DOC(aq) -> HgS(s)
        REAL Vhga
!       REAL, PARAMETER :: ksa = 0.0
        REAL, PARAMETER :: ksa = 9.50E-6        ! [s-1]       !Pian Jiang, 2016/Slowey 2010

C...... (R9) HgS oxidation
C...... kso : HgS(s) + 2(O2) -> Hg2+ + SO4(2-)
        REAL Vhgo
!       REAL, PARAMETER :: kso = 0.0
        REAL, PARAMETER :: kso = 1.00E-4        ! [l/ml(O2) l/ng(Hg) s-1] !Evgenyi/Petrochowa 2019

C...... (R10) Constant methylation Hg2+(aq) -> MeHg+
        REAL Vmx
!       REAL, PARAMETER :: kmx = 0.0            ! default
!       REAL, PARAMETER :: kmx = 4.05E-09       ![s-1]        !Olsen et al., 2018
!       REAL, PARAMETER :: kmx = 1.74E-08       !]s-1]        !
        REAL, PARAMETER :: kmx = 3.47E-08       ![s-1]        !Duran, 2008
!       REAL, PARAMETER :: kmx = 2.21E-07       ![s-1]        !Monperuss, 2007 avg (*)

C...... (R11) mono-methylation (abiotic anearob/hypoxic waters and sediments only)
C...... km1 : Hg2+ --> MMHg @(O2 < 0)
        REAL  Vm1
!       REAL, PARAMETER :: km1 = 0.0
!       REAL, PARAMETER :: km1 = 2.30E-09       ![s-1]        !Monperrus, 2007 min
!       REAL, PARAMETER :: km1 = 4.05E-09       ![s-1]        !Olsen et al., 2018
!       REAL, PARAMETER :: km1 = 7.50E-08       ![s-1]        !Lehnherr, 2011 in SCM
!       REAL, PARAMETER :: km1 = 7.87E-08       ![s-1]        !Lehnherr, 2011 average
        REAL, PARAMETER :: km1 = 2.21E-07       ![s-1]        !Monperuss, 2007 avg (*)
!       REAL, PARAMETER :: km1 = 4.40E-07       ![s-1]        !Monperrus, 2007 max

C...... (R12) kmb : Hg2+ --> MMHg (biologically induced)
        REAL  Vmb
!       REAL, PARAMETER :: kmb = 0.0
!       REAL, PARAMETER :: kmb = 6.00E-10 !1.00E-09       ![mgC s-1]
        REAL, PARAMETER :: kmb = 4.05E-09       ![s-1]        !Olsen et al., 2018
!       REAL, PARAMETER :: kmb = 6.05E-09       ![s-1]        !Olsen/Duran average
!       REAL, PARAMETER :: kmb = 3.47E-08       ![s-1]        !Duran, 2008
!       REAL, PARAMETER :: km1 = 2.21E-07       ![s-1]        !Monperuss, 2007 avg (*)
!       REAL, PARAMETER :: kmb = 4.40E-07       ![s-1]        !Monperrus, 2007 max


C...... (R13) bi-methylation
C...... km2 : Hg2+ --> DMHg
        REAL  Vm2
!       REAL, PARAMETER :: km2 = 0.0
!       REAL, PARAMETER :: km2 = 1.40E-10       ![s-1]        !Lehnherr, 2011 in SCM
        REAL, PARAMETER :: km2 = 4.63E-10       ![s-1]        !Lehnherr, 2011 average
!       REAL, PARAMETER :: km2 = 5.71E-10       ![s-1]        !Mason, 1995 (see also Mason and Fitzgerald, 1993)

C...... (R14) secondary methylation
C...... km3 : MMHg --> DMHg
        REAL  Vm3
        REAL, PARAMETER :: km3 = 1.51E-08       ![s-1]        !Lehnherr, 2011 average
!       REAL, PARAMETER :: km3 = 1.85E-08       ![s-1]        !Lehnherr, 2011 in SCM

C...... (R15) demethylation
C...... kdm : MMHg --> Hg2+
        REAL  Vdm
        REAL  kdm_T
!       REAL, PARAMETER :: kdm = 5.71E-09       ![s-1]        !Mason, 1996
        REAL, PARAMETER :: kdm = 6.94E-07       ![s-1]        !Monperrus, 2007
!       REAL, PARAMETER :: kdm = 1.0648E-06     ![s-1]        !Olsen et al., 2018
!       REAL, PARAMETER :: kdm = 3.24E-06       ![s-1]        !Lehnherr, 2011 min
!       REAL, PARAMETER :: kdm = 4.05E-06       ![s-1]        !Lehnherr, 2011 avg (*)
!       REAL, PARAMETER :: kdm = 4.86E-06       ![s-1]        !Lehnherr, 2011 max

C...... (R16 & R17 & R18) photo demethylation
C...... kpd : MMHg + phot --> Hg2+  &&  DMHg + phot --> MMHg && DMHg + phot --> Hg2+ (experimental)
        REAL  Vpd
        REAL, PARAMETER :: kpd = 4.57E-09       ![m2 W-1 s-1] = 1.00E-03 [m2 E-1 s-1]   !Lehnherr, 2011
!       REAL, PARAMETER :: kpd = 3.47E-08       ![s-1]        !Monperrus, 2007 min
!       REAL, PARAMETER :: kpd = 1.04E-06       ![s-1]        !Whalin, 2007 min
!       REAL, PARAMETER :: kpd = 1.62E-06       ![s-1]!       !Monperrus, 2007 max
!       REAL, PARAMETER :: kpd = 4.98E-06       ![s-1]        !Whalin, 2007 max

C...... (R19) de-dimethylation
C...... kdd : DMHg --> MMHg
        REAL  Vdd
!       REAL, PARAMETER :: kdd = 0.0
        REAL, PARAMETER :: kdd = 2.22E-09       ![s-1]        !Mason, 1996
!       REAL, PARAMETER :: kdd = 2.30E-06       ![s-1]        !Mason, 1999

C...... (R20) reductive de-methylation
C...... krm : MMHg --> Hg0
        REAL  Vrm
!       REAL, PARAMETER :: krm = 0.0
        REAL, PARAMETER :: krm = 2.22E-09       ![s-1]        !

C...... (R21) atmospheric de-dimethylation
C...... kda : DMHg(air) --> MMHg(air)
        REAL  Vda
!       REAL, PARAMETER :: kda = 1.16E-06       ![s-1]        !10% per day estimate
        REAL, PARAMETER :: kda = 2.78E-05       ![s-1]        !10% per hour estimate

C...... POC/SPM ration
!       REAL, PARAMETER :: PSR = 0.16           ![]           !Wozniak, 2010 (2.3%-42.0%)
!       REAL, PARAMETER :: PSR = 0.10           ![]           !Sharif, 2014  (9.3%-10.8%)
        REAL, PARAMETER :: PSR = 1.00           ![]           !No SPM case

        REAL, PARAMETER :: FRED = 0.4           ![]           !reducible fraction of Hg2+

        REAL  Vdoc, Vdet
        REAL  dDOC, dDET
        REAL TDREM, TPREM

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        INTEGER                          , INTENT(IN)     :: ii         ! cell number           []
        REAL, DIMENSION( ndrei,nchem )   , INTENT(INOUT)  :: Achem      ! chemistry array       [ng/L]
        REAL, DIMENSION( ndrei,3:ninbio ), INTENT(IN)     :: Abio       ! biology array         [mgC/m^3]
        REAL, DIMENSION( ndrei )         , INTENT(IN)     :: SALT       ! salinity in PSU       [PSU]
        REAL, DIMENSION( ndrei )         , INTENT(IN)     :: TEMP       ! temperature in °C     [°C]
        REAL, DIMENSION( ndrei )         , INTENT(IN)     :: RADI       ! radiation in W/m**2   [W/m^2]

        INTEGER, INTENT(IN) :: VERBOSE  ! print debug information				[0,1,2]

C...........   LOCAL VARIABLES and their descriptions:
        INTEGER NNC     ! species counter for loops
        REAL    HGTOT, FNORM, HGSTOT

        REAL    CORG    ! total mass of organic carbon CORG = DOC + POC + BPC + BZC
        REAL    TPHYT   ! total phytoplankton mass [mgC/m**3]
        REAL    FDOC    ! fraction of DOC / POM
        REAL    FPOC    ! fraction of POC / POM
        REAL    FBPC    ! fraction of BPC / POM
        REAL    FBZC    ! fraction of BZC / POM

        REAL    HG2TOT  ! total mass of all mercury oxidized species DOM + POM
        REAL    HG2RED  ! reducible fraction of dissolved mercury DOM * FDIST * FFORT

        REAL :: OXI_D   ! amount oxidized mercury in current timestep
        REAL :: OXI_P   ! amount oxidized mercury in current timestep

        REAL :: REDD_DIS
        REAL :: REDD_DOC
        REAL :: REDB_DIS
        REAL :: REDB_DOC        ! amount of reduced mercury in current timestep
        REAL :: REDB_POC        ! amount of reduced mercury in current timestep
        REAL :: REDP_DIS        ! amount of reduced mercury in current timestep
        REAL :: REDP_DOC       
        REAL :: REDP_POC

        REAL    MET_RED1, MET_RED2   ! MeHg-DOM reduction
        REAL    MET1_DIS,MET1_DOC,MET1_POC,MET1_HGS
        REAL    METB_DIS,METB_DOC,METB_POC,METB_HGS
        REAL    MET2_DIS,MET2_DOC,MET2_POC,MET2_HGS
        REAL                      MET3_DIS,MET3_DOC,MET3_POC
        REAL                      DEMP_DIS,DEMP_DOC
        REAL                      DEM_DIS,DEM_DOC,DEM_POC
        REAL    METX_D  ,METX_S
        REAL    DED, DEP, DEP2
        REAL    LABFRAC   ! non-labile fraction of DOM associated Hg

        REAL    C_O2
        REAL    C_H2S
        REAL    DIS_FORM,DOC_FORM 
        REAL    FORM_OXI,FORM_DISS
        REAL    DISS_ADS,DISS_OXI,DISS_MET1,DISS_METB,DISS_MET2
        REAL    DISS_REDB,DISS_REDP

        REAL REDDOC,REDBPC,REDBZC,REDPOC

        REAL RDET, RDOC ! reminerilazation rate of detritus ond DOM
        REAL FREM       ! total reminerilazation
        REAL TTOM       ! total terrestrial organic matter
        REAL DMOM

        CHARACTER( 265 ) MESG
        
C***********************************************************************
C   begin body of subroutine CHEM
!        print*, "START OF SUBROUTINE CHEM", ii, ndrei

        C_O2 = Abio( ii,O2 ) - 2.5

!       total mercury
        HGTOT = 0.
        DO NNC = 1,NHGTRANS
          HGTOT = HGTOT + Achem( ii,NNC )
          IF( Achem( ii,NNC ) < 0. ) THEN
            print*, 'X9', NNC, Achem( ii,NNC )
            Achem( ii,NNC ) = 0.
          END IF
        END DO

!...... Chemical reactions
!       C.1 red-ox chemisty

!	C. dark oxidation (Hg0 --> Hg+2)
        Vod = 1. - EXP( -ITS * kod )
        IF( Vod .GT. 0.9 ) print*, 'ERROR 725: Vod = ',Vod, ii
        OXI_D = Achem( ii,HGDEM3D ) * Vod

!       C. photolytic oxidation (Hg0 --> Hg*) now Hg0 --> Hg2+ direct
        Vpo = 1. - EXP( -ITS * kpo * RADI( ii ) * RADI_TO_PAR )
        IF( Vpo .GT. 0.9 ) print*, 'ERROR 726: Vor = ',Vpo, ii
        OXI_P = Achem( ii,HGDEM3D ) * Vpo

!       C. biogenic reduction (Hg2+ --> Hg0)
        TPHYT = Abio( ii,CYA )

        Vbr = 1. - EXP( -ITS * krb * TPHYT )  
        IF( Vbr .GT. 0.9 ) print*, 'ERROR 729: Vbr = ',Vbr, ii

        REDB_DIS = Achem( ii,HGDIS3D ) * FRED * Vbr   
        REDB_DOC = 0.
        REDB_POC = 0.   
        DISS_REDB = 0.

!       C. Photolytic reduction (Hg2+ --> Hg0)
        Vpr = 1. - EXP( -ITS * kpr * RADI( ii ) * RADI_TO_PAR )
        IF( Vpr .GT. 0.9 ) print*, 'ERROR 730: Vpr = ',Vpr, ii
        REDP_DIS  = Achem( ii,HGDIS3D ) * FRED * Vpr    ! reducible fraction
        REDP_DOC  = 0.
        DISS_REDP = 0.
        REDP_POC  = 0.

!       C. dark reduction (Hg2+ --> Hg0)
!       krd_T = krd
        IF( krd .EQ. 0. ) THEN
            krd_T = 0.
        ELSE
            krd_T = 2.92E-07 * EXP( 0.045 * TEMP( ii ) )
        END IF
        Vrd = 1. - EXP( -ITS * krd_T )
        IF( Vrd .GT. 0.9 ) print*, 'ERROR 770: Vrd = ',Vrd, ii
        REDD_DIS = Achem( ii,HGDIS3D ) * Vrd * FRED
        REDD_DOC = 0.

!       C.2 methylation / demethylation
!       C. reductive de-methylation
        Vrm = 1. - EXP( -ITS * krm )
        MET_RED1 = Achem( ii,MEHGDIS3D ) * Vrm
        MET_RED2 = 0.

!       C. oxic methylation
        Vmx = 1. - EXP( -ITS * kmx )
        METX_D = Achem( ii,HGDIS3D ) * Vmx !* FRED
        METX_S = Achem( ii,HGS_P3D ) * Vmx
        IF( Vmx .GT. 0.9 ) print*, 'ERROR 991: Vmx = ',Vmx ,ii

!       C. anaerob methylation
        IF( C_O2 .LE. 0.0 ) THEN
          Vm1 = 1. - EXP( -ITS * km1 ) 
        ELSE
          Vm1 = 0.
        END IF
        IF( Vm1 .GT. 0.9 ) print*, 'ERROR 731: Vm1 = ',Vm1, ii

        MET1_DIS = Achem( ii,HGDIS3D )  * Vm1
        MET1_HGS = Achem( ii,HGS_P3D )  * Vm1
        DISS_MET1 = Achem( ii,HGS_DOC3D )  * Vm1
        MET1_DOC = 0.
        MET1_POC = 0.

        IF( MET1_HGS < 0. ) THEN
          print*, 'Z7', MET1_HGS, 3.0-O2
          MET1_HGS = 0.
          IF( Achem( ii,HGS_P3D ) < 0. ) THEN
            print*, 'Y8', Achem( ii,HGS_P3D )
            Achem( ii,HGS_P3D ) = 0.
          END IF
        END IF

!       C. biogenic methylation
        RDET = 0.006 
     &       * (1 + 20 * ( TEMP( ii )**2 / (13**2 + TEMP( ii )**2) ) )  ! [day-1]
        RDET = RDET / SEC_PER_DAY                                       ! [s-1]
        RDOC = 10. * RDET                                               ! [s-1]

        Vdoc = 1. - EXP( -ITS * RDOC )
        TDREM = Achem( ii,DTOM3D ) * Vdoc
        Achem( ii,DTOM3D ) = Achem( ii,DTOM3D ) - TDREM
        dDOC = Vdoc * Abio( ii,DOM )
        Vdet = 1. - EXP( -ITS * RDET )
        TPREM = Achem( ii,PTOM3D ) * Vdet
        Achem( ii,PTOM3D ) = Achem( ii,PTOM3D ) - TPREM
        dDET = Vdet * Abio( ii,DET )
        FREM = dDOC + dDET

        TTOM = Achem( ii,PTOM3D ) + Achem( ii,DTOM3D )                 ! [mgC m-3]

        IF( Achem( ii,PTOM3D ) .LT. 0. ) THEN
          print*, 'ERROR 081: PTOM ', Achem( ii,PTOM3D )
          Achem( ii,PTOM3D ) = CMIN
        END IF
        IF( Achem( ii,DTOM3D ) .LT. 0. ) THEN
          print*, 'ERROR 082: DTOM ', Achem( ii,DTOM3D )
          Achem( ii,DTOM3D ) = CMIN
        END IF

        DMOM = 0.                             ! DOM fraction available and not inhibited by TOM
        IF( Abio( ii,DOM ) .GT. 0. ) THEN
          DMOM = MAX( 0., Abio( ii,DOM ) - 5.0 * TTOM ) / Abio( ii,DOM )! Soerensen et al., 2017    ! [] [0-1]
        END IF

!       disable methylation inhibition
        DMOM = 1.0

!       print*, 'Vmb = ', kmb, FREM, dDET, dDOC, Abio( ii,DET ), Abio( ii,DOM), TEMP(ii), RDET
        Vmb = 1. - EXP( -ITS * kmb * FREM )
        IF( Vmb .GT. 0.9 ) print*, 'ERROR 732: Vmb = ',Vmb, ii
        METB_DIS = Achem( ii,HGDIS3D ) * Vmb * DMOM ! * FRED
        METB_HGS = Achem( ii,HGS_P3D ) * Vmb * DMOM
        METB_POC = 0.
        METB_DOC = Achem( ii,HGDOC3D ) * Vmb * DMOM
        DISS_METB = 0.

!       C. di-methylation
        Vm2 = 1. - EXP( -ITS * km2 )
        IF( Vm2 .GT. 0.9 ) print*, 'ERROR 735: Vm2 = ',Vm2, ii
        MET2_DIS = Achem( ii,HGDIS3D ) * Vm2
        MET2_HGS = Achem( ii,HGS_P3D ) * Vm2
        MET2_DOC = 0.
        MET2_POC = 0.
        DISS_MET2 = 0.

!       C. bi-methylation
        Vm3 = 1. - EXP( -ITS * km3 )
        IF( Vm3 .GT. 0.9 ) print*, 'ERROR 736: Vm3 = ',Vm3, ii
        MET3_DIS = Achem( ii,MEHGDIS3D ) * Vm3
        MET3_DOC = 0.
        MET3_POC = 0.

!       C. de-methylation
        kdm_T = kdm
        Vdm = 1. - EXP( -ITS * kdm_T )
        IF( Vdm .GT. 0.9 ) print*, 'ERROR 737: Vdm = ',Vdm, ii
        DEM_DIS = Achem( ii,MEHGDIS3D ) * Vdm
        DEM_DOC = 0.
        DEM_POC = 0.

!       C. de-di-methylation
        Vdd = 1. - EXP( -ITS * kdd )
        IF( Vdd .GT. 0.9 ) print*, 'ERROR 738: Vdd = ',Vdd, ii
        DED = Achem( ii,DMHG3D ) * Vdd

!       C. photo de-methylation
!       Vpd = 1. - EXP( -ITS * kpd * RADI( ii ) * RADI_TO_PAR )
        Vpd = 1. - EXP( -ITS * kpd * RADI( ii ) )
        IF( Vpd .GT. 0.9 ) print*, 'ERROR 739: Vpd = ',Vpd, ii
        DEMP_DIS = Achem( ii,MEHGDIS3D ) * Vpd
        DEMP_DOC = 0.
        DEP      = Achem( ii,DMHG3D )    * Vpd
        DEP2     = DEP

!       C.3 sulfate chemistry
!       Set fraction of non-labile Hg-DOC
        IF( Abio( ii,DOM ) .LE. 1.0 ) THEN
          LABFRAC = 1.0 - 0.04
        ELSE
          LABFRAC = 1.0 - 0.20
        END IF

!       Distinguish reductive and oxidative regimes
!       C. cinnabar oxidation
        IF( C_O2 .GT. 0.0 ) THEN
          DIS_FORM = 0.0
          DOC_FORM = 0.0

          IF( C_O2 .GT. 1.0 ) C_O2 = 1.0

          Vhgo = 1. - EXP( -ITS * kso * C_O2 )
          IF( Vhgo .GT. 0.9 ) print*, 'ERROR 740: Vhgo = ',Vhgo, ii
          FORM_OXI = Achem( ii,HGS_P3D ) * Vhgo

          Vhgo = 1. - EXP( -ITS * kso * C_O2 ) 
          IF( Vhgo .GT. 0.9 ) print*, 'ERROR 741: Vhgo = ',Vhgo, ii
          DISS_OXI = Achem( ii,HGS_DOC3D ) * Vhgo * LABFRAC
!       C. cinnabar formation
        ELSE
          C_H2S = 0.5 * ABS( C_O2 )            
          IF( C_H2S .GT. 1.0 ) C_H2S = 1.0

          Vhgf = 1. - EXP( -ITS * ksf * C_H2S )
          IF( Vhgf .GT. 0.9 ) print*, 'ERROR 742: Vhgf = ',Vhgf, ii
          DIS_FORM = Achem( ii,HGDIS3D )  * Vhgf
          DOC_FORM = Achem( ii,HGDOC3D )  * Vhgf * LABFRAC

          FORM_OXI = 0.0
          DISS_OXI = 0.0
        END IF

!       C. cinnabar dissolution
        Vhgd = 1. - EXP( -ITS * ksd )
        IF( Vhgd .GT. 0.9 ) print*, 'ERROR 743: Vhgd = ',Vhgd, ii
        FORM_DISS = Achem( ii,HGS_P3D ) * Vhgd

!       C. re-crystallisation
        Vhga = 1. - EXP( -ITS * ksa )
        IF( Vhga .GT. 0.9 ) print*, 'ERROR 744: Vhga = ',Vhga, ii
        DISS_ADS = Achem( ii,HGS_DOC3D ) * Vhga * LABFRAC

C...... Apply Tendencies to chemistry fields
!       Inorganic chemistry
        Achem( ii,HGDEM3D ) =
     &  Achem( ii,HGDEM3D ) + REDB_DIS + REDB_DOC + REDB_POC + DISS_REDB
     &                      + REDP_DIS + REDP_DOC + REDP_POC + DISS_REDP
     &                      + REDD_DIS + REDD_DOC
     &                      - OXI_D    - OXI_P

        Achem( ii,HGDIS3D ) =
     &  Achem( ii,HGDIS3D ) + OXI_D    + OXI_P    + FORM_OXI
     &                      - REDB_DIS - REDP_DIS - DIS_FORM - REDD_DIS
     &                      - METX_D

        Achem( ii,HGDOC3D ) =
     &  Achem( ii,HGDOC3D ) + DISS_OXI - DOC_FORM
     &                      - REDB_DOC - REDP_DOC - REDD_DOC

        Achem( ii,HGPOC3D ) =
     &  Achem( ii,HGPOC3D ) - REDB_POC - REDP_POC

        Achem( ii,HGS_P3D ) =
     &  Achem( ii,HGS_P3D ) + DIS_FORM + DOC_FORM + DISS_ADS
     &                      - FORM_OXI - FORM_DISS
     &                      - METX_S

        Achem(ii,HGS_DOC3D) =
     &  Achem(ii,HGS_DOC3D) + FORM_DISS
     &                      - DISS_OXI  - DISS_ADS
     &                      - DISS_REDB - DISS_REDP

!       MET1_HGS = 0.

!       Organic chemistry
        Achem( ii,HGDEM3D ) =
     &  Achem( ii,HGDEM3D ) + MET_RED1 + MET_RED2

        Achem( ii,HGDIS3D ) =
     &  Achem( ii,HGDIS3D ) + DEM_DIS  + DEM_DOC  + DEMP_DIS + DEMP_DOC
     &                      + DEP2
     &                      - MET1_DIS - MET2_DIS - METB_DIS - MET_RED1

        Achem( ii,HGDOC3D ) = 
     &  Achem( ii,HGDOC3D ) - MET1_DOC - MET2_DOC - METB_DOC

        Achem( ii,HGPOC3D ) =
     &  Achem( ii,HGPOC3D ) + DEM_POC
     &                      - MET1_POC - MET2_POC - METB_POC

        Achem(ii,MEHGDIS3D) =
     &  Achem(ii,MEHGDIS3D) + MET1_DIS + DISS_MET1 + METB_HGS
     &                      + METB_DIS + DISS_METB + DED
     &                      - MET3_DIS - DEM_DIS   - DEMP_DIS
     &                      + METX_D   + METX_S    + DEP

        Achem(ii,HGS_DOC3D) =
     &  Achem(ii,HGS_DOC3D) - DISS_MET1 - DISS_METB - DISS_MET2

        Achem( ii,HGS_P3D ) =
     &  Achem( ii,HGS_P3D ) - MET1_HGS - METB_HGS - MET2_HGS

        Achem(ii,MEHGDOC3D) =
     &  Achem(ii,MEHGDOC3D) + MET1_DOC + METB_DOC
     &                      - MET3_DOC - DEM_DOC  - DEMP_DOC - MET_RED2

        Achem(ii,MEHGPOC3D) =
     &  Achem(ii,MEHGPOC3D) + MET1_POC + METB_POC + MET1_HGS 
     &                      - MET3_POC - DEM_POC 

        Achem( ii,DMHG3D )  =
     &  Achem( ii,DMHG3D )  + MET2_DIS + MET2_DOC + MET2_POC
     &                      + MET3_DIS + MET3_DOC + MET3_POC + DISS_MET2
     &                      + MET2_HGS
     &                      - DED      - DEP      - DEP2

!       terrestrial organic carbon deterioration

!       C. biogenic methylation
        RDET = 0.006
     &       * (1 + 20 * ( TEMP( ii )**2 / (13**2 + TEMP( ii )**2) ) ) ! [day-1]
        RDET = 0.1 * RDET / SEC_PER_DAY                                ! [s-1]
        RDOC = 10. * RDET                                              ! [s-1]

        Vpam = 1. - EXP( -ITS * RDET )
        Achem( ii,PTOM3D ) = Achem( ii,PTOM3D )
     &                     - Achem( ii,PTOM3D ) * Vpam
        Vpam = 1. - EXP( -ITS * RDOC )
        Achem( ii,DTOM3D ) = Achem( ii,DTOM3D )
     &                     - Achem( ii,DTOM3D ) * Vpam

!...... check minimum concentration and mass conservation
        FNORM = 0.
        DO NNC = 1,NHGTRANS
          IF( Achem( ii,NNC ) .LT. 0.0 ) THEN
            print*, 'ERROR 801: ', VNAMOUT( NNC ), Achem( ii,NNC )
          END IF
          Achem( ii,NNC ) = MAX( CMIN,Achem( ii,NNC ) )
          FNORM = FNORM + Achem( ii,NNC )
        END DO

        IF( FNORM .GT. 0 ) THEN
            FNORM = HGTOT / FNORM
        ELSE
            print*, 'ERROR 808: negative FNORM ', FNORM
        END IF

        IF( FNORM .GT. 1.01 .OR. FNORM .LT. 0.99 ) THEN
            print*, 'WARNING 123: chem.f(624) FNORM=',FNORM
        END IF

        DO NNC = 1, nchem
          Achem( ii,NNC ) = Achem( ii,NNC ) * FNORM
        END DO

C.........  End of subroutine: 
        END SUBROUTINE CHEM
