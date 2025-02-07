	MODULE COUT

	USE CPARAM
	USE CUTIL

        IMPLICIT NONE

C...........
!
! Global arrays and subroutines for chemistry output
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........


!        include 'C_model_inc'

C...........   Mercury species in water
        INTEGER, PARAMETER :: HGDEM3D   =  1      ! Dissolved Elemental Mercury          (Hg0)
        INTEGER, PARAMETER :: HGDIS3D   =  2      ! Dissolved oxidized mercury           (Hg2+)
        INTEGER, PARAMETER :: HGSTAR3D  =  3      ! Dissolved fortified oxidized mercury (Hg+)
        INTEGER, PARAMETER :: HGDOC3D   =  4      ! Oxidized Mercury associated to dissolved organic matter
        INTEGER, PARAMETER :: HGPOC3D   =  5      ! Particulate Mercury associated to Organic Carbon
        INTEGER, PARAMETER :: HGS_P3D   =  6      ! Cinnabar HgS-alpha & Metacinnabar Hgs-beta particulate
        INTEGER, PARAMETER :: HGS_DOC3D =  7      ! HgS associated to dossolved oganic matter

        INTEGER, PARAMETER :: MEHGDIS3D =  8      ! Dissolved methylmercury (CH3Hg)
        INTEGER, PARAMETER :: MEHGDOC3D =  9      ! Dissolved Methylmercury associated to Organic Carbon
        INTEGER, PARAMETER :: MEHGPOC3D =  10     ! Particulate Methylmercury associated to Organic Carbon
        INTEGER, PARAMETER :: DMHG3D    =  11     ! Di-methyl-mercury ((CH3)2Hg

C...........   Mercury species in biota
        INTEGER, PARAMETER :: HGDET3D   = 12      ! Mercury accumulated in Detritus
        INTEGER, PARAMETER :: HGFLA3D   = 13      ! Mercury bio-accumlation in flaggelates
        INTEGER, PARAMETER :: HGDIA3D   = 14      ! Mercury bio-accumlation in diatomes
        INTEGER, PARAMETER :: HGCYA3D   = 15      ! Mercury bio-accumlation in cyanobacteria
        INTEGER, PARAMETER :: HGZOSe3D  = 16      ! Mercury bio-accumlation on small Zooplankton
        INTEGER, PARAMETER :: HGZOLe3D  = 17      ! Mercury bio-accumlation on large Zooplankton
        INTEGER, PARAMETER :: HGZOSi3D  = 18      ! Mercury bio-accumlation in small Zooplankton
        INTEGER, PARAMETER :: HGZOLi3D  = 19      ! Mercury bio-accumlation in large Zooplankton

        INTEGER, PARAMETER :: MEHGDET3D = 20      ! Methylmercury accumulated in Detritus
        INTEGER, PARAMETER :: MEHGFLA3D = 21      ! Methylmercury bio-accumlation in flaggelates
        INTEGER, PARAMETER :: MEHGDIA3D = 22      ! Methylmercury bio-accumlation in diatomes
        INTEGER, PARAMETER :: MEHGCYA3D = 23      ! Methylmercury bio-accumlation in cyanobacteria
        INTEGER, PARAMETER :: MEHGZOSe3D = 24      ! Methylmercury bio-accumlation on small Zooplankton
        INTEGER, PARAMETER :: MEHGZOLe3D = 25      ! Methylmercury bio-accumlation on large Zooplankton
        INTEGER, PARAMETER :: MEHGZOSi3D = 26      ! Methylmercury bio-accumlation in small Zooplankton
        INTEGER, PARAMETER :: MEHGZOLi3D = 27      ! Methylmercury bio-accumlation in large Zooplankton

        INTEGER, PARAMETER :: NHGTRANS  = 27      ! C_model_inc IHGT = NHGTRANS

        INTEGER, PARAMETER :: PTOM3D       = NHGTRANS + 1
        INTEGER, PARAMETER :: DTOM3D       = NHGTRANS + 2 
        INTEGER, PARAMETER :: HGMACe3D     = NHGTRANS + 3      ! Mercury bio-accumulation on macro benthos
        INTEGER, PARAMETER :: HGMACi3D     = NHGTRANS + 4      ! Mercury bio-accumulation in macro benthos
        INTEGER, PARAMETER :: MEHGMACe3D   = NHGTRANS + 5      ! Methylmercury bio-accumulation on macro benthos
        INTEGER, PARAMETER :: MEHGMACi3D   = NHGTRANS + 6      ! Methylmercury bio-accumulation in macro benthos
        INTEGER, PARAMETER :: SEDHG3D      = NHGTRANS + 7      ! ng/m²
        INTEGER, PARAMETER :: SEDMEHG3D    = NHGTRANS + 8      ! ng/m²
        INTEGER, PARAMETER :: NON_HG_TRANS            = 8 ! C_model_inc NCHEM = NHGTRANS + NON_HG_TRANS

!nchem and imehg in C_model_inc
!        INTEGER, PARAMETER :: NHG = 23         ! Last index for Hg total
!        INTEGER, PARAMETER :: IMEHG = 13       ! starting index for MeHg total

!       non advected variables
        INTEGER, PARAMETER :: HGFISHi3D   = NHGTRANS + NON_HG_TRANS + 1
        INTEGER, PARAMETER :: MEHGFISHi3D = NHGTRANS + NON_HG_TRANS + 2 
        INTEGER, PARAMETER :: HGFISHe3D   = NHGTRANS + NON_HG_TRANS + 3
        INTEGER, PARAMETER :: MEHGFISHe3D = NHGTRANS + NON_HG_TRANS + 4
        INTEGER, PARAMETER :: MEHGTOT3D  = NHGTRANS + NON_HG_TRANS + 5 ! Total methylated mercury
        INTEGER, PARAMETER :: HGTOT3D    = NHGTRANS + NON_HG_TRANS + 6 ! Total mercury
        INTEGER, PARAMETER :: MEHGBIO3D  = NHGTRANS + NON_HG_TRANS + 7 ! Total MeHg in biota


        INTEGER, PARAMETER :: DOC3D = NHGTRANS + NON_HG_TRANS +  8
        INTEGER, PARAMETER :: POC3D = NHGTRANS + NON_HG_TRANS +  9
        INTEGER, PARAMETER :: BPC3D = NHGTRANS + NON_HG_TRANS + 10
        INTEGER, PARAMETER :: BZC3D = NHGTRANS + NON_HG_TRANS + 11
        INTEGER, PARAMETER :: FLA3D = NHGTRANS + NON_HG_TRANS + 12
        INTEGER, PARAMETER :: DIA3D = NHGTRANS + NON_HG_TRANS + 13
        INTEGER, PARAMETER :: CYA3D = NHGTRANS + NON_HG_TRANS + 14 
        INTEGER, PARAMETER :: ZOS3D = NHGTRANS + NON_HG_TRANS + 15
        INTEGER, PARAMETER :: ZOL3D = NHGTRANS + NON_HG_TRANS + 16
        INTEGER, PARAMETER :: DET3D = NHGTRANS + NON_HG_TRANS + 17
        INTEGER, PARAMETER :: FSH3D = NHGTRANS + NON_HG_TRANS + 18
        INTEGER, PARAMETER :: SAL3D = NHGTRANS + NON_HG_TRANS + 19
        INTEGER, PARAMETER :: OXY3D = NHGTRANS + NON_HG_TRANS + 20
!       INTEGER, PARAMETER :: PHO3D = NHGTRANS + NON_HG_TRANS + 19
        INTEGER, PARAMETER :: RAD3D = NHGTRANS + NON_HG_TRANS + 21
        INTEGER, PARAMETER :: TEM3D = NHGTRANS + NON_HG_TRANS + 22
        INTEGER, PARAMETER :: NON_TRANS                       = 22 !Number of variables that are not transported

        INTEGER, PARAMETER :: NHGOUT3D  = NHGTRANS + NON_HG_TRANS + 22 !number of 3D variables for output file (last one is TEM3D)

	INTEGER, PARAMETER :: NHGOUT2D  = 22	! number of 2D variables for output file
	INTEGER, PARAMETER :: HGSAX2D   =  1	! mercury flux water to air
	INTEGER, PARAMETER :: HGDEP2D   =  2	! atmospheric mercury deposition
	INTEGER, PARAMETER :: HGICE2D   =  3	! mercury in ice
	INTEGER, PARAMETER :: HGSED2D   =  4	! mercury in sediment
	INTEGER, PARAMETER :: HGBUR2D   =  5	! mercury burried in deep sediment
	INTEGER, PARAMETER :: HGRES2D   =  6	! mercury resuspended from sediment
	INTEGER, PARAMETER :: USTAR2D   =  7	! u* at ocean floor
	INTEGER, PARAMETER :: MEHGSED2D =  8	! methylmercury in sediment
	INTEGER, PARAMETER :: MEHGBUR2D =  9	! methylmercury burried in deep sediment
	INTEGER, PARAMETER :: MEHGRES2D = 10	! methylmercury resuspended from sediment
	INTEGER, PARAMETER :: LAYERS2D  = 11	! 
	INTEGER, PARAMETER :: U2D       = 12	! 
	INTEGER, PARAMETER :: V2D       = 13	! 
	INTEGER, PARAMETER :: FRICE2D   = 14	!
        INTEGER, PARAMETER :: STOT2D    = 15    ! net methylation in sediment
        INTEGER, PARAMETER :: MAC2D     = 16    ! macrobenthos
        INTEGER, PARAMETER :: HGMACe2D   = 17    ! Hg on macrobenthos
        INTEGER, PARAMETER :: HGMACi2D   = 18    ! Hg in macrobenthos
        INTEGER, PARAMETER :: MEHGMACe2D = 19    ! MeHg on macrobenthos
        INTEGER, PARAMETER :: MEHGMACi2D = 20    ! MeHg in macrobenthos
        INTEGER, PARAMETER :: HGFLU2D   = 21    ! Hg flux from lowest wet grid cell to sediment
        INTEGER, PARAMETER :: MEHGFLU2D = 22    ! MeHg flux from lowest wet grid cell to sediment

        CHARACTER( 16 ), PARAMETER :: VNAMOUT( NHGOUT3D+NHGOUT2D ) =
     &          ( / 'HGDEM'  ,'HGDIS'  ,'HGSTAR'  ,'HGDOC'  ,'HGPOC'   ,
     &              'HGS_P'  ,'HGS_DOC',
     &              'MEHGDIS','MEHGDOC','MEHGPOC' ,'DMHG'   ,           
     &    'HGDET'  ,'HGFLA'  ,'HGDIA'   ,'HGCYA'  ,'HGZOSe' ,'HGZOLe'  ,
     &                                             'HGZOSi' ,'HGZOLi'  ,
     &    'MEHGDET','MEHGFLA','MEHGDIA','MEHGCYA','MEHGZOSe','MEHGZOLe',
     &                                            'MEHGZOSi','MEHGZOLi',
     &    'PTOM'   ,'DTOM'    ,'HGMACe' ,'HGMACi','MEHGMACe','MEHGMACi',
     &                                   'SEDHG'  ,'SEDMEHG',
     &    'HGFISHi','MEHGFISHi'        ,'HGFISHe' ,'MEHGFISHe'         ,
     &              'MEHGTOT','HGTOT'  ,'MEHGBIO' ,
     &              'DOC'    ,'POC'    ,'BPC'     ,'BZC'    ,'FLA'     ,
     &              'DIA'    ,'CYA'    ,'ZOS'     ,'ZOL'    ,'DET'     ,
     &              'FSH'    ,'SALT'   ,'OXYGEN'  ,'RADI','TEMPERATURE',
     &              'HGSAX'  ,'HGDEP'  ,'HGICE'   ,'HGSED'  ,'HGBUR'   ,
     &              'HGRES'  ,'USTAR'  ,'MEHGSED' ,'MEHGBUR','MEHGRES' ,
     &              'LAYERS' ,'U'      ,'V'       ,'FRICE'  ,'STOT'    ,
     &    'MAC'    ,'HGMACe' ,'HGMACi' ,'MEHGMACe','MEHGMACi',
     &              'HGFLU2D','MEHGFLU2D' 
     &           / )

C...........   LOCAL VARIABLES and their descriptions:
!	REAL, ALLOCATABLE, DIMENSION( :,:,:,: ) :: HGOUT3D	! 3D mercury fields
!	REAL, ALLOCATABLE, DIMENSION( :,:,:,: ) :: HGOUT2D	! 2D mercury fields
       REAL, ALLOCATABLE, TARGET, DIMENSION( :,:,:,: ) :: HGOUT3D      ! 3D mercury fields
       REAL, ALLOCATABLE, TARGET, DIMENSION( :,:,:,: ) :: HGOUT2D      ! 2D mercury fields

	CONTAINS
	
	SUBROUTINE WRITE_CHEM( INITIAL )

	USE CPARAM

	IMPLICIT NONE

C...........
!
! Global arrays and subroutines for chemistry output
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)
	INCLUDE 'C_model_inc'	!  Model domain dimension size

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(16)      PROMPTMFILE
        LOGICAL            SETENVVAR

	EXTERNAL  PROMPTMFILE, SETENVVAR

C...........   PARAMETERS and their descriptions:
        CHARACTER( 16 ), PARAMETER  :: PROGNAME = 'WRITE_CHEM' ! program name

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
	LOGICAL, INTENT(IN) ::	INITIAL		! initial call open file new

C...........   LOCAL VARIABLES and their descriptions:
	INTEGER			I,V,X,R,C,L
	REAL, DIMENSION( n,m,ilo ) :: TMP	! helper array

	REAL, PARAMETER	     :: UFAC = 1.	! 1000 to write ng/m³ instead of ng/L

	REAL		     :: COUNTER = 0	! time step counter to calculate average
	REAL		     :: VGLVS = 0	! layer information
	INTEGER			WDATE		! write date
	INTEGER			WTIME           ! write time
	INTEGER			WSTEP2D, WSTEP3D, XSTEP	! write time step
	LOGICAL		     :: WFLAG = .FALSE.	! write flag (.FALSE. add data, .TRUE. write data)

	INTEGER			WYEAR		! Year for zearly output files
	CHARACTER( 16 )		O3NAME 		! 3D output file name
	CHARACTER( 16 )		O2NAME 		! 2D output file name
	CHARACTER( 256 )	CHEMOUT3D	! logical name for O3NAME
	CHARACTER( 256 )	CHEMOUT2D	! locical name for O2NAME
	CHARACTER( 256 )	MESG		! message buffer
        INTEGER       		IOS		! i/o status
        REAL TOT

C***********************************************************************
C   begin body of subroutine WRITE_CHEM
	WYEAR = CDATE / 1000
	WDATE = CDATE
	WTIME = CTIME
!	WSTEP3D = 120000
        WSTEP3D = 240000

	IF( WTIME .EQ. 0 ) THEN	! write daily averages after the 24th time step has been added
	    WFLAG = .TRUE.
	ELSE
	    WFLAG = .FALSE.
	END IF

!for hourly output
        WSTEP3D = 240000
!        WSTEP3D =  10000
        IF( MOD( WTIME,WSTEP3D ) .EQ. 0 ) THEN
            WFLAG = .TRUE.
        ELSE
            WFLAG = .FALSE.
        END IF

! 25.03.2016 temporary change to hourly output
!        WSTEP2D =  10000
       WSTEP2D = 240000
!
        IF( WFLAG ) THEN
            COUNTER = 1.0 / 24. !!! for daily average output
        END IF

!for hourly output              !!!! IMPORTANT !!!!
!!!     COUNTER = 1.0

!        WRITE( MESG,* ) 'COUNTER CHECK: ', COUNTER, WDATE, WTIME, WFLAG
!        CALL M3MSG2( MESG )

        IF( WFLAG .XOR. INITIAL ) THEN
          WRITE( MESG,* ) 'WRITE_CHEM: ', WDATE, WTIME, WSTEP3D, WSTEP2D, WFLAG, INITIAL
          CALL M3MSG2( MESG )
        END IF

	IF( INITIAL ) THEN
	    WFLAG = .FALSE.
	    COUNTER = 0.

	    CALL M3MSG2( DASHLINE )
	    WRITE( MESG,* ) 'Creating ocean chemistry output file for ', WYEAR
	    CALL M3MSG2( MESG )
	    CALL M3MSG2( '' )

	    FTYPE3D = GRDDED3
	    GDTYP3D = 1
	    GDNAM3D = 'ECOSMO_NSBS'
	    FDESC3D = 'ECOSMO_Hg 3D chemistry output file'

	    SDATE3D = WDATE
	    STIME3D = WTIME
	    TSTEP3D = WSTEP3D

	    NCOLS3D = n
	    NROWS3D = m
	    NLAYS3D = ilo
	    NTHIK3D = 1

	    VGTYP3D = IMISS3
	    VGTOP3D = 0
	    DO I = 1,ilo+1
	        VGLVS3D( I ) = VGLVS
		VGLVS = VGLVS + 0.01
	    END DO

	    NVARS3D = NHGOUT3D
	    VNAME3D( 1 ) = 'HGDEM'
	    UNITS3D( 1 ) = 'ng/L'
	    VDESC3D( 1 ) = 'elemental mercury (Hg0)'
	    VTYPE3D( 1 ) = M3REAL
	    VNAME3D( 2 ) = 'HGDIS'
	    UNITS3D( 2 ) = 'ng/L'
	    VDESC3D( 2 ) = 'dissolved oxidized mercury (Hg2+)'
	    VTYPE3D( 2 ) = M3REAL
            VNAME3D( 3 ) = 'HGSTAR'
            UNITS3D( 3 ) = 'ng/L'
            VDESC3D( 3 ) = 'dissolved non-reducible oxidized (Hg+)'
            VTYPE3D( 3 ) = M3REAL
            VNAME3D( 4 ) = 'HGDOC'
            UNITS3D( 4 ) = 'ng/L'
            VDESC3D( 4 ) = 'Hg2+ complexed to DOM'
            VTYPE3D( 4 ) = M3REAL
	    VNAME3D( 5 ) = 'HGPOC'
	    UNITS3D( 5 ) = 'ng/L'
	    VDESC3D( 5 ) = 'particulate oxidized mercury'
	    VTYPE3D( 5 ) = M3REAL
            VNAME3D( 6 ) = 'HGS_P'
            UNITS3D( 6 ) = 'ng/L'
            VTYPE3D( 6 ) = M3REAL
            VDESC3D( 6 ) =
     &             'Particulate (meta-)cinnabar HgS-alpha & Hg-Sbeta'
            VNAME3D( 7 ) = 'HGS_DOC'
            UNITS3D( 7 ) = 'ng/L'
            VTYPE3D( 7 ) = M3REAL
            VDESC3D( 7 ) = 'HgS-DOM complex dissolved'
            VNAME3D( 8 ) = 'MEHGDIS'
            UNITS3D( 8 ) = 'ng/L'
            VDESC3D( 8 ) = 'dissolved methylmercury (CH3Hg)'
            VTYPE3D( 8 ) = M3REAL
	    VNAME3D( 9 ) = 'MEHGDOC'
	    UNITS3D( 9 ) = 'ng/L'
	    VDESC3D( 9 ) = 'methylmercury associated to DOM'
	    VTYPE3D( 9 ) = M3REAL
	    VNAME3D( 10 ) = 'MEHGPOC'
	    UNITS3D( 10 ) = 'ng/L'
	    VDESC3D( 10 ) = 'particulate methylmercury'
	    VTYPE3D( 10 ) = M3REAL
            VNAME3D( 11 ) = 'DMHG'
            UNITS3D( 11 ) = 'ng/L'
            VDESC3D( 11 ) = 'dissolved dimethylmercury ((CH3)2Hg)'
            VTYPE3D( 11 ) = M3REAL
            VNAME3D( 12 ) = 'HGDET'
            UNITS3D( 12 ) = 'ng/L'
            VDESC3D( 12 ) = 'mercury accumulated in detritus'
            VTYPE3D( 12 ) = M3REAL
            VNAME3D( 13 ) = 'HGFLA'
            UNITS3D( 13 ) = 'ng/L'
            VDESC3D( 13 ) = 'mercury accumulated in flaggelates'
            VTYPE3D( 13 ) = M3REAL
            VNAME3D( 14 ) = 'HGDIA'
            UNITS3D( 14 ) = 'ng/L'
            VDESC3D( 14 ) = 'mercury accumulated in diatomes'
            VTYPE3D( 14 ) = M3REAL
            VNAME3D( 15 ) = 'HGCYA'
            UNITS3D( 15 ) = 'ng/L'
            VDESC3D( 15 ) = 'mercury accumulated in cyanobacteria'
            VTYPE3D( 15 ) = M3REAL
            VNAME3D( 16 ) = 'HGZOSe'
            UNITS3D( 16 ) = 'ng/L'
            VDESC3D( 16 ) = 'mercury accumulated on small zoop.'
            VTYPE3D( 16 ) = M3REAL
            VNAME3D( 17 ) = 'HGZOLe'
            UNITS3D( 17 ) = 'ng/L'
            VDESC3D( 17 ) = 'mercury accumulated on large zoop.'
            VTYPE3D( 17 ) = M3REAL
            VNAME3D( 18 ) = 'HGZOSi'
            UNITS3D( 18 ) = 'ng/L'
            VDESC3D( 18 ) = 'mercury accumulated in small zoop.'
            VTYPE3D( 18 ) = M3REAL
            VNAME3D( 19 ) = 'HGZOLi'
            UNITS3D( 19 ) = 'ng/L'
            VDESC3D( 19 ) = 'mercury accumulated in large zoop.'
            VTYPE3D( 19 ) = M3REAL
            VNAME3D( 20 ) = 'MEHGDET'
            UNITS3D( 20 ) = 'ng/L'
            VDESC3D( 20 ) = 'methylmercury accumulated in detritus'
            VTYPE3D( 20 ) = M3REAL
            VNAME3D( 21 ) = 'MEHGFLA'
            UNITS3D( 21 ) = 'ng/L'
            VDESC3D( 21 ) = 'methylmercury accumulated in flaggelates'
            VTYPE3D( 21 ) = M3REAL
            VNAME3D( 22 ) = 'MEHGDIA'
            UNITS3D( 22 ) = 'ng/L'
            VDESC3D( 22 ) = 'methylmercury accumulated in diatomes'
            VTYPE3D( 22 ) = M3REAL
            VNAME3D( 23 ) = 'MEHGCYA'
            UNITS3D( 23 ) = 'ng/L'
            VDESC3D( 23 ) = 'methylmercury accumulated in cyanobacteria'
            VTYPE3D( 23 ) = M3REAL
            VNAME3D( 24 ) = 'MEHGZOSe'
            UNITS3D( 24 ) = 'ng/L'
            VDESC3D( 24 ) = 'methylmercury accumulated on small zoop.'
            VTYPE3D( 24 ) = M3REAL
            VNAME3D( 25 ) = 'MEHGZOLe'
            UNITS3D( 25 ) = 'ng/L'
            VDESC3D( 25 ) = 'methylmercury accumulated on large zoop.'
            VTYPE3D( 25 ) = M3REAL
            VNAME3D( 26 ) = 'MEHGZOSi'
            UNITS3D( 26 ) = 'ng/L'
            VDESC3D( 26 ) = 'methylmercury accumulated in small zoop.'
            VTYPE3D( 26 ) = M3REAL
            VNAME3D( 27 ) = 'MEHGZOLi'
            UNITS3D( 27 ) = 'ng/L'
            VDESC3D( 27 ) = 'methylmercury accumulated in large zoop.'
            VTYPE3D( 27 ) = M3REAL

!       non Hg advected variables
            VNAME3D( NHGTRANS +  1 ) = 'PTOM'
            UNITS3D( NHGTRANS +  1 ) = 'mgC/L'
            VDESC3D( NHGTRANS +  1 ) = 'terrestrial part. org. matter'
            VTYPE3D( NHGTRANS +  1 ) = M3REAL
            VNAME3D( NHGTRANS +  2 ) = 'DTOM'
            UNITS3D( NHGTRANS +  2 ) = 'mgC/L'
            VDESC3D( NHGTRANS +  2 ) = 'terrestrial diss. org. matter'
            VTYPE3D( NHGTRANS +  2 ) = M3REAL
            VNAME3D( NHGTRANS +  3 ) = 'HGMACe'
            UNITS3D( NHGTRANS +  3 ) = 'ng/m**2'
            VDESC3D( NHGTRANS +  3 ) = 'mercury on macrobenthos'
            VTYPE3D( NHGTRANS +  3 ) = M3REAL
            VNAME3D( NHGTRANS +  4 ) = 'HGMACi'
            UNITS3D( NHGTRANS +  4 ) = 'ng/m**2'
            VDESC3D( NHGTRANS +  4 ) = 'mercury in macrobenthos'
            VTYPE3D( NHGTRANS +  4 ) = M3REAL
            VNAME3D( NHGTRANS +  5 ) = 'MEHGMACe'
            UNITS3D( NHGTRANS +  5 ) = 'ng/m**2'
            VDESC3D( NHGTRANS +  5 ) = 'methylmercury on macrobenthos'
            VTYPE3D( NHGTRANS +  5 ) = M3REAL
            VNAME3D( NHGTRANS +  6 ) = 'MEHGMACi'
            UNITS3D( NHGTRANS +  6 ) = 'ng/m**2'
            VDESC3D( NHGTRANS +  6 ) = 'methylmercury in macrobenthos'
            VTYPE3D( NHGTRANS +  6 ) = M3REAL
            VNAME3D( NHGTRANS +  7 ) = 'SEDHG'
            UNITS3D( NHGTRANS +  7 ) = 'ng/mgC'
            VDESC3D( NHGTRANS +  7 ) = 'Mercury in sediment'
            VTYPE3D( NHGTRANS +  7 ) = M3REAL
            VNAME3D( NHGTRANS +  8 ) = 'SEDMEHG'
            UNITS3D( NHGTRANS +  8 ) = 'ng/mgC'
            VDESC3D( NHGTRANS +  8 ) = 'methylmercury in sediment'
            VTYPE3D( NHGTRANS +  8 ) = M3REAL

!       non advected totals
            VNAME3D( NHGTRANS + NON_HG_TRANS +  1 ) = 'HGFISHi'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  1 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  1 ) = 'inorganic mercury in fish'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  1 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS +  2 ) = 'MEHGFISHi'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  2 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  2 ) = 'methyl mercury in fish'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  2 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS +  3 ) = 'HGFISHe'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  3 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  3 ) = 'inorganic mercury on fish'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  3 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS +  4 ) = 'MEHGFISHe'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  4 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  4 ) = 'methyl mercury on fish'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  4 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS +  5 ) = 'MEHGTOT'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  5 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  5 ) = 'particulate methylmercury'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  5 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS +  6 ) = 'HGTOT'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  6 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  6 ) = 'total mercury'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  6 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS +  7 ) = 'MEHGBIO'
            UNITS3D( NHGTRANS + NON_HG_TRANS +  7 ) = 'ng/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS +  7 ) = 'total methylmercury in biota'
            VTYPE3D( NHGTRANS + NON_HG_TRANS +  7 ) = M3REAL

!       non advected variables
	    VNAME3D( NHGTRANS + NON_HG_TRANS +  8 ) = 'DOC'
	    UNITS3D( NHGTRANS + NON_HG_TRANS +  8 ) = 'mg/L'
	    VDESC3D( NHGTRANS + NON_HG_TRANS +  8 ) = 'dissolved organic matter'
	    VTYPE3D( NHGTRANS + NON_HG_TRANS +  8 ) = M3REAL
	    VNAME3D( NHGTRANS + NON_HG_TRANS +  9 ) = 'POC'
	    UNITS3D( NHGTRANS + NON_HG_TRANS +  9 ) = 'mg/L'
	    VDESC3D( NHGTRANS + NON_HG_TRANS +  9 ) = 'particulate organic matter'
	    VTYPE3D( NHGTRANS + NON_HG_TRANS +  9 ) = M3REAL
	    VNAME3D( NHGTRANS + NON_HG_TRANS + 10 ) = 'BPC'
	    UNITS3D( NHGTRANS + NON_HG_TRANS + 10 ) = 'mg/L'
	    VDESC3D( NHGTRANS + NON_HG_TRANS + 10 ) = 'phytoplankton'
	    VTYPE3D( NHGTRANS + NON_HG_TRANS + 10 ) = M3REAL
	    VNAME3D( NHGTRANS + NON_HG_TRANS + 11 ) = 'BZC'
	    UNITS3D( NHGTRANS + NON_HG_TRANS + 11 ) = 'mg/L'
	    VDESC3D( NHGTRANS + NON_HG_TRANS + 11 ) = 'zooplankton'
	    VTYPE3D( NHGTRANS + NON_HG_TRANS + 11 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 12 ) = 'FLA'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 12 ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 12 ) = 'flaggelates'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 12 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 13 ) = 'DIA'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 13 ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 13 ) = 'diatomes'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 13 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 14 ) = 'CYA'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 14 ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 14 ) = 'cyanobacteria'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 14 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 15 ) = 'ZOS'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 15 ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 15 ) = 'small zooplankton'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 15 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 16 ) = 'ZOL'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 16 ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 16 ) = 'large zooplankton'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 16 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 17  ) = 'DET'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 17  ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 17  ) = 'detritus'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 17  ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 18 ) = 'FSH'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 18 ) = 'mgC/m**3'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 18 ) = 'fish biomass'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 18 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 19 ) = 'SALT'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 19 ) = 'PSU'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 19 ) = 'salinity'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 19 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 20 ) = 'OXYGEN'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 20 ) = 'mL/L'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 20 ) = 'oxygen'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 20 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 21 ) = 'RADI'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 21 ) = 'W/m**2'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 21 ) = 'radiation'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 21 ) = M3REAL
            VNAME3D( NHGTRANS + NON_HG_TRANS + 22 ) = 'TEMPERATURE'
            UNITS3D( NHGTRANS + NON_HG_TRANS + 22 ) = '°C'
            VDESC3D( NHGTRANS + NON_HG_TRANS + 22 ) = 'water temp'
            VTYPE3D( NHGTRANS + NON_HG_TRANS + 22 ) = M3REAL

            VNAME3D( NHGOUT3D+ 1 ) = 'HGSAX'
	    UNITS3D( NHGOUT3D+ 1 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 1 ) = 'mercury flux water to air'
	    VTYPE3D( NHGOUT3D+ 1 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 2 ) = 'HGDEP'
	    UNITS3D( NHGOUT3D+ 2 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 2 ) = 'mercury flux water to air'
	    VTYPE3D( NHGOUT3D+ 2 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 3 ) = 'HGICE'
	    UNITS3D( NHGOUT3D+ 3 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 3 ) = 'mercury in ice'
	    VTYPE3D( NHGOUT3D+ 3 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 4 ) = 'HGSED'
	    UNITS3D( NHGOUT3D+ 4 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 4 ) = 'mercury sedimentation'
	    VTYPE3D( NHGOUT3D+ 4 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 5 ) = 'HGBUR'
	    UNITS3D( NHGOUT3D+ 5 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 5 ) = 'mercury in sediment'
	    VTYPE3D( NHGOUT3D+ 5 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 6 ) = 'HGRES'
	    UNITS3D( NHGOUT3D+ 6 ) = 'ng/L'
	    VDESC3D( NHGOUT3D+ 6 ) = 'mercury resupension'
	    VTYPE3D( NHGOUT3D+ 6 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 7 ) = 'USTAR'
	    UNITS3D( NHGOUT3D+ 7 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 7 ) = 'ustar at ocean floor'
	    VTYPE3D( NHGOUT3D+ 7 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 8 ) = 'MEHGSED'
	    UNITS3D( NHGOUT3D+ 8 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 8 ) = 'methylmercury in sediment'
	    VTYPE3D( NHGOUT3D+ 8 ) = M3REAL
	    VNAME3D( NHGOUT3D+ 9 ) = 'MEHGBUR'
	    UNITS3D( NHGOUT3D+ 9 ) = 'ng/m**2'
	    VDESC3D( NHGOUT3D+ 9 ) = 'methylmercury in sediment'
	    VTYPE3D( NHGOUT3D+ 9 ) = M3REAL
	    VNAME3D( NHGOUT3D+10 ) = 'MEHGRES'
	    UNITS3D( NHGOUT3D+10 ) = 'ng/L'
	    VDESC3D( NHGOUT3D+10 ) = 'methylmercury resuspended'
	    VTYPE3D( NHGOUT3D+10 ) = M3REAL
	    VNAME3D( NHGOUT3D+11 ) = 'LAYERS'
	    UNITS3D( NHGOUT3D+11 ) = 'none'
	    VDESC3D( NHGOUT3D+11 ) = 'number of vertical layers'
	    VTYPE3D( NHGOUT3D+11 ) = M3REAL
	    VNAME3D( NHGOUT3D+12 ) = 'U'
	    UNITS3D( NHGOUT3D+12 ) = 'm/s'
	    VDESC3D( NHGOUT3D+12 ) = 'west-east velocity'
	    VTYPE3D( NHGOUT3D+12 ) = M3REAL
	    VNAME3D( NHGOUT3D+13 ) = 'V'
	    UNITS3D( NHGOUT3D+13 ) = 'm/s'
	    VDESC3D( NHGOUT3D+13 ) = 'north-south velocity'
	    VTYPE3D( NHGOUT3D+13 ) = M3REAL
	    VNAME3D( NHGOUT3D+14 ) = 'FRICE'
	    UNITS3D( NHGOUT3D+14 ) = '%/100'
	    VDESC3D( NHGOUT3D+14 ) = 'fraction covered by ice'
	    VTYPE3D( NHGOUT3D+14 ) = M3REAL
	    VNAME3D( NHGOUT3D+15 ) = 'STOT'
	    UNITS3D( NHGOUT3D+15 ) = 'mgC/m**2'
	    VDESC3D( NHGOUT3D+15 ) = 'total sediment C+N'
	    VTYPE3D( NHGOUT3D+15 ) = M3REAL
            VNAME3D( NHGOUT3D+16 ) = 'MAC'
            UNITS3D( NHGOUT3D+16 ) = 'mgC/m**2'
            VDESC3D( NHGOUT3D+16 ) = 'Macrobenthos biomass'
            VTYPE3D( NHGOUT3D+16 ) = M3REAL
            VNAME3D( NHGOUT3D+17 ) = 'HGMACe2D'
            UNITS3D( NHGOUT3D+17 ) = 'ng/m**2'
            VDESC3D( NHGOUT3D+17 ) = 'Hg on Macrobenthos'
            VTYPE3D( NHGOUT3D+17 ) = M3REAL
            VNAME3D( NHGOUT3D+18 ) = 'HGMACi2D'
            UNITS3D( NHGOUT3D+18 ) = 'ng/m**2'
            VDESC3D( NHGOUT3D+18 ) = 'Hg in Macrobenthos'
            VTYPE3D( NHGOUT3D+18 ) = M3REAL
            VNAME3D( NHGOUT3D+19 ) = 'MEHGMACe2D'
            UNITS3D( NHGOUT3D+19 ) = 'ng/m**2'
            VDESC3D( NHGOUT3D+19 ) = 'MeHg on Macrobenthos'
            VTYPE3D( NHGOUT3D+19 ) = M3REAL
            VNAME3D( NHGOUT3D+20 ) = 'MEHGMACi2D'
            UNITS3D( NHGOUT3D+20 ) = 'ng/m**2'
            VDESC3D( NHGOUT3D+20 ) = 'MeHg in Macrobenthos'
            VTYPE3D( NHGOUT3D+20 ) = M3REAL
            VNAME3D( NHGOUT3D+21 ) = 'HGFLUX'
            UNITS3D( NHGOUT3D+21 ) = 'ng/L'
            VDESC3D( NHGOUT3D+21 ) = 'Hg sedimentation flux'
            VTYPE3D( NHGOUT3D+21 ) = M3REAL
            VNAME3D( NHGOUT3D+22 ) = 'MEHGFLUX'
            UNITS3D( NHGOUT3D+22 ) = 'ng/L'
            VDESC3D( NHGOUT3D+22 ) = 'MeHg sedimentation flux'
            VTYPE3D( NHGOUT3D+22 ) = M3REAL

            IF( ALLOCATED( HGOUT3D ) ) THEN
!	        WRITE( MESG,* ) 'HGOUT3D IS ALLOCATED'
!	        CALL M3MSG2( MESG )
            ELSE
                WRITE( MESG,* ) 'HGOUT3D IS NOT ALLOCATED'
                CALL M3MSG2( MESG )
                ALLOCATE( HGOUT3D( NROWS3D,NCOLS3D,NLAYS3D,NHGOUT3D ),STAT=IOS)
            END IF

            IF( ALLOCATED( HGOUT2D ) ) THEN
!                WRITE( MESG,* ) 'HGOUT2D IS ALLOCATED'
!                CALL M3MSG2( MESG )
            ELSE
                WRITE( MESG,* ) 'HGOUT2D IS NOT ALLOCATED'
                CALL M3MSG2( MESG )
                ALLOCATE( HGOUT2D( NROWS3D,NCOLS3D,1,NHGOUT2D ), STAT=IOS )
            END IF

C........  Create a new file for every year
           CALL ENVSTR( 'CHEMOUT3D','CHEMOUT3D file','',CHEMOUT3D,IOS )
            IF( IOS .GT. 0 ) THEN
                MESG = 'Unknown Value for CHEMOUT3D'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            WRITE( CHEMOUT3D, '(A,I4)' ), TRIM( CHEMOUT3D ), WYEAR

            IF( .NOT. SETENVVAR( 'CCHEMOUT3D', CHEMOUT3D ) ) THEN
                MESG =
     &			'Could not set environment variable CCHEMOUT3D = '
     &		   				          // CHEMOUT3D
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Open 3D output file
            MESG = 'Enter logical name for the CHEMISTRY OUTPUT file'
            O3NAME = 'CCHEMOUT3D'
            O3NAME = PROMPTMFILE( MESG, FSNEW3, O3NAME, PROGNAME )

C.........  The inital time step creates the file but does not write
	    IF( .NOT. CLOSE3( O3NAME ) ) THEN
	        MESG = 'Could not close file ' // O3NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	    END IF

C.........  Adjust global attributes for 2D file
	    FDESC3D = 'ECOSMO_Hg 2D chemistry output file'
	    NLAYS3D = 1
	    NVARS3D = NHGOUT3D + NHGOUT2D
        
C........  Create a new file for every year
           CALL ENVSTR( 'CHEMOUT2D', 'CHEMOUT2D file', '',
     &			 CHEMOUT2D, IOS )
            IF( IOS .GT. 0 ) THEN
                MESG = 'Unknown Value for CHEMOUT2D'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

	    WRITE( CHEMOUT2D, '(A,I4)' ), TRIM( CHEMOUT2D ), WYEAR

            IF( .NOT. SETENVVAR( 'CCHEMOUT2D', CHEMOUT2D ) ) THEN
                MESG = 
     &			'Could not set environment variable CCHEMOUT2D = '
     &		   				          // CHEMOUT2D
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        WRITE( MESG,* ) 'TEST', NVARS3D
        CALL M3MSG2( MESG )

        DO V = 1, NVARS3D
                WRITE( MESG,* ) V, VNAME3D( V ), UNITS3D( V ), VDESC3D( V )
                CALL M3MSG2( MESG )
        END DO

C.........  Open 2D output file
            TSTEP3D = WSTEP2D !Different time step for 2D file (1h) 25.03.2016

            MESG = 'Enter logical name for the CHEMISTRY OUTPUT file'
            O2NAME = 'CCHEMOUT2D'
            O2NAME = PROMPTMFILE( MESG, FSNEW3, O2NAME, PROGNAME )

C.........  The inital time step creates the file but does not write
	    IF( .NOT. CLOSE3( O2NAME ) ) THEN
	        MESG = 'Could not close file ' // O2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	    END IF

	ELSE IF ( WFLAG ) THEN
	    CALL M3MSG2( DASHLINE )
            XSTEP = -1 * WSTEP3D
            CALL NEXTIME( WDATE,WTIME,XSTEP ) ! this is now more flexible 25.03.2016
	    WYEAR = WDATE / 1000
	    WRITE( MESG,* ) 'Opening ocean chemistry output file for',
     &							 WDATE, WTIME
	    CALL M3MSG2( MESG )

C.........  Counter for number of timesteps aggregated in HGOUT?D fields
        WRITE( MESG,* ) 'Temporal averaging and unit conversion factor:', COUNTER, UFAC
        CALL M3MSG2( MESG )


C........  Check for 3D output timestep
!           IF( WTIME .EQ. 0 .OR. WTIME .EQ. 120000) THEN
C........  Create a new file for every year
           CALL ENVSTR( 'CHEMOUT3D', 'CHEMOUT3D file', '',
     &			 CHEMOUT3D, IOS )
            IF( IOS .GT. 0 ) THEN
                MESG = 'Unknown Value for CHEMOUT3D'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

	    WRITE( CHEMOUT3D, '(A,I4)' ), TRIM( CHEMOUT3D ), WYEAR

            IF( .NOT. SETENVVAR( 'CCHEMOUT3D', CHEMOUT3D ) ) THEN
                MESG = 
     &		       'Could not set environment variable CCHEMOUT3D = '
     &		   				         // CHEMOUT3D
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Open 3D output file
            MESG = 'Enter logical name for the CHEMISTRY OUTPUT file'
            O3NAME = 'CCHEMOUT3D'
            O3NAME = PROMPTMFILE( MESG, FSRDWR3, O3NAME, PROGNAME )

C.........  Get header from file and store necessary info
            IF ( .NOT. DESC3( O3NAME ) ) THEN
                MESG = 'Could not get description of file ' // O3NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	    END IF

	    DO V = 1, NVARS3D

		CALL DIMTURN3D( HGOUT3D( :,:,:,V ), TMP, COUNTER*UFAC )
		IF( .NOT. WRITE3( O3NAME, VNAME3D( V ),WDATE,WTIME,TMP ) ) THEN
		    MESG = 'Could not write timestep to file' // O3NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
		END IF

		WRITE( MESG,* ) 'MIN:', MINVAL( TMP ), 'MAX:', MAXVAL( TMP ),
     &				'SUM:', SUMVAL( TMP ), 'AVG:', SUMAVG( TMP )
		CALL M3MSG2( MESG )
	    END DO

	    IF( .NOT. CLOSE3( O3NAME ) ) THEN
	        MESG = 'Could not close file ' // O3NAME
	        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	    END IF
!            END IF !for check 3D output time step

C.........  Write 2D (surface) output file

C........  Create a new file for every year
           CALL ENVSTR( 'CHEMOUT2D', 'CHEMOUT2D file', '',
     &			 CHEMOUT2D, IOS )
            IF( IOS .GT. 0 ) THEN
                MESG = 'Unknown Value for CHEMOUT2D'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

	    WRITE( CHEMOUT2D, '(A,I4)' ), TRIM( CHEMOUT2D ), WYEAR

            IF( .NOT. SETENVVAR( 'CCHEMOUT2D', CHEMOUT2D ) ) THEN
                MESG = 
     &			'Could not set environment variable CCHEMOUT2D = '
     &							  // CHEMOUT2D
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.........  Open 2D output file
            MESG = 'Enter logical name for the CHEMISTRY OUTPUT file'
            O2NAME = 'CCHEMOUT2D'
            O2NAME = PROMPTMFILE( MESG, FSRDWR3, O2NAME, PROGNAME )

C.........  Get header from file and store necessary info
            IF ( .NOT. DESC3( O2NAME ) ) THEN
                MESG = 'Could not get description of file ' // O2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO V = 1, NHGOUT3D
                CALL DIMTURN2D( HGOUT3D( :,:,1,V ), TMP( :,:,1 ), COUNTER*UFAC )
                IF( .NOT. WRITE3( O2NAME, VNAME3D( V ), WDATE, WTIME,
     &                                                  TMP( :,:,1 ) ) ) THEN
                MESG = 'Could not write timestep to file' // O2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                WRITE( MESG,* ) 'MIN:', MINVAL( TMP( :,:,1 ) ),
     &                          'MAX:', MAXVAL( TMP( :,:,1 ) ),
     &                          'SUM:', SUMVAL( TMP( :,:,1:1 ) ),
     &                          'AVG:', SUMAVG( TMP( :,:,1:1 ) )
                CALL M3MSG2( MESG )
            END DO

            DO V = 1, NHGOUT2D

C..... output sum (cummulativ daily fluxes) for air-sea exchange and deposition, etc.
!      but average concentrations for sediment species, layer number etc.
              IF( V .EQ. HGSED2D .OR. V .EQ. MEHGSED2D .OR. V .EQ. LAYERS2D ) THEN
      CALL DIMTURN2D( HGOUT2D( :,:,1,V ), TMP( :,:,1 ), COUNTER*UFAC )
              ELSE
      CALL DIMTURN2D( HGOUT2D( :,:,1,V ), TMP( :,:,1 ), UFAC )
              END IF

                IF( .NOT. WRITE3( O2NAME, VNAME3D( V + NHGOUT3D ),
     &                             WDATE, WTIME, TMP( :,:,1 ) ) ) THEN
                    MESG = 'Could not write timestep to file' // O2NAME
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                WRITE( MESG,* ) 'MIN:', MINVAL( TMP( :,:,1 ) ),
     &                          'MAX:', MAXVAL( TMP( :,:,1 ) ),
     &                          'SUM:', SUMVAL( TMP( :,:,1:1 ) ),
     &                          'AVG:', SUMAVG( TMP( :,:,1:1 ) )
                CALL M3MSG2( MESG )
            END DO

        IF( .NOT. CLOSE3( O2NAME ) ) THEN
            MESG = 'Could not close file ' // O2NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        WRITE( MESG,* ) 'COUNTER CHECK: ', COUNTER, WDATE, WTIME, WFLAG, UFAC
     &                                   , HGOUT3D( 100,40,1,HGTOT3D )
     &                                   , HGOUT3D( 100,40,1,1 )
     &                                   , HGOUT3D( 100,135,1,HGTOT3D )
     &                                   , HGOUT3D( 100,135,1,1 )

        CALL M3MSG2( MESG )

	COUNTER = 0.

C...... Reset output arrays
        WRITE( MESG,* ) 'CHEM resetting data arrays'
        CALL M3MSG2( MESG )
        HGOUT3D = 0.
        HGOUT2D = 0.

	END IF

C.........  End of subroutine
	END SUBROUTINE WRITE_CHEM




	SUBROUTINE READIC( ) 

	USE CPARAM
	USE CUTIL

	IMPLICIT NONE

C...........
!
!	Reads oceanic 3D chemistry fields from outputfile as initial conditions
!
! HISTORY:
!
! 30.05.2013	Creation
!
C...........


C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)
	INCLUDE 'C_model_inc'	!  Model domain dimension size

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(16)      PROMPTMFILE
        LOGICAL            SETENVVAR

	EXTERNAL  PROMPTMFILE, SETENVVAR

C...........   INPUT VARIABLES and their descriptions:
!	INTEGER, INTENT(IN) :: CDATE		! Date to read from file

C..........	LOCAL FIELDS
	INTEGER		    I,R,C,L,V,T		  ! counters and indices
	CHARACTER( 265 )    INAME, BNAME	  ! IC 3d file name
        CHARACTER( 265 )    I2NAME                ! IC 2d file name
	CHARACTER( 265 )    DNAME		  ! dummy output file name

!	REAL, DIMENSION( n,m,ilo ) :: TMP3D	  ! data array
	REAL, DIMENSION( n,m )     :: TMP2D	  ! data array
	REAL, DIMENSION( m,n )     :: TMPSWAP	  ! data array
!	REAL, DIMENSION( m,n )     :: TMPTURN	  ! data array

	INTEGER		    RDATE		  ! read startin timestep - 1
	INTEGER		    RTIME		  ! read startin timestep - 1
        INTEGER             ICDATE                ! helper variable
        INTEGER             MXREC

C...........	LOCAL PARAMETERS
	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'READIC'	  ! program name
        CHARACTER( 16 )     ICVARNAME
        REAL SUM

C***********************************************************************
C   begin body of subroutine READIC
	CALL M3MSG2( DASHLINE )
	WRITE( MESG,* ) 'Reading initial conditions'
	CALL M3MSG2( MESG )
	CALL M3MSG2( '' )

C.........  Reset variables
	HGOUT3D = NAN
	HGOUT2D = 0.

C.........  Open 3D output file
        MESG = 'Enter logical name for the CHEMISTRY IC file'
        INAME = 'CHEM_IC'
        INAME = PROMPTMFILE( MESG, FSRDWR3, INAME, PROGNAME )

C.........  Get header from file and store necessary info
        IF ( .NOT. DESC3( INAME ) ) THEN
            MESG = 'Could not get description of file ' // INAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

	RDATE = SDATE3D
	RTIME = STIME3D

        DO T = 1,MXREC3D - 1
          CALL NEXTIME( RDATE,RTIME,TSTEP3D )
        END DO

C........  Check if dimensions match model domain
	IF( NROWS3D .NE. m ) THEN
            WRITE( MESG,* ) 
     &		'ERROR: Wrong dimensions in IC file NROWS != ', m, NROWS3D
 !     	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF
	IF( NCOLS3D .NE. n ) THEN
            WRITE( MESG,* ) 
     &		'ERROR: Wrong dimensions in IC file NCOLS != ' , n, NCOLS3D
 !     	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

C........  Read time independent IC data from 2D profile file
	IF( NVARS3D .EQ. 1 .AND. VNAME3D( 1 ) .EQ. 'HG' ) THEN
	    CALL M3MSG2( '' )
	    WRITE( MESG,* ) 'Reading data from IC profile file'
	    CALL M3MSG2( MESG )

C........  Read mercury concnetrations [ng/L]
            IF ( .NOT. READ3( INAME, 'HG', 1, 0, 0, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // INAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C........  Because x and y axis are swapped in netcdf files.
            CALL DIMSWAP2D( TMP2D,TMPSWAP )
!            HGOUT2D( :,:,1,HGSED2D )   = 100.
!            HGOUT2D( :,:,1,MEHGSED2D ) =  10.

            DO R = 1,m
                DO C = 1,n

                HGOUT2D( R,C,1,HGSED2D )   = 950.
                HGOUT2D( R,C,1,MEHGSED2D ) =  50.

                    DO L = 1,ilo

                IF( TMPSWAP( R,C ) .GT. 0. ) THEN
                    HGOUT3D( R,C,L,PTOM3D ) = 1E-10
                    HGOUT3D( R,C,L,DTOM3D ) = 1E-10
                END IF
!                HGOUT3D( R,C,L,HGDIS3D )   = 1.0
!                HGOUT3D( R,C,L,HGDEM3D )   = 1.0
!                HGOUT3D( R,C,L,HGDET3D )   = 1.0

!       Populate only HGDEM HGDIS MEHGDIS DMHG - the others will follow ! (e.g. pertpart.f)
! 
! 1000 to convert to ng/L
                        IF( L .LT. 14 ) THEN
               HGOUT3D( R,C,L,HGDEM3D )   = 0.04 * TMPSWAP( R,C ) * 5.0!* 0.6!* 1000
               HGOUT3D( R,C,L,HGDIS3D )   = 0.90 * TMPSWAP( R,C ) * 5.0!* 0.6!* 1000
               HGOUT3D( R,C,L,DMHG3D  )   = 0.01 * TMPSWAP( R,C ) * 5.0!* 0.6!* 1000
               HGOUT3D( R,C,L,MEHGDIS3D ) = 0.05 * TMPSWAP( R,C ) * 5.0!* 0.6!* 1000

                        ELSE
               HGOUT3D( R,C,L,HGDEM3D )   = 0.04 * TMPSWAP( R,C ) * 10.0!* 0.8!* 1000
               HGOUT3D( R,C,L,HGDIS3D )   = 0.90 * TMPSWAP( R,C ) * 10.0!* 0.8!* 1000
               HGOUT3D( R,C,L,DMHG3D  )   = 0.01 * TMPSWAP( R,C ) * 10.0!* 0.8!* 1000
               HGOUT3D( R,C,L,MEHGDIS3D ) = 0.05 * TMPSWAP( R,C ) * 10.0!* 0.8!* 1000
                        END IF

                    END DO
                END DO
            END DO

!           C....... Close file
            IF( .NOT. CLOSE3( INAME ) ) THEN
                MESG = 'Could not close file ' // INAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C........  Read Hg in biota IC file -> Use last timestep of any CHEM3D file for now
C                                      Using zero values breaks inorganic Hg cycling

C.........  Open 3D output file
        MESG = 'Enter logical name for the BIO IC file'
        BNAME = 'BIO_IC'
        BNAME = PROMPTMFILE( MESG, FSRDWR3, BNAME, PROGNAME )

C.........  Get header from file and store necessary info
        IF ( .NOT. DESC3( BNAME ) ) THEN
            MESG = 'Could no get description of file ' // BNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RDATE = SDATE3D
        RTIME = STIME3D
        DO T = 2,MXREC3D
            CALL NEXTIME( RDATE,RTIME,TSTEP3D )
        END DO

C........  Read mercury concnetrations [ng/L]
        DO L = 1,ilo
            IF ( .NOT. READ3( BNAME, 'MEHGFLA', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // BNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGFLA3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF
            IF ( .NOT. READ3( BNAME, 'MEHGDIA', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // BNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGDIA3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF
            IF ( .NOT. READ3( BNAME, 'MEHGCYA', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // BNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGCYA3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF
            IF ( .NOT. READ3( BNAME, 'MEHGZOSe', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
!               MESG = 'ERROR: Reading timestep from file ' // BNAME
!               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGZOSe3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF
            IF ( .NOT. READ3( BNAME, 'MEHGZOLe', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
!               MESG = 'ERROR: Reading timestep from file ' // BNAME
!               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGZOLe3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF

            IF ( .NOT. READ3( BNAME, 'MEHGZOSi', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
!               MESG = 'ERROR: Reading timestep from file ' // BNAME
!               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGZOSi3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF
            IF ( .NOT. READ3( BNAME, 'MEHGZOLi', L, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
!               MESG = 'ERROR: Reading timestep from file ' // BNAME
!               CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                CALL DIMSWAP2D( TMP2D,TMPSWAP )
                DO R = 1,m
                    DO C = 1,n
                        HGOUT3D( R,C,L,MEHGZOLi3D ) = TMPSWAP( M-R+1,C )
                    END DO
                END DO
            END IF

        END DO

C........  Read initial conditions from output file from previous run
        ELSE !916
            CALL M3MSG2( '' )
            WRITE( MESG,* ) 'Reading data from ECOSMO CHEM3D file'
            CALL M3MSG2( MESG )

            IF (RDATE / 1000 .EQ. 2000) THEN
                ICDATE = RDATE-1
            ELSE   
               ICDATE = RDATE
            END IF

C........  Read mercury concnetrations [ng/L]
            DO V = 1, NHGOUT3D
              ICVARNAME = VNAMOUT( V )
!             IF( ICVARNAME == 'HGFISHi' ) ICVARNAME = 'MEHGFISH'
!             IF( ICVARNAME == 'HGFISHe' ) ICVARNAME = 'MEHGFISH'
!             IF( ICVARNAME == 'MEHGFISHi' ) ICVARNAME = 'MEHGFISH'
!             IF( ICVARNAME == 'MEHGFISHe' ) ICVARNAME = 'MEHGFISH'
!             IF( ICVARNAME == 'HGMACi' ) ICVARNAME = 'MEHGMACe'
!             IF( ICVARNAME == 'HGMCAe' ) ICVARNAME = 'MEHGMACe'
!             IF( ICVARNAME == 'MEHGMACi' ) ICVARNAME = 'MEHGMACe'
!             IF( ICVARNAME == 'MEHGMACe' ) ICVARNAME = 'MEHGMACe'
                WRITE( MESG,* ) '    ',ICVARNAME, RDATE, RTIME, VNAMOUT( V )
                CALL M3MSG2( MESG )

                DO L = 1,ilo
                    IF ( .NOT. READ3( INAME, ICVARNAME, L, ICDATE, RTIME, TMP2D( :,: ) ) ) THEN
                        MESG = 'ERROR: Reading timestep from file ' // INAME
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                    CALL DIMSWAP2D( TMP2D,TMPSWAP )

                    DO R = 1,m
                        DO C = 1,n
                            HGOUT3D( R,C,L,V ) = MIN( 1.0,TMPSWAP( M-R+1,C ) )
                            IF( HGOUT3D( R,C,L,V ) < 0. ) THEN
                                HGOUT3D( R,C,L,V ) = 0.
                            END IF

                        END DO
                    END DO

!		C. Get rid of those NaNs from the NaNs'n'Centre
                    DO R = 1,m
                        DO C = 1,n
                            IF( TMPSWAP( R,C ) .GE. 1E20 ) THEN
                                HGOUT3D( R,C,L,V ) = 0.
                            ELSE IF( ISNAN( TMPSWAP( R,C ) ) ) THEN
                                HGOUT3D( R,C,L,V ) = 0.
                            END IF
                        END DO
                    END DO

                END DO
                HGOUT3D(  M- 93-1, 94,:,V ) = 0.
                HGOUT3D(  M-117-1,153,:,V ) = 0.
                HGOUT3D(  M- 83-1, 85,:,V ) = 0.
!              END IF
            END DO

!           C....... Close file
            IF( .NOT. CLOSE3( INAME ) ) THEN
                MESG = 'Could not close file ' // INAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!               READ SEDIMENT RESTART DATA
            CALL M3MSG2( '' )
            WRITE( MESG,* ) 'Reading data from ECOSMO CHEM2D file'
            CALL M3MSG2( MESG )

C.........  Open 2D output file
            MESG = 'Enter logical name for the CHEMISTRY IC file'
            I2NAME = 'CHEM_IC2D'
            I2NAME = PROMPTMFILE( MESG, FSRDWR3, I2NAME, PROGNAME )

            IF ( .NOT. READ3( I2NAME, 'MEHGSED', 1, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // I2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL DIMSWAP2D( TMP2D,TMPSWAP )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,MEHGSED2D ) =  TMPSWAP( M-R+1,C )
              END DO
            END DO

            IF ( .NOT. READ3( I2NAME, 'HGSED', 1, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // I2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL DIMSWAP2D( TMP2D,TMPSWAP )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,HGSED2D ) = TMPSWAP( M-R+1,C ) 
              END DO
            END DO

            IF ( .NOT. READ3( I2NAME, 'MEHGMACe2D', 1, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // I2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL DIMSWAP2D( TMP2D,TMPSWAP )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,MEHGMACe2D ) = TMPSWAP( M-R+1,C )
              END DO
            END DO

            IF ( .NOT. READ3( I2NAME, 'MEHGMACi2D', 1, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // I2NAME
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL DIMSWAP2D( TMP2D,TMPSWAP )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,MEHGMACi2D ) = TMPSWAP( M-R+1,C )
              END DO
            END DO

            IF ( .NOT. READ3( I2NAME, 'HGMACe2D', 1, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // I2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL DIMSWAP2D( TMP2D,TMPSWAP )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,HGMACe2D ) = TMPSWAP( M-R+1,C )
              END DO
            END DO

            IF ( .NOT. READ3( I2NAME, 'HGMACi2D', 1, RDATE, RTIME, TMP2D( :,: ) ) ) THEN
                MESG = 'ERROR: Reading timestep from file ' // I2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL DIMSWAP2D( TMP2D,TMPSWAP )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,HGMACi2D ) = TMPSWAP( M-R+1,C )
              END DO
            END DO

            IF( .NOT. CLOSE3( I2NAME ) ) THEN
                MESG = 'Could not close file ' // I2NAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO V = 1, NHGOUT2D
                HGOUT2D(  M- 93-1, 94,:,V ) = 0.
                HGOUT2D(  M-117-1,153,:,V ) = 0. 
                HGOUT2D(  M- 83-1, 85,:,V ) = 0. 
            END DO

        END IF
        CALL M3MSG2( '' )

            DO R = 1,m
              DO C = 1,n
                HGOUT2D( R,C,1,HGSED2D ) = 
     &          MAX( 0.,HGOUT2D( R,C,1,HGSED2D ) )
                HGOUT2D( R,C,1,MEHGSED2D ) =
     &          MAX( 0.,HGOUT2D( R,C,1,MEHGSED2D ) )
              END DO
            END DO

C........  End of initial conditions
        DO I = 1,NHGOUT2D
            WRITE( MESG,* ) VNAMOUT( NHGOUT3D + I ), 
     &                      MINVAL( HGOUT2D( :,:,1,I ) ),
     &                      MAXVAL( HGOUT2D( :,:,1,I ) )
            CALL M3MSG2( MESG )
        END DO

C........  Give some statistics output
	DO I = 1,NHGOUT3D
	    WRITE( MESG,* ) VNAMOUT( I ), SUMVAL( HGOUT3D( :,:,:,I ) )/I_NDREI,
     &					  MINVAL( HGOUT3D( :,:,:,I ) ),
     &			    		  MAXVAL( HGOUT3D( :,:,:,I ) )
	    CALL M3MSG2( MESG )
	END DO
        SUM = 0.
        DO I = 1,19
          DO L = 1,ilo
            DO R = 1,m
              DO C = 1,n
                SUM = SUM + HGOUT3D( R,C,L,I )
              END DO
            END DO
          END DO
        END DO
        print*, 'Hg in water = ', SUM
        SUM = 0.
        DO I = 20,27
          DO L = 1,ilo
            DO R = 1,m
              DO C = 1,n
                SUM = SUM + HGOUT3D( R,C,L,I )
              END DO
            END DO
          END DO
        END DO
        print*, 'Hg in biota = ', SUM
        SUM = 0.
        DO L = 1,ilo
          DO R = 1,m
            DO C = 1,n
              SUM = SUM + HGOUT3D( R,C,L,HGMACe3D )
     &                  + HGOUT3D( R,C,L,HGMACi3D )
     &                  + HGOUT3D( R,C,L,MEHGMACe3D )
     &                  + HGOUT3D( R,C,L,MEHGMACi3D )
            END DO
          END DO
        END DO
        print*, 'Hg in macro benthos = ', SUM
        SUM = 0.
        DO L = 1,ilo
          DO R = 1,m
            DO C = 1,n
              SUM = SUM + HGOUT3D( R,C,L,HGFISHe3D )
     &                  + HGOUT3D( R,C,L,HGFISHi3D )
     &                  + HGOUT3D( R,C,L,MEHGFISHe3D )
     &                  + HGOUT3D( R,C,L,MEHGFISHi3D )
            END DO
          END DO
        END DO
        print*, 'Hg in fish = ', SUM
        SUM = 0.
        DO L = 1,ilo
          DO R = 1,m
            DO C = 1,n
              SUM = SUM + HGOUT3D( R,C,L,SEDHG3D )
     &                  + HGOUT3D( R,C,L,SEDMEHG3D ) 
            END DO
          END DO
        END DO
        print*, 'Hg in sediment = ', SUM


!	C....... Close file
!        IF( .NOT. CLOSE3( INAME ) ) THEN
!            MESG = 'Could not close file ' // INAME
!            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!        END IF

C........  Write dummy output to visually check IC fields
!	NVARS3D = NHGOUT3D

!        MESG = 'Enter logical name for the CHEMISTRY IC file'
!        DNAME = 'CHEM_CHECK'
!        DNAME = PROMPTMFILE( MESG, FSUNKN3, DNAME, PROGNAME )

!	DO V = 1, NHGOUT3D

!		CALL DIMTURN3D( HGOUT3D( :,:,:,V ), TMP3D( :,:,: ), 1.0 )
!	    IF ( .NOT. WRITE3( DNAME, VNAME3D( V ), SDATE3D, STIME3D,
!     &						TMP3D( :,:,1 ) ) ) THEN
!        	MESG = 'ERROR: Reading timestep from file ' // DNAME
!        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!            END IF

!	END DO

!        IF( .NOT. CLOSE3( DNAME ) ) THEN
!            MESG = 'Could not close file ' // DNAME
!            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
!        END IF

	CALL M3MSG2( '' )

C.......... End of subroutine
	END SUBROUTINE READIC


	REAL FUNCTION SUMVAL( ARRAY )

	IMPLICIT NONE

	REAL ARRAY( :,:,: )
	INTEGER D1
	INTEGER D2
	INTEGER D3
	INTEGER X,Y,Z
	REAL NAN

	D1 = SIZE( ARRAY,DIM=1 )
	D2 = SIZE( ARRAY,DIM=2 )
	D3 = SIZE( ARRAY,DIM=3 )
	SUMVAL = 0; 
	NAN = 0.
	NAN = NAN / 0.

	DO X = 1, D1
	    DO Y = 1, D2
		DO Z = 1, D3
		    IF( ARRAY( X,Y,Z ) .LT. 1E20 ) THEN
			SUMVAL = SUMVAL + ARRAY( X,Y,Z )
		    END IF
		END DO
	    END DO
	END DO

	END FUNCTION SUMVAL


	REAL FUNCTION SUMAVG( ARRAY )

	IMPLICIT NONE

	REAL ARRAY( :,:,: )
	INTEGER D1
	INTEGER D2
	INTEGER D3
	INTEGER X,Y,Z
	REAL NAN
	REAL C

	D1 = SIZE( ARRAY,DIM=1 )
	D2 = SIZE( ARRAY,DIM=2 )
	D3 = SIZE( ARRAY,DIM=3 )
	SUMAVG = 0; 
	NAN = 0.
	NAN = NAN / 0.
	C = 0.

	DO X = 1, D1
	    DO Y = 1, D2
		DO Z = 1, D3
		    IF( ARRAY( X,Y,Z ) .LT. 1E20 ) THEN
			SUMAVG = SUMAVG + ARRAY( X,Y,Z )
			C = C + 1.
		    END IF
		END DO
	    END DO
	END DO

	SUMAVG = SUMAVG / C

	END FUNCTION SUMAVG

C.........  End of module
	END MODULE COUT
