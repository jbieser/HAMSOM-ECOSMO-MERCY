	MODULE CIN

	USE CPARAM
	USE CUTIL

	IMPLICIT NONE
C...........
!
! DESC:	This module handles all atmospheric input fields
!	This includes met fields, concentration fields, and deposition fields
!
!	It contains a global data array for atmospheric met and chem data
!		 called ATMIN( m,n,HOUR,MFNUM+MCNUM ) ATMIN = ATMospheric INput
!	
!	It contains all subroutines to read this data from CMAQ and MCIP netcdf files
!		and subroutines to interpolate them to the ECOSMO grid
!
!	Basially all one needs to do is to
!		1) USE CIN
!		2) CALL CIN( )where CDATE is a Julian date <YYYYDDD>
!
!		CIN calls READCHEM and READMET which use INTERPOL to interpolate
!			the input data and write it into the global array ATMIN
!			Optionally WRITEIPOL can write ATMIN into a netcdf output file
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........

C...........   Indices for meteorological parameters
	INTEGER,	PARAMETER :: MFNUM = 9	! Number of meteorological variables
	INTEGER,	PARAMETER :: MTEMP2 = 1	  ! Temperature at 2m 		      	[K]    2D
	INTEGER,	PARAMETER :: MPRSFC = 2	  ! Surface pressure 		      	[Pa]   2D
	INTEGER,	PARAMETER :: MUSTAR = 3	  ! atmospheric friction velocity     	[m/s]  2D
	INTEGER,	PARAMETER :: MRADYN = 4   ! invers of aerodynamic resistance  	[m/s]  2D
	INTEGER,	PARAMETER :: MRGRND = 5	  ! Solar Radiation reaching ground   	[W/m²] 2D
	INTEGER,	PARAMETER :: MWSDP10 = 6  ! Wind speed at 10m		      	[m/s]  2D
	INTEGER,	PARAMETER :: MTA   = 7	  ! Temperature (average of grid cell)	[K]    3D
	INTEGER,	PARAMETER :: MPRES = 8	  ! Pressure   (average of grid cell) 	[Pa]   3D
	INTEGER,	PARAMETER :: MZF   = 9	  ! Grid cell height		     	[m]    3D

C...........   Indices for chemical parameters (atmosphere)
	INTEGER,	PARAMETER :: MCNUM = 7	! Number of chemical variables
	INTEGER,	PARAMETER :: MCHG0 = 10	  ! Atmospheric concentration of Hg0 [ppm]
	INTEGER,	PARAMETER :: MDHG0 = 11	  ! Dry deposition of Hg0 [kg/hectare]
	INTEGER,	PARAMETER :: MDHG2 = 12	  ! Dry deposition of Hg2 [kg/hectare]
	INTEGER,	PARAMETER :: MDHGP = 13	  ! Dry deposition of HgP [kg/hectare]
	INTEGER,	PARAMETER :: MWHG0 = 14	  ! Wet deposition of Hg0 [kg/hectare]
	INTEGER,	PARAMETER :: MWHG2 = 15	  ! Wet deposition of Hg2 [kg/hectare]
	INTEGER,	PARAMETER :: MWHGP = 16	  ! Wet deposition of HgP [kg/hectare]

!       Flux array indices
        INTEGER, PARAMETER :: PROD_FLA   =  1 ! Flaggelate primary  production  [mgC/m**3]
        INTEGER, PARAMETER :: PROD_DIA   =  2 ! Diatom primary production       [mgC/m**3]
        INTEGER, PARAMETER :: ZOS_ON_FLA =  3 ! Feeding of small zooplankton on flaggelates           [mgC/m**3]
        INTEGER, PARAMETER :: ZOS_ON_DIA =  4 ! Feeding of small zooplankton on diatoms               [mgC/m**3]
        INTEGER, PARAMETER :: ZOS_ON_DET =  5 ! Feeding of small zooplankton on detritus              [mgC/m**3]
        INTEGER, PARAMETER :: ZOL_ON_FLA =  6 ! Feeding of large zooplankton on flaggelates           [mgC/m**3]
        INTEGER, PARAMETER :: ZOL_ON_DIA =  7 ! Feeding of large zooplankton on diatoms               [mgC/m**3]
        INTEGER, PARAMETER :: ZOL_ON_ZOS =  8 ! Feeding of large zooplankton on small zooplankton     [mgC/m**3]
        INTEGER, PARAMETER :: ZOL_ON_DET =  9 ! Feeding of large zooplankton on detritus              [mgC/m**3]
        INTEGER, PARAMETER :: UNKN_FLUX  = 10
        INTEGER, PARAMETER :: FSH_ON_ZOS = 11
        INTEGER, PARAMETER :: FSH_ON_ZOL = 12
        INTEGER, PARAMETER :: FSH_ON_DET = 13
        INTEGER, PARAMETER :: FSH_ON_MAC = 14
        INTEGER, PARAMETER :: MAC_ON_SED = 15 ! [mgC/m**2 per day]
        INTEGER, PARAMETER :: MAC_ON_DET = 16 ! Detritus and DOM combined
        INTEGER, PARAMETER :: MAC_ON_ZOO = 17 ! total zooplankton
        INTEGER, PARAMETER :: MAC_ON_PHY = 18 ! total phytoplankton
        INTEGER, PARAMETER :: UP_N       = 19
        INTEGER, PARAMETER :: UP_P       = 20
        INTEGER, PARAMETER :: UP_Si      = 21
        INTEGER, PARAMETER :: B_light    = 22
        INTEGER, PARAMETER :: PROD_CYA   = 23 ! Cyanobacteria primary production [mgC/m**3]
        INTEGER, PARAMETER :: ZOS_ON_CYA = 24
        INTEGER, PARAMETER :: ZOL_ON_CYA = 25
        INTEGER, PARAMETER :: SEDI       = 26 ! [mgC/m**3]
        INTEGER, PARAMETER :: RESU       = 27 ! [mgC/m**2] should be per m²
        INTEGER, PARAMETER :: DENIT      = 28
        INTEGER, PARAMETER :: BURIAL     = 29
        INTEGER, PARAMETER :: TAUBOT     = 30
        INTEGER, PARAMETER :: MORT_ZOS   = 31
        INTEGER, PARAMETER :: MORT_ZOL   = 32
        INTEGER, PARAMETER :: MORT_FLA   = 33
        INTEGER, PARAMETER :: MORT_DIA   = 34
        INTEGER, PARAMETER :: MORT_CYA   = 35
        INTEGER, PARAMETER :: MORT_FSH   = 36
        INTEGER, PARAMETER :: MORT_MAC   = 37

!       Bio array indices
        INTEGER, PARAMETER :: FLA  =  3 ! Flaggelates   [mgC/m**3]
        INTEGER, PARAMETER :: DIA  =  4 ! Diatoms       [mgC/m**3]
        INTEGER, PARAMETER :: ZOS  =  5 ! Small Zooplankton     [mgC/m**3]  
        INTEGER, PARAMETER :: ZOL  =  6 ! Large Zooplankton     [mgC/m**3]
        INTEGER, PARAMETER :: DET  =  7 ! Detritus
        INTEGER, PARAMETER :: NH4  =  8
        INTEGER, PARAMETER :: DOM  =  9
        INTEGER, PARAMETER :: NO3  = 10
        INTEGER, PARAMETER :: PO4  = 11
        INTEGER, PARAMETER :: SiO2 = 12
        INTEGER, PARAMETER :: O2   = 13
        INTEGER, PARAMETER :: OPAL = 14
        INTEGER, PARAMETER :: CYA  = 15 ! Cyanobacteria  [mgC/m**3]
        INTEGER, PARAMETER :: FSH  = 16 ! Fish           [mgC/m**3]
        INTEGER, PARAMETER :: STOT = 17 ! Sediment C & N [mgC/m**2]
        INTEGER, PARAMETER :: PSO4 = 18 ! Sediment P     [mgC/m**2]
        INTEGER, PARAMETER :: SSiO2 = 19! Sediment Si    [mgC/m**2]
        INTEGER, PARAMETER :: MAC   = 20 ! Macro benthos [mgC/m**2]

C............	Data array
	REAL, ALLOCATABLE :: ATMIN( :,:,:,: )	! dimensions: rows,cols,hour,MFNUM+MCNUM

C........... Global variables for riverine mercury input
	REAL, ALLOCATABLE :: RIVER( :,:,: )	! dimensions: rows,cols,month	unit: [kg Hg/month]

C........... Global attributes (used for interpolation)
	INTEGER INUM		! interpolation cell count
	INTEGER YORIG		! northing	[km] (positive)
	INTEGER XORIG		! easting 	[km] (positive)
	INTEGER XCELL		! cell size	[km]
	INTEGER NROWS		! number of rows
	INTEGER NCOLS		! number of columns
	REAL	P_ALP		! 1st standard parallel
	REAL	P_BET		! 2nd standard parallel
	REAL	YCENT		! central latitude
	REAL	XCENT		! central longitude
	LOGICAL :: IREADY = .FALSE.
	CHARACTER( 265 )    GRID		  ! for QA (not used atm)
	CHARACTER( 265 )    PROJECTION		  ! for QA (not used atm)

C........... interpolation matrix
	INTEGER, ALLOCATABLE :: IROW( : )
	INTEGER, ALLOCATABLE :: ICOL( : )
	INTEGER, ALLOCATABLE :: OROW( : )
	INTEGER, ALLOCATABLE :: OCOL( : )
	REAL, ALLOCATABLE    :: DIST( : )


C...........

	CONTAINS

	SUBROUTINE CIN_START( )

	USE CPARAM

	IMPLICIT NONE

C...........
!
!	This is the central subroutine that invokes all other subroutines
!		reading atmospheric input files
!
! INPUT:
!	INTEGER 	     CDATE	Expected start date of input files
!
! OUTPUT:
!
! HISTORY: 26.12.2012 Centralized funcionality of module CIN
!
C...........

C...........   INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C...........   EXTERNAL FUNCTIONS and their descriptions:
 
C...........   INPUT VARIABLES and their descriptions:
!        INTEGER, INTENT(IN) :: CDATE     ! start Julian date (YYYYDDD)

C...........   LOCAL VARIABLES
	INTEGER		    YEAR		  ! current year

	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer
        CHARACTER( 16 )  :: PROGNAME = 'CIN'	  ! program name


C***********************************************************************
C   begin body of subroutine READMET

C........... initiate IO API
        LDEV = INIT3()

	! Read MCIP files for met input
	CALL READMET( CDATE )
	! Read CMAQ files for chemistry input

	YEAR = CDATE / 1000
	YEAR = CDATE - 1000 * YEAR + 2009000
	IF( YEAR .EQ. 2009366 ) THEN
	    YEAR = 2009365
	END IF
	WRITE( MESG,* ) 'Year for CMAQ input is fixed to 2009 '
	CALL READCHEM( YEAR )

	END SUBROUTINE CIN_START



	SUBROUTINE INTERPOL( IN, HOUR, VAR ) 

	IMPLICIT NONE

C...........
!
!	Interpolates Lambert Conformal atmospheric fields to model domain
!	When invoked for the first time it automatically calls READIPOLMAT
!
! HISTORY:  24.05.2013	Creation
!
C...........

C...........   INCLUDES:
	INCLUDE 'C_model_inc'

C...........   EXTERNAL FUNCTIONS and their descriptions:

C...........   INPUT VARIABLES and their descriptions:
	REAL,    INTENT(IN) :: IN( :,: )	! input array of size NCOLS,NROWS,1,1
	INTEGER, INTENT(IN) :: HOUR		! 3rd dimension of output array
	INTEGER, INTENT(IN) :: VAR		! 4th dimension of output array

C..........	LOCAL FIELDS
	INTEGER		    I,J,X,Y,C,R		  ! counters and indices
	INTEGER 	    NLINES		  ! number of lines in ASCII file

	REAL		    IVAL		  ! interpolated value
	REAL		    NORM		  ! normalization factor

C...........	LOCAL PARAMETERS
	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'INTERPOL' ! program name

C***********************************************************************
C   begin body of subroutine INTERPOL

!....... Read file name for interpolation matrix
	IF( IREADY == .FALSE. ) THEN
	    WRITE( MESG,* ) ''
	    CALL M3MSG2( MESG )
	    WRITE( MESG,* ) 'Inizializing interpolation matrix', NROWS, NCOLS,
     &		m,n,INUM
	    CALL M3MSG2( MESG )
	    CALL READIPOLMAT( )
	END IF

!	WRITE( MESG,* ) SIZE( MF ), SIZE( MF,DIM=1 ), SIZE( MF,DIM=2 ),
!     &			SIZE( MF,DIM=3 ), SIZE( MF,DIM=4 )
!	CALL M3MSG2( MESG )

	DO I = 1, m*n*INUM-1 , INUM
	    IVAL = 0
	    NORM = 0
	    Y = OROW( I ) + 1
	    X = OCOL( I ) + 1

	    DO J = 1, INUM
		NORM = NORM + 1/DIST( I+J-1 )
	    END DO

	    DO J = 1, INUM
		R = IROW( I+J-1 ) + 1
		C = ICOL( I+J-1 ) + 1
		IVAL = IVAL + IN( C,R ) / DIST( I+J-1 ) / NORM
	    END DO

	    ATMIN( Y, X, HOUR, VAR ) = IVAL

	END DO

C.......... End of subroutine
	END SUBROUTINE INTERPOL


	SUBROUTINE READIPOLMAT( ) 

	IMPLICIT NONE

C...........
!
!      Reads interpolation matrix used by subroutine INTERPOL
!
! HISTORY:
!
! 24.05.2013	Creation
!
C...........

C...........   INCLUDES:
	INCLUDE 'C_model_inc'

C...........   EXTERNAL FUNCTIONS and their descriptions:

C...........   INPUT VARIABLES and their descriptions:

C.......... 	LOCAL FIELDS
	INTEGER		    I,X,Y

	REAL		    X1,X2,X3,X4		  ! dummy variables
        CHARACTER( 1 )      H                     ! dummy variable
        CHARACTER( 256 )    HEAD                  ! file header
	INTEGER 	    NLINES		  ! number of lines in ASCII file
	CHARACTER( 265 )    FNAME		  ! file name of interpolation matrix
	INTEGER		    U, V		  ! unit number for input file
	LOGICAL		    EXS			  ! unit exists?
	LOGICAL		    OPN			  ! unit is open?

C...........   LOCAL PARAMETERS
	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'READIPOLMAT' ! program name

C***********************************************************************
C   begin body of subroutine INTERPOL

!....... This subroutine my only be run once - it is invoked automatically by INTERPOL
	IF( IREADY .EQ. .FALSE. ) THEN

	    WRITE( MESG,* ) 'Allocating memory'
	    CALL M3MSG2( MESG )
	    WRITE( MESG,* ) ''
	    CALL M3MSG2( MESG )

!	    C....... Allocate variables
	    ALLOCATE( ATMIN( m,n,24,MFNUM+MCNUM ), STAT=IOS )
!	    CALL CHECKMEM( IOS, 'MF', PROGNAME )
	    ATMIN = 0

	    WRITE( MESG,* ) 'Reading interpolation matrix file'
	    CALL M3MSG2( MESG )

!....... Read file name for interpolation matrix
           CALL ENVSTR( 'IPOLMAT', 'Interpoaltion matrix',
     &					 '', FNAME, IOS )
            IF( IOS .GT. 0 ) THEN
                MESG = 'Unknown Value for IPOLMAT'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!....... Check for unused unit to open file
	    U = 1
	    DO V = 7,99
	        INQUIRE( UNIT=V,EXIST=EXS,OPENED=OPN )
	        IF( EXS .AND. .NOT. OPN ) THEN
		    U = V
	        END IF
	    END DO

!....... Open interpolation matrix file
	    OPEN( U,FILE=FNAME,STATUS='OLD',ACTION='READ' )

!....... Read header
	    READ( U,'(A1,A)' ) H, HEAD
	    CALL M3MSG2( HEAD )
	    READ( U,'(A1,A)' ) H, GRID
	    CALL M3MSG2( GRID )
	    READ( U,'(A1,A)' ) H, PROJECTION
	    CALL M3MSG2( PROJECTION )
	    READ( U,'(A1,I1)' ) H, INUM
	    WRITE( MESG,* ) INUM
	    CALL M3MSG2( MESG )
	    READ( U,'(A1,I)' ) H, X
	    READ( U,'(A1,I)' ) H, Y
	!...... check if model domain matches interpolation domain
	    IF( X .NE. M .OR. Y .NE. N ) THEN
		WRITE( MESG,* ) X, ' == ', M
	 	CALL M3MSG2( MESG )
		WRITE( MESG,* ) Y, ' == ', N
	 	CALL M3MSG2( MESG )
		WRITE( MESG,* ) 'Dimension mismatch - wrong matrix file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	    END IF
	    READ( U,'(A1,A)' ) H, HEAD
	    READ( U,'(A1,A)' ) H, HEAD
	    CALL M3MSG2( HEAD )

!....... Allocate interpolation arrays
	    NLINES = m*n*INUM
	    ALLOCATE( IROW( NLINES ), STAT=IOS )
	    ALLOCATE( ICOL( NLINES ), STAT=IOS )
	    ALLOCATE( OROW( NLINES ), STAT=IOS )
	    ALLOCATE( OCOL( NLINES ), STAT=IOS )
	    ALLOCATE( DIST( NLINES ), STAT=IOS )

!....... Read interpolaton data and store in global fields
	    DO I = 1, NLINES
	        READ( U,* ) IROW( I ), ICOL( I ), X1, X2,
     &			    OROW( I ), OCOL( I ), X3, X4, DIST( I )
		!........ print first data line for debugging
	        IF( I .EQ. 1 ) THEN
	            WRITE( MESG,* ) IROW( I ), ICOL( I ),
     &			            OROW( I ), OCOL( I ), DIST( I )
	            CALL M3MSG2( MESG )
		    WRITE( MESG, * ) ''
		    CALL M3MSG2( MESG )
	        END IF
	    END DO

!....... Close interpolation matrix file
	    CLOSE( U )

	    IREADY = .TRUE.
	END IF

C.......... End of subroutine
	END SUBROUTINE READIPOLMAT



	SUBROUTINE WRITEIPOL( CDATE )

	IMPLICIT NONE

C...........
!
!	Writes the interpolated data to a netcdf output file
!
! INPUT:
!	INTEGER 	     CDATE	Expected start date of input files
!
! OUTPUT:
!
! HISTORY:
!
! 26.12.2012	Extracted from READMET
!
C...........

C...........   INCLUDES:
	INCLUDE 'C_model_inc'	!  Global variables
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(16)      PROMPTMFILE
        LOGICAL            SETENVVAR

        EXTERNAL  PROMPTMFILE, SETENVVAR

C...........   INPUT VARIABLES and their descriptions:
        INTEGER, INTENT(IN) :: CDATE     ! start Julian date (YYYYDDD)

C...........   LOCAL VARIABLES and their descriptions:
	INTEGER		H,V	! HOUR, ROW, COL, VAR
        INTEGER         SDATE  	! start Julian date (YYYYDDD)
        INTEGER         STIME   ! start time (HHMMSS)

        REAL   , ALLOCATABLE :: TMP( :,: )   ! temporal data field

	INTEGER		     YEAR	 ! year (CDATE / 1000)
        INTEGER              JDAY
        CHARACTER( 256 )     METCRO2D    ! location of Meteorological Fields 2D
        CHARACTER( NAMLEN3 ) M2NAME      ! logical name for Meteorological Fields 2D
        CHARACTER( NAMLEN3 ) ONAME       ! logical name for interpolated output file

C...........   LOCAL PARAMETERS
	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'WRITEIPOL'! program name


C***********************************************************************
C   begin body of subroutine WRITEIPOL

C........... initiate IO API
        LDEV = INIT3()

C....... Read CMAQ METCRO2D file in order to retrieve global attributes
       CALL ENVSTR( 'METCRO2D', 'METCRO2D file',
     &					 '', METCRO2D, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for METCRO2D'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        YEAR = CDATE / 1000
        JDAY = CDATE - 1000 * YEAR
        IF( JDAY .EQ. 366 ) JDAY = 365
        YEAR = 1990
        YEAR = 2006
        SDATE = 1000 * YEAR + JDAY
        STIME = STIME3D

	WRITE( METCRO2D, '(A,I4,A15,I7)' ), TRIM( METCRO2D ),
     &			YEAR, '/METCRO2D_cd72_', SDATE

        IF( .NOT. SETENVVAR( 'CMETCRO2D', METCRO2D ) ) THEN
            MESG = 'Could not set environment variable CMETCRO2D = '
     &						    // METCRO2D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	C....... Open file
        MESG = 'Enter logical name for the METCRO3D file'
        M2NAME = 'CMETCRO2D'
        M2NAME = PROMPTMFILE( MESG, FSREAD3, M2NAME, PROGNAME )

!	C....... Read global attributes
        IF ( .NOT. DESC3( M2NAME ) ) THEN
            MESG = 'Could not get description of file ' // M2NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!       SDATE = SDATE3D
!       STIME = STIME3D

!	C....... Close file
        IF( .NOT. CLOSE3( M2NAME ) ) THEN
            MESG = 'Could not close file ' // M2NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C..........	Adjust global attributes
	NVARS3D = MFNUM+MCNUM
	VNAME3D( 1 ) = 'TEMP2'
	UNITS3D( 1 ) = 'K'
	VNAME3D( 2 ) = 'PRSFC'
	UNITS3D( 2 ) = 'Pa'
	VNAME3D( 3 ) = 'USTAR'
	UNITS3D( 3 ) = 'm/s'
	VNAME3D( 4 ) = 'RADYNI'
	UNITS3D( 4 ) = 'm/s'
	VNAME3D( 5 ) = 'RGRND'
	UNITS3D( 5 ) = 'W/m**2'
	VNAME3D( 6 ) = 'WSDP10'
	UNITS3D( 6 ) = 'm/s'
	VNAME3D( 7 ) = 'TA'
	UNITS3D( 7 ) = 'K'
	VNAME3D( 8 ) = 'PRES'
	UNITS3D( 8 ) = 'Pa'
	VNAME3D( 9 ) = 'ZF'
	UNITS3D( 9 ) = 'm'
	VNAME3D( 10 ) = 'HG'
	UNITS3D( 10 ) = 'ng/m**3'
	VNAME3D( 11 ) = 'DRYHG'
	UNITS3D( 11 ) = 'ng/m**2'
	VNAME3D( 12 ) = 'DRYHGIIGAS'
	UNITS3D( 12 ) = 'ng/m**2'
	VNAME3D( 13 ) = 'DRYHGP'
	UNITS3D( 13 ) = 'ng/m**2'
	VNAME3D( 14 ) = 'WETHG'
	UNITS3D( 14 ) = 'ng/m**2'
	VNAME3D( 15 ) = 'WETHGIIGAS'
	UNITS3D( 15 ) = 'ng/m**2'
	VNAME3D( 16 ) = 'WETHGP'
	UNITS3D( 16 ) = 'ng/m**2'

	NROWS3D = m
	NCOLS3D = n
	NLAYS3D = 1

!	C....... Allocate tmp array with swapped dimensions for IO API output
	ALLOCATE( TMP( n,m ), STAT=IOS )
	TMP = 0

!	C....... Open file
        MESG = 'Enter logical name for the TESTOUT file'
        ONAME = 'TESTOUT'
	ONAME = PROMPTMFILE( MESG, FSUNKN3, ONAME, PROGNAME )

	DO H = 1, 24
	    DO V = 1, MFNUM+MCNUM
	        CALL DIMTURN2D( ATMIN( :,:,H,V ), TMP, 1. )
                IF( .NOT. WRITE3( ONAME, VNAME3D( V ),
     &          	      	  SDATE, STIME, TMP ) ) THEN
                    MESG = 'Could not write timestep to file' // ONAME
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

	    END DO
	    CALL NEXTIME( SDATE,STIME,TSTEP3D )
	END DO

!	C....... Close file
        IF( .NOT. CLOSE3( ONAME ) ) THEN
            MESG = 'Could not close file ' // ONAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	END SUBROUTINE WRITEIPOL



	SUBROUTINE READMET( CDATE )

	IMPLICIT NONE

C...........
!
! Meteorological data is read for 24 time steps once each day
!
! INPUT:
!	INTEGER 	     CDATE	Expected start date of input files
!	Reads file names from environmental variables METCRO2D and METCRO3D
!
! OUTPUT:
!	REAL MF( ROW,COL,HOUR,VAR ) 	Meteorological fields in lowest model layer
!
! HISTORY:
!
! 24.12.2012	Read daily 2D and 3D met fields
!
! 24.05.2013	Moved to module cin.f
!		Added checksum output (MINVAL, MAXVAL)
!		Added grid and projection parameters for interpolation
!
C...........

C...........   INCLUDES:
	INCLUDE 'C_model_inc'	!  Global variables
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(16)      PROMPTMFILE
        LOGICAL            SETENVVAR

        EXTERNAL  PROMPTMFILE, SETENVVAR

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        INTEGER, INTENT(IN) :: CDATE     ! start Julian date (YYYYDDD)

C...........   LOCAL VARIABLES and their descriptions:
	INTEGER		H,R,C,V	! HOUR, ROW, COL, VAR
        INTEGER         SDATE  	! start Julian date (YYYYDDD)
        INTEGER         STIME   ! start time (HHMMSS)

        REAL   , ALLOCATABLE :: TMP( :,: )   ! temporal data field

	INTEGER		 YEAR	     ! Year
        INTEGER          JDAY
        CHARACTER( 256 ) METCRO2D    ! location of Meteorological Fields 2D
        CHARACTER( 256 ) METCRO3D    ! location of Meteorological Fields 3D
        CHARACTER( NAMLEN3 ) M2NAME  ! logical name for Meteorological Fields 2D
        CHARACTER( NAMLEN3 ) M3NAME  ! logical name for Meteorological Fields 3D

C...........   LOCAL PARAMETERS
	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'READMET'  ! program name


C***********************************************************************
C   begin body of subroutine READMET

C........... initiate IO API
        LDEV = INIT3()

!##########################################################################################
	CALL M3MSG2( DASHLINE )
	WRITE( MESG,* ) 'Reading meteo fields for day ', CDATE
	CALL M3MSG2( MESG )
	CALL M3MSG2( '' )

!....... Read CMAQ METCRO3D file
       CALL ENVSTR( 'METCRO3D', 'METCRO3D file',
     &					 '', METCRO3D, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for METCRO3D'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	WRITE( METCRO3D, '(A,I7)' ), TRIM( METCRO3D ), CDATE

        YEAR = CDATE / 1000
        JDAY = CDATE - 1000 * YEAR
        IF( JDAY .EQ. 366 ) JDAY = 365
        YEAR = 1990
        YEAR = 2006
        SDATE = 1000 * YEAR + JDAY
        STIME = STIME3D
!        IF( SDATE .EQ. 1990366 ) SDATE = 1990365

	WRITE( METCRO3D, '(A,I4,A15,I7)' ), TRIM( METCRO3D ),
     &			YEAR, '/METCRO3D_cd72_', SDATE

        IF( .NOT. SETENVVAR( 'CMETCRO3D', METCRO3D ) ) THEN
            MESG = 'Could not set environment variable CMETCRO3D = '
     &						    // METCRO3D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	C....... Open file
        MESG = 'Enter logical name for the METCRO3D file'
        M3NAME = 'CMETCRO3D'
        M3NAME = PROMPTMFILE( MESG, FSREAD3, M3NAME, PROGNAME )

!	C....... Read global attributes
        IF ( .NOT. DESC3( M3NAME ) ) THEN
            MESG = 'Could not get description of file ' // M3NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!       SDATE = CDATE
!       SDATE = 1000 * YEAR + JDAY
!       STIME = STIME3D
!	C....... This needs to be cross checked with chemistry input dimensions
	NCOLS = NCOLS3D
	NROWS = NROWS3D

!	C....... Check consistency
	IF( NCOLS3D .NE. NCOLS ) THEN
            MESG = 'ERROR: Dimension NCOLS in file ' // METCRO3D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	ELSE IF( NROWS3D .NE. NROWS ) THEN
            MESG = 'ERROR: Dimension NROWS in file ' // METCRO3D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

!	C....... Allocate variables!
	ALLOCATE( TMP( NCOLS,NROWS ), STAT=IOS )
!	CALL CHECKMEM( IOS, 'TMP', PROGNAME )

	DO H = 1, 24

!		C....... Read data and call INTERPOL to store interpolated fields in ATMIN array
        	IF ( .NOT. READ3( M3NAME, 'TA', 1, SDATE, STIME,
     &						 TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M3NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MTA )

        	IF ( .NOT. READ3( M3NAME, 'PRES', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M3NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MPRES )

        	IF ( .NOT. READ3( M3NAME, 'ZF', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M3NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MZF )

		CALL NEXTIME( SDATE,STIME,TSTEP3D )

	END DO !H = 1, 24

!	C....... Close file
        IF( .NOT. CLOSE3( M3NAME ) ) THEN
            MESG = 'Could not close file ' // M3NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!##########################################################################################


C....... Read CMAQ METCRO2D file
       CALL ENVSTR( 'METCRO2D', 'METCRO2D file',
     &					 '', METCRO2D, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for METCRO2D'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	WRITE( METCRO2D, '(A,I7)' ), TRIM( METCRO2D ), CDATE

	YEAR = CDATE / 1000
        JDAY = CDATE - 1000 * YEAR
        YEAR = 1990
        YEAR = 2006
        SDATE = 1000 * YEAR + JDAY
        IF( SDATE .EQ. 1990366 ) SDATE = 1990365

	WRITE( METCRO2D, '(A,I4,A15,I7)' ), TRIM( METCRO2D ),
     &			YEAR, '/METCRO2D_cd72_', SDATE

        IF( .NOT. SETENVVAR( 'CMETCRO2D', METCRO2D ) ) THEN
            MESG = 'Could not set environment variable CMETCRO2D = '
     &						    // METCRO2D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	C....... Open file
        MESG = 'Enter logical name for the METCRO3D file'
        M2NAME = 'CMETCRO2D'
        M2NAME = PROMPTMFILE( MESG, FSREAD3, M2NAME, PROGNAME )

!	C....... Read global attributes
        IF ( .NOT. DESC3( M2NAME ) ) THEN
            MESG = 'Could not get description of file ' // M2NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	SDATE = CDATE
!       SDATE = 1000 * YEAR + JDAY

!	C....... Check consistency
	IF( NCOLS3D .NE. NCOLS ) THEN
            MESG = 'ERROR: Dimension NCOLS in file ' // METCRO2D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	ELSE IF( NROWS3D .NE. NROWS ) THEN
            MESG = 'ERROR: Dimension NROWS in file ' // METCRO2D
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

	DO H = 1, 24

!		C....... Read data and call INTERPOL to store interpolated fields in ATMIN array
        	IF ( .NOT. READ3( M2NAME, 'TEMP2', 1, SDATE, STIME, 
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M2NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MTEMP2 )

        	IF ( .NOT. READ3( M2NAME, 'PRSFC', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M2NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MPRSFC )

        	IF ( .NOT. READ3( M2NAME, 'USTAR', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M2NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MUSTAR )

        	IF ( .NOT. READ3( M2NAME, 'RADYNI', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M2NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MRADYN )

        	IF ( .NOT. READ3( M2NAME, 'RGRND', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M2NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MRGRND )

        	IF ( .NOT. READ3( M2NAME, 'WSPD10', 1, SDATE, STIME,
     &						TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // M2NAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		CALL INTERPOL( TMP, H, MWSDP10 )


		CALL NEXTIME( SDATE,STIME,TSTEP3D )

	END DO !H = 1, 24

!	C....... Close file
        IF( .NOT. CLOSE3( M2NAME ) ) THEN
            MESG = 'Could not close file ' // M2NAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
	
!............	Print out min and max values for atmospheric onput fields
	WRITE( MESG,* ) ''
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'Successfully read meteo fields for day ', CDATE
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'TEMP2  ',MINVAL(ATMIN( :,:,:,MTEMP2 ) ),
     &				  MAXVAL(ATMIN( :,:,:,MTEMP2 ) ), ' K'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'PRSFC  ',MINVAL(ATMIN( :,:,:,MPRSFC ) ),
     &			 	  MAXVAL(ATMIN( :,:,:,MPRSFC ) ), ' Pa'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'USTAR  ',MINVAL(ATMIN( :,:,:,MUSTAR ) ),
     &				  MAXVAL(ATMIN( :,:,:,MUSTAR ) ), ' m/s'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'RADYNI ',MINVAL(ATMIN( :,:,:,MRADYN ) ),
     &				  MAXVAL(ATMIN( :,:,:,MRADYN ) ), ' m/s'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'RGRND  ',MINVAL(ATMIN( :,:,:,MRGRND ) ),
     &				  MAXVAL(ATMIN( :,:,:,MRGRND ) ), ' W/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'TA     ',MINVAL(ATMIN( :,:,:,MTA ) ),
     &				  MAXVAL(ATMIN( :,:,:,MTA ) ), ' K'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'PRES   ',MINVAL(ATMIN( :,:,:,MPRES ) ),
     &				  MAXVAL(ATMIN( :,:,:,MPRES ) ), ' Pa'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'ZF     ',MINVAL(ATMIN( :,:,:,MZF ) ),
     &				  MAXVAL(ATMIN( :,:,:,MZF ) ), ' m'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) ''
	CALL M3MSG2( MESG )

!##########################################################################################

C.........  End of subroutine
	END SUBROUTINE READMET



	SUBROUTINE READCHEM( CDATE )

	IMPLICIT NONE

C...........
!
! Chemistry data is read for 24 time steps once each day
!
! INPUT:
!	INTEGER 	     SDATE	Expected start date of input files
!	Reads file names from environmental variables CONCFILE, DRYDEP, WETDEP
!
! OUTPUT:
!	REAL CA( ROW,COL,HOUR,VAR ) 	Chemical Species in upper most layer
!
! HISTORY:
!
! 24.12.2012	Read concentration and deposition fields
!
! 25.12.2012	Convert units to ng/m**3 and ng/m**2
!
! 24.05.2013	Moved to module cin.f
!		Added checksum output (MINVAL, MAXVAL)
!
C...........

C...........   INCLUDES:

	INCLUDE 'C_model_inc'	!  Global variables
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(16)      PROMPTMFILE
        LOGICAL            SETENVVAR

        EXTERNAL  PROMPTMFILE, SETENVVAR

C...........   INPUT/OUTPUT VARIABLES and their descriptions:
        INTEGER, INTENT(IN) :: CDATE     ! start Julian date (YYYYDDD)

C...........   LOCAL VARIABLES and their descriptions:
	INTEGER		H, R, C	! HOUR, ROW, COL
        INTEGER         SDATE   ! start Julian date (YYYYDDD)
        INTEGER         STIME   ! start time (HHMMSS)

        CHARACTER( 256 ) CONCFILE  !  location of CMAQ concentration file
        CHARACTER( 256 ) DRYDEP    !  location of CMAQ dry deposition file
        CHARACTER( 256 ) WETDEP    !  location of CMAQ wet deposition file
        CHARACTER( NAMLEN3 ) CNAME     !  logical name for CONCFILE
        CHARACTER( NAMLEN3 ) DNAME     !  logical name for DRYDEP
        CHARACTER( NAMLEN3 ) WNAME     !  logical name for WETDEP

        REAL   , ALLOCATABLE :: TMP( :,: )  	 ! temporal data field
        REAL   , ALLOCATABLE :: ADD( :,: )  	 ! data array for aggregation of particle bins

C...........   Local parameters
	INTEGER        	    IOS 	          ! temporay IO status
        INTEGER             LDEV    	 	  ! log-device
        CHARACTER( 256 )    MESG                  ! message buffer 
        CHARACTER( 16 )  :: PROGNAME = 'READCHEM' ! program name


C***********************************************************************
C   begin body of subroutine READCHEM

C........... initiate IO API
        LDEV = INIT3()

!##########################################################################################

	CALL M3MSG2( DASHLINE )
	WRITE( MESG,* ) 'Reading atmospheric chemistry for day ',
     &			 CDATE
	CALL M3MSG2( MESG )
	CALL M3MSG2( '' )

       CALL ENVSTR( 'CONCFILE', 'CONCFILE file',
     &					 '', CONCFILE, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for CONCFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	WRITE( CONCFILE, '(A,I7)' ), TRIM( CONCFILE ), CDATE

        IF( .NOT. SETENVVAR( 'CCONCFILE', CONCFILE ) ) THEN
            MESG = 'Could not set environment variable CCONCFILE = '
     &						    // CONCFILE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


       CALL ENVSTR( 'DRYDEP', 'DRYDEP file',
     &					 '', DRYDEP, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for DRYDEP'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	WRITE( DRYDEP, '(A,I7)' ), TRIM( DRYDEP ), CDATE

        IF( .NOT. SETENVVAR( 'CDRYDEP', DRYDEP ) ) THEN
            MESG = 'Could not set environment variable CDRYDEP = '
     &						    // DRYDEP
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


       CALL ENVSTR( 'WETDEP', 'WETDEP file',
     &					 '', WETDEP, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for WETDEP'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	WRITE( WETDEP, '(A,I7)' ), TRIM( WETDEP ), CDATE

        IF( .NOT. SETENVVAR( 'CWETDEP', WETDEP ) ) THEN
            MESG = 'Could not set environment variable CWETDEP = '
     &						    // WETDEP
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	C....... Open concentration file
        MESG = 'Enter logical name for the CONC file'
        CNAME = 'CCONCFILE'
        CNAME = PROMPTMFILE( MESG, FSREAD3, CNAME, PROGNAME )

!	C....... Read global attributes
        IF ( .NOT. DESC3( CNAME ) ) THEN
            MESG = 'Could not get description of file ' // CNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	C....... Save dimensions in local variables
	SDATE = CDATE
	STIME = STIME3D
	NROWS = NROWS3D
	NCOLS = NCOLS3D

!	C....... Allocate variables
	ALLOCATE( TMP( NROWS,NCOLS ), STAT=IOS )
!	CALL CHECKMEM( IOS, 'TMP', PROGNAME )
	ALLOCATE( ADD( NROWS,NCOLS ), STAT=IOS )
!	CALL CHECKMEM( IOS, 'ADD', PROGNAME )

	DO H = 1, 24

!		C....... Read data and save interpolated fields in ATMIN array
        	IF ( .NOT. READ3( CNAME, 'HG', 1, SDATE, STIME,
     &				 	TMP ) ) THEN
	            MESG = 'ERROR: Reading timestep from file ' // CNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	        END IF

!		C....... convert units to ng/m**3
		DO R = 1, NROWS
		    DO C = 1, NCOLS
			TMP( R,C ) = TMP( R,C ) * 8E6	! ppm to ng/m**3
		    END DO
		END DO

		CALL INTERPOL( TMP, H, MCHG0 )

		CALL NEXTIME( SDATE,STIME,TSTEP3D )

	END DO !H = 1, 24

!	C....... Close file
        IF( .NOT. CLOSE3( CNAME ) ) THEN
            MESG = 'Could not close file ' // CNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!##########################################################################################

!	C....... Open dry deposition file
        MESG = 'Enter logical name for the DRYDEP file'
        DNAME = 'CDRYDEP'
        DNAME = PROMPTMFILE( MESG, FSREAD3, DNAME, PROGNAME )

!	C....... Read global attributes
        IF ( .NOT. DESC3( DNAME ) ) THEN
            MESG = 'Could not get description of file ' // DNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!	C....... DRYDEP file begins at 010000 and only contains 24 time steps
	SDATE = CDATE
	STIME = STIME3D

!	C....... Check consistency
	IF( NCOLS3D .NE. NCOLS ) THEN
            MESG = 'ERROR: Dimension NCOLS in file ' // DNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	ELSE IF( NROWS3D .NE. NROWS ) THEN
            MESG = 'ERROR: Dimension NROWS in file ' // DNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

	DO H = 1, 24

!		C....... Read data and save interpolated fields in ATMIN array
        	IF ( .NOT. READ3( DNAME, 'HG', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // DNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			TMP( R,C ) = TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

		CALL INTERPOL( TMP, H, MDHG0 )

        	IF ( .NOT. READ3( DNAME, 'HGIIGAS', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // DNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			TMP( R,C ) = TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

		CALL INTERPOL( TMP, H, MDHG2 )

        	IF ( .NOT. READ3( DNAME, 'APHGI', 1, SDATE, STIME,
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // DNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			ADD( R,C ) = TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

        	IF ( .NOT. READ3( DNAME, 'APHGJ', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // DNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			ADD( R,C ) = ADD( R,C ) + TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

        	IF ( .NOT. READ3( DNAME, 'APHGK', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // DNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			ADD( R,C ) = ADD( R,C ) + TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

		CALL INTERPOL( ADD, H, MDHGP )

		CALL NEXTIME( SDATE,STIME,TSTEP3D )

	END DO !H = 1, 24

!	C....... Close file
        IF( .NOT. CLOSE3( DNAME ) ) THEN
            MESG = 'Could not close file ' // DNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!##########################################################################################

!	C....... Open file
        MESG = 'Enter logical name for the WETDEP file'
        WNAME = 'CWETDEP'
        WNAME = PROMPTMFILE( MESG, FSREAD3, WNAME, PROGNAME )

!	C....... Read global attributes
        IF ( .NOT. DESC3( WNAME ) ) THEN
            MESG = 'Could not get description of file ' // WNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	SDATE = CDATE
	STIME = STIME3D

!	C....... Check consistency
	IF( NCOLS3D .NE. NCOLS ) THEN
            MESG = 'ERROR: Dimension NCOLS in file ' // WNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	ELSE IF( NROWS3D .NE. NROWS ) THEN
            MESG = 'ERROR: Dimension NROWS in file ' // WNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
	END IF

	DO H = 1, 24
!		C....... Read data
        	IF ( .NOT. READ3( WNAME, 'HG', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // WNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			TMP( R,C ) = TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

		CALL INTERPOL( TMP, H, MWHG0 )

        	IF ( .NOT. READ3( WNAME, 'HGIIGAS', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // WNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			TMP( R,C ) = TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

		CALL INTERPOL( TMP, H, MWHG2 )

        	IF ( .NOT. READ3( WNAME, 'APHGI', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // WNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			ADD( R,C ) = TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

        	IF ( .NOT. READ3( WNAME, 'APHGJ', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // WNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			ADD( R,C ) = ADD( R,C ) + TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

        	IF ( .NOT. READ3( WNAME, 'APHGK', 1, SDATE, STIME, 
     &					TMP ) ) THEN
        	    MESG = 'ERROR: Reading timestep from file ' // WNAME
        	    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        	END IF

		DO R = 1, NROWS
		    DO C = 1, NCOLS
			ADD( R,C ) = ADD( R,C ) + TMP( R,C ) * 1E8	! kg/hectar -> ng/m**2
		    END DO
		END DO

		CALL INTERPOL( ADD, H, MWHGP )

		CALL NEXTIME( SDATE,STIME,TSTEP3D )

	END DO !H = 1, 24

!	C....... Close file
        IF( .NOT. CLOSE3( WNAME ) ) THEN
            MESG = 'Could not close file ' // WNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	WRITE( MESG,* ) ''
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'Successfully read atmospheric chemistry for day ',
     &			 CDATE
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'C(Hg0) ',MINVAL(ATMIN( :,:,:,MCHG0 ) ),
     &			 	  MAXVAL(ATMIN( :,:,:,MCHG0 ) ), ' ng/m**3'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'D(Hg0) ',MINVAL(ATMIN( :,:,:,MDHG0 ) ),
     &				  MAXVAL(ATMIN( :,:,:,MDHG0 ) ), ' ng/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'D(Hg2) ',MINVAL(ATMIN( :,:,:,MDHG2 ) ),
     &				  MAXVAL(ATMIN( :,:,:,MDHG2 ) ), ' ng/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'D(HgP) ',MINVAL(ATMIN( :,:,:,MDHGP ) ),
     &				  MAXVAL(ATMIN( :,:,:,MDHGP ) ), ' ng/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'W(Hg0) ',MINVAL(ATMIN( :,:,:,MWHG0 ) ),
     &				  MAXVAL(ATMIN( :,:,:,MWHG0 ) ), ' ng/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'W(Hg2) ',MINVAL(ATMIN( :,:,:,MWHG2 ) ),
     &				  MAXVAL(ATMIN( :,:,:,MWHG2 ) ), ' ng/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'W(HgP) ',MINVAL(ATMIN( :,:,:,MWHGP ) ),
     &				  MAXVAL(ATMIN( :,:,:,MWHGP ) ), ' ng/m**2'
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) ''
	CALL M3MSG2( MESG )
!##########################################################################################

C.........  End of subroutine:
	END SUBROUTINE READCHEM



	SUBROUTINE RIVERIN(  )

	USE CPARAM

	IMPLICIT NONE

C...........
!
! Reads ASCII file with monthly river inputs
!
! INPUT:
!
! OUTPUT:
!	REAL RIVER( ROW,COL,MONTH )	! monthly total emissions
!
! HISTORY:
!
! 01.06.2013	Implementation of river input subroutine
!
C...........

C...........   INCLUDES:

	INCLUDE 'C_model_inc'	!  Global variables
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C...........   LOCAL VARIABLES and their descriptions:
	INTEGER		I, ID, R, C, U, V	 ! indices
        INTEGER         SDATE   		 ! start Julian date (YYYYDDD)
        INTEGER         STIME   		 ! start time (HHMMSS)

        CHARACTER( 265 )     RNAME      	 ! logical name for CONCFILE
	CHARACTER( 265 )     HEAD		 ! file header information
	LOGICAL		     EXS		 ! unit exists
	LOGICAL		     OPN		 ! unit already opened
	LOGICAL		 ::  EOL = .FALSE.	 ! end of line flag

	REAL		 ::	 TOTAL = 0.
	REAL		 ::	 MAXIMUM = -999.
	REAL			 YEAR		 ! data report year
	REAL			 PERC		 ! percentage of annual
	REAL		         ANNUAL		 ! annual input
	REAL, DIMENSION( 12 ) :: MONTHLY	 ! monthly input

C...........   Local parameters
	INTEGER        	    IOS 	         ! temporay IO status
        INTEGER             LDEV    	 	 ! log-device
        CHARACTER( 256 )    MESG                 ! message buffer 
        CHARACTER( 265 )    DESC                 ! river description
        CHARACTER( 16 )  :: PROGNAME = 'RIVERIN' ! program name


C***********************************************************************
C   begin body of subroutine RIVERIN

	CALL M3MSG2( DASHLINE )
        WRITE( MESG,* ) 'Reading river input'
        CALL M3MSG2( MESG )
        WRITE( MESG,* ) ''
        CALL M3MSG2( MESG )

!       C....... Allocate variables
	IF( ALLOCATED( RIVER ) ) THEN
	    DEALLOCATE( RIVER )
	END IF
		
        ALLOCATE( RIVER( m,n,12 ), STAT=IOS )
!       CALL CHECKMEM( IOS, 'MF', PROGNAME )
        RIVER = 0

!....... Read file name for North Sea river data
        CALL ENVSTR( 'NSRIVERS', 'River Hg inflow from North Sea',
     &					'', RNAME, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for RIVERIN'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!....... Check for unused unit to open file
        U = 1
        DO V = 7,99
            INQUIRE( UNIT=V,EXIST=EXS,OPENED=OPN )
            IF( EXS .AND. .NOT. OPN ) THEN
                U = V
            END IF
        END DO

!....... Open file
        OPEN( U,FILE=RNAME,STATUS='OLD',ACTION='READ' )

!....... Read monthly values from file
	READ( U,'(A)' ) HEAD

	DO WHILE( .NOT. EOL )
	    READ( U,*, IOSTAT=IOS ) ID, R, C, ANNUAL, MONTHLY, YEAR ! ANNUAL [kg/a] MONTHLY [kg/month]
                IF( IOS .LT. 0 ) THEN
                    EOL = .TRUE.
                    CYCLE
                END IF

		IF( YEAR .EQ. 2004 ) THEN
		    IF( ANNUAL .GT. MAXIMUM ) THEN
			MAXIMUM = ANNUAL
		    END IF
		    TOTAL = TOTAL + ANNUAL

		    DO I = 1, 12
			RIVER( R, C, I ) = RIVER( R, C, I ) + MONTHLY( I )	! [kg/month]
		    END DO

                    WRITE( MESG,* ) ID, R, C, ANNUAL
                    CALL M3MSG2( MESG )
		END IF
	END DO

!....... Print info
	WRITE( MESG,* ) TOTAL, 'kg/a  max:', MAXIMUM
	CALL M3MSG2( MESG )

!....... Close file
        CLOSE( U )
	CALL M3MSG2( '' )

!....... Read Baltic Sea river input 
	EOL = .FALSE.
	TOTAL = 0
	MAXIMUM = 0

!....... Read file name for Baltic Sea river data
        CALL ENVSTR( 'BSRIVERS', 'River Hg inflow from Baltic Sea',
     &					'', RNAME, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for RIVERIN'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!....... Open file
        OPEN( U,FILE=RNAME,STATUS='OLD',ACTION='READ' )

!...... Read monthly values from file
	READ( U,'(A)' ) HEAD

	DO WHILE( .NOT. EOL )
	    READ( U,*, IOSTAT=IOS ) R, C, PERC, ANNUAL, DESC	! [t/a]

            IF( IOS .LT. 0 ) THEN
                EOL = .TRUE.
                CYCLE
            END IF

	    ANNUAL = ANNUAL * 1000. * PERC / 100. ! [kg/a]
	    WRITE( MESG,* ) R, C, ANNUAL, TRIM( DESC )
	    CALL M3MSG2( MESG )
	    IF( ANNUAL .GT. MAXIMUM ) THEN
		MAXIMUM = ANNUAL
	    END IF
	    TOTAL = TOTAL + ANNUAL

	    DO I = 1, 12
		RIVER( R, C, I ) = RIVER( R, C, I ) + ANNUAL / 12. ! [kg/month]
	    END DO
	END DO

!....... Print info
	WRITE( MESG,* ) TOTAL, 'kg/a max:', MAXIMUM
	CALL M3MSG2( MESG )

!....... Close file
        CLOSE( U )
	CALL M3MSG2( '' )

	END SUBROUTINE RIVERIN


C.........  End of module
	END MODULE CIN
