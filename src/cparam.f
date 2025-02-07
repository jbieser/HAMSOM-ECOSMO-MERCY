	MODULE CPARAM

	IMPLICIT NONE

C...........
!
! DESC:         Global parameters
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........


C........  GLOBAL PARAMETERS
	REAL, PARAMETER :: ITS = 600.	! INTERNAL TIME STEP for airsea and chem [s]
	REAL, PARAMETER	:: F_GRAVITY = 9.80665	! Earth's acceleration 			 [m/s**2]
	REAL, PARAMETER :: RHO_0 = 1030		! Standard density of ocean water	 [kg/m**3]

C........  GLOBAL INDICES
	INTEGER, PARAMETER	:: M_ROWS = 177
	INTEGER, PARAMETER	:: N_COLS = 207
	INTEGER, PARAMETER	:: I_LAYS = 20
	INTEGER, PARAMETER	:: I_KHOR = 8216
	INTEGER, PARAMETER	:: I_KASOR = 8206
	INTEGER, PARAMETER	:: I_NDREI = 82108

C........  GLOBAL VARIABLES
	REAL		   MSTEP	! ECOMSO timestep
	REAL		:: NAN = 0.	! NaN value - created in subroutine to trick the compiler


!	C. Variables containing infroamtion about date and time
	INTEGER		   CDATE		! Date <YYYYDDD>
	INTEGER		   CYEAR		! Year <YYYY>
	INTEGER		   CJDAY		! Julian day <DDD>
	INTEGER		   CDECA		! Decade <YY>
	INTEGER		:: CTIME = 000000	! Time <HHMMSS>
	INTEGER		:: CHOUR = 1		! Hour 1-24
	INTEGER		:: TSTEP = 002000  	! Time Step <HHMMSS> = deltat*10/6

	INTEGER		   EYEAR		! Last Year to process <YYYY>
!	INTEGER		   EDATE		! End date <YYYYDDD>
!	LOGICAL		:: VERBOSE = .FALSE.

!	C. Output
	CHARACTER( 60 ) :: DASHLINE =
     &	'############################################################'

	CONTAINS

	SUBROUTINE INIT_CPARAM

	IMPLICIT NONE

	print*,'Initializing chemistry module'

!	NAN = NAN / 0.

!	MSTEP = 3600 * 24 / ( IPER * IVIERT * ITS )

	END SUBROUTINE INIT_CPARAM

	END MODULE CPARAM
