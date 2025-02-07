	MODULE CDATA

	USE CPARAM
        USE COUT, ONLY : NHGOUT3D, NHGOUT2D

	IMPLICIT NONE

C...........
!
! DESC:         Global arrays used by the chemistry routines
!               Contains subroutine CDATA_INIT to allocate the chemistry fields
!
! HISTORY:      13.10.2021    Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........


C...........   PARAMETERS and their descriptions
	INTEGER, PARAMETER	:: NCVARS = 4

C...........   LOCAL VARIABLES and their descriptions:
!	C. dimensions: rows,cols,lays,vars
!TODO:	Rescue this poor array, which is held captive by the evil common blocks!
!		(maybe start a petition to the EU Parliament)
!		Tw <-- Tc( :,9 ) + Tc( :,1:4 )
!	REAL, ALLOCATABLE, DIMENSION( :,:,:,: ) :: Tw		! mercury in water		[ng/L]

!	C. dimensions: rows,cols,top/bottom,vars
!		Ti( :,:,1,: ) = Mercury inside the ice
	REAL, ALLOCATABLE, DIMENSION( :,:,:,: ) :: Ti	 	! mercury concentration in ice 	[ng/L]
	REAL, ALLOCATABLE, DIMENSION( : )       :: Ticevol	! ice volume of last time step 	[m**3]
	REAL, ALLOCATABLE, DIMENSION( :,: )     :: Tdep		! mercury deposited onto ice	[ng/m**2]
 
!	C. dimensions: rows,cols,layer,vars (burial layer = 2, exchange layer = 1)
	REAL, ALLOCATABLE, DIMENSION( :,:,:,: ) :: Ts		! mercury in sediment		[ng/m**2]

!	C.
	REAL, ALLOCATABLE, DIMENSION( : )   :: Trho		! salinity and temperatur dependent density [psu]
	REAL, ALLOCATABLE, DIMENSION( : )   :: Tpres		! water pressure			    [kg/m**3]
	REAL, ALLOCATABLE, DIMENSION( : )   :: Tarea		! grid cell area			    [m**2]
        REAL, ALLOCATABLE, DIMENSION( : )   :: Tdept            ! grid cell depth                           [m]    
	REAL, ALLOCATABLE, DIMENSION( :,: ) :: Thn		! heigt of lowest grid cell		    [m]
	INTEGER, ALLOCATABLE, DIMENSION( :,: ) :: Tdz		! number of vertical layers		    [layers]
	REAL, ALLOCATABLE, DIMENSION( :,:,: ) :: Tun		! vertically integrated u component in lowest grid cell	    [m**2/s]
	REAL, ALLOCATABLE, DIMENSION( :,:,: ) :: Tvn		! vertically integrated v component in lowest grid cell	    [m**2/s]


	CONTAINS

	SUBROUTINE INIT_CDATA( )

	INCLUDE 'C_model_inc'

	INTEGER			IOS	! status flag for ALLOCATE
	CHARACTER( 265 )	MESG	! output buffer
	INTEGER D
	REAL Kd, DOC

	CALL M3MSG2( DASHLINE )
	WRITE( MESG,* ) 'Initializing chemistry fields'
	CALL M3MSG2( MESG )

	ALLOCATE( Tpres( ndrei ), STAT=IOS )
        Tpres = 0.
	ALLOCATE( Trho( ndrei ), STAT=IOS )
        Trho = 0.
	ALLOCATE( Tarea( ndrei ), STAT=IOS )
        Tarea = 0.
        ALLOCATE( Tdept( ndrei ), STAT=IOS )
        Tdept = 0.
	ALLOCATE( Thn( m,n ), STAT=IOS )
	Thn = 0
	ALLOCATE( Tdz( m,n ), STAT=IOS )
	Tdz = 0
	ALLOCATE( Tun( m,n,ilo ), STAT=IOS )
	Tun = 0.
	ALLOCATE( Tvn( m,n,ilo ), STAT=IOS )
	Tvn = 0.
!	ALLOCATE( Tw( M,N,ILO,NCVARS ), STAT=IOS )
!	ALLOCATE( Ti( M,N,3,NCVARS ), STAT=IOS )
!	ALLOCATE( Ticevol( M,N ), STAT=IOS )
!	ALLOCATE( Ts( M,N,NCVARS ), STAT=IOS )

        ALLOCATE( Ti( M,N,1,NHGOUT3D ), STAT=IOS )
        Ti = 0.
        ALLOCATE( Ticevol( NDREI ), STAT=IOS )
        Ticevol = -999.
        ALLOCATE( Tdep( M,N ), STAT=IOS )
        Tdep = 0.

        ALLOCATE( Ts( M,N,2,NHGOUT2D ), STAT=IOS )
        Ts = 0. ! [ng/^2]

	CALL M3MSG2( MESG )

	END SUBROUTINE INIT_CDATA

	END MODULE CDATA
