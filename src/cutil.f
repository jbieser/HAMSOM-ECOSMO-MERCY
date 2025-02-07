	MODULE CUTIL

	IMPLICIT NONE

	CONTAINS


!--------------------------------------------------------------------------
! density and enthalpy, based on the 48-term expression for density
!--------------------------------------------------------------------------

!==========================================================================
	function gsw_rho(sa,ct,p)
!==========================================================================
!
!	This subroutine was extraqcted from the gsw_oceanographic_toolbox
!
!  Calculates in-situ density from Absolute Salinity and Conservative 
!  Temperature, using the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (McDougall et al., 2011).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho  : in-situ density (48 term equation)

	implicit none

	real, parameter :: v01 =  9.998420897506056d2
	real, parameter :: v02 =  2.839940833161907d0
	real, parameter :: v03 = -3.147759265588511d-2
	real, parameter :: v04 =  1.181805545074306d-3
	real, parameter :: v05 = -6.698001071123802d0
	real, parameter :: v06 = -2.986498947203215d-2
	real, parameter :: v07 =  2.327859407479162d-4
	real, parameter :: v08 = -3.988822378968490d-2
	real, parameter :: v09 =  5.095422573880500d-4
	real, parameter :: v10 = -1.426984671633621d-5
	real, parameter :: v11 =  1.645039373682922d-7
	real, parameter :: v12 = -2.233269627352527d-2
	real, parameter :: v13 = -3.436090079851880d-4
	real, parameter :: v14 =  3.726050720345733d-6
	real, parameter :: v15 = -1.806789763745328d-4
	real, parameter :: v16 =  6.876837219536232d-7
	real, parameter :: v17 = -3.087032500374211d-7
	real, parameter :: v18 = -1.988366587925593d-8
	real, parameter :: v19 = -1.061519070296458d-11
	real, parameter :: v20 =  1.550932729220080d-10
	real, parameter :: v21 =  1.0d0
	real, parameter :: v22 =  2.775927747785646d-3
	real, parameter :: v23 = -2.349607444135925d-5
	real, parameter :: v24 =  1.119513357486743d-6
	real, parameter :: v25 =  6.743689325042773d-10
	real, parameter :: v26 = -7.521448093615448d-3
	real, parameter :: v27 = -2.764306979894411d-5
	real, parameter :: v28 =  1.262937315098546d-7
	real, parameter :: v29 =  9.527875081696435d-10
	real, parameter :: v30 = -1.811147201949891d-11
	real, parameter :: v31 = -3.303308871386421d-5
	real, parameter :: v32 =  3.801564588876298d-7
	real, parameter :: v33 = -7.672876869259043d-9
	real, parameter :: v34 = -4.634182341116144d-11
	real, parameter :: v35 =  2.681097235569143d-12
	real, parameter :: v36 =  5.419326551148740d-6
	real, parameter :: v37 = -2.742185394906099d-5
	real, parameter :: v38 = -3.212746477974189d-7
	real, parameter :: v39 =  3.191413910561627d-9
	real, parameter :: v40 = -1.931012931541776d-12
	real, parameter :: v41 = -1.105097577149576d-7
	real, parameter :: v42 =  6.211426728363857d-10
	real, parameter :: v43 = -1.119011592875110d-10
	real, parameter :: v44 = -1.941660213148725d-11
	real, parameter :: v45 = -1.864826425365600d-14
	real, parameter :: v46 =  1.119522344879478d-14
	real, parameter :: v47 = -1.200507748551599d-15
	real, parameter :: v48 =  6.057902487546866d-17 

	real :: sa, ct, p, sqrtsa, v_hat_denominator, v_hat_numerator,
     &		      gsw_rho

	sqrtsa = sqrt(sa)

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
     &         + sa*(v05 + ct*(v06 + v07*ct)
     &         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
     &         + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
     &         + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

	v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
     &     + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
     &     + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
     &     + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
     &     + sa*(v41 + v42*ct)
     &     + p*(v43 + ct*(v44 + v45*ct + v46*sa)
     &     + p*(v47 + v48*ct)))

	gsw_rho = v_hat_denominator/v_hat_numerator

	return
	END function gsw_rho

!--------------------------------------------------------------------------





!	Calculates area of rectangle on a sphere
	SUBROUTINE RECTAREA( AREA, EAST, NORTH, WEST, SOUTH, VERBOSE )

	IMPLICIT NONE

	REAL, INTENT(OUT)	:: AREA		! area of polygon		[km]
	REAL, INTENT(IN)	:: EAST		! estern side of rectangle	[km]
	REAL, INTENT(IN)	:: NORTH	! northern side of rectangle	[km]
	REAL, INTENT(IN)	:: WEST		! western side of rectangle	[km]
	REAL, INTENT(IN)	:: SOUTH	! southern side of rectangle	[km]
	INTEGER, INTENT(IN)	:: VERBOSE

	REAL, PARAMETER		:: R = 6371			! radius of spheric earth	 [km]
	REAL, PARAMETER		:: R2 = 40589641		! square of radius		 [km**2]
	REAL, PARAMETER		:: PI = 3.14159265358979323846	! The infamous PI
	REAL, PARAMETER		:: CIRCUM = 40015.086796021	! circumference (2*PI*R)	 [km]
	REAL, PARAMETER		:: DPR = 57.29577951308232	! degrees per radian (180/PI)	 [deg/rad]
	REAL, PARAMETER		:: DPK = 0.008993216059187	! degree per km	(360/U)		 [deg/km]

	REAL	A,B,C		! sides of triangle
	REAL	S		! semiperimeter
	REAL	T1,T2,T3,T4	! helper variables
	REAL	E		! spherical excess
	REAL	ATRI1		! area of first triangel
	REAL	ATRI2		! area of second triangel

!	C. Calculate area of first triangle
	A = EAST / R					! [rad]
	B = NORTH / R					! [rad]
	C = SQRT( EAST * EAST + NORTH * NORTH ) / R	! [rad]
	if (verbose .eq. 1) print*,A,B,C
	S = 0.5 * ( A + B + C )
	T1 = TAN( 0.5 * S )
	T2 = TAN( 0.5 * ( S - A ) )
	T3 = TAN( 0.5 * ( S - B ) )
	T4 = TAN( 0.5 * ( S - C ) )
	E = 4 * ABS( ATAN( SQRT( ABS( T1 * T2 * T3 * T4 ) ) ) )
	if (verbose .eq. 1)print*,S,T1,T2,T3,T4,E
	ATRI1 = E * R2

!	C. Calculate area of second triangle
	A = WEST / R					! [rad]
	B = SOUTH / R					! [rad]
	C = SQRT( WEST * WEST + SOUTH * SOUTH ) / R	! [rad]
	if (verbose .eq. 1)print*,A,B,C
	S = 0.5 * ( A + B + C )
	T1 = TAN( 0.5 * S )
	T2 = TAN( 0.5 * ( S - A ) )
	T3 = TAN( 0.5 * ( S - B ) )
	T4 = TAN( 0.5 * ( S - C ) )
	E = 4 * ABS( ATAN( SQRT( ABS( T1 * T2 * T3 * T4 ) ) ) )
	if (verbose .eq. 1)print*,S,T1,T2,T3,T4,E
	ATRI2 = E * R2

!	C. Return area of polygon
	AREA = ATRI1 + ATRI2
	if (verbose .eq. 1)print*,AREA,ATRI1,ATRI2
	END SUBROUTINE RECTAREA



	CHARACTER( 265 ) FUNCTION GETENVVAR( VNAME )

	IMPLICIT NONE

        LOGICAL            ENVYN
        INTEGER            ENVINT

      	EXTERNAL  ENVYN, ENVINT

	CHARACTER( 265 ), INTENT(IN) :: VNAME

        INTEGER       	   IOS                     ! i/o status
        CHARACTER(256)     MESG              	   ! message field
        CHARACTER(256)  :: PROGNAME  = 'GETENVVAR' ! function name

!       CALL ENVSTR( 'TOPOFILE', 'bathymetry file',
!     &					 '', topofile, IOS )

        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for ' // TRIM( VNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
	
	RETURN

	END FUNCTION GETENVVAR

c-----------------------------------------------------------------------
      subroutine Cdeco1d3d (u,ucomp,ntot,xdef,IWET,LDEP,LAZC,INDEND) 
c-----------------------------------------------------------------------
c      1-d-field ->  3-d-field 
c-----------------------------------------------------------------------
	IMPLICIT NONE

      include 'C_model_inc'
!      INCLUDE 'C_declaration_inc'

	INTEGER I,J,K,NWET,LW,LWE,LWA,LUMP

	REAL, INTENT(OUT), dimension(n,m,ilo)  :: U	! 3D output array
	REAL, INTENT(IN), dimension(ntot) :: ucomp	! 1D input array
	INTEGER NTOT					! dimension of UCOMP (usually NDREI)
	REAL XDEF						! value to fill in dry grid cells

	INTEGER IWET( : )
	INTEGER LDEP( : )
	INTEGER LAZC( : )
	INTEGER INDEND( : )

      lwe = 0 
      nwet = 0 
      do 2 k=1,n 
        do i=1,m 
        do j=1,ilo 
         u(i,k,j) = xdef 
        enddo
        enddo
        lwa = lwe+1 
        lwe = indend(k) 
        do 2 lw=lwa,lwe 
          lump = lazc(lw) 
          i = iwet(lw) 
          do j=1,lump 
             nwet = nwet+1 
             u(k,m-i+1,j) = ucomp(nwet) 
          enddo
    2 continue 


      return 
      end SUBROUTINE CDECO1D3D



	SUBROUTINE DIMSWAP2D( IN, OUT )

	IMPLICIT NONE

	REAL, INTENT(IN) :: IN( :,: )
	REAL, INTENT(INOUT):: OUT( :,: )

	INTEGER R,C,NR,NC

	NR = SIZE( IN,DIM=1 )
	NC = SIZE( IN,DIM=2 )

	DO R = 1, NR
	    DO C = 1, NC
		OUT( C,R ) = IN( R,C )
	    END DO
	END DO

	END SUBROUTINE DIMSWAP2D

	SUBROUTINE DIMSWAP3D( IN, OUT )

	IMPLICIT NONE

	REAL, INTENT(IN) :: IN( :,:,: )
	REAL, INTENT(INOUT):: OUT( :,:,: )

	INTEGER R,C,NR,NC,T,NT

	NR = SIZE( IN,DIM=1 )
	NC = SIZE( IN,DIM=2 )
	NT = SIZE( IN,DIM=3 )

	DO R = 1, NR
	    DO C = 1, NC
		DO T = 1,NT
		    OUT( C,R,T ) = IN( R,C,T )
		END DO
	    END DO
	END DO

	END SUBROUTINE DIMSWAP3D


	SUBROUTINE DIMTURN3D( IN, OUT, FAC )

	IMPLICIT NONE

	REAL, INTENT(IN) :: IN( :,:,: )
	REAL, INTENT(INOUT):: OUT( :,:,: )
	REAL, INTENT(IN) :: FAC

	INTEGER R,C,NR,NC,T,NT

	NR = SIZE( IN,DIM=1 )
	NC = SIZE( IN,DIM=2 )
	NT = SIZE( IN,DIM=3 )

	DO R = 1, NR
	    DO C = 1, NC
		DO T = 1,NT
		    OUT( C,NR-R+1,T ) = IN( R,C,T ) * FAC
		END DO
	    END DO
	END DO

	END SUBROUTINE DIMTURN3D


	SUBROUTINE DIMTURN2D( IN, OUT, FAC )

	IMPLICIT NONE

	REAL, INTENT(IN) :: IN( :,: )
	REAL, INTENT(INOUT):: OUT( :,: )
	REAL, INTENT(IN) :: FAC

	INTEGER R,C,NR,NC

	NR = SIZE( IN,DIM=1 )
	NC = SIZE( IN,DIM=2 )

	DO R = 1, NR
	    DO C = 1, NC
		OUT( C,NR-R+1 ) = IN( R,C ) * FAC
	    END DO
	END DO

	END SUBROUTINE DIMTURN2D

	END MODULE CUTIL
