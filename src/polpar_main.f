         program main
!c------------------------------------------------------------------------------------
!c      North Sea Eulerian Advection model   
!c      version: 04.2012 /  North Sea & Baltic Sea  
!c      setup according Schrum & Bachhaus, 1999 & Schrum et al., 2006,
!c      advektion algorithmn Barthel et al. (2012)
!c      input data: POM, DOM, Plankton from ECOSMO II Daewel & Schrum, 2012, submitted
!c      modelling pollution in the sea
!c
!c      Fr 13.04.2018 Mercury --> Includes complete marine mercury cycling
!c
!c--------------------------------------------------------------------------------------
        USE CPARAM      ! parameters for chemistry, partitioning, and air-sea exchange
        USE CDATA       ! fields (Tc),Ts,Ti for chemistry
        USE CIN         ! atmospheric input module
        USE COUT        ! chemistry output module

!      IMPLICIT NONE    !One day baby.....
       include 'C_model_inc'
       include 'C_loads_inc'
       include 'C_declaration_inc'
       include 'C_biomod_inc'

c------------------------------------------------------------------------------------------------
c    Variables needed for chem module and input file names - added by Johannes Bieser 24.12.2012
c-----------------------------------------------------------------------------------------------

C..........  INCLUDES:
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters (in modfileset)
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations 
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures (in modfileset)

C..........  EXTERNAL FUNCTIONS
        LOGICAL            ENVYN
        INTEGER            ENVINT

        EXTERNAL  ENVYN, ENVINT

C..........  openMP VARIABLES
        INTEGER         :: TID = 0
        INTEGER         BATCH_S, BATCH_E
        REAL            T_START, T_FINISH
        INTEGER         :: NTHREADS = 1
        INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

!     INTEGER :: TID = -1

        REAL, PARAMETER :: CMIN = 1.0E-06
C..........  VARIABLES FOR STROM

C..........  LOCAL VARIABLES:
        INTEGER            SN
	INTEGER		   CSTEP		! counter for chemistry loop

	INTEGER		:: VERBOSE = 0		! debug flag

!	REAL :: CDEPTH		! depth of grid cell	[m]
!	REAL :: CAREA		! area of grid cell 	[m**2]
	REAL		            DEP1,DEP2		! deposition		[ng/(L * ITS)]
        REAL                        ICE_FRAC, ICE_THIC
        REAL                        RIVER_LOAD

!       transform 1d array to 3d coordinates
        INTEGER, DIMENSION( NDREI ) :: ii_k
        INTEGER, DIMENSION( NDREI ) :: ii_i
        INTEGER, DIMENSION( NDREI ) :: ii_j
        INTEGER, DIMENSION( NDREI ) :: ii_j_max
        INTEGER, DIMENSION( NDREI ) :: ii_ldep
        INTEGER, DIMENSION( NDREI ) :: ii_lzac

        REAL            Dmin
        REAL            :: Faw   = 0.           ! helper for air-sea exchange
        REAL            :: Faw_SUM  = 0.        ! helper for hourly air-sea flux sum written to 2D output

        REAL            :: Kbur = 1.1575E-10    ! Burial rate [s-1]

!       These are variables for debugging purposes (the good stuff)
        REAL, DIMENSION( 1:IHGT+3 )      :: DBG_CHG
        REAL            :: DBG_DEM = 0.
        REAL            :: DBG_FAW = 0.
        REAL            :: DBG_SAX = 0.
        REAL            :: DBG_GEM = 0.
        REAL            :: DBG_SED = 0.
        REAL            :: DBG_RES = 0.
        REAL            :: DBG_PAM = 0.
        REAL            :: DBG_RADI = 0.
        REAL            :: DBG_FISH = 0.
        REAL, DIMENSION( 1:2 )       :: DBG_DEP
        REAL            :: DBG_NUM = 0.
!..........................................................................

        INTEGER            TMON

        INTEGER       	   IOS               	! i/o status
        CHARACTER(1024)     MESG              	! message field
	LOGICAL		:: FIRSTTIME = .TRUE.	! first timestep flag
        LOGICAL         :: READ_DAILY = .FALSE. ! input time step
        LOGICAL         :: READ_HOURLY = .TRUE. ! input time step

        INTEGER         :: igstep = 0
c-----------------------------------------------------------------------
C    End - Variables needed for chem module and input file names
c-----------------------------------------------------------------------


       common /atmosphere/ atmflux(m,n,nbio,12,5)

        include 'C_info_inc'  
        include 'C_files_inc'
        include 'C_logcontr_inc'
cc      include '/usr/include/netcdf.inc'
        include 'C_verticalgrid_inc'

c-----------------------------------------------------------------------
       print*,' '
       print*,'STARTING MECOSMO v2.1'
       print*,' '
c-----------------------------------------------------------------------
       include 'C_timesettings_inc'
       include 'C_initialize_inc'

c    Read file names from environment variables - added by Johannes Bieser 24.12.2012
c-----------------------------------------------------------------------
        DBG_CHG = 0.
	CALL INIT_CPARAM( )
	CALL INIT_CDATA( )
	MSTEP = 3600 * 24 / ( IPER * IVIERT * ITS )

	WRITE( MESG,* ) ''
	CALL M3MSG2( MESG )
	WRITE( MESG,* ) 'READING INPUT FILES'
	CALL M3MSG2( MESG )

        CDATE = ENVINT( 'SDATE', 'Start date <YYYYDDD>', -1, IOS )
	CYEAR = CDATE / 1000
	CDECA = CDATE - CDATE/100*100
	CJDAY = CDATE - CDATE/1000*1000

	EYEAR = ENVINT( 'EYEAR', 'Last year <YYYY>', -1, IOS )

        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for SDATE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'BIOFILE', 'bio input file', '', biofile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for BIOFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'FLUFILE', 'flux input file', '', flufile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for FLUFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'PHYFILE', 'phyto input file', '', phyfile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for PHYFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'OUTFILE', 'output file', '', outfile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for OUTFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'CONTROL', 'control file', '', control, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for CONTROL'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'GRIDINF', 'grid input file', '', gridinf, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for GRIDINF'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'TOPOFILE', 'bathymetry file', '', topofile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for TOPOFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

      CALL ENVSTR( 'CHEMOUT', 'chemistry output file', '', chemout, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for CHEMOUT'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        WRITE( MESG,* ) ''
        CALL M3MSG2( MESG )
c-----------------------------------------------------------------------
c    End - Read file names from environment variables
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c       print*,'  open output files'
c-----------------------------------------------------------------------

        nfugri = 88
        open(unit=nfugri,file=gridinf,form='unformatted')

        nin     = 23
        nbio_in = 24
        nflu_in = 25
        nout    = 33

        open(unit=99,file=home//'topinfo',form='formatted',
     &       status='unknown')
        open(unit=98,file=home//'seamask',form='formatted')

!c-----------------------------------------------------------------------
!c       print*,'  open input files'
!c-----------------------------------------------------------------------

      nfutop = 10
      open(unit=nfutop,file= topofile,form='formatted')

c-----------------------------------------------------------------------
      print*,'  read modelltopo             file : ',topofile
c-----------------------------------------------------------------------

      read(nfutop,2000) izet
 2000 format(18i4)
      close(nfutop)

c-----------------------------------------------------------------------
      print *,'compute topo-dependent grid values'
c      compute topo-dependent grid values and
c      time-dependent grid values with SR kotief and SR gatmit
c-----------------------------------------------------------------------

      call kotief (dz,izet,nhor,ntot)
      call gatmit (deltat)
      call gittop (ltief,dzbod)
      print *,'grid initialization done:kotief,gatmit,gittop'

c----------------------------------------------------------------------
c      set grid resolution in m
c----------------------------------------------------------------------

      do i = 1,m
        do j = 1,n
          dx(i,j)=dln(i)
          dy(i,j)=dl
        enddo
      enddo

!	set grid cell area using spherical trigonometry (added by Johannes Bieser 10.06.2013)
!       This should be used but needs checking 15.04.2018
	DO I = 1,M-1
	    DO J = 1,N
!		CALL RECTAREA( Tarea( : ), DL, DX( I,J ), DL, DX( I+1,J ), 0 )
	    END DO
	END DO

	DO J = 1,N
!	    CALL RECTAREA( Tarea( : ), DL, DX( I,J ), DL, DX( I,J ), 0 )
	END DO

c
c-----------------------------------------------------------------------
c      set dummy arrays, compute index fields 
c-----------------------------------------------------------------------

      lwe = 0
      nwet=0
      iidrei=0
    
      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i = iwet(lw)
          lump = lazc(lw)
          jjc(i,k) = 1
          iindex(i,k) = nwet
          id3sur(lw) = nwet+1
          ibotlay(i,k)=lump
          Tdz( I,K ) = LAZC( LW )
          Thn( I,K ) = LDEP( LW )
           do jj=1,lump
             nwet = nwet+1
             iidrei(jj,i,k)=nwet
             ii_k( nwet ) = k
             ii_i( nwet ) = i
             ii_j( nwet ) = jj
             ii_j_max( nwet ) = lump            ! number of vertical grid cells
             ii_ldep( nwet ) = ldep( lw )       ! depth of lowest grid cell
             ii_lzac( nwet ) = zac( lw )        ! Oberflächenauslenkung

!        IF( k .gt. 204 ) then
!          print*, "LARGER 205: ",i, k, ii_i(nwet-1), ii_k(nwet-1),jj,nwet
!82106 82107 82108
!        END IF

             !         C. Calculate grid cell area (spherical
             !         calculation seems unneccesary and expensive)
              IF( I .LT. M ) THEN
                 Tarea( nwet ) = ( DLN( I ) + DLN( I+1 ) ) / 2 * 11120     ! grid cell area [m**2]
!                 CALL RECTAREA( Tarea( nwet ), 11.120,DLN(I)/1000,11.120,DLN(I+1)/1000,0)
!                 Tarea( nwet ) = Tarea( nwet ) * 1E6 !rectarea gives km²
              ELSE
                 Tarea( nwet ) = DLN( I ) * 11120                          ! grid cell area [m**2]
!                 CALL RECTAREA( Tarea( nwet ), 11.120,DLN(I)/1000,11.120,DLN(I)/1000,0 )
!                 Tarea( nwet ) = Tarea( nwet ) * 1E6
              END IF

          enddo
           zic(lw) = izet(i,k)
          izet(i,k) = lw
        enddo
      enddo

!        do ii = 1, ndrei
!          print*, ii, ii_i(ii),ii_k(ii),ii_j(ii),
!     &            ii_j_max(ii),ii_ldep(ii)!,ii_lzac(ii)
!        end do

      write(99,*) iwet,indend,lazc,ldep
      write (99,*) ltief,dzbod
      write(99,*) jjc
      close(99)
      
      do i=1,m
        write(98,'(i3)') (jjc(i,j),j=1,n)
      enddo
      close(98)

      call setice(iindex,izet,m,n)

c-----------------------------------------------------------------------
c    river loads/ discharge
c-----------------------------------------------------------------------
 1001    continue
              
      include 'C_discharges_inc'
c-----------------------------------------------------------------------
      print*,'  write grid info             file : ',gridinf
c-----------------------------------------------------------------------

      write (nfugri) dt,m,n,ilo
      write (nfugri) dz,nhor,ntot,iwet,ldep,lazc,indend,islab
      write (nfugri) dlr,rdln,dlvo,dlvu
      close (nfugri)

c***********************************************************************
      IJAHR = CDATE/1000
      IYEARS_START = IJAHR - IJAHR/100*100
      MJAR = IYEARS_START
	print*, CDATE,IJAHR,IYEARS_START,MJAR,EYEAR
      MJAR = MJAR - 1
      IF( MJAR .EQ. 100 ) THEN
          MJAR = 01
      END IF
      CYEAR = CYEAR -1
      DO 5099  IIIYEAR = IJAHR,EYEAR

	CYEAR = CYEAR + 1

	IF( MJAR .EQ. 99 ) THEN
	    MJAR = 00
	ELSE
	    MJAR = MJAR + 1
	END IF
	CALL M3MSG2( '' )
	CALL M3MSG2( '++++++++++++++++++++++++++++++++++++++++++++++++++' )
	WRITE( MESG,'(A20,I4,A2,I2,A2,I7,A2,I4)' ) 'Now starting year:  '
     &				,IIIYEAR,'  ',MJAR,'  ',CDATE,'  ',CYEAR
	CALL M3MSG2( MESG )
	CALL M3MSG2( '++++++++++++++++++++++++++++++++++++++++++++++++++')
	CALL M3MSG2( '' )

	IF( .NOT. IIIYEAR .EQ. CYEAR ) THEN
	    CYEAR = IIIYEAR
	    CDATE = IIIYEAR * 1000 + 1
	END IF

!	C. Open new file for each year and allocate output arrays
        print*, 'Create new annual output files'
	CALL WRITE_CHEM( .TRUE. )

! TODO:	Move this out of the main code
!	C. Read IC before first timestep
	IF( FIRSTTIME ) THEN
	    CALL RIVERIN( )
	    CALL READIC( )
            Ts = 0.
            Tc = 0.
            Tfish = 0.

!     C. Print some info to log file
      WRITE( MESG,* ) 'AFTER INITIAL CONDITION'
      CALL M3MSG2( MESG )
      CALL M3MSG2( '3D variables' )

      lwe = 0
      nwet = 0
      do k = 1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw = lwa,lwe
          i = iwet(lw)
          lump = lazc(lw)

!         2D fields
          Ts( i,k,1,HGSED2D ) = HGOUT2D( I,K,1,HGSED2D )
          HGOUT2D( I,K,1,HGSED2D ) = 0.
          Ts( i,k,1,MEHGSED2D ) = HGOUT2D( I,K,1,MEHGSED2D )
          HGOUT2D( I,K,1,MEHGSED2D ) = 0.
          Ts( i,k,1,HGMACe2D ) = HGOUT2D( I,K,1,HGMACe2D )
          HGOUT2D( I,K,1,HGMACe2D ) = 0.
          Ts( i,k,1,HGMACi2D ) = HGOUT2D( I,K,1,HGMACi2D )
          HGOUT2D( I,K,1,HGMACeiD ) = 0.
          Ts( i,k,1,MEHGMACe2D ) = HGOUT2D( I,K,1,MEHGMACe2D )
          HGOUT2D( I,K,1,MEHGMACe2D ) = 0.
          Ts( i,k,1,MEHGMACi2D ) = HGOUT2D( I,K,1,MEHGMACi2D )
          HGOUT2D( I,K,1,MEHGMACi2D ) = 0.
          Ts( i,k,1,LAYERS2D ) = lump
          HGOUT2D( I,K,1,LAYERS2D ) = 0.

!         Default initial conditions average from Leipe et al., 2013
          IF( lump .LE. 5 ) THEN
!             Ts( i,k,1,HGSED2D )   = 20. * 0.75
!             Ts( i,k,1,MEHGSED2D ) = 20. * 0.25
          ELSE
!             Ts( i,k,1,HGSED2D )   = 100. * 0.75
!             Ts( i,k,1,MEHGSED2D ) = 100. * 0.25
          END IF

!         3d fields
          do jj = 1,lump
             nwet = nwet+1
             itit = nwet

             DO NNC = 1,NCHEM
                 Tc( itit,NNC ) = HGOUT3D( I, K, JJ, NNC )
                 HGOUT3D( I, K, JJ, NNC ) = 0.
             END DO

!            Bugfix
!            IF( NNC .EQ. MEHGMAC3D .AND. jj .NE. lump ) THEN
!                Tc( itit,NNC ) = 0.
!            END IF

          enddo
        enddo
      enddo

      CALL M3MSG2( '' )
      CALL M3MSG2( '3D variables' )
      DO NNC = 1,NCHEM
          WRITE( MESG,* ) VNAMOUT( NNC ), SUM( TC( :,NNC ) ) / NDREI,
     &                    MINVAL( TC( :,NNC ) ), MAXVAL( TC(:,NNC ) )
          CALL M3MSG2( MESG )
      END DO

      CALL M3MSG2( '' )
      CALL M3MSG2( '2D variables' )
      DO NNC =  1, NHGOUT2D
        WRITE( MESG,* ) VNAMOUT(  NHGTRANS + NON_HG_TRANS + NON_TRANS + NNC ),
     &        SUM( TS( :,:,1,NNC ) ) / NDREI,
     &        MINVAL( TS( :,:,1,NNC ) ), MAXVAL( TS(:,:,1,NNC ) )
        CALL M3MSG2( MESG )
      END DO

	    FIRSTTIME = .FALSE.
	    CALL M3MSG2( '' )
	END IF

c-----------------------------------------------------------------------
c      month loop
c-----------------------------------------------------------------------

      WRITE( MESG,* ) 'STARTING TIME LOOP'
      CALL M3MSG2( MESG )

      istart=1
      do 2111 lmon=mona,mone
        monat = lmon

	print*,'mona lisa', lmon,mona,mone
c-----------------------------------------------------------------------
c    get file names from environment variables - added by Johannes Bieser 20.05.2013
c	phyfile and biofile need to be reset each month
c-----------------------------------------------------------------------
	CALL M3MSG2( '' )
	CALL M3MSG2( DASHLINE )
	WRITE( MESG,'(A45,I2,A2,I2)' )
     &		'Reading ECOSMO physics and biomass data for Y',MJAR,' M',LMON
	CALL M3MSG2( MESG )
	CALL M3MSG2( '' )
	CALL ENVSTR( 'BIOFILE', 'bio input file', '', biofile, IOS )
	IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for BIOFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL ENVSTR( 'FLUFILE', 'flux input file', '', flufile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for FLUFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

       CALL ENVSTR( 'PHYFILE', 'phyto input file', '', phyfile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for PHYFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

	IF( LMON .LT. 10 .AND. MJAR .LT. 10 ) THEN
	    WRITE( phyfile,'(A,4I1)' )
     &			 TRIM( phyfile ), 0, MJAR, 0, LMON
	    WRITE( biofile,'(A,4I1,A1)' )
     &			 TRIM( biofile ), 0, MJAR, 0, LMON, 'b'
            WRITE( flufile,'(A,4I1,A1)' )
     &                   TRIM( flufile ), 0, MJAR, 0, LMON, 'f'
	ELSE IF( LMON .LT. 10 .AND. MJAR .GE. 10 ) THEN
	    WRITE( phyfile,'(A,I2,2I1)' )
     &			 TRIM( phyfile ), MJAR, 0, LMON
	    WRITE( biofile,'(A,I2,2I1,A1)' )
     &			 TRIM( biofile ), MJAR, 0, LMON, 'b'
            WRITE( flufile,'(A,I2,2I1,A1)' )
     &                   TRIM( flufile ), MJAR, 0, LMON, 'f'
	ELSE IF( LMON .GE. 10 .AND. MJAR .LT. 10 ) THEN
	    WRITE( phyfile,'(A,2I1,I2)' )
     &			 TRIM( phyfile ), 0, MJAR, LMON
	    WRITE( biofile,'(A,2I1,I2,A1)' )
     &			 TRIM( biofile ), 0, MJAR, LMON, 'b'
            WRITE( flufile,'(A,2I1,I2,A1)' )
     &                   TRIM( flufile ), 0, MJAR, LMON, 'f'
	ELSE IF( LMON .GE. 10 .AND. MJAR .GE. 10 ) THEN
	    WRITE( phyfile,'(A,2I2)' )
     &			 TRIM( phyfile ), MJAR, LMON
	    WRITE( biofile,'(A,2I2,A1)' )
     &			 TRIM( biofile ), MJAR, LMON, 'b'
            WRITE( flufile,'(A,2I2,A1)' )
     &                   TRIM( flufile ), MJAR, LMON, 'f'
	ELSE
	    CALL M3MSG2( 'ERROR: Unexpected month xor year' )
	END IF

	CALL M3MSG2( phyfile )
	CALL M3MSG2( biofile )
        CALL M3MSG2( flufile )

	open(unit=nin,file=phyfile,form='unformatted')
	open(unit=nbio_in,file=biofile,form='unformatted')
        open(unit=nflu_in,file=flufile,form='unformatted')

!	CALL M3MSG2( '' )

c-----------------------------------------------------------------------
c      day loop (days per month)
c-----------------------------------------------------------------------
c      leap-year check
       iss = 0
       if((mod(iiiyear,4).eq.0).and.mod(iiiyear,100).ne.0)iss = 1
       if(mod(iiiyear,400).eq.0)iss = 1
       itagan(2) = 28+iss

      if(lmon.eq.mone)then
        itagend = itage
      else
        itagend = itagan(monat)
      endif

	AJDAY = CJDAY
	DO TMON = 1, MONAT-1
	    AJDAY = AJDAY - ITAGAN( TMON )
	END DO

!	THIS IS USED TO TEST LEAP DAY AND NEW YEAR TO GO SMOOTHLY
!	do nday=1,itagend
!     	    CALL NEXTIME( CDATE,CTIME,240000 )
!	enddo
!	itagend = 0

c----------------------------------------------------------------------- 


      do 1111 nday=1,itagend

        CALL CPU_TIME( T_START )
        print*, "STARTING DAY ", nday

!      do 1111 nday=AJDAY,itagend
      ltag = nday 


       CALL ENVSTR( 'OUTFILE', 'output file', '', outfile, IOS )
        IF( IOS .GT. 0 ) THEN
            MESG = 'Unknown Value for OUTFILE'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( MJAR .LT. 10 ) THEN
            WRITE( OUTFILE( 4:5 ),'(A1,I1)' ) '0', MJAR
        ELSE
            WRITE( OUTFILE( 4:5 ),'(I2)' ) MJAR
        END IF

        IF( LMON .LT. 10 ) THEN
            WRITE( OUTFILE( 6:7 ),'(A1,I1)' ) '0', LMON
        ELSE
            WRITE( OUTFILE( 6:7 ),'(I2)' ) LMON
        END IF

        IF( NDAY .LT. 10 ) THEN
            WRITE( OUTFILE( 8:9 ),'(A1,I1)' ) '0', NDAY
        ELSE
            WRITE( OUTFILE( 8:9 ),'(I2)' ) NDAY
        END IF
        WRITE( MESG,*) outfile, nday
        CALL M3MSG2( MESG )

c-----------------------------------------------------------------------
c    read daily model output from netcdf files, transports etc
c-----------------------------------------------------------------------
        IF( READ_DAILY ) THEN
        print*, 'READ_DAILY .EQ. .TRUE.'
          CALL M3MSG2( 'Read daily model output...' )
          iwrite=2 ! iwrite=1 writing, iwrite=2 reading
          print*,nin,mjar,lmon,nday,ivier,imal,iwrite
          call IO_daily_physics(nin,mjar,lmon,nday,ivier,imal,iwrite)! nnn write=1 or read=2
          read (nbio_in) Tcin
          read (nflu_in) Tflu
          istart=0

          uc=umit;vc=vmit;szahl=szmit;avc=acmit
          sac=scmit;tec=tcmit

          call uvrand (izet,iindex,lvrp,m,n,ndrei,szahl)
          large=0;gross=0
        END IF

C..... biofile
c      ibio=3 - flaggelates  [mgC/m**3]
c      ibio=4 - diatomes     [mgC/m**3]
c      ibio=5 - small zoo    [mgC/m**3]
c      ibio=6 - large zoo    [mgC/m**3]
c      ibio=7 - Detritus     [mgC/m**3]
c      ibio=8 - NH4          [mmolN/m**3]
c      ibio=9 - DOM          [mgC/m**3]
c      ibio=10- NO3          [mmolN/m**3]
c      ibio=11- PO4          [mmolP/m**3]
c      ibio=12- SiO2         [mmolSi/m**3]
c      ibio=13- O2           [ml/l]
c      ibio=14- Opal         [mmolSi/m**3]
c      ibio=15 (nxbio)- cyanopbacteria [mgC/m**3]
c      ibio=16 fish
c      ibio=nxbio+1:nbio   !possible other variables including carbon chemistry and Fish
c      ibio=nbio-2- sediment 1 detritus [mgC/m**3]
c      ibio=nbio-1- sediment 1 detritus phosphate [mgC/m**3]
c      ibio=nbio- sediment 2 silicate [mgC/m**3]
c      ibio=17 sediment total carbon and nitrogen          [mgC/m**2]
c      ibio=18 sediment phosphate       [mgC/m**2]
c      ibio=19 sediment silicate        [mgC/m**2]
c      ibio=20 MB                       [mgC/m**2]

C..... flufile
c      l=1,nflu: biological fluxes of matter (units for output)  
c      l=1 - flaggelates production           [mgC/m**3 per day]
c      l=2 - diatomes  production             [mgC/m**3 per day]
c      l=3 - small zoo grazing on flaggelates [mgC/m**3 per day]
c      l=4 - small zoo grazing on diatomes    [mgC/m**3 per day]
c      l=5 - small zoo grazing on Detritus    [mgC/m**3 per day]
c      l=6 - large zoo grazing on flaggelates [mgC/m**3 per day]
c      l=7 - large zoo grazing on diatomes    [mgC/m**3 per day]
c      l=8 - large zoo grazing on small zoo   [mgC/m**3 per day]
c      l=9 - large zoo grazing on Detritus    [mgC/m**3 per day]
c      l=10 -   [mgC/m**3 per day]
c      l=11 - fish grazing on small zooplankton  [mgC/m**3 per day]
c      l=12 - fish grazing on large zooplankton  [mgC/m**3 per day]
c      l=13 - fish grazing on Detritus           [mgC/m**3 per day]
c      l=14 - fish grazing on marco benthos      [mgC/m**3 per day]
c      l=15 - macro benthos grazing on sediment  [mgC/m**2 per day]
c      l=16 - macro benthos grazing on det+dom   [mgC/m**? per day]
c      l=17 - macro benthos grazing on total zooplankton   [mgC/m**2 per day]
c      l=18 - macro benthos grazing on total phytoplankton [mgC/m**2 per day]
c      l=19 - ??? UP_N       [mgC/m**3 per day]
c      l=20 - ??? UP_P       [mgC/m**3 per day]
c      l=21 - ??? UP_Si      [mgC/m**3 per day]
c      l=22 - ??? Blight     [mgC/m**3 per day]
c      l=23 - cyanobacteria production                   [mgC/m**3 per day]
c      l=24 - small zooplankton grazing on cyanobacteria [mgC/m**3 per day]
c      l=25 - large zooplankton grazing on cyanobacteria [mgC/m**3 per day]
c      l=26 - ??? sediment resuspension [mgC/m**3 per day]
c      l=27 - ??? sediment resuspension [mgC/m**3 per day]
c      l=28 - ??? denitrification       [mgC/m**3 per day]
c      l=29 - ??? burial                [mgC/m**3 per day]
c      l=30 - ??? taubot/72.
c      l=31 - mortality of small zooplankton                    [mgC/m**3 per day]
c      l=32 - mortality of large zooplankton                    [mgC/m**3 per day]
c      l=33 - mortality of flaggelates                          [mgC/m**3 per day]
c      l=34 - mortality of diatomes                             [mgC/m**3 per day]
c      l=35 - mortality of cyanobacteria                        [mgC/m**3 per day]
c      l=36 - mortality of fish                                 [mgC/m**3 per day]
c      l=37 - mortality + sedimentation of macro benthos        [mgC/m**2 per day]

c-----------------------------------------------------------------------
c    calculate additional values (density, area) - added by Johannes Bieser 07.06.2013
c-----------------------------------------------------------------------
!	PAC = surface pressure (KHOR NOT NDREI )
!	CTEST_CT = gsw_ct_from_t( SAC( 20 ),TEC( 20 ),PAC( 20 ) )
!	XRHOTEST = gsw_rho( SAC( 20 ), TEC( 20 ), PAC( 20 ) )
!	print*, 'RHO=',XRHOTEST

c-----------------------------------------------------------------------
c    read daily model input - added by Johannes Bieser 24.12.2012
c-----------------------------------------------------------------------
	CALL M3MSG2( '' )
	CALL CIN_START( )

!	C. calculate density using hourly PSRFC values from met fields
      Tun = 0.
      Tvn = 0.

      lwe = 0
      nwet=0

      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i = iwet(lw)
          lump = lazc(lw)

		PRES_DU = 0.0
		PRES_PS = ATMIN( I,K,CHOUR,MPRSFC )
!		Trho( ii ) = gsw_rho( SAC( ii ), TEC( ii ), ATMIN( I,K,CHOUR,MPRSFC ) 

          do jj=1,lump
             nwet = nwet+1
             ii=nwet

!		PRES_DO = RHO_DU
                PRES_DO = RHO_0
		PRES_DU = DZ( JJ )/2 * F_GRAVITY * RHO_0		! [Pa] [kg/(m * s**2)] based on default density
		RHO_CORR = gsw_rho( SAC( ii ), TEC( ii ), PRES_DU )	! [kg/m**3] corrected density
		PRES_DU = DZ( JJ )/2 * F_GRAVITY * RHO_CORR		! [Pa] corrected pressure

		PRES_PS = PRES_DO + PRES_DU + PRES_PS
		Tpres( ii ) = PRES_PS
		Trho( ii ) = gsw_rho( SAC( ii ), TEC( ii ), Tpres( ii ) )
                Tdept( ii ) = dz( jj )
                IF( jj .EQ. lump ) THEN
                    Tdept( ii ) = ldep( lw )          ! depth of lowest layer
                END IF
                IF( jj .EQ. 1 ) THEN
                  Tdept( ii ) = Tdept( ii ) + zac( lw )
                  IF( Tdept( ii ) .LT. 0.1 ) THEN
                    Tdept( ii ) = 0.1
                  END IF
                END IF

!	Save integrated velocities for lowest grid cell in 2D fields
!		IF( JJ .EQ. LUMP ) THEN
		    Tun( I,K,JJ ) = UC( ii )
		    Tvn( I,K,JJ ) = VC( ii )
!		END IF
      enddo
      enddo
      enddo

	DO I = 1,M-1
	    DO K = 1,N-1
		LUMP = Tdz( I,K )
		DO JJ = 1, LUMP

!		    IF( JJ .EQ. 1 ) THEN
!			Tun( I,K,JJ ) = Tun( I,K,JJ ) / DZ( JJ )
!			Tvn( I,K,JJ ) = Tvn( I,K,JJ ) / DZ( JJ )
		    IF( JJ .EQ. LUMP ) THEN

			IF( JJ .LT. Tdz( I,K+1 ) ) THEN
			    Tun( I,K,JJ ) = 2. * Tun( I,K,JJ )
     &				          / ( Thn( I,K ) + DZ( JJ ) )
			ELSE IF( JJ .EQ. Tdz( I,K+1 ) ) THEN
			    Tun( I,K,JJ ) = 2. * Tun( I,K,JJ )
     &				          / ( Thn( I,K ) + Thn( I,K+1 ) )
			ELSE
			    Tun( I,K,JJ ) = 0.
			END IF

			IF( JJ .LT. Tdz( I+1,K ) ) THEN
			    Tvn( I,K,JJ ) = 2. * Tvn( I,K,JJ )
     &					  / ( Thn( I,K ) + DZ( JJ ) )
			ELSE IF( JJ .EQ. Tdz( I+1,K ) ) THEN
			    Tvn( I,K,JJ ) = 2. * Tvn( I,K,JJ )
     &					  / ( Thn( I,K ) + Thn( I+1,K ) )
			ELSE
			    Tvn( I,K,JJ ) = 0.
!			    Tvn( I,K,JJ ) = Tvn( I,K,JJ ) / Thn( I,K )
			END IF
		    ELSE
			IF( JJ .EQ. Tdz( I,K+1 ) ) THEN
			    Tun( I,K,JJ ) =  2. * Tun( I,K,JJ )
     &					  / ( DZ( JJ ) + Thn( I,K+1 ) )
			ELSE IF ( JJ .LT. Tdz( I,K+1 ) ) THEN
			    Tun( I,K,JJ ) = Tun( I,K,JJ ) / DZ( JJ )
			ELSE 
			    Tun( I,K,JJ ) = 0.
			END IF

			IF( JJ .EQ. Tdz( I+1,K ) ) THEN
			    Tvn( I,K,JJ ) = 2. * Tvn( I,K,JJ )
     &					  / ( DZ( JJ ) + Thn( I+1,K ) )
			ELSE IF( JJ .LT. Tdz( I+1,K ) ) THEN
			    Tvn( I,K,JJ ) = Tvn( I,K,JJ ) / DZ( JJ )
			ELSE
			    Tvn( I,K,JJ ) = 0.
			END IF

		    END IF

		END DO
	    END DO
	END DO

	print*,'U = ', MINVAL( Tun ), MAXVAL( Tun )
	print*,'V = ', MINVAL( Tvn ), MAXVAL( Tvn )
c-----------------------------------------------------------------------

c-------------------------------------------------------------------
c    solve continuity equation, ensure divergence free transport
c------------------------------------------------------------------- 

          call konti(ht, dzbod,wtest,wsurf,ndrei)

          lwe=0
            do k=1,n
              lwa = lwe+1
              lwe = indend(k)
              do lw=lwa,lwe
                i = iwet(lw)
                zinc(i,k) =wsurf(lw)*dble(deltat)
                zac(lw)=zinc(i,k)
              enddo
           enddo
          dd(1)=dz(1)
          do kk=2,ilo
           dd(kk)=dz(kk)-dz(kk-1)
          enddo

           print*, istart
           istart=0
           
c----------------------------------------------------------------
c    ice parameter, convert
c---------------------------------------------------------------

      icemold=frice*his      
      frice=frimit;his=hismit;hisr=hisrmit
   
      include 'C_pol_ice_inc'

        print*, '    #DBG        Hg(T) Hg0 Hg* Hg2+(aq) Hg-DOM HgP | '
     &  //'DET PHYT ZOOP | MeHg(T)  MeHg(aq) MeHg-DOM  MeHg-DET MeHg(P)'

c----------------------------------------------------------------------- 
c      date loop (dates per day)
c-----------------------------------------------------------------------

      do 1010 ivier=1,iviert ! half-day loop 2 x 12 hours

      gross=0.;large=0

c----------------------------------------------------------------------
c      dt-loop (dt per half day)
c-----------------------------------------------------------------------

      do 1000 ip=1,iper ! 36 x 20min

c-----------------------------------------------------------------------
c    read hourly model output from netcdf files, transports etc
c-----------------------------------------------------------------------
!       ip == 1 full hour in 20min time step loop
        IF (READ_HOURLY .AND. MOD(ip-1,3) .EQ. 0) THEN
!!!          CALL M3MSG2( 'Read hourly model output...' )
          iwrite=2 ! iwrite=1 writing, iwrite=2 reading
!!!          print*,nin,mjar,lmon,nday,ivier,imal,iwrite
          call IO_daily_physics(nin,mjar,lmon,nday,ivier,imal,iwrite)! nnn write=1 or read=2
          read (nbio_in) Tcin
          read (nflu_in) Tflu

          istart=0

          uc=umit;vc=vmit;szahl=szmit;avc=acmit
          sac=scmit;tec=tcmit

          call uvrand (izet,iindex,lvrp,m,n,ndrei,szahl)
          large=0;gross=0.

C..... biofile
c      ibio=3 - flaggelates  [mgC/m**3]
c      ibio=4 - diatomes     [mgC/m**3]
c      ibio=5 - small zoo    [mgC/m**3]
c      ibio=6 - large zoo    [mgC/m**3]
c      ibio=7 - Detritus     [mgC/m**3]
c      ibio=8 - NH4          [mmolN/m**3]
c      ibio=9 - DOM          [mgC/m**3]
c      ibio=10- NO3          [mmolN/m**3]
c      ibio=11- PO4          [mmolP/m**3]
c      ibio=12- SiO2         [mmolSi/m**3]
c      ibio=13- O2           [ml/l]
c      ibio=14- Opal         [mmolSi/m**3]
c      ibio=15 (nxbio)- cyanopbacteria [mgC/m**3]i
c      ibio=nxbio+1:nbio   !possible other variables including carbon chemistry and Fish
c      ibio=nbio-2- sediment 1 detritus [mgC/m**3]
c      ibio=nbio-1- sediment 1 detritus phosphate [mgC/m**3]
c      ibio=nbio- sediment 2 silicate [mgC/m**3]

C..... flufile
c      l=1,nflu: biological fluxes of matter (units for output)
c      l=1 - flaggelates production           [mgC/m**3 per day]
c      l=2 - diatomes  production             [mgC/m**3 per day]
c      l=3 - small zoo grazing on flaggelates [mgC/m**3 per day]
c      l=4 - small zoo grazing on diatomes    [mgC/m**3 per day]
c      l=5 - small zoo grazing on Detritus    [mgC/m**3 per day]
c      l=6 - large zoo grazing on flaggelates [mgC/m**3 per day]
c      l=7 - large zoo grazing on diatomes    [mgC/m**3 per day]
c      l=8 - large zoo grazing on small zoo   [mgC/m**3 per day]
c      l=9 - large zoo grazing on Detritus    [mgC/m**3 per day]

!       C. calculate density using hourly PSRFC values from met fields
      Tun = 0.
      Tvn = 0.

      lwe = 0
      nwet=0

      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i = iwet(lw)
          lump = lazc(lw)

                PRES_DU = 0.0
                PRES_PS = ATMIN( I,K,CHOUR,MPRSFC )
!               Trho( ii ) = gsw_rho( SAC( ii ), TEC( ii ), ATMIN( I,K,CHOUR,MPRSFC )
          do jj=1,lump
             nwet = nwet+1
             ii=nwet

!                PRES_DO = RHO_DU
                PRES_DO = RHO_0
                PRES_DU = DZ( JJ )/2 * F_GRAVITY * RHO_0                ! [Pa] [kg/(m * s**2)] based on default density
                RHO_CORR = gsw_rho( SAC( ii ), TEC( ii ), PRES_DU )     ! [kg/m**3] corrected density
                PRES_DU = DZ( JJ )/2 * F_GRAVITY * RHO_CORR             ! [Pa] corrected pressure

                PRES_PS = PRES_DO + PRES_DU + PRES_PS
                Tpres( ii ) = PRES_PS
                Trho( ii ) = gsw_rho( SAC( ii ), TEC( ii ), Tpres( ii ) )
                Tdept( ii ) = dz( jj )
                IF( jj .EQ. lump ) THEN
                    Tdept( ii ) = ldep( lw )          ! depth of lowest layer
                END IF
                IF( jj .EQ. 1 ) THEN
                  Tdept( ii ) = Tdept( ii ) + zac( lw )
                  IF( Tdept( ii ) .LT. 0.1 ) THEN
                    Tdept( ii ) = 0.1
                  END IF
                END IF

!       Save integrated velocities for lowest grid cell in 2D fields
!               IF( JJ .EQ. LUMP ) THEN
                    Tun( I,K,JJ ) = UC( ii )
                    Tvn( I,K,JJ ) = VC( ii )
!               END IF
          enddo
        enddo
      enddo

!        DO I = 1,ndrei
!          print*, 'DEPTH: ', ii_j(i), Tdept(i)
!        END DO

        DO I = 1,M-1
            DO K = 1,N-1
                LUMP = Tdz( I,K )
                DO JJ = 1, LUMP

!                IF ( JJ .EQ. LUMP .AND. I .EQ. 82 ) THEN
!                        print*, Tun( I,K,JJ ),Tvn(I,K,JJ),Tdz(I,K),Thn(I,K),K
!                END IF

!                   IF( JJ .EQ. 1 ) THEN
!                       Tun( I,K,JJ ) = Tun( I,K,JJ ) / DZ( JJ )
!                       Tvn( I,K,JJ ) = Tvn( I,K,JJ ) / DZ( JJ )
                    IF( JJ .EQ. LUMP ) THEN

                        IF( JJ .LT. Tdz( I,K+1 ) ) THEN
                            Tun( I,K,JJ ) = 2. * Tun( I,K,JJ )
     &                                    / ( Thn( I,K ) + DZ( JJ ) )
                        ELSE IF( JJ .EQ. Tdz( I,K+1 ) ) THEN
                            Tun( I,K,JJ ) = 2. * Tun( I,K,JJ )
     &                                    / ( Thn( I,K ) + Thn( I,K+1 ) )
                        ELSE
                            Tun( I,K,JJ ) = 0.
                        END IF

                        IF( JJ .LT. Tdz( I+1,K ) ) THEN
                            Tvn( I,K,JJ ) = 2. * Tvn( I,K,JJ )
     &                                    / ( Thn( I,K ) + DZ( JJ ) )
                        ELSE IF( JJ .EQ. Tdz( I+1,K ) ) THEN
                            Tvn( I,K,JJ ) = 2. * Tvn( I,K,JJ )
     &                                    / ( Thn( I,K ) + Thn( I+1,K ) )
                        ELSE
                            Tvn( I,K,JJ ) = 0.
!                           Tvn( I,K,JJ ) = Tvn( I,K,JJ ) / Thn( I,K )
                        END IF
                    ELSE
                        IF( JJ .EQ. Tdz( I,K+1 ) ) THEN
                            Tun( I,K,JJ ) =  2. * Tun( I,K,JJ )
     &                                    / ( DZ( JJ ) + Thn( I,K+1 ) )
                        ELSE IF ( JJ .LT. Tdz( I,K+1 ) ) THEN
                            Tun( I,K,JJ ) = Tun( I,K,JJ ) / DZ( JJ )
                        ELSE
                            Tun( I,K,JJ ) = 0.
                        END IF

                        IF( JJ .EQ. Tdz( I+1,K ) ) THEN
                            Tvn( I,K,JJ ) = 2. * Tvn( I,K,JJ )
     &                                    / ( DZ( JJ ) + Thn( I+1,K ) )
                        ELSE IF( JJ .LT. Tdz( I+1,K ) ) THEN
                            Tvn( I,K,JJ ) = Tvn( I,K,JJ ) / DZ( JJ )
                        ELSE
                            Tvn( I,K,JJ ) = 0.
                        END IF

                    END IF

                END DO
            END DO
        END DO

!        print*,'U = ', MINVAL( Tun ), MAXVAL( Tun )
!        print*,'V = ', MINVAL( Tvn ), MAXVAL( Tvn )
!        print*,'H = ', MINVAL( Thn ), MAXVAL( Thn )

c-----------------------------------------------------------------------


c-------------------------------------------------------------------
c    solve continuity equation, ensure divergence free transport
c------------------------------------------------------------------- 
          
          call konti(ht, dzbod,wtest,wsurf,ndrei)

          lwe=0
            do k=1,n
              lwa = lwe+1
              lwe = indend(k)
              do lw=lwa,lwe
                i = iwet(lw)
                zinc(i,k) =wsurf(lw)*dble(deltat)
                zac(lw)=zinc(i,k)
              enddo
            enddo
          dd(1)=dz(1)
          do kk=2,ilo
            dd(kk)=dz(kk)-dz(kk-1)
          enddo

!          print*, istart
           istart=0
           
c----------------------------------------------------------------
c    ice parameter, convert
c---------------------------------------------------------------

      icemold=frice*his
      frice=frimit;his=hismit;hisr=hisrmit

!      include 'C_pol_ice_inc'  !already in outer loop


        END IF !End read hourly input
c-----------------------------------------------------------------------
c    Julian time and date counter for netcdf files I/O - added by Johannes Bieser 24.12.2012
c-----------------------------------------------------------------------
!.... Edited by Johannes Bieser
!........ increases CDATE and CTIM by 20minutes (internal timestep)
	CALL NEXTIME( CDATE,CTIME,TSTEP )
	CJDAY = CDATE - CDATE/1000*1000
!........ calculates hour a index for the ATMIN array in CIN
	CHOUR = CTIME / 10000 + 1
c-----------------------------------------------------------------------

      imal = imal+1
!     print*,'outfile=',outfile

      call settime(time,ijulu,ihouu,iminu,isecu)

!	WRITE( MESG,2829 ) ijahr,lmon,nday,ihouu,iminu,isecu,ijulu
!	CALL M3MSG2( MESG )
!      write(*,2829)ijahr,lmon,nday,ihouu,iminu,isecu,ijulu
! 2829 format('y m d h m s jul ',i5,5i3,i4)
       julianday=ijulu
!       print*, 'relative julianday' ,   ijulu

c-----------------------------------------------------------------------
c      check parameters (min,max,mean,sum) every timestep
c-----------------------------------------------------------------------

      if(loxstatinput)then !      if(.true.)then
        call xstat(vc,ndrei,vcmax,vcmin,vcmean,vcsum)
        call xstat(sac,ndrei,sacmax,sacmin,sacmean,sacsum)
        call xstat(tec,ndrei,tecmax,tecmin,tecmean,tecsum)
        call xstat(zac,khor,zacmax,zacmin,zacmean,zacsum)
        call xstat(wc,ndrei,ucmax,ucmin,ucmean,ucsum)
        call xstat(szahl,ndrei,vcmax,vcmin,vcmean,vcsum)
        call xstat(avc,ndrei,wcmax,wcmin,wcmean,wcsum)
!       print*,'vc min,max,mean,sum ',vcmin,vcmax,vcmean,vcsum
!       print*,'sac min,max,mean,sum ',sacmin,sacmax,sacmean,sacsum
!        print*,'tec min,max,mean,sum ',tecmin,tecmax,tecmean,tecsum
!        print*,'zac min,max,mean,sum ',zacmin,zacmax,zacmean,zacsum
!        print*,' wc min,max,mean,sum ',ucmin,ucmax,ucmean,ucsum
!        print*,' szahl min,max,mean,sum ',vcmin,vcmax,vcmean,vcsum
!        print*,' avc min,max,mean,sum ',wcmin,wcmax,wcmean,wcsum
       end if

c-----------------------------------------------------------------------
c      dummy arrays set to zero
c      call advection routine
c-----------------------------------------------------------------------
       large=0;gross=0.

        DO NNC = 1, nchem
          Tc( 82106,NNC ) = Tc( 82107,NNC )
          Tc( 82108,NNC ) = Tc( 82107,NNC )
        END DO

!       Apply sinking velocity at full hour time step for consistency with input data fields (MOD(ip-1,3) == 0.)
        CALL STROM( zinc,vtmit,vtmic,avmax,szahl,MOD(ip-1,3) )

        DO NNC = 1, nchem
          Tc( 82106,NNC ) = Tc( 82107,NNC )
          Tc( 82108,NNC ) = Tc( 82107,NNC )
        END DO

        DO ii = 1, ndrei
          DO nnc = 1, nchem
            IF( Tc( ii,nnc ) .LE. 0. ) THEN
                Tc( ii,nnc ) = 1E-6
            END IF
          END DO
        END DO
c-----------------------------------------------------------------------
c        set boundary values
c-----------------------------------------------------------------------

          call tsrneu3s(izet,iindex,lvrp,m,n,ndrei,ijahr,lmon,
     &              rsteps,julianday)

c-----------------------------------------------------------------------
c    read hourly radiation fields from MCIP input - edited by Johannes Bieser 15.04.2018
c-----------------------------------------------------------------------

        CALL LIGHT( ATMIN( :,:,CHOUR,MRGRND ), Tc, Tcin, XLIGHT, DD )  ! apply ice fraction FRICE
        DO ii = 1, ndrei
            IF( HISMIT( ii_i( ii ), ii_k( ii ) ) .GE. 1. ) THEN
                XLIGHT( ii ) = XLIGHT( ii )
     &                       * ( 1 - FRIMIT( ii_i( ii ),ii_k( ii ) ) )
            ELSE
                XLIGHT( ii ) = XLIGHT( ii )
     &                       * ( 1 - FRIMIT( ii_i( ii ),ii_k( ii ) ) )
     &                       * ( 1 - HISMIT( ii_i( ii ),ii_k( ii ) ) )
            END IF
        END DO

c-----------------------------------------------------------------------
!     C. Begin of mercury main loop
!     C. 1D loop for deposition, chemistry, partitioning, air-sea exchange, river input, and ice chemistry/exchange/transport

        DBG_CHG      = 0.
        DBG_FAW      = 0.
        DBG_SAX      = 0.
        DBG_DEP( 1 ) = 0.
        DBG_DEP( 2 ) = 0.
        DBG_GEM      = 0.
        DBG_SED      = 0.
        DBG_RES      = 0.
        DBG_PAM      = 0.
        DBG_RADI     = 0.
        DBG_FISH     = 0.
        DBG_NUM      = 0.

        DO CSTEP = 1, MSTEP
          do ii = 1, ndrei

           DO NNC = 1, NHGTRANS
             IF( Tc( ii,NNC ) .LT. 0. ) THEN
                print*, 'ERROR 000: ', VNAMOUT( NNC ), Tc( ii,NNC )
                Tc( ii,NNC ) = CMIN
             END IF
           END DO

!           C.1 Surface processes
            IF( ii_j( ii ) .EQ. 1 ) THEN
                ICE_THIC = HISMIT( ii_i( ii ), ii_k( ii ) )
                IF( ICE_THIC .GT. 0.1 ) THEN ! Ice thickness > 10cm
                    ICE_FRAC = FRIMIT( ii_i( ii ), ii_k( ii ) )
                ELSE
                    ICE_FRAC = 0.
                END IF
                HGOUT2D( ii_i(ii),ii_k(ii),1,FRICE2D ) = ICE_FRAC

!               C.1.1 Riverine inflow
                Dmin = MAX( 1.0,Tdept( ii ) )
                RIVER_LOAD = RIVER( ii_i( ii ), ii_k( ii ), LMON )
                RIVER_LOAD = RIVER_LOAD * 1000.

                IF ( RIVER_LOAD .LT. 0.0 ) THEN
                  print*, 'ERROR 404: River load ', RIVER_LOAD
                END IF
                IF ( Dmin .LT. 0.0 ) THEN
                  print*, 'ERROR 405: River depth ', Dmin
                  d0 = 1.0
                END IF
                IF( RIVER_LOAD .GT. 100000. ) THEN
                  print*, 'ERROR 097: HIGH RIVER INFLOW ',
     &                     RIVER_LOAD, ii_i( ii ), ii_k( ii )
                END IF

                IF( RIVER_LOAD .GT. 0 ) THEN

                  Tc( ii,HGDEM3D ) =
     &            Tc( ii,HGDEM3D ) + RIVER_LOAD * 0.10                   ! g/month
     &                             / ITAGAN( LMON ) / 24. / 3600.        ! g/s
     &                             * 1E09 * ITS                          ! ng
     &                             / Tarea( ii ) / Dmin / 1000.          ! ng/l

                  Tc( ii,HGPOC3D ) =
     &            Tc( ii,HGPOC3D ) + RIVER_LOAD * 0.80                   ! g/month
     &                             / ITAGAN( LMON ) / 24. / 3600.        ! g/s
     &                             * 1E09 * ITS                          ! ng
     &                             / Tarea( ii ) / Dmin / 1000.          ! ng/l

                  Tc( ii,MEHGPOC3D ) =
     &            Tc( ii,MEHGPOC3D ) + RIVER_LOAD * 0.10                 ! g/month
     &                               / ITAGAN( LMON ) / 24. / 3600.      ! g/s
     &                               * 1E09 * ITS                        ! ng
     &                               / Tarea( ii ) / Dmin / 1000.        ! ng/l

                  Tc( ii,PTOM3D ) = 
     &            Tc( ii,PTOM3D ) + RIVER_LOAD * 1000.                   ! kg/month
     &                            / ITAGAN( LMON ) / 24. / 3600.         ! kg/s
     &                            * 1E09 * ITS                           ! ng
     &                            / Tarea( ii ) / Tdept( ii ) / 1000.    ! ng/l

                  Tc( ii,DTOM3D ) = 
     &            Tc( ii,DTOM3D ) + RIVER_LOAD * 1000.            ! kg/month
     &                           / ITAGAN( LMON ) / 24. / 3600.          ! kg/s
     &                           * 1E09 * ITS                            ! ng
     &                           / Tarea( ii ) / Tdept( ii ) / 1000.     ! ng/l

                END IF

                IF( Tc( ii,HGDEM3D ) .LT. 0.0 ) THEN
                  print*, 'ERROR 222 Rdem: ',  Tc( ii,HGDEM3D ),
     &                    ITAGAN( LMON ), MSTEP, Tarea(ii), Dmin,
     &                    RIVER_LOAD, ITS
                  Tc( ii,HGDEM3D ) = CMIN
                END IF

!               C.1.2 Atmospheric deposition
!                DEP1 = ( ATMIN( ii_i( ii ), ii_k( ii ), CHOUR, MDHG0 ) +      ! [ng/L]
!     &                   ATMIN( ii_i( ii ), ii_k( ii ), CHOUR, MWHG0 ) )
!     &               / ( Dmin * 3600. * 1000. ) * ITS ! [ng/m**2*h]/[m]*[m**3/L] = [ng/L] per time step ITS [s]
                DEP1 = 0. !Use Hg0 exchange from MECOSMO

                DEP2 = ( ATMIN( ii_i( ii ), ii_k( ii ), CHOUR, MDHG2 ) +      ! [ng/L]
     &                   ATMIN( ii_i( ii ), ii_k( ii ), CHOUR, MWHG2 ) +
     &                   ATMIN( ii_i( ii ), ii_k( ii ), CHOUR, MDHGP ) +
     &                   ATMIN( ii_i( ii ), ii_k( ii ), CHOUR, MWHGP ) )
     &               / ( Dmin * 3600. * 1000. ) * ITS ! [ng/m**2*h]/[m]*[m**3/L] = [ng/L] per time step ITS [s]

!               Include ice coverage 
                IF( ICE_THIC .GT. 0.1 ) THEN
                    DEP1 = DEP1 * ( 1 - ICE_FRAC )
                    DEP2 = DEP2 * ( 1 - ICE_FRAC )
                END IF

                !Write to output and apply to tracer array
                Tc( ii,HGDEM3D ) = Tc( ii,HGDEM3D ) + DEP1 ! [ng/L]
                Tc( ii,HGDIS3D ) = Tc( ii,HGDIS3D ) + DEP2 ! [ng/L]

                HGOUT2D( ii_i(ii),ii_k(ii),1,HGDEP2D ) =
     &          HGOUT2D( ii_i(ii),ii_k(ii),1,HGDEP2D ) + ( DEP1 + DEP2 )
     &                   * Tarea(ii) * 1.E-12 * Dmin * 1000.        ! [kg]

!	        Not implemented
!               C.1.3 Ice chemistry
!                CALL ICECHEM( Tc, ii, ICE_FRAC, ICE_THIC, UMIT, VMIT,
!     &                  DX(ii_i(ii),ii_k(ii)), DY(ii_i(ii),ii_k(ii)),
!     &                        ii_i( ii ), ii_k( ii ), sac, tec, 0, 0 )


!               C.0.5 Surface statistics
                DO NNC = 1, IHGT
                    DBG_CHG( NNC ) = DBG_CHG( NNC ) + Tc( ii,NNC )
     &                     + Tc( ii,HGMACe3D )   + TC( ii,HGMACi3D )
     &                     + Tc( ii,MEHGMACe3D ) + Tc( ii,MEHGMACi3D )
                END DO

                DBG_SAX      = DBG_SAX + Faw / ( Tdept( ii ) * 1000 )
                DBG_FAW      = DBG_FAW + Faw
                DBG_DEP( 1 ) = DBG_DEP( 1 ) + DEP1
                DBG_DEP( 2 ) = DBG_DEP( 2 ) + DEP2
                DBG_DEM      = DBG_DEM + Tc( ii,HGDEM3D )
                DBG_GEM      = DBG_GEM + ATMIN( ii_i(ii),ii_k(ii),CHOUR,MCHG0 )
                DBG_RADI     = DBG_RADI + XLIGHT( ii )
                DBG_PAM      = DBG_PAM + Tc( ii,PTOM3D )
                DBG_FISH     = DBG_FISH + Tfish( ii,2 )
                DBG_NUM      = DBG_NUM + 1.
!                print*, ii, ii_i(ii), ii_k(ii), ii_j(ii), DBG_NUM

            END IF

!           C.3.1 Particle partitioning
            CALL PARTPART( Tc, Tcin, VERBOSE, ii )

!           C.4 Bioaccumulation
            CALL BIOACC_HG( Tc, Tfish, Tcin, Tflu, tec, sac, VERBOSE,
     &                      ii, ii_i, ii_k, ii_j, ii_j_max, ii_ldep )

            CALL BIOACC_MEHG( Tc, Tfish, Tcin, Tflu, tec, sac, VERBOSE, 
     &                        ii, ii_i, ii_k, ii_j, ii_j_max, ii_ldep )

!           C.3.2 Particle partitioning
            CALL PARTPART( Tc, Tcin, VERBOSE, ii )

!           C.5 Chemistry
            CALL CHEM( Tc, Tcin, tec, sac, xlight, VERBOSE, ii )

!           C.1 Surface processes
            IF( ii_j( ii ) .EQ. 1 ) THEN

!               C.1.4 Air-sea exchange

                DO NNC = 1, 10  !smaller timestep for asx
                    CALL AIRSEA( Faw, Tc, ii, ii_i(ii), ii_k(ii), CHOUR,
     &                 sac, tec, Tdept, Tarea, ICE_THIC, ICE_FRAC,VERBOSE )
                END DO
           
!           C.2 Benthic processes
            ELSE IF( ii_j( ii ) .EQ. ii_j_max( ii ) ) THEN
                CALL SEDIMENT( ii, ii_i, ii_k, ii_j, Tc, Tcin, VERBOSE )
            END IF
          end do ! loop over all wet grid cells
        END DO ! Internal time step loop

!       !HGTOT
        DBG_CHG( IHGT+1 ) = 0.
        DO NNC = 1, IBIOHG-1 ! -1 excludes MEHGDET3D
            DBG_CHG( IHGT+1 ) = DBG_CHG( IHGT+1 ) + DBG_CHG( NNC )
        END DO

!       !MEHGTOT
        DBG_CHG( IHGT+2 ) = 0.
        DO NNC = IMEHG, IBIOHG-1 ! -1 excludes MEHGDET3D
            DBG_CHG( IHGT+2 ) = DBG_CHG( IHGT+2 ) + DBG_CHG( NNC )
        END DO
        DBG_CHG( IHGT+2 ) = DBG_CHG( IHGT+2 ) / DBG_CHG( IHGT+1 ) * 100

!       !MEHGBIO
        DBG_CHG( IHGT+3 ) = 0.
        DO NNC = IBIOHG, IHGT
            DBG_CHG( IHGT+3 ) = DBG_CHG( IHGT+3 ) + DBG_CHG( NNC )
        END DO
        DBG_CHG( IHGT+3 ) = DBG_CHG( IHGT+3 ) / DBG_CHG( IHGT+1 ) * 100

        DO NNC = 1, IHGT
          DBG_CHG( NNC ) = DBG_CHG( NNC ) / DBG_CHG( IHGT+1 ) * 100
        END DO

        DBG_FISH = DBG_FISH / DBG_CHG( IHGT+1 ) * 100

!       Normalization to ng/L
        DBG_CHG( IHGT+1 ) = DBG_CHG( IHGT+1 ) / DBG_NUM

        DBG_DEM = DBG_DEM / DBG_NUM * 1000.
        DBG_GEM = DBG_GEM / DBG_NUM 
        DBG_PAM = DBG_PAM / DBG_NUM
        DBG_RADI = DBG_RADI / DBG_NUM
        DBG_SAX = DBG_SAX / DBG_NUM / DBG_DEM * 1000. * 100.
        DBG_NUM = DBG_NUM / MSTEP

        WRITE( MESG,* )
     &         'DBG ', DBG_CHG( IHGT+1 ),
     &         DBG_CHG( HGDEM3D ), DBG_CHG( HGSTAR3D ), DBG_CHG( HGDIS3D ), DBG_CHG( HGDOC3D ), DBG_CHG( HGPOC3D ), ' | ',
     &         DBG_CHG( IHGT+2 ), DBG_CHG( MEHGDIS3D ), DBG_CHG( MEHGDOC3D ), DBG_CHG( MEHGPOC3D ), DBG_CHG( DMHG3D ), ' | ',
     &         DBG_CHG( IHGT+3 ), DBG_FISH, ' | ',
     &         DBG_SAX, DBG_FAW, DBG_DEM, DBG_GEM, ' | ',
     &         DBG_DEP( 1 ), DBG_DEP( 2 ), DBG_SED, DBG_RES, DBG_PAM,
     &         DBG_NUM, DBG_RADI, CDATE, CTIME
        CALL M3MSG2( MESG )

        DBG_CHG = 0.
        DBG_DEM = 0.
        DBG_FAW = 0.
        DBG_SAX = 0.
        DBG_GEM = 0.
        DBG_DEP = 0.
        DBG_SED = 0.
        DBG_RES = 0.
        DBG_PAM = 0.
        DBG_NUM = 0.
        DBG_RADI = 0.

99	FORMAT( A20,4E14.4E3 )
!	C. Write hourly values to output field and daily means to output file
C...... This adds data for every hour - data is then divided by 24. in cout.f WRITE_CHEM
        IF ( MODULO( CTIME ,10000 ) .EQ. 0 ) THEN
            DO SN = 1, IHGT
                CALL ADD1DTO3D( HGOUT3D( :,:,:,SN ), Tc( :,SN ), NDREI )

                IF( SN .LT. IBIOHG ) THEN ! .LT. excludes MEHGDET3D
!                print*, "HGTOT: ", VNAMOUT( SN ), SN
           CALL ADD1DTO3D( HGOUT3D( :,:,:,HGTOT3D ), Tc( :,SN ), NDREI )

                  IF( SN .GE. IMEHG ) THEN
!                  print*, "MEHGTOT: ", VNAMOUT( SN ), SN
         CALL ADD1DTO3D( HGOUT3D( :,:,:,MEHGTOT3D ), Tc( :,SN ), NDREI )
                  END IF

                ELSE
!                print*, "MEHGBIO: ", VNAMOUT( SN ), SN
         CALL ADD1DTO3D( HGOUT3D( :,:,:,MEHGBIO3D ), Tc( :,SN ), NDREI )
                END IF

            END DO

            DO SN = IHGT+1, NCHEM
!              print*, "NON-HG: ", VNAMOUT( SN ), SN
                CALL ADD1DTO3D( HGOUT3D( :,:,:,SN ), Tc( :,SN ), NDREI )
            END DO

      lwe = 0
      nwet=0
      do k=1,n
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i = iwet(lw)
          lump = lazc(lw)
          do jj=1,lump
             nwet = nwet+1
             itit=nwet

          enddo
        enddo
      enddo
      CALL ADD1DTO3D( HGOUT3D( :,:,:,FLA3D ), TCin( :,FLA ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,DIA3D ), TCin( :,DIA ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,CYA3D ), TCin( :,CYA ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,ZOS3D ), TCin( :,ZOS ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,ZOL3D ), TCin( :,ZOL ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,DET3D ), TCin( :,DET ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,DOC3D ), TCin( :,DOM ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,OXY3D ), TCin( :,O2  ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,RAD3D ), xlight( : ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,SAL3D ), sac( : ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,TEM3D ), tec( : ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,FSH3D ), TCin( :,FSH ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,POC3D ), TCin( :,DET ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,BPC3D ), TCin( :,FLA ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,BPC3D ), TCin( :,DIA ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,BPC3D ), TCin( :,CYA ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,BZC3D ), TCin( :,ZOS ), NDREI )
      CALL ADD1DTO3D( HGOUT3D( :,:,:,BZC3D ), TCin( :,ZOL ), NDREI )

!       2D Vars are now not devided by count and thus written each
!       timestep in the respective routines
      do  k=1,n
        do i=1,m
          DO NNC = 1, NHGOUT2D
            HGOUT2D( i,k,1,NNC ) = HGOUT2D( i,k,1,NNC ) + Ts( i,k,1,NNC )
          END DO
        end do
      end do

      CALL WRITE_CHEM( .FALSE. )	! .FALSE. = write data / .TRUE. = create new file 

        END IF

c-----------------------------------------------------------------------
c        set boundary values
c-----------------------------------------------------------------------

          call tsrneu3s(izet,iindex,lvrp,m,n,ndrei,ijahr,lmon,
     &              rsteps,julianday)

c-----------------------------------------------------------------------
c      check end of the day & compute daily mean
c-----------------------------------------------------------------------

      igstep = igstep+1
      if(igstep.eq.nsteps)then
!       print*,'End daily mean',igstep,nsteps,rsteps,iper
        igstep = 0

      end if
      time = time+dt

 1000 continue

c----------------------------------------------------------------------- 
c      end of dt loop
c-----------------------------------------------------------------------

      write(7,606) ivier,mjar,lmon,nday
 606  format (2x,'termin',i2,3x,'j,m,d, /',3i4,3x) 


 1010 continue  
c-----------------------------------------------------------------------
c      end of half-day loop
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c          write daily output in netcdf format

c-----------------------------------------------------------------------

!        if(nc_out) then
!           print*,'outfile',outfile,'vor ncwriter'
!!          write(infile(6:7),'(i2.2)') monat
!CCC          write(infile(8:9),'(i2.2)') nday
!             call ncwriter(pathnc_out, outfile)
!            print*,'ncwriter erledigt'
!         end if
!
      Tc_mit=0. ! set daily mean to zero

        CALL CPU_TIME( T_FINISH )
        print*, 'CPU_TIME', T_FINISH - T_START, nday

 1111 continue

c-----------------------------------------------------------------------
c      end of day loop
c-----------------------------------------------------------------------


        print*, 'End of daily loop ', mjar, monat, itaga, ltag
      write(7,6999) mjar,monat,itaga,ltag
 6999 format (/3x,'jahr, monat /',2i5,'  von tag',i4, 
     * '  bis tag',i4,'  simuliert') 

 2111 continue 

c-----------------------------------------------------------------------
c      end of month loop
c-----------------------------------------------------------------------
c

c-----------------------------------------------------------------------
c      system dependent stuff to change backup filename
c-----------------------------------------------------------------------

c      write(setupcep(4:5),'(i2.2)') mjar
c      write(setupcep(6:7),'(i2.2)') lmon-1
c      write(setupcep(8:9),'(i2.2)') nday-1

!      call system('cp '//pathout//fbackup//' '//pathset//setupcep)


 5099    continue !end of year loop  
        goto 1002
        goto 1001
c-----------------------------------------------------------------------
 1002      stop'normal program **** s t o p ****'
c-----------------------------------------------------------------------
      end 

c-----------------------------------------------------------------------
      block data 
c----------------------------------------------------------------------- 
      real i0,i1,i2,j0,m0,m1,m2,k0,k1,k2 
      common /coef/b0,b1,b2,b3,b4,c0,c1,c2,d0,a0,a1,a2,a3,a4,a5, 
     * f0,f1,f2,f3,g0,g1,g2,i0,i1,i2,j0,m0,m1,m2,e0,e1,e2,e3,e4, 
     * h0,h1,h2,h3,k0,k1,k2 

c      coefficients  
      data b0/+8.24493e-1/, 
     * b1/-4.0899 e-3 /, 
     * b2/ 7.6438 e-5 /, 
     * b3/-8.2467 e-7 /, 
     * b4/ 5.3875 e-9 /, 
     * c0/-5.72466e-3 /, 
     * c1/ 1.0227 e-4 /, 
     * c2/-1.6546 e-6 /, 
     * d0/ 4.8314 e-4 /, 
     * a0/999.842594 /, 
     * a1/ 6.793952e-2 /, 
     * a2/-9.095290e-3 /, 
     * a3/ 1.001685e-4 /, 
     * a4/-1.120083e-6 /, 
     * a5/ 6.536332e-9 / 

      data f0/54.6746 /, 
     * f1/-0.603459 /, 
     * f2/ 1.09987 e-2 /, 
     * f3/-6.1670  e-5 /,  
     * g0/ 7.944  e-2 /, 
     * g1/ 1.6483 e-2 /, 
     * g2/-5.3009 e-4 /, 
     * i0/ 2.2838 e-3 /, 
     * i1/-1.0981 e-5 /, 
     * i2/-1.6078 e-6 /, 
     * j0/ 1.91075 e-4 /, 
     * m0/-9.9348 e-7 /, 
     * m1/ 2.0816 e-8 /, 
     * m2/ 9.1697 e-10 / 

      data e0/19652.21  /, 
     * e1/ 148.4206 /, 
     * e2/-2.327105 /, 
     * e3/ 1.360477e-2 /, 
     * e4/-5.155288e-5 /, 
     * h0/ 3.239908 /, 
     * h1/ 1.43713 e-3 /, 
     * h2/ 1.16092 e-4 /, 
     * h3/-5.77905 e-7 /, 
     * k0/ 8.50935 e-5 /, 
     * k1/-6.12293 e-6 /, 
     * k2/ 5.2787  e-8 /  


      end


c-----------------------------------------------------------------------
      subroutine deco1d2d (z,zc,iwet,indend,m,n,fac,xdef) 
c-----------------------------------------------------------------------
c      creates 2-d fields from 1-d fields 
c-----------------------------------------------------------------------
      
      include 'C_khor_inc'

      dimension iwet(khor1),indend(n) 
      dimension z(m,n),zc(khor) 

      do k=1,n 
      do i=1,m 
         z(i,k) = xdef 
      enddo
      enddo

      lwe = 0 
      do k=1,n 
        lwa = lwe+1 
        lwe = indend(k) 
        do lw=lwa,lwe 
          i  = iwet(lw) 
          zz = zc(lw) 
          z(i,k) = zz*fac 
        enddo
      enddo


      return 
      end 

c-----------------------------------------------------------------------
      subroutine deco1d2di (iz,izc,iwet,indend,m,n,ifac,ixdef)
c-----------------------------------------------------------------------
c      creates 2-d fields from 1-d fields
c-----------------------------------------------------------------------
    
      include 'C_khor_inc'

      dimension iwet(khor1),indend(n)
      dimension iz(m,n),izc(khor)  

      do k=1,n 
      do i=1,m 
         iz(i,k) = ixdef 
      enddo
      enddo

      lwe = 0
      do k=1,n 
        lwa = lwe+1
        lwe = indend(k)
        do lw=lwa,lwe
          i  = iwet(lw)
          izz = izc(lw)
          iz(i,k) = izz*ifac
        enddo
      enddo


      return
      end

c-----------------------------------------------------------------------
      function bomgrd(x) 
c-----------------------------------------------------------------------
c      wandelt bogenmass in altgrad 
c----------------------------------------------------------------------- 

      pi = atan(1.) * 4. 
      r180pi = 180. / pi 
      bomgrd = x * r180pi 
      return 
      end 


c-----------------------------------------------------------------------
      subroutine comp3d1d (d3,d1,mwet) 
c-----------------------------------------------------------------------
c      compress 3-d to 1d 
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter (khor1=khor+1) 

      dimension d3(m,n,ilo),d1(mwet) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *  indend(n),isornr(n),isorsr(n),islab(n)

      nwet = 0 
      lwe  = 0 
      do k=1,n 
        lwa = lwe+1 
        lwe = indend(k) 
        do lw=lwa,lwe 
          i = iwet(lw) 
          ldown = lazc(lw) 
          do j=1,ldown 
             d1(nwet+j) = d3(i,k,j) 
          enddo
          nwet = nwet+ldown 
       enddo
      enddo 

      return 
      end


c-----------------------------------------------------------------------
      subroutine comp2d1d (zic,z) 
c-----------------------------------------------------------------------
c      compress 2-d to 1d 
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 

      dimension zic(khor),z(m,n) 

      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *  indend(n),isornr(n),isorsr(n),islab(n) 

      lwe = 0 
      nwet = 0 
      do k=1,n 
        lwa = lwe+1 
        lwe = indend(k) 
        do lw=lwa,lwe 
          i = iwet(lw) 
          zic(lw) = z(i,k) 
        enddo
      enddo


      return 
      end 

c-----------------------------------------------------------------------
      subroutine gittop(ltief,dzbod)
c-----------------------------------------------------------------------

      include 'C_model_inc'

      parameter(khor1=khor+1)

      real dzbod(m,n)
      integer ltief(m,n)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *  indend(n),isornr(n),isorsr(n),islab(n)

      do i=1,m
      do j=1,n
        dzbod(i,j) = 0.
        ltief(i,j) = 0
      enddo
      enddo

      lwe = 0
      do j=1,n
        lwa = lwe+1
        lwe = indend(j)
        do lw=lwa,lwe
            i = iwet(lw)
            ltief(i,j) = lazc(lw)
            dzbod(i,j) = real(ldep(lw))
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      function grdbom(x) 
c-----------------------------------------------------------------------
c      altgrad in radian 
c-----------------------------------------------------------------------

      pi = atan(1.) * 4. 
      pi180 = pi / 180. 
      grdbom = x * pi180 

      return 
      end 

c-----------------------------------------------------------------------
      subroutine stromTVD(zinc,vtmit,vtmic,avmax,ftac,szahl,
     & nday,Tc_flux_all, EPSIL, ISPEC) 
c-----------------------------------------------------------------------
c      transport of tracers, advection & turbuleny diffusion
c      TVD total variation diminishing scheme with superbee
c      Advection scheme (Barthel et al., 2012)
c      
c-----------------------------------------------------------------------

      USE COUT

      include 'C_model_inc'
      parameter(ilop1=ilo+1)  
      parameter(khor1=khor+1)  

      common/T_4d/ Tc(ndrei,nchem),Tc_mit(ndrei,nchem) 
      common/Tin_4d/Tcin(ndrei,3:ninbio),Tflu(ndrei,nflu),Tfis(ndrei,2)  
      dimension Tc_flux_all(khor,nbio)  !river and P-E ![m**s/sec*timestep]

c parameters used for the TVD procedure:
       integer ISPEC ! index for chemistry array for opemMP
      integer NSIG
       double precision COURUL, COURUR, COURVO, COURVU, COURWO, COURWU
       double precision D1, D2, DDI, DDIZ, DENUM
       double precision EPSIL
       double precision PHIR, PHIR1
       double precision RATIO, RSUL, RSUR, RSVO, RSVU, RSWO, RSWU
      double precision SADV

      double precision dwdz, divcour
       double precision zzz
c#######
      double precision  s,sn,sd,sbet,tsa,tsb,tsc,tsal 
      double precision tsbb

      double precision  XYadv,svert
      double precision D_o2,D_u2,SinkD_o(nchem),SinkD_u(nchem)

c       double precision uot,uwt,vnt,vst
c       double precision  wo, wu
      real*8 zinc
      dimension XYadv(ndrei),svert(ndrei)
!      dimension Tl (m,ilo,nbio,3),Tc_k(ilo,nbio),Tl_k(ilo,nbio)
!      dimension Tl_bio(nbio)
      dimension s(ndrei),sn(ndrei),sd(ndrei)
      dimension sbet(ndrei),ftac(khor),szahl(ndrei) 
      dimension tsa(ndrei),tsb(ndrei),tsc(ndrei) 
      dimension tsal(ndrei)
      dimension zinc(m,n)
      logical obwas,advekt,rfals 
c 
       dimension e(ilo),taa(ilo),tac(ilo),tab(ilo),sad(ilo),tad(ilo)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 
c 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 

c 
      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m), 
     * sincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m) 
c 
 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,dtrdln(m),dtdlr

c 
c      ------------ laenge von common /aux/ = m*(ilo*3*10+12) 
      common /aux/ u(ilo,m,2),v(ilo,m,2),sa(ilo,m,3),te(ilo,m,3) 
     * ,za(m,2),su(ilo,m,2),sv(ilo,m),wd(ilo,m,2)
      dimension avd(ndrei) 
c      ---------------------- laenge von common /auxint/ = m*(ilo*6+3) 
      common /auxint/ lay(ilo,m,2),laz(m,2) 
c      new array
      integer jc5(m,n,ilo)
c       common /rivflow/ rivflux(m,n)

      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo), 
     * prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),d_(ilo), 
     * qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
      logical init
      common /vec44/ init,ljumm(khor),vol0,dtdx(m),dvo(m),dvu(m),
     1 lww(-3:3,-3:3,ndrei),llay(0:khor,ilo)

      common /vecind/ indwet(khor),lb(n),le(n),indver(ndrei),
     * jwet(khor),llw(ndrei),Dvoltj(m,ilo,n),
     * indi(ndrei),irbsor(kasor,2,2),nrbsor(2,2)


      common /julian/ julianday

      data init /.true./
      dimension d(ilo,m),Tc0(0:ndrei),dd(0:ndrei),r_dt_dd(ndrei)
      real*8 dd_dd(ndrei), dt_lr
c--------  set model variables to be advected -----------

!      ibio0=1  !
!      print*,'ibio0 nbio',ibio0,nbio
c       call Rzero(Bmass,nbio+1) !integrated biomass, output each time step

      lob = 1

      mz = m-1 
      nz = n-1 
ccc  epsilon minimum level thicknes        
      epsilon=0.5
c--------------------------------
      nwet = 0 
      lwe = 0
      do j = 1,n 
         lwa = lwe+1 
         lwe = indend(j) 
         lb(j) = lwa
         le(j) = lwe
         do lw = lwa,lwe 
             indwet(lw) = nwet
             jwet(lw) = j
             do k = 1, lazc(lw)
                !indi(nwet+k) = iwet(lw) 
                indver(nwet+k) = k
                llw(nwet+k)=lw
             end do
             nwet = nwet+lazc(lw) 
         end do
      end do

c----------------------------
      if (init) then
         do lw = 1,khor
             ljumm(lw) = 0
         end do
         do j = 2,nz 
             do lw = lb(j),le(j)
                i = iwet(lw) 
                ljum = min0(i-2,1) * min0(j-2,1) 
                ljum = max0(0,ljum) 
                ljumm(lw) = ljum
             end do
         end do

         dd(0)=0.
         vol0=0.
           nwet=0
         do j = 1,n 
             do lw = lb(j),le(j) 
                i = iwet(lw) 
                ldown = lazc(lw) 
                do k = 1,ldown 
             nwet=nwet+1
                    if (k.ne.1) then
                        if(k.eq.ldown) then
c                            zzz=float(ldep(lw))
                        zzz=dfloat(ldep(lw))
                        else
                           zzz=dble(pd(k)) ! +zac(lw))
cc                           if(k.eq.ldown)then
c                           zzz=dfloat(ldep(lw))+zac(lw)
c                           end if
                        end if
c                         Dvoltj(i,k,j)=zzz*dl*dln(i) ![m**3]
                        Dvoltj(i,k,j)=zzz*dprod(dl,dln(i))
                        vol0=vol0+Dvoltj(i,k,j)            
                        dd(nwet)=zzz
                     dd_dd(nwet)=dprod(dd(nwet),dd(nwet))
                    endif
                end do
             end do
         end do



         do i = 1, m
             dtdx(i) = dt/  dln(i)    !dln(j)  = dlam*cosfiu 
             dvu(i)=dln(i)*dlvu(i)
             dvo(i)=dln(i)*dlvo(i)
         end do
!	  print*, dt, dlr 
             dt_lr=dt*dlr
!	  print*, dt_lr, dt, dlr 
         do i=1,m
          do j=1,n
            do k=1,ilo
         jc5(i,j,k)=0
            enddo
          enddo
         enddo

      nwet = 0 
      lwe = 0

      do j = 1,n 
         lwa = lwe+1 
         lwe = indend(j) 
         do lw = lwa,lwe 
          i = iwet(lw)
          j5= jwet(lw)
             do k = 1, lazc(lw)
             jc5(i,j5,k)=nwet+k
             end do
             nwet = nwet+lazc(lw) 
         end do
      end do

cMM: changed from -2,2
       do ii = -3,3
       do jj = -3,3
            do nwet = 1,ndrei
             lw=llw(nwet)
c        since there are 4 additional land lines at west and north, no limits
c        for i=iwet(..)-3 and j=jwet(..)-3 :
        i=max(iwet(lw)+ii,1)       !!!CS NEW 
        j = max(jwet(lw)+jj,1)     !!!CS NEW
c
        i = min(m,i)
        j = min(n,j)

        k=indver(nwet)
        lww(ii,jj,nwet) = jc5(i,j,k)
        end do
        enddo
        enddo



         do k = 1,ilo
             llay(0,k) = 0
         end do

         do j=2,nz
          do lw = lb(j),le(j)
            if (ljumm(lw).eq.1) then
             nwet = indwet(lw)
             ldown = lazc(lw) 
             ld  = ldown-1 
             lein= min0(1,ld) 
             if (lein.eq.1) then
                do k = ld,lob,-1 
                    llay(0,k) = llay(0,k)+1
                    llay(llay(0,k),k) = nwet+k
                end do
             end if
            end if
          end do
         end do

cc         init = .false.
      end if
c





*------------Volume of the 1st layer---------------
      volbio=0.
      nwet=0
      do j = 1,n 
      do lw = lb(j),le(j) 
             i = iwet(lw) 
                           zzz=dble(pd(1))+dble(zac(lw))
       if(lazc(lw).eq.1) zzz=dfloat(ldep(lw))+dble(zac(lw))

            Dvoltj(i,1,j)=zzz*dprod(dl,dln(i))
            volbio=volbio+Dvoltj(i,1,j)            
            dd(nwet+1) =  max(epsilon,zzz)
            dd_dd(nwet+1)=dprod(dd(nwet+1),dd(nwet+1))

       nwet=nwet+lazc(lw)
      end do
      end do
      volbio = volbio+vol0

       do ll=1,ndrei
       r_dt_dd(ll)=dt/dd(ll)
       enddo

c 
       wmax=0.
       umax=0.
       admax=avmax 
       dtdt = 2.*dt

CCC  EPSIL is a minimum value used by the TVD procedure
      EPSIL=1.0d-5
c 
      avmin = 0.134e-6 
      stfak = 1./1.35 

*-----------------------------
       do ll = 1,ndrei
           Tc0(ll) = Tc(ll,ISPEC)
       end do
       Tc0(0) = 0.       
c--------------------------------------------

c------------------- transport cycle -------------------------------------

      do j = 2,nz
      lwa = lb(j)
      lwe = le(j)
c#else
c       lwa = lb(2)
c       lwe = le(nz)
c#endif

      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
c#ifndef MPI
c              j = jwet(llw(ll))
c#endif
             i = iwet(llw(ll))
             k = ll-indwet(llw(ll))

             if (ljumm(llw(ll)).eq.1) then
       if (j.lt.n-2 .and. i.lt.m-2) then  
*------------ horizontal Advektion terms, TVD -----------

        uo = uc(ll)*min(lww(0,1,ll),1)
     &  /( 0.5* ( dd(ll)+dd(lww(0,1,ll)) ) )
        uot=uc(ll)* min(lww(0,1,ll),1)

       COURUR=abs(dtdx(i)*uo)

       lu = lww(0,-1,ll)
       if (lu.eq.0) then
          uw = 0.
          uwt= 0.0d+0
       else
          uw = uc(lu)/( 0.5* ( dd(ll)+dd(lww(0,-1,ll)) ) )
          uwt= uc(lu)
       end if

        COURUL=abs(dtdx(i)*uw)

        lv = lww(-1,0,ll)
        if (lv.eq.0) then
            vn = 0.
            vnt= 0.0d+0
        else
            vn=vc(lv)/( 0.5* ( dd(ll)+dd(lww(-1,0,ll)) ) )

           vnt= vc(lv)
        end if

         COURVO=abs(dt_lr*vn)

        vs=vc(ll)* min(lww(1,0,ll),1)
     &  /( 0.5* ( dd(ll)+dd(lww(1,0,ll)) ) )
        vst=vc(ll)*min(lww(1,0,ll),1)

         COURVU=abs(dt_lr*vs)
*------------ vertical advection terms ------------------------
!       SinkD_o =-BioC(23)*r_dt_dd(ll)    !BioC(23) [m/s]
!       SinkD_u = SinkD_o
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 8.644 m/day ~ 0.000116 m/s  SINKING velosity
        SinkD_o=0.
        SinkD_o( HGDET3D ) =-0.000058*r_dt_dd(ll) !    5 m/day
        SinkD_o( HGPOC3D ) =-0.000058*r_dt_dd(ll) !    5 m/day
        SinkD_o( HGS_P3D ) =-0.000058*r_dt_dd(ll) !    5 m/day
        SinkD_o( HGCYA3D ) = 0.000001157*r_dt_dd(ll)    ! -0.1 m/day upwelling
        SinkD_o( MEHGDET3D ) =-0.000058*r_dt_dd(ll) !    5 m/day
        SinkD_o( MEHGPOC3D ) =-0.000058*r_dt_dd(ll) !    5 m/day
        SinkD_o( MEHGCYA3D ) = 0.000001157*r_dt_dd(ll)    ! -0.1 m/day upwelling
        SinkD_u = SinkD_o             
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
      
      if (k.eq.1) then
         llo = ll
      else
         llo = ll-1
cc          llo=ll
       endif

      if (k.eq.lazc(llw(ll))) then
         llu = ll
      else
         llu = ll+1
ccc       llu=ll
      endif
      wo = wc(ll)         !	k
      wu = wc(llu)        !	k+1
      
      if (k.eq.1) then
        wo = 0.
      endif

      if (k.eq.lazc(llw(ll))) then
          wu=0.
      endif

      wo=wo
      Wu=wu

      dwdz=(wo-wu)*r_dt_dd(ll)

!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM cycle for ibio
!      do ibio= ibio0,nbio-nsed
!     do ibio = 1,nchem

c west
      RSUL=0.0d+0
      do irsu=1,min0(jc5(i,j-1,k),1)
      DENUM=Tc0(ll)-Tc0(lww(0,-1,ll))
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.0d+0
      else
       if (uwt.ge.0.0d+0) NSIG=-1
       if (uwt.lt.0.0d+0) NSIG=1
       if ((NSIG.eq.-1 .and. jc5(i,j-2,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i,j+1,k).eq.0)) then
       PHIR=0.

       else
       RATIO=(Tc0(lww(0,NSIG,ll))-Tc0(lww(0,NSIG-1,ll)))/DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(uwt)*(1.-COURUL)*D1
       RSUL=D2+0.5*(dprod((uwt+abs(uwt)),Tc0(lww(0,-1,ll)))+
     &                dprod((uwt-abs(uwt)),Tc0(ll)))
       RSUL=RSUL*dble(dl)
       enddo

c east
      RSUR=0.0d+0
      do irsu=1,min0(jc5(i,j+1,k),1)
      DENUM=(Tc0(lww(0,1,ll))-Tc0(ll) )
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.0d+0
      else
       if (uot.ge.0.0d+0) NSIG=-1
       if (uot.lt.0.0d+0) NSIG=1
       if ((NSIG.eq.-1 .and. jc5(i,j-1,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i,j+2,k).eq.0)) then

       PHIR=0.0d+0
       else
c
       RATIO=(Tc0(lww(0,NSIG+1,ll))-Tc0(lww(0,NSIG,ll)))
     & /DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(uot)*(1.-COURUR)*D1
       RSUR=D2+0.5*(dprod((uot+abs(uot)),Tc0(ll))+
     &                dprod((uot-abs(uot)),Tc0(lww(0,1,ll))))
       RSUR=RSUR*dble(dl)
       enddo
c north
      RSVO=0.0d+0
      do irsu=1,min0(jc5(i-1,j,k),1)
      DENUM=(Tc0(lww(-1,0,ll))-Tc0(ll) )
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
       if (vnt.ge.0.0d+0) NSIG=1
       if (vnt.lt.0.0d+0) NSIG=-1
       if ((NSIG.eq.-1 .and. jc5(i-2,j,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i+1,j,k).eq.0)) then
       PHIR=0.
       else
       RATIO=(Tc0(lww(NSIG-1,0,ll))-Tc0(lww(NSIG,0,ll)))
     & /DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(vnt)*(1.-COURVO)*D1
       RSVO=D2+0.5*(dprod((vnt+abs(vnt)),Tc0(ll))+
     &                dprod((vnt-abs(vnt)),Tc0(lww(-1,0,ll))))
       RSVO=RSVO*dvo(i)
       enddo
c south
      RSVU=0.0d+0
       do irsu=1,min0(jc5(i+1,j,k),1)
      DENUM=(Tc0(ll)-Tc0(lww(1,0,ll)))
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
       if (vst.ge.0.0d+0) NSIG=1
       if (vst.lt.0.0d+0) NSIG=-1
      if ((NSIG.eq.-1 .and. jc5(i-1,j,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i+2,j,k).eq.0)) then

       PHIR=0.
       else
       RATIO=(Tc0(lww(NSIG,0,ll))-Tc0(lww(NSIG+1,0,ll)))/DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(vst)*(1.-COURVU)*D1
       RSVU=D2+0.5*(dprod((vst-abs(vst)),Tc0(ll))+
     &                dprod((vst+abs(vst)),Tc0(lww(1,0,ll))))
       RSVU=RSVU*dvu(i)
       enddo

*-------- the horizontal advection part to the new value -------------------
      XYadv(ll)=(RSUL-RSUR+RSVU-RSVO)*dt/Dvoltj(i,k,j)

! additional sinking velosity

cc      if (ibio.eq.4.or.ibio.eq.8) then  ! Sinking for POM
         D_o2=SinkD_o(ISPEC)      ! in the "upper" layer
         D_u2=SinkD_u(ISPEC)      ! in the "bottom" layer

        D_o2 = 0.
        D_u2 = 0.

         if(k.eq.1) D_o2=0.
         if(k.eq.lazc(llw(ll))) D_u2=0.   
c      TVD scheme
      SADV=0.0d+0
c      from above
      RSWO=0.0d+0
       if (k.gt.1) then
       DDI=(dd(ll)+dd(llo))*0.5
       COURWO=abs(wo)*dt/DDI
      DENUM=( Tc0(llo)-Tc0(ll))/DDI
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
        if (wo.ge.0.) NSIG=1
        if (wo.lt.0.) NSIG=-1
        if ((k.eq.lazc(llw(ll)) .and. NSIG.eq.1) .or. (k.eq.2 .and. NSIG
     &          .eq.-1)) then
          PHIR=0.
        else
        DDIZ=(dd(ll+NSIG)+dd(llo+NSIG))*0.5
        RATIO=( Tc0(llo+NSIG)-Tc0(ll+NSIG))/(DDIZ*DENUM)
        PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
        PHIR=max(0.,PHIR1)
        endif
      endif
      D1=PHIR*DENUM*DDI
      D2=0.5*abs(wo)*(1.-COURWO)*D1
      RSWO=D2+0.5*(dprod((wo+abs(wo)),Tc0(ll))+
     &               dprod((wo-abs(wo)),Tc0(llo)))
      endif
c      from below
      RSWU=0.0d+0
       if (k.lt.lazc(llw(ll))) then !!!!
       DDI=(dd(ll)+dd(llu))*0.5 
      COURWU=abs(wu)*(dt)/DDI
      DENUM=( Tc0(ll)-Tc0(llu))/DDI
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
      if (wu.ge.0.) NSIG=1
      if (wu.lt.0.) NSIG=-1
      if ((k.eq.lazc(llw(ll))-1 .and. NSIG.eq.1) .or. (k.eq.1 .and. NSIG
     &          .eq.-1)) then
      PHIR=0.
      else
      DDIZ=(dd(ll+NSIG)+dd(llu+NSIG))*0.5
       RATIO=( Tc0(ll+NSIG)-Tc0(llu+NSIG))/(DDIZ*DENUM)
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
      endif
      endif
      D1=PHIR*DENUM*DDI
      D2=0.5*abs(wu)*(1.-COURWU)*D1
      RSWU=D2+0.5*(dprod((wu+abs(wu)),Tc0(llu))+
     &               dprod((wu-abs(wu)),Tc0(ll)))
      endif

c vertical advection contribution

      SADV=(RSWU-RSWO)*r_dt_dd(ll)


cc----sinking vel is added by upstream---  NEW
        w_up   = max(0.,D_o2*Tc0(ll))  +min(0.,D_o2*Tc0(llo))
        w_down = max(0.,D_u2*Tc0(llu)) +min(0.,D_u2*Tc0(ll))
cc---------------------------------------
cic       w_up   = max(0.,0.01*Tc0(ll))  +min(0.,0.01*Tc0(llo))
cc       w_down = max(0.,0.01*Tc0(llu)) +min(0.,0.01*Tc0(ll))
      svert(ll)= SADV
     &        -w_up+w_down             ! ----sinking---
     
      if (k.eq.1) then
      svert(ll) = svert(ll) - zinc(i,j)*Tc0(ll)
     & /dd(ll)
      endif ! if (k.eq.1) then

!      end do       ! ib = 1, nbio  

      endif !if (ljumm(llw(ll)).eq.1) then

             end if !if (j.lt.n-2 .and. i.gt.5 .and. i.lt.m-2) then  
         end do  !i
      end if
c#ifdef MPI
      end do      !j
c#endif

222   format(1x,3(e9.3,2x))
c 
c -----------vertikal diffusion ---------------- 
c 
!                do ibio=ibio0,nbio
!        do ibio = 1,nchem
c#ifdef MPI
      do j = 2,nz 
      lwa = lb(j)
      lwe = le(j)
c#else
c       lwa = lb(2)
c       lwe = le(nz)
c#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
             if (ljumm(llw(ll)).eq.1) then
                    salref = Tc0(ll)
                    s (ll)  = Tc0(ll)
                    sn(ll)  = Tc0(ll)
                    stfak=szahl(ll)
!HHHHHHHHHHH avc - viscosity, avd - diffusivuty
        avd(ll) = amin1(stfak*avc(ll),admax)
        if(avd(ll).le.0.) avd(ll)=admax
              avd(ll)=max(avd(ll),avmin)
cccc         avd(ll)=0.
 
             end if
         end do
      end if
c#ifdef MPI
      end do
c#endif
!                end do ! ibio = 1,nchem
!      do ibio=ibio0,nbio
!      do ibio = 1,nchem
c#ifdef MPI
      do j = 2,nz 
      lwa = lb(j)
      lwe = le(j)
c#else
c       lwa = lb(2)
c       lwe = le(nz)
c#endif
      if (lwa.le.lwe) then
         llb = indwet(lwa)+1
         lle = indwet(lwe)+lazc(lwe)
         do ll = llb,lle
             k = ll-indwet(llw(ll))
             if (ljumm(llw(ll)).eq.1) then
                ldown = lazc(llw(ll))
                if (ldown.gt.1.and.k.gt.1.and.k.le.ldown-1) then
         tsa(ll) = dprod(dtdt,avd(ll))  
     1    /(dprod(dd(ll-1),dd(ll))+dd_dd(ll)) 
         tsc(ll) = dprod(dtdt,avd(ll+1)) 
     1    /(dprod(dd(ll+1),dd(ll))+dd_dd(ll)) 
         tsb(ll) = 1.d0+vtmit*tsa(ll)+vtmit*tsc(ll) 
         tsbb   = 1.d0-vtmic*tsa(ll)-vtmic*tsc(ll) 
         sd(ll) = vtmic*tsa(ll)*s(ll-1)+tsbb
     1    *s(ll)+vtmic*tsc(ll)*s(ll+1) 
                end if
             end if
         end do
      end if
c#ifdef MPI
      end do
c#endif
!      end do   ! ibio = 1,nchem
c
c --------- ende vertikal diffusion 
c
!             do ibio=ibio0,nbio
!        do ibio = 1,nchem
      do j = 2,nz 
      do lw = lb(j),le(j)
         nwet = indwet(lw)
         ldown = lazc(lw) 
         if (ljumm(lw).eq.1) then
                ld  = ldown-1 
                lein= min0(1,ld) 
                if (lein.eq.1) then

      tsa(nwet+1) = 0.d0 
      tsc(nwet+1) = dprod(dtdt,avd(nwet+2)) 
     1 /(dprod(dd(nwet+2),dd(nwet+1))+dd_dd(nwet+1) ) 
      tsb(nwet+1) = 1.d0+vtmit*tsc(nwet+1) 
      tsbb   = 1.d0-vtmic*tsc(nwet+1) 
      sd(nwet+1) = tsbb*s(nwet+1)+vtmic*tsc(nwet+1)
     1 *s(nwet+2) 

      tsa(nwet+ldown) = dprod(dtdt,avd(nwet+ldown)) 
     1 /(dprod(dd(nwet+ld),dd(nwet+ldown))+dd_dd(nwet+ldown) ) 
      tsc(nwet+ldown) = 0.d0 
      tsb(nwet+ldown) = 1.d0+vtmit*tsa(nwet+ldown) 
      tsbb        = 1.d0-vtmic*tsa(nwet+ldown) 
      sd(nwet+ldown)  = vtmic*tsa(nwet+ldown)*s(nwet+ld)
     1 +tsbb*s(nwet+ldown) 

                end if
         end if
      end do
      end do
!             end do    ! ibio = 1,nchem

!                do ibio=ibio0,nbio
!        do ibio = 1,nchem
c#ifdef MPI
      do j = 2,nz 
      do lw = lb(j),le(j)
c#else
c       do lw = lb(2),le(nz)
c#endif
         if (ljumm(lw).eq.1) then
             nwet = indwet(lw)
             ldown = lazc(lw) 
             ld  = ldown-1 
             lein= min0(1,ld) 
             if (lein.eq.1) then
                    tsal(nwet+lob) =tsb(nwet+lob) ! lob = 1 (or 2)
                    sbet(nwet+lob) = sd(nwet+lob) 
             end if
         end if
      end do
c#ifdef MPI
      end do
c#endif
!                end do ! ibio = 1,nchem
c
!                do ibio=ibio0,nbio
!        do ibio = 1,nchem
      do k = 1,ilo
!CDIR VECTOR NODEP
         do lh = 1,llay(0,k)
             ll = llay(lh,k)
                        tsal(ll+1) = tsb(ll+1)-(vtmit
     1                   *tsa(ll+1)/tsal(ll))*vtmit
     2                   *tsc(ll) 
                        sbet(ll+1) = sd(ll+1)+(vtmit
     1                   *tsa(ll+1)
     1                   /tsal(ll))*sbet(ll) 
         end do
      end do
!                end do ! ibio = 1,nchem
c
!                do ibio=ibio0,nbio
!        do ibio = 1,nchem
c#ifdef MPI
      do j = 2,nz 
      do lw = lb(j),le(j)
c#else
c       do lw = lb(2),le(nz)
c#endif
         if (ljumm(lw).eq.1) then
             nwet = indwet(lw)
             ldown = lazc(lw) 
             ld  = ldown-1 
             lein= min0(1,ld) 
             if (lein.eq.1) then
                    sn(nwet+ldown) = sbet(nwet+ldown)
     1               /tsal(nwet+ldown) 
             end if
         end if
      end do
c#ifdef MPI
      end do
c#endif
!                end do ! ibio = 1,nchem

!             do ibio=ibio0,nbio
!        do ibio = 1,nchem
      do k = ilo,1,-1
!CDIR VECTOR NODEP
         do lh = 1,llay(0,k)
             ll = llay(lh,k)
                sn(ll) = (sbet(ll)+vtmit
     1            *tsc(ll)*sn(ll+1))
     2            /tsal(ll) 
         end do
      end do
!             end do    ! ibio = 1,nchem
!      do ibio = ibio0,nbio-nsed
!        do ibio = 1,nchem
c#ifdef MPI
      do j = 2,nz
         lwa = lb(j)
         lwe = le(j)
c#else
c          lwa = lb(2)
c          lwe = le(nz)
c#endif
         if (lwa.le.lwe) then
             llb = indwet(lwa)+1
             lle = indwet(lwe)+lazc(lwe)
             do ll = llb,lle
                if (ljumm(llw(ll)).eq.1) then

cc        Tc(ll,ISPEC)= XYadv(ll)+svert(ll)
        Tc(ll,ISPEC) =  
     + sn(ll)
     +  + XYadv(ll)
     +   +svert(ll)
                end if
             end do
         end if
c#ifdef MPI
      end do
c#endif
!      end do   ! ibio = 1,nchem
c       if(nbio.gt.2)then
c       call bio(Tc,dd,dz,sh_wave,sh_depth) ! biological source	DXbio(k,nbio)
c       end if

*________ end of   find biology  _______________________

      return 
      end 
c-----------------------------------------------------------------------
      subroutine deco1d3d (u,ucomp,ntot,xdef) 
c-----------------------------------------------------------------------
c      1-d-field ->  3-d-field 
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter (khor1=khor+1) 
      dimension u(m,n,ilo) 
      dimension ucomp(ntot) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 

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
             u(i,k,j) = ucomp(nwet) 
          enddo
    2 continue 


      return 
      end

c-----------------------------------------------------------------------
c adds TC( :,X ) values to chem output 3D fields - added by Johannes Bieser 27.05.2012
c-----------------------------------------------------------------------
      subroutine add1dto3d (u,ucomp,ntot) 
c-----------------------------------------------------------------------
c      1-d-field ->  3-d-field 
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter (khor1=khor+1) 
      dimension u(m,n,ilo) 
      dimension ucomp(ntot) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 

      lwe = 0 
      nwet = 0 
      do 2 k=1,n 
        do i=1,m 
        do j=1,ilo 
!         u(i,k,j) = xdef 
        enddo
        enddo
        lwa = lwe+1 
        lwe = indend(k) 
        do 2 lw=lwa,lwe 
          lump = lazc(lw) 
          i = iwet(lw) 
          do j=1,lump 
             nwet = nwet+1 
             u(i,k,j) = u(i,k,j) + ucomp(nwet) 
          enddo
    2 continue 


      return 
      end

c-----------------------------------------------------------------------
c clears chem fields and marks dry grid cells - added by Johannes Bieser 27.05.2012
c-----------------------------------------------------------------------
      subroutine clear3d (u,ntot,xdef) 
c-----------------------------------------------------------------------
c      1-d-field ->  3-d-field 
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter (khor1=khor+1) 
      dimension u(m,n,ilo)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     * indend(n),isornr(n),isorsr(n),islab(n) 

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
             u(i,k,j) = 0.
          enddo
    2 continue 


      return 
      end


c-----------------------------------------------------------------------
      subroutine settime(time,ijulu,ihouu,iminu,isecu)
c-----------------------------------------------------------------------
c      time set to rhour,rminu,rsecu
c-----------------------------------------------------------------------

      itime = int(time)
      isecu = mod(itime,60)
      iminu = mod((itime-isecu),3600)/60
      ihouu = mod((itime-(iminu*60)-isecu),86400)/3600
      ijulu = (itime-(ihouu*3600)-(iminu*60)-isecu)/86400+1

      return
      end



c-----------------------------------------------------------------------
      subroutine xstat(para,ianz,pmax,pmin,pmean,psum)
c-----------------------------------------------------------------------

      dimension para(ianz)

      pmax  = -9999999.9
      pmin  =  9999999.9
      psum  = 0.0
      pmean = 0.0

      do i=1,ianz
        pmax = amax1(para(i),pmax)
        pmin = amin1(para(i),pmin)
        psum = psum+para(i)
      enddo
      pmean = psum/real(ianz)

      return
      end

c-----------------------------------------------------------------------
      subroutine xstat2(para,m,n,jjc,pmax,pmin,pmean,psum)
c-----------------------------------------------------------------------

      dimension para(m,n),jjc(m,n)

      pmax  = -9999999.9
      pmin  =  9999999.9
      psum  = 0.0
      pmean = 0.0
      ianz  = 0

      do i=1,m
      do k=1,n
      do j=1,jjc(i,k)
        ianz = ianz+1
        pmax = amax1(para(i,k),pmax)
        pmin = amin1(para(i,k),pmin)
        psum = psum+para(i,k)
      enddo
      enddo
      enddo
      pmean = psum/ianz

      return
      end


c-----------------------------------------------------------------------
      subroutine means(m,n,ilo,ntot2d,ntot3d,dz,izet,iindex,ltief,dzbod,
     &   jjc,rval,zac,janf,jend,ianf,iend,rvalmean)
c-----------------------------------------------------------------------

      integer jjc(m,n),ltief(m,n),iindex(m,n),izet(m,n)
      real dz(ilo),dzbod(m,n),rval(ntot3d),zac(ntot2d)
      double precision rvalsum,ddzsum

      do j=janf,jend
      do i=ianf,iend
        if(jjc(i,j).eq.1)then
          do k=1,ltief(i,j)
             if((k.eq.1).and.(ltief(i,j).gt.1)) ddz=dz(1)
             if(k.eq.ltief(i,j))                   ddz=dzbod(i,j)
             if((k.gt.1).and.(k.lt.ltief(i,j))) ddz=(dz(k)-dz(k-1))
             if(k.eq.1) ddz=ddz+zac(izet(i,j))
             rvalsum=rvalsum+dble(ddz*rval(iindex(i,j)+k))
             ddzsum=ddzsum+dble(ddz)
c             write(*,*) rval(iindex(i,j)+k),rvalsum,ddz,ddzsum,i,j,k,
c    &                    iindex(i,j)+k,izet(i,j)
          enddo
        endif
      enddo
      enddo

      rvalmean=real(rvalsum/ddzsum)
      write(*,*) rvalmean

      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      function sigma (s,t,p) 
c-----------------------------------------------------------------------
c     testrechnung mit autodbl(dblpad) ist erfolgreich gewesen 
c     testdaten:  s=0. t=25. p=1000.  rho=1037.90204 
c                 s=35.t=5.  p=0.     rho=1027.67547 
c                 s=35.t=25. p=0.     rho=1023.34306 
c                 s=35.t=25. p=1000.  rho=1062.53817 
c                 kai jancke13.01.83 
c-----------------------------------------------------------------------
c 
c                       i n p u t 
c 
c **********   p must be given as depth in meters * 0.1  ********** 
c **********   s must be given as salinity in promille   ********** 
c **********   t must be given as degree centigrade      ********** 
c 
c-----------------------------------------------------------------------

      real i0,i1,i2,j0,m0,m1,m2,k0,k1,k2,ksbm
      common /coef/b0,b1,b2,b3,b4,c0,c1,c2,d0,a0,a1,a2,a3,a4,a5
     ,,f0,f1,f2,f3,g0,g1,g2,i0,i1,i2,j0,m0,m1,m2,e0,e1,e2,e3,e4
     ,,h0,h1,h2,h3,k0,k1,k2
      s2=s*s 
      s3=s2*s 
      t2=t*t 
      t3=t2*t 
      t4=t2*t2 
     
*************** density of reference pure water (snow) 
***************  density at standard atmosphere (p=0) 
      rhonul=   a1*t+a2*t2+a3*t3+a4*t4+a5*t4*t 
     +           +(b0+b1*t+b2*t2+b3*t3+b4*t4)*s 
     +           +(c0+c1*t+c2*t2) *sqrt(s3)+d0*s2 
***************  secant bulk modulus 
      ksbm =  e0+e1*t+e2*t2+e3*t3+e4*t4 
     +          + (f0+f1*t+f2*t2+f3*t3)*s + (g0+g1*t+g2*t2)*sqrt(s3) 
     +     + (h0+h1*t+h2*t2+h3*t3 
     +          + (i0+i1*t+i2*t2)*s       + j0*sqrt(s3)) * p 
     +     + (k0+k1*t+k2*t2 
     +          + (m0+m1*t+m2*t2)*s) * p * p 
      
      sigma  = 1.e-3*(rhonul+a0*p/ksbm)/(1.-p/ksbm) 


      return 
      end 

c-----------------------------------------------------------------------
      subroutine setice(iindex,izet,mm,nn)
c-----------------------------------------------------------------------
c     initialisation for ice model and sea ice parameter
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter(nx=n-1,ny=m-1)
      parameter(khor1=khor+1)   

      dimension iindex(mm,nn),izet(mm,nn)

      include 'somice' 
      
      dtsec=dt

c     drag coefficient for ice water stress 
      cwa = 0.0025

c     volumetric specific heat of sea-water cw (j/m**3 k)) 
c      ccw = 4.07e6         
 
c     latent heat of fusion of ice (j/kg) 
      hlatis = 3.36e5 
 
c     density of sea ice (kg/m**3) 
      rois = 930. 
 
c     thermal conductivity of ice (j/(m s k)) 
      tci = 2.033 
 
c     ice-salinity (psu) 
c      sice = 5. 
 
c     minimum ice-thickness 
      epsis = 1.e-7 
 
c     solar constant 
      solcon  = 1365.0 
 
c     cloud cover between 1(min) and 8(max) 
c      icloud = 5 
 
c     air temperature (if no data) 
       airtem = -15.0 
 
c     start (in hours) for shortwave radiation 
c      start = 30.*5.*24.0 * 3600. 


c-----derived constants  for ice-model 
c----- no changes here
 
c     freezing temperature of seaice 
c      tf = tfreez(sice,0.0) 

c     volumetric latent heat of fusion of ice , (j/m**3) 
      roil = rois*hlatis 
 
c     dt divided by latent heat of fusion of ice 
      dtroil = dtsec/roil 


      return 
      end 

c-----------------------------------------------------------------------
       real function tfreez (s,p)
c-----------------------------------------------------------------------
c     freezing point of sea water (millero 1978)
c-----------------------------------------------------------------------

       data ca/-0.0575/,cb/1.710523e-3/,cc/-2.154996e-4/,cd/-7.53e-3/

       tfreez = s*(ca+cb*sqrt(s)+cc*s)+cd*p
       return
       end
c------------------------------------------------------------------------
      subroutine konti(ht,dzbod,wtest,wsurf,kkhor)
c-----------------------------------------------------------------------
c     solving of continuity equation, calculating w-component  
c-----------------------------------------------------------------------
      include 'C_model_inc'
c Attention!!! ngro must be recalculated depending on model domain
c              ngro = (max(m*(ilo*20+10),(kasor*8)))
      parameter (ngro=m*(ilo*40+10)) ! north_b

      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)

      double precision wo, wsurf(kkhor)
      real ht(m,n),dzbod(m,n),wtest(m,n)

      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     cindend(n),isornr(n),isorsr(n),islab(n)

      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)

      common /dreh/ sinfu(m),sinfv(m),cosfu(m),cosfv(m),sincx(m),
     ssincy(m),bx(m),by(m),pxu(m),pyv(m),pyu(m),pxv(m)

      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     ddlvo(m),dlvu(m),gh,rdt,dt2,r4,dtrdln(m),dtdlr


      common/num/dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo),
     pprd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo),
     qqa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1)


      common /aux/ u(ilo,m,2),v(ilo,m)
      mz = m-1
      nz = n-1
c 
c     prepare moving  arrays for  kk = 1 
c 
      nwet = 0
      lwa = 1
      lwe = indend(1)
      do 4 lw = lwa,lwe
      i = iwet(lw)
      ldown = lazc(lw)
      do 3 j = 1,ldown
      u(j,i,2) = uc(nwet+j)
   3  continue
      nwet = nwet+ldown
   4  continue
c 
c 
      do 76 kk = 2,nz
c 
c 
c     --------------------- shifting slabs west
c 
      do 16 i = 1,m
      do 16 j = 1,ilo
      u(j,i,1) = u(j,i,2)
      u(j,i,2) = 0.0
      v(j,i) = 0.0
   16 continue
c 
c     --------------------- fill w for  slab kk 
      lwa = indend(kk-1)+1
      lwe = indend(kk)
      nwet = islab(kk-1)
      do 76 lw = lwa,lwe
      i = iwet(lw)
      im1 = max0(i-1,1)
      ldown = lazc(lw)
      dddvo = dlvo(i)*dlr
      dddvu = dlvu(i)*dlr
      wo = 0.d0

c     ---------------- up-down-loop ----------------------- 
c     ---------------------------- (u,v,w  up,  rest down)- 
      do 97 jj = ldown,1,-1
      u(jj,i,2) = uc (nwet+jj)
      v(jj,i) = vc (nwet+jj)
c     ------------------ continuity equation for  k + 1 
      wo = wo+dprod(rdln(i),u(jj,i,1))-dprod(rdln(i),u(jj,i,2))
     ++dprod(dddvu,v(jj,i))-dprod(dddvo,v(jj,im1))
   97 wc(nwet+jj) = wo
c     ------- w at sea surface (compressed) -------- 
      wobc(lw) = wo
      wsurf(lw)=wo
c 
      nwet = nwet+ldown
   76 continue


      return
      end


c------------------------------------------------------------------------
