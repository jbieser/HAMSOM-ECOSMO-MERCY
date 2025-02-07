c-------------------------------------------------------------------------
c       calculating under water light
c---------------------------------------------------------------
        SUBROUTINE LIGHT( RGRND, Achem, Abio, RADI, df )

        USE CPARAM
        USE COUT
        USE CIN

C...........
!
! DESC:         Light extinction
!
! HISTORY:      13.10.2021 Final Version MERCY v2.0
!
! AUTHOR:       johannes.bieser@hereon.de
!
! LICENSE:      GPL Version 3
!
C...........


!        IMPLICIT NONE

        include 'C_loads_inc'
        include 'C_model_inc'
        include 'C_declaration_inc'

!       Input/Output variables
        REAL, DIMENSION( m,n )           , INTENT(IN)    :: RGRND
        REAL, DIMENSION( ndrei,nchem )   , INTENT(IN)    :: Achem
        REAL, DIMENSION( ndrei,3:ninbio ), INTENT(IN)    :: Abio
        REAL, DIMENSION( ndrei )         , INTENT(INOUT) :: RADI
        REAL, DIMENSION( ilo )           , INTENT(IN)    :: df

!       Local variables
        REAL, DIMENSION( 0:ilo ) :: EXwater
        REAL, DIMENSION( 0:ilo ) :: EXphyto
        REAL, DIMENSION( 0:ilo ) :: EXdom
        REAL, DIMENSION( 0:ilo ) :: EXdet
        REAL, DIMENSION( 0:ilo ) :: EXtot
        REAL                        TPHYT
        REAL                        TPART
        REAL                        TDISS

!        INTEGER lwe, nwet, k, lwa, lwe, lw, i, lump, jj, ii
        
!       Local parameters
        REAL, PARAMETER :: EX_WATER = 0.05     !light extinction [1/m]
        REAL, PARAMETER :: EX_PHYTO = 0.000377 ! 0.03 / ( 6.625 * 12.01 ) !phyto self shading [m**2/(mgC)]
        REAL, PARAMETER :: EX_DOM   = 0.00029  ![m**2/(mgC)]
        REAL, PARAMETER :: EX_DET   = 0.00020  ![m**2/(mgC)]

!       REAL, PARAMETER :: PSR = 0.16   ![]      !Wozniak, 2010
!       (2.3%-42.0%)
        REAL, PARAMETER :: PSR = 0.10   ![]      !Sharif, 2014 (9.3%-10.8%)
!       REAL, PARAMETER :: PSR = 1.00   ![]      !No SPM case (POC/SPM
!       ratio)

!       Begin of program

        !1. Set surface values (k=0) to zero
        EXwater( 0 ) = 0.0
        EXphyto( 0 ) = 0.0
        EXdom( 0 )   = 0.0
        EXdet( 0 )   = 0.0
        Extot( 0 )   = 0.0

        !2. Calculate cummulative extinction coefficient at layer bottom
        DO jj = 1,ilo
           EXwater(jj) = EXwater(jj-1) + EX_WATER * df(jj)
        END DO

        lwe = 0
        nwet = 0

        do k = 1,n
          lwa = lwe+1
          lwe = indend(k)

          do lw = lwa,lwe
            i = iwet(lw)
            lump = lazc(lw)

            do jj = 1,lump
               nwet = nwet+1
               ii = nwet
               
               TPHYT = Abio( ii,FLA ) + Abio( ii,DIA ) + Abio( ii,CYA )
               EXphyto(jj) = EXphyto(jj-1) + TPHYT * EX_PHYTO * df(jj)

               TDISS = Abio( ii,DOM ) + Achem( ii,DTOM3D )
               EXdom(jj) = EXdom(jj-1) + TDISS * EX_DOM * df(jj)

!!!     19.02.2021 added additional particle load based on Sharif et !al., 2014 (PSR)
               TPART = Abio( ii,DET ) / PSR + Achem( ii,PTOM3D )
               EXdet(jj) = EXdet(jj-1) + TPART * EX_DET * df(jj)

               EXtot(jj)    = EXwater(jj) + EXphyto(jj) + EXdom(jj) + EXdet(jj)

!               RADI(ii) = RGRND(i,k) * EXP( -EXtot(jj-1) ) ! Radiation at cell surface (old version)
!               RADI(ii) = RGRND(i,k) * EXP( -EXtot(jj) )   ! Radiation at cell bottom

!              Average extinction
               RADI( ii ) = RGRND(i,k) * EXP( -( EXtot(jj-1) + (EXtot(jj) - EXtot(jj-1) ) / 2  ) )

            enddo
          enddo
        enddo

        return

        END SUBROUTINE LIGHT
