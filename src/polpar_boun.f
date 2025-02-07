!====================================================================================
      subroutine tsrneu3s(izet,iindex,llvrp,mm,nn,n3d,j_year,lmon,
     &                    rsteps,julday)
c-----------------------------------------------------------------------
c     boundary conditions, tracer
c-----------------------------------------------------------------------

      USE COUT

      include 'C_model_inc'

      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 
      dimension trum(llvrp),srum(llvrp)
cc     
      integer izet(m,n)   
      dimension iindex(mm,nn)
      dimension f(m,n,ilo)
c      
      parameter (nr=97) ! Number of rivers

      REAL  nut1sumN(ilo),nut1winN(ilo),nut1sumE(ilo),nut1winE(ilo)
 
      common/T_4d/ Tc(ndrei,nchem),Tc_mit(ndrei,nchem) 
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     cindend(n),isornr(n),isorsr(n),islab(n) 
c 
      integer nord,sued,ost,west 
      integer*2 isor,ksor 
c 
      common /intit/ itel(lrp),jc(m,3),ikk(kasor),isor(kasor), 
     cksor(kasor),nofrez,lone,itende,nfreez 
c 
      common /iterat/ omega(kasor),za(m,3),cx(m,3),cy(m,3),rhoz(m,3), 
     ctcc(8) 
c 
c      common /uirand/ ilwrnd(lrp),islrnd(lrp),ldownr(lrp), 
c     ,              nord,sued,ost,west,kna,kne,ksa,kse,iea,iee,iwa,iwe 
c 
c 
ccc   common /tide/ costom,sintom,cosinc,sininc,hcosg(lrp),hsing(lrp) 
c 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 
c 
c 
       common / iceconc/ icemold(m,n) , delconc (m,n, nbio)
          real concice (nbio) 

!     atmospheric flux
      common /atmosphere/ atmflux(m,n,nbio,12,5)
c     LOCAL for rivers
 !     common /rivers/ nvassdr(nr),ix_riv(nr),iy_riv(nr),dis_av(nr,nbio)
       common /rivers/ nvassdr(nr),ix_riv(nr),iy_riv(nr),
     &                 dis_yr(nr,4),dis_mn(nr,12,4)


!      print *,'nbio',nbio
!----------------------------------------------------------------------
!                 open ocean boundaries 
!----------------------------------------------------------------------
!	C. Divisor to raise or lower Hg BC
	XBC_FAC = 1.0

! 1st order cond.: 3 gridpoints across the boundary should be equal
c....Western boundary: La Manch: ~48-51 N along ~-2 W ................
        do i=156,164       ! "i" - index from N to S
          l=iindex(i,16)   ! here we list 3 poins from W to E ("j")
          ll=iindex(i,17)
          lll=iindex(i,18)
          lh=izet(i,16)
          ltt=lazc(lh)
         do lt=1,ltt   ! depth-index
!             do iii=1,nbio-nsed
             do iii=1,nchem

!               C. Fixed boundaries
                 IF (iii .EQ. HGDIS3D) THEN
                     Tc(lll+lt,iii) = 0.8
                 ELSE IF (iii .EQ. MEHGDIS3D) THEN
                     Tc(lll+lt,iii)  = 0.06
                 ELSE IF (iii .EQ. HGDEM3D) THEN
                     Tc(lll+lt,iii) = 0.05
                 ELSE
                     Tc(lll+lt,iii) = 0.0
                 END IF
!               C.

                 Tc(l+lt,iii) = Tc(lll+lt,iii) 
                 Tc(ll+lt,iii) = Tc(lll+lt,iii)
!                 Tc(lll+lt,iii) = 1

             enddo

        enddo
        enddo

c.... Northern boundary along ~60N ..................................................

      kna1=12
      kne1=57

      ll = 0 
      lwe = indend(kna1-1)
      ll  = islab(kna1-1)
      do k = kna1,kne1 
        lw = lwe+1 
        lwe = indend(k) 
        ldown = lazc(lw) 
        lll = ll+ldown 
        llll = ll+ldown*2 

        do l = 1,ldown 
!            do iii=1,nbio-nsed 
            do iii=1,nchem

!               C. Fixed boundaries
                 IF (iii .EQ. HGDIS3D) THEN
                     Tc(iindex(68,k)+l,iii) = 0.2
                 ELSE IF (iii .EQ. MEHGDIS3D) THEN
                      Tc(iindex(68,k)+l,iii) = 0.03 
                 ELSE IF (iii .EQ. HGDEM3D) THEN
                      Tc(iindex(68,k)+l,iii) = 0.02
                 ELSE
                      Tc(iindex(68,k)+l,iii) = 0.0
                 END IF
!               C.

                Tc(iindex(67,k)+l,iii) = Tc(iindex(68,k)+l,iii)
                Tc(iindex(66,k)+l,iii) = Tc(iindex(68,k)+l,iii)
!                Tc(iindex(68,k)+l,iii) = 1
            enddo

!	    Tc( iindex( 66,k )+l, 9 ) = Tc( iindex( 68,k )+l, 9 )
!	    Tc( iindex( 66,k )+l, 1 ) = Tc( iindex( 68,k )+l, 1 )
!	    Tc( iindex( 66,k )+l, 2 ) = Tc( iindex( 68,k )+l, 2 )
!	    Tc( iindex( 66,k )+l, 3 ) = Tc( iindex( 68,k )+l, 3 )
!	    Tc( iindex( 66,k )+l, 4 ) = Tc( iindex( 68,k )+l, 4 )

!	    Tc( iindex( 67,k )+l, 9 ) = Tc( iindex( 68,k )+l, 9 )
!	    Tc( iindex( 67,k )+l, 1 ) = Tc( iindex( 68,k )+l, 1 )
!	    Tc( iindex( 67,k )+l, 2 ) = Tc( iindex( 68,k )+l, 2 )
!	    Tc( iindex( 67,k )+l, 3 ) = Tc( iindex( 68,k )+l, 3 )
!	    Tc( iindex( 67,k )+l, 4 ) = Tc( iindex( 68,k )+l, 4 )

!	    Tc( iindex( 68,k )+l, 9 ) = 0.4*0.0566 / XBC_FAC
!	    Tc( iindex( 68,k )+l, 1 ) = 0.4*0.566 / XBC_FAC
!	    Tc( iindex( 68,k )+l, 2 ) = 0.4*0.126 / XBC_FAC
!	    Tc( iindex( 68,k )+l, 3 ) = 0.4*0.126 / XBC_FAC
!	    Tc( iindex( 68,k )+l, 4 ) = 0.4*0.126 / XBC_FAC

        enddo 
       ll = islab(k) 
      enddo

      return 
      end 

c-----------------------------------------------------------------------
      subroutine kotief (dz,tz,nhor,ktot) 

c-----------------------------------------------------------------------
c     estimation of topography dependent arrays 
c-----------------------------------------------------------------------

      include 'C_model_inc'
c Attention!!! ngro must be recalculated depending on model domain geometry by
c              ngro = (max(m*(ilo*20+10),(kasor*8)))
      parameter (ngro=m*(ilo*40+10)) ! north_b
c______   north_b _________________

      parameter(khor1=khor+1) 
      integer tz(m,n) 
      dimension dz(ilo) 

       common /vecind/ indwet(khor),lb(n),le(n),indver(ndrei), 
     1 jwet(khor), 
     2 llw(ndrei),Dvoltj(m,ilo,n), 
     3 indi(ndrei),irbsor(kasor,2,2),nrbsor(2,2) 

       common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     *  indend(n),isornr(n),isorsr(n),islab(n) 

      idm = 1 
      nhor = 0 
      ktot = 0 
      do 6 k = 1,n 
      do 5 i = 1,m 
      it = tz(i,k) 
      itt = it 
      if (it.le.0) goto 5 
      do 2 j = 1,ilo 
      jj = j 
      if (it.le.dz(j)) goto 3 
    2 continue 
      it = dz(ilo) 
      tz(i,k) = it 
cc      IF(myid .EQ. 0)THEN
      write (7,603) tz(i,k),itt,i,k
  603 format (2x,'gitief, depth reduction to  dz(ilo),', 
     2'  tneu, talt, (i,k)/',2x,2i8,4x,2i4) 
cc      endif 
    3 iz = 0 
      j = jj 
      if (j.gt.1) iz = dz(j-1) 
      idep = it-iz 
      if(idep-idm+1) 45,44,4 
   44 tz(i,k) = tz(i,k)+1 
      it = tz(i,k) 
ccc      IF(myid .EQ. 0)THEN 
      write (6,602) tz(i,k),itt,i,k
ccccc      endif 
      goto 3 
   45 j = j-1 
      jj = max  (j,0) 
      tz(i,k) = iz 
      it = iz 
cc      IF(myid .EQ. 0)THEN 
      write (7,602) tz(i,k),itt,i,k 
  602 format (2x,'gitief, depth changed, tneu, talt, (i,k)/', 
     22x,2i8,4x,2i3) 
cc      endif 
      if (jj-1) 5,3,3 
    4 nhor = nhor+1 
      ktot = ktot+jj 
      ldep(nhor) = idep 
      iwet(nhor) = i 
      lazc(nhor) = jj 
    5 continue 
      indend(k) = nhor 
      islab (k) = ktot 
    6 continue 
c     -----------------beginning and end- indicees for sor-slabs 
      lwe = 0 
      isor=3 
      do 9 k = 1,n 
*     if(k.gt.100) isor=2 
      lwa = lwe+1 
      lwe = indend(k) 
      lxan=-2 
      do 7 lw = lwa,lwe 
      i = iwet(lw) 
      lxan=lw 
c--------------------------------------------------------------------
c      Attention!!!!  this has to be adapted to new topography
c      exclude wet open boundary lines from sor iteration
c  N and W boundaries:inner First and Second lines to be excluded
c  S and E boundaries:inner First lines to be excluded
c--------------------------------------------------------------------
      if((k.le.17).and.(i.gt.150)) goto 7 ! has to be changed 
      if((i.le.67).and.(k.le.57 )) goto 7 ! has to be changed

      if (i.ge.isor) goto 88 
    7 continue 
   88 isornr(k) = lxan 
      lwp = lwe+lwa 
      lxen=-4 
      do 8 ll = lwa,lwe 
      lxen=lwp-ll 
      lw = lwp-ll 
      i = iwet(lw) 
c--------------------------------------------------------------
c      Attention!!!!  this has to be adapted to new topography
c      exclude wet open boundary lines from sor iteration
c      the same as above 
c--------------------------------------------------------------

      if((k.le.17).and.(i.gt.150)) goto 8  !has to be changed
      if((i.le.67).and.(k.le.57 )) goto 8  !has to be changed

      if (i.le.m-1) goto 99 

    8 continue 
   99 isorsr(k) = lxen 
    9 continue 
c

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
               indi(nwet+k) = iwet(lw) 
               indver(nwet+k) = k
               llw(nwet+k)=lw
            end do

            nwet = nwet+lazc(lw) 
         end do
      end do
cc      IF(myid .EQ. 0)THEN
      write (7,600) nhor,ktot 
      write (7,'(a2,15(1x,i5))')'nr',(isornr(j),j=1,n) 
      write (7,'(a2,15(1x,i5))')'sr',(isorsr(j),j=1,n) 
 600  format (/1x,'number of wet grid points/cells (2d,3d) /',2i8) 
cc      endif
        k2=0
        do k=1,n
        do i=isornr(k),isorsr(k)
        k2=k2+1
        enddo
        enddo

cc      IF(myid .EQ. 0)THEN ^M
      write (7,*)'//////' 
      write (7,*)' KSOR  ',k2 
      write (*,*)' KSOR  ',k2 
cc      endif ^M
      return 
      end

c-----------------------------------------------------------------------
      subroutine gatmit (deltat) 
c-----------------------------------------------------------------------
c     time step dependent and grid dependent variables  
c----------------------------------------------------------------------- 

      include 'C_model_inc'
      parameter(ilop1=ilo+1) 
      parameter(khor1=khor+1) 
 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n) 
      common zac(khor),wobc(khor),stuvc(khor)
      common fricu(khor),cxc(khor) 
      common cyc(khor),pac(khor),txc(khor),tyc(khor) 
      common stpc(ndrei),sac(ndrei),tec(ndrei) 
      common pres(ilo),wc(ndrei),fricv(khor) 

      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth, 
     *  dlvo(m),dlvu(m),gh,rdt,dt2,r4,dtrdln(m),dtdlr
    
      common /gitter2/ fianf,yambdanf,dphi,dlambda,rearth

      common /num/ dc(ilo),av(ilo),ad(ilo),dh(ilo),pd(ilo), 
     *  prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),dd(ilo), 
     *  qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 

c-----------------------------------------------------------------------
c     x,y e.g. phi lam coordinates in grad
c-----------------------------------------------------------------------

      common /xy_coor/  xt(m), yt(n),cosxt(m)    ! coordinates of T-points
      common /xy_cooru/ xv(m), yu(n)           ! coordinates of uv-points
      common /part_dens/ part_d(m,n,ilo)
      common /begin_xy/ phigrad,yamgrad

      g = 9.81   ! defined in program main
      rearth=6350000.

*---------------------------------------------------------------------
*		set XY grid for your model domain           
*         print out in control file #7:  pathpc//pathout//control
*---------------------------------------------------------------------

c---- set  first point  i=1, j=1 --------------
c     northern boundary  phi= 65 grad 56. min  N 
       f = 65.+59./60.  
c     western boundary  lam= -4 grad 10. min (W)
       y = -4.0-10./60.  ! 


c---- set  XY step im minutes
       phimin= 6. ! phi
       yammin=10. ! lam
*--------------------------------------------------------------

               fianf = f 
             yambdanf= y

*		XY step in degree 
             phigrad  = phimin/60.
             yamgrad  = yammin/60.

          do i=1,m
          xt(i) = fianf-real(i-1)*phigrad  !phi = 1 min
         cosxt(i)=cos(xt(i)*rad)
          enddo

        do j=1,n
        yt(j) = yambdanf+real(j-1)*yamgrad            !lam = 2 min
        enddo

        do j=1,n
        yu(j) = yambdanf+real(j-1+0.5)*yamgrad            !lam = 2 min
        enddo

        do i=1,m
          xv(i) = fianf-real(i-1+0.5)*phigrad  !phi = 1 min
        enddo

c     gezeit - und gitterkonstanten 
      dt = deltat 
      dth = dt/3600. 
      rdt = r*dt 
      r4 = 0.25*r 
      dt2 = dt*0.5 
      gh = 0.5*g 
      gdt = g*dt 
      pi = 4.*atan(1.) 
      rad = pi/180. 
      aeqgrd = 2.*pi*6.37104e6/360.   ! 111200 [m/grad]

          do i=1,m
         cosxt(i)=cos(xt(i)*rad)
          enddo


      dl     = phigrad *aeqgrd   !step phi in [m] 
      dlam   = yamgrad *aeqgrd   !step lam in [m]
      dphi   =phigrad * rad   !step phi in rad
      dlambda=yamgrad * rad   !step lam in rad


      dlr = 1./dl 
      dtdlr = dt*dlr 
      dlrh = 0.5*dlr 
      omto = 4.*pi/86164. 
 
 
      fiu = fianf * rad + dphi 
      fiv = fiu-.5*dphi 
        write(7,*)' j,     xv(j),     dlvo(j),     dlvu(j)'

      do j=1,m 
        fiu = fiu-dphi 		        ! 1=fianf,2,3,4   u-point
        fivo = fiv 					! u(i-1)=fianf-1/2 North
        fiv = fiv-dphi                ! u (i) =fianf+1/2 South
        cosfiu = cos(fiu) 
        dlvu(j) = cos(fiv)/cosfiu  !u (i)
        dlvo(j) = cos(fivo)/cosfiu !u(i-1)
        write(7,*)j,xv(j),dlvo(j),dlvu(j)
        dln(j) = dlam*cosfiu            ! at T(j) point
        rdln(j) = 1./dln(j) 
        dtrdln(j) = dt*rdln(j) 
        dlnv = dlam*cosfiu*dlvu(j) 
      enddo 

      return 
      end 

         subroutine IO_daily_physics(nout,mjar,lmon,nday,ivier,ip,nnn) 
ccc !nnn-->       write=1 or read=2^M
      include 'C_model_inc'

      common /daily_phy/ zmit(khor),
     &  umit(ndrei),vmit(ndrei),wcmit(ndrei),
     &  acmit(ndrei),szmit(ndrei),
     &  tcmit(ndrei),scmit(ndrei),
     &  frimit(m,n),hismit(m,n),hisrmit(m,n),
     &  tismit(m,n),uimit(m,n),vimit(m,n),
     &  qois(m,n),qiis(m,n)

      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)

         if(nnn.eq.1)then
          write (nout) mjar,lmon,nday,ivier,ip
          write (nout) zmit
          write (nout) umit
          write (nout) vmit
          write (nout) wcmit
          write (nout) acmit
          write (nout) szmit
            write (nout) tcmit
            write (nout) scmit
          write(nout) frimit
          write(nout) hismit
          write(nout) hisrmit
          write(nout) tismit
          write(nout) uimit
          write(nout) vimit
            write(nout) fqgmit
            write(nout) fqrmit
            write(nout) fqsmit
            write(nout) fqlmit
          write(nout) qois
          write(nout) qiis
         else !nnn=2

          read (nout) mjarcc ,lmoncc,ndaycc ,iviercc,ipcc
!          write (*,*) mjarcc ,lmoncc,ndaycc,iviercc,ipcc

          read (nout) zmit
          read (nout) umit
          read (nout) vmit
          read (nout) wcmit
          read (nout) acmit
          read (nout) szmit
            read (nout) tcmit
            read (nout) scmit
          read(nout) frimit
          read(nout) hismit
          read(nout) hisrmit
          read(nout) tismit
          read(nout) uimit
          read(nout) vimit
            read(nout) fqgmit
            read(nout) fqrmit
            read(nout) fqsmit
            read(nout) fqlmit
          read(nout) qois
          read(nout) qiis
        endif      

        return
        end

C--------------------------------------------------------------------------
       subroutine change_units(f,n_choice)
c-------------------------------------------------------------------------
      include 'C_model_inc'
      common /bi/  BioC(422),BioOM(4,m,ilo),GI(2,6),
     &               REDF(20),BWup(ilo),
     &               sedy,QNH4,QBup

         common/Tin_4d/Tcin(ndrei,ninbio),Tflu(ndrei,nflu)

        common /bio_out/bio_out(ndrei,3:ninbio)
        Double Precision Dvoltj
        Integer biolag

      dimension f(ndrei,ninbio)
      if(n_choice.eq.1)then
*------change units
         do l=1,ndrei
         do ibio=3,7  !Ps, Pl, Zs, Zl, De^M
         bio_out(l,ibio)=f(l,ibio)  !*REDF(6)  !input in   mgC
                                    !  /m**3 to [mmolC/m**3]^M
         enddo
         do ibio=8,10 ! NH4, NO2, NO3^M
         bio_out(l,ibio)=f(l,ibio)*REDF(1)*REDF(6)  !input in mmolN
                                                    !/m**3 to [mmolC/m**3]^M
         enddo
         bio_out(l,11)  =f(l,11)  *REDF(2)*REDF(6)  !input in mmolP
                                                    !/m**3 to [mmolC/m**3]^M
         bio_out(l,12)  =f(l,12)  *REDF(3)*REDF(6)  !input in
                                                    !mmolSi/m**3 to [mmolC/m**3]^M
         bio_out(l,14)  =f(l,14)  *REDF(3)*REDF(6)  !input in
                                                    !mmolSi/m**3 to [mmolC/m**3]^M
         bio_out(l,13)  =f(l,13)  *REDF(7)          !input in
                                                 !milliliters O2/liter to  mmolO2/m**3^M
         bio_out(l,15)  =f(l,15)                    !bg^M
         bio_out(l,16)  =f(l,16)                    !sed1^M
         bio_out(l,17)  =f(l,17)                    !sed2^M
         enddo !l=1,nrei
      else
*------change units^M
         do l=1,ndrei
         do ibio=3,7  !Ps, Pl, Zs, Zl, De
         bio_out(l,ibio)=f(l,ibio)         !*REDF(16)  !output from
                                           ! [mmolC/m**3] to mgC /m**3^M
         enddo
        do ibio=8,10 ! NH4, NO2, NO3
         bio_out(l,ibio)=f(l,ibio)*REDF(11)*REDF(16)  !output from
                                                     ![mmolC/m**3] to mmolN /m**3^M
         enddo
         bio_out(l,11)  =f(l,11)  *REDF(12)*REDF(16)  !output from
                                                   !  [ mmolC/m**3] to mmolP /m**3^M
         bio_out(l,12)  =f(l,12)  *REDF(13)*REDF(16)  !output from
                                                      ![mmolC/m**3] to mmolSi/m**3^M
         bio_out(l,14)  =f(l,14)  *REDF(13)*REDF(16)  !output from
                                                      ! [mmolC/m**3] to mmolSi/m**3^M
         bio_out(l,13)  =f(l,13)  *REDF(17)           !output from
                                                      ! mmolO2/m**3 to mlO2/l^M
         bio_out(l,15)  =f(l,15)                    !bg^M
         bio_out(l,16)  =f(l,16)                    !sed1^M
         bio_out(l,17)  =f(l,17)                    !sed2^M
         enddo !l=1,nrei
      endif
      return
      end
c-----------------------------------------------------------------------
       subroutine uvrand (izet,iindex,llvrp,mm,nn,n3d,szahl)

c-----------------------------------------------------------------------
c      transports, boundary conditions
c-----------------------------------------------------------------------

      include 'C_model_inc'
      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)
      dimension trum(llvrp),srum(llvrp)
      dimension szahl(n3d)
cc     
      integer izet(m,n)
      dimension iindex(mm,nn)
      dimension f(m,n,ilo)
c      
      common/T_4d/ Tc(ndrei,nchem),Tc_mit(ndrei,nchem)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor),
     cindend(n),isornr(n),isorsr(n),islab(n)
c 
      integer nord,sued,ost,west
      integer*2 isor,ksor
c 
      common /intit/ itel(lrp),jc(m,3),ikk(kasor),isor(kasor),
     cksor(kasor),nofrez,lone,itende,nfreez
c 
      common /iterat/ omega(kasor),za(m,3),cx(m,3),cy(m,3),rhoz(m,3),
     ctcc(8)
c 
      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)
c 



ccc set u, v, w,avc, szahl boundary conditions
c... Western boundary  
c....................................................

        do i=156,164
        l   = iindex(i,16) !First line
        ll  = iindex(i,17) !Second line
        lll = iindex(i,18) !Third line

        lw  = izet(i,16)
                     lump=0
        if(lw.gt.0)lump = lazc(lw)
        do k=1,lump
          vc (l +k)  = vc (lll+k)
          vc (ll+k)  = vc (lll+k)
          uc (l +k)  = uc (ll +k)
 
          wc (l +k)  = wc (lll+k)
          wc (ll+k)  = wc (lll+k)
          avc(l +k)  = avc(lll+k)
          avc(ll+k)  = avc(lll+k)
          szahl(l +k)  = szahl(lll+k)
          szahl(ll+k)  = szahl(lll+k)

        enddo
      enddo

c.... Northern boundary
c..................................................

      kna1=12     
      kne1=57      
       do j=kna1,kne1 
       l=iindex(66,j)
      ll=iindex(67,j)
      lll=iindex(68,j)
      lw=izet(66,j)
       lump=0
       if(lw.gt.0) lump=lazc(lw)
 
      do k = 1,lump 
         uc(l+k)=uc(lll+k)
         uc(ll+k)=uc(lll+k)
          vc(l+k)=vc(ll+k)
          wc(l+k)=wc(lll+k)
          wc(ll+k)=wc(lll+k)
          avc(l+k)=avc(lll+k)
          avc(ll+k)=avc(lll+k)
          szahl(l+k)=szahl(lll+k)
          szahl(ll+k)=szahl(lll+k)
        enddo
      enddo

        return
        end

