
      subroutine stromTcin(zinc,vtmit,vtmic,avmax,szahl)

c-----------------------------------------------------------------------
c      transport of tracers, advection & turbuleny diffusion
c      TVD total variation diminishing scheme with superbee
c      Advection scheme (Barthel et al., 2012)
c-----------------------------------------------------------------------

!        IMPLICIT NONE

        USE CIN

        INCLUDE 'C_model_inc'

!        INTEGER TID, NTHREADS
!        INTEGER OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS, OMP_GET_WTIME

        REAL, DIMENSION(khor) :: ljumm
        REAL vol0
        REAL, DIMENSION(m)    :: dtdx
        REAL, DIMENSION(m)    :: dvo
        REAL, DIMENSION(m)    :: dvu
        INTEGER, DIMENSION(-3:3,-3:3,ndrei)        :: lww
        INTEGER, DIMENSION(0:khor,ilo)             :: llay



      parameter(ilop1=ilo+1)
      parameter(khor1=khor+1)

      common/T_4d/ Tc(ndrei,nchem),Tc_mit(ndrei,nchem) 
      common/Tin_4d/Tcin(ndrei,3:nbio),Tflu(ndrei,nflu),Tfis(ndrei,2)  

c parameters used for the TVD procedure:
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

      real*8 zinc
      dimension XYadv(ndrei,nchem),svert(ndrei,nchem)
      dimension s(ndrei,nchem),sn(ndrei,nchem),sd(ndrei,nchem)
      dimension sbet(ndrei,nchem),szahl(ndrei)
      dimension tsa(ndrei,nchem),tsb(ndrei,nchem),tsc(ndrei,nchem)
      dimension tsal(ndrei,nchem)
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
      dimension avd(ndrei,nchem) 
c      ---------------------- laenge von common /auxint/ = m*(ilo*6+3) 
      common /auxint/ lay(ilo,m,2),laz(m,2) 
c      new array
      integer jc5(m,n,ilo)
c       common /rivflow/ rivflux(m,n)

      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo), 
     * prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),d_(ilo), 
     * qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
      logical init

!      common /vec444/ init,ljumm(khor),vol0,dtdx(m),dvo(m),dvu(m),
!     1 lww(-3:3,-3:3,ndrei),llay(0:khor,ilo)

      common /vecind/ indwet(khor),lb(n),le(n),indver(ndrei),
     * jwet(khor),llw(ndrei),Dvoltj(m,ilo,n),
     * indi(ndrei),irbsor(kasor,2,2),nrbsor(2,2)


      common /julian/ julianday

      data init /.true./
      dimension d(ilo,m),Tc0(0:ndrei,nchem),dd(0:ndrei),r_dt_dd(ndrei)
      real*8 dd_dd(ndrei), dt_lr


!        print*, vtmit, vtmic, avmax


c--------  set model variables to be advected -----------
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
        dt_lr=dt*dlr
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
       do ibio = 3,nbio

!        print*, "STROM ",VNAMOUT( ibio )

         do ll = 1,ndrei
             Tc0(ll,ibio) = Tcin(ll,ibio)
         end do
       end do

       do ibio = 3,nbio
         Tc0(0,ibio) =0.       
       end do                  
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
! 5m/day = 0.00005787 m/s
! 1m/day = 0.000116 m/s
!        print*, "SINKING SHIP ",r_dt_dd(ll), ll
         SinkD_o=0.
         SinkD_o( DET )   = -0.00005787 * r_dt_dd(ll)!/2. !    5 m/day
         SinkD_o( CYA )   =  0.000001157* r_dt_dd(ll)!/2. ! -0.1 m/day upwelling
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
!        XTIME1 = OMP_GET_WTIME()
!        !$OMP PARALLEL PRIVATE (TID,NTHREADS)
!        TID = OMP_GET_THREAD_NUM()
!        NTHREADS = OMP_GET_NUM_THREADS()
!        !$OMP DO
      do ibio = 3,nbio

c west
      RSUL=0.0d+0
      do irsu=1,min0(jc5(i,j-1,k),1)
      DENUM=Tc0(ll,ibio)-Tc0(lww(0,-1,ll),ibio)
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.0d+0
      else
       if (uwt.ge.0.0d+0) NSIG=-1
       if (uwt.lt.0.0d+0) NSIG=1
       if ((NSIG.eq.-1 .and. jc5(i,j-2,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i,j+1,k).eq.0)) then
       PHIR=0.

       else
       RATIO=(Tc0(lww(0,NSIG,ll),ibio)-Tc0(lww(0,NSIG-1,ll),ibio))/DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(uwt)*(1.-COURUL)*D1
       RSUL=D2+0.5*(dprod((uwt+abs(uwt)),Tc0(lww(0,-1,ll),ibio))+
     &                dprod((uwt-abs(uwt)),Tc0(ll,ibio)))
       RSUL=RSUL*dble(dl)
       enddo

c east
      RSUR=0.0d+0
      do irsu=1,min0(jc5(i,j+1,k),1)
      DENUM=(Tc0(lww(0,1,ll),ibio)-Tc0(ll,ibio) )
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
       RATIO=(Tc0(lww(0,NSIG+1,ll),ibio)-Tc0(lww(0,NSIG,ll),ibio))
     & /DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(uot)*(1.-COURUR)*D1
       RSUR=D2+0.5*(dprod((uot+abs(uot)),Tc0(ll,ibio))+
     &                dprod((uot-abs(uot)),Tc0(lww(0,1,ll),ibio)))
       RSUR=RSUR*dble(dl)
       enddo
c north
      RSVO=0.0d+0
      do irsu=1,min0(jc5(i-1,j,k),1)
      DENUM=(Tc0(lww(-1,0,ll),ibio)-Tc0(ll,ibio) )
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
       if (vnt.ge.0.0d+0) NSIG=1
       if (vnt.lt.0.0d+0) NSIG=-1
       if ((NSIG.eq.-1 .and. jc5(i-2,j,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i+1,j,k).eq.0)) then
       PHIR=0.
       else
       RATIO=(Tc0(lww(NSIG-1,0,ll),ibio)-Tc0(lww(NSIG,0,ll),ibio))
     & /DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(vnt)*(1.-COURVO)*D1
       RSVO=D2+0.5*(dprod((vnt+abs(vnt)),Tc0(ll,ibio))+
     &                dprod((vnt-abs(vnt)),Tc0(lww(-1,0,ll),ibio)))
       RSVO=RSVO*dvo(i)
       enddo
c south
      RSVU=0.0d+0
       do irsu=1,min0(jc5(i+1,j,k),1)
      DENUM=(Tc0(ll,ibio)-Tc0(lww(1,0,ll),ibio))
      if (abs(DENUM).lt.EPSIL) THEN
       PHIR=0.
      else
       if (vst.ge.0.0d+0) NSIG=1
       if (vst.lt.0.0d+0) NSIG=-1
      if ((NSIG.eq.-1 .and. jc5(i-1,j,k).eq.0) .or. (NSIG.eq.1 .and.
     & jc5(i+2,j,k).eq.0)) then

       PHIR=0.
       else
       RATIO=(Tc0(lww(NSIG,0,ll),ibio)-Tc0(lww(NSIG+1,0,ll),ibio))/DENUM
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
       endif
      endif
       D1=PHIR*DENUM
       D2=0.5*abs(vst)*(1.-COURVU)*D1
       RSVU=D2+0.5*(dprod((vst-abs(vst)),Tc0(ll,ibio))+
     &                dprod((vst+abs(vst)),Tc0(lww(1,0,ll),ibio)))
       RSVU=RSVU*dvu(i)
       enddo

*-------- the horizontal advection part to the new value -------------------
      XYadv(ll,ibio)=(RSUL-RSUR+RSVU-RSVO)*dt/Dvoltj(i,k,j)

! additional sinking velosity

!cc      if (ibio.eq.4.or.ibio.eq.8) then  ! Sinking for POM
         D_o2=SinkD_o(ibio)      ! in the "upper" layer
         D_u2=SinkD_u(ibio)      ! in the "bottom" layer	
         if(k.eq.1) D_o2=0.
         if(k.eq.lazc(llw(ll))) D_u2=0.   
c      TVD scheme
      SADV=0.0d+0
c      from above
      RSWO=0.0d+0
       if (k.gt.1) then
       DDI=(dd(ll)+dd(llo))*0.5
       COURWO=abs(wo)*dt/DDI
      DENUM=( Tc0(llo,ibio)-Tc0(ll,ibio))/DDI
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
        RATIO=( Tc0(llo+NSIG,ibio)-Tc0(ll+NSIG,ibio))/(DDIZ*DENUM)
        PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
        PHIR=max(0.,PHIR1)
        endif
      endif
      D1=PHIR*DENUM*DDI
      D2=0.5*abs(wo)*(1.-COURWO)*D1
      RSWO=D2+0.5*(dprod((wo+abs(wo)),Tc0(ll,ibio))+
     &               dprod((wo-abs(wo)),Tc0(llo,ibio)))
      endif
c      from below
      RSWU=0.0d+0
       if (k.lt.lazc(llw(ll))) then !!!!
       DDI=(dd(ll)+dd(llu))*0.5 
      COURWU=abs(wu)*(dt)/DDI
      DENUM=( Tc0(ll,ibio)-Tc0(llu,ibio))/DDI
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
       RATIO=( Tc0(ll+NSIG,ibio)-Tc0(llu+NSIG,ibio))/(DDIZ*DENUM)
       PHIR1=max(min(1.,2.*RATIO),min(2.,RATIO))
       PHIR=max(0.,PHIR1)
      endif
      endif
      D1=PHIR*DENUM*DDI
      D2=0.5*abs(wu)*(1.-COURWU)*D1
      RSWU=D2+0.5*(dprod((wu+abs(wu)),Tc0(llu,ibio))+
     &               dprod((wu-abs(wu)),Tc0(ll,ibio)))
      endif

c vertical advection contribution

      SADV=(RSWU-RSWO)*r_dt_dd(ll)


cc----sinking vel is added by upstream---  NEW
        w_up   = max(0.,D_o2*Tc0(ll,ibio))  +min(0.,D_o2*Tc0(llo,ibio))
        w_down = max(0.,D_u2*Tc0(llu,ibio)) +min(0.,D_u2*Tc0(ll,ibio))
cc---------------------------------------
cic       w_up   = max(0.,0.01*Tc0(ll,ibio))  +min(0.,0.01*Tc0(llo,ibio))
cc       w_down = max(0.,0.01*Tc0(llu,ibio)) +min(0.,0.01*Tc0(ll,ibio))
      svert(ll,ibio)= SADV
     &        -w_up+w_down             ! ----sinking---
     
      if (k.eq.1) then
      svert(ll,ibio) = svert(ll,ibio) - zinc(i,j)*Tc0(ll,ibio)
     & /dd(ll)
      endif ! if (k.eq.1) then

      end do       ! ib = 1, nbio  
!        !$OMP END DO
!        XTIME = OMP_GET_WTIME() - XTIME1
!        print*, TID, XTIME, 'time', XTIME1
!        !$OMP END PARALLEL
      endif !if (ljumm(llw(ll)).eq.1) then

             end if !if (j.lt.n-2 .and. i.gt.5 .and. i.lt.m-2) then  
         end do  !i loop
      end if
c#ifdef MPI
      end do      !j loop

c#endif

222   format(1x,3(e9.3,2x))
c 
c -----------vertikal diffusion ---------------- 
c 
        DO ibio = 3,nbio
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
                    salref = Tc0(ll,ibio)
                    s (ll,ibio)  = Tc0(ll,ibio)
                    sn(ll,ibio)  = Tc0(ll,ibio)
                    stfak=szahl(ll)
!HHHHHHHHHHH avc - viscosity, avd - diffusivuty
        avd(ll,ibio) = amin1(stfak*avc(ll),admax)
        if(avd(ll,ibio).le.0.) avd(ll,ibio)=admax
              avd(ll,ibio)=max(avd(ll,ibio),avmin)
cccc         avd(ll, ibio)=0.
 
             end if
         end do
      end if
c#ifdef MPI
      end do
c#endif
      END DO

      DO ibio = 3,nbio
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
         tsa(ll,ibio) = dprod(dtdt,avd(ll,ibio))  
     1    /(dprod(dd(ll-1),dd(ll))+dd_dd(ll)) 
         tsc(ll,ibio) = dprod(dtdt,avd(ll+1,ibio)) 
     1    /(dprod(dd(ll+1),dd(ll))+dd_dd(ll)) 
         tsb(ll,ibio) = 1.d0+vtmit*tsa(ll,ibio)+vtmit*tsc(ll,ibio) 
         tsbb   = 1.d0-vtmic*tsa(ll,ibio)-vtmic*tsc(ll,ibio) 
         sd(ll,ibio) = vtmic*tsa(ll,ibio)*s(ll-1,ibio)+tsbb
     1    *s(ll,ibio)+vtmic*tsc(ll,ibio)*s(ll+1,ibio) 
                end if
             end if
         end do
      end if
c#ifdef MPI
      end do
c#endif
      END DO
c
c --------- ende vertikal diffusion 
c
!             do ibio=ibio0,nbio
        DO ibio = 3,nbio
      do j = 2,nz 
      do lw = lb(j),le(j)
         nwet = indwet(lw)
         ldown = lazc(lw) 
         if (ljumm(lw).eq.1) then
                ld  = ldown-1 
                lein= min0(1,ld) 
                if (lein.eq.1) then

      tsa(nwet+1,ibio) = 0.d0 
      tsc(nwet+1,ibio) = dprod(dtdt,avd(nwet+2,ibio)) 
     1 /(dprod(dd(nwet+2),dd(nwet+1))+dd_dd(nwet+1) ) 
      tsb(nwet+1,ibio) = 1.d0+vtmit*tsc(nwet+1,ibio) 
      tsbb   = 1.d0-vtmic*tsc(nwet+1,ibio) 
      sd(nwet+1,ibio) = tsbb*s(nwet+1,ibio)+vtmic*tsc(nwet+1,ibio)
     1 *s(nwet+2,ibio) 

      tsa(nwet+ldown,ibio) = dprod(dtdt,avd(nwet+ldown,ibio)) 
     1 /(dprod(dd(nwet+ld),dd(nwet+ldown))+dd_dd(nwet+ldown) ) 
      tsc(nwet+ldown,ibio) = 0.d0 
      tsb(nwet+ldown,ibio) = 1.d0+vtmit*tsa(nwet+ldown,ibio) 
      tsbb        = 1.d0-vtmic*tsa(nwet+ldown,ibio) 
      sd(nwet+ldown,ibio)  = vtmic*tsa(nwet+ldown,ibio)*s(nwet+ld,ibio)
     1 +tsbb*s(nwet+ldown,ibio) 

                end if
         end if
      end do
      end do
             END DO

        DO ibio = 3,nbio
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
                    tsal(nwet+lob,ibio) =tsb(nwet+lob,ibio) ! lob = 1 (or 2)
                    sbet(nwet+lob,ibio) = sd(nwet+lob,ibio) 
             end if
         end if
      end do
c#ifdef MPI
      end do
c#endif
        END DO

        DO ibio = 3,nbio
      do k = 1,ilo
!CDIR VECTOR NODEP
         do lh = 1,llay(0,k)
             ll = llay(lh,k)
                        tsal(ll+1,ibio) = tsb(ll+1,ibio)-(vtmit
     1                   *tsa(ll+1,ibio)/tsal(ll,ibio))*vtmit
     2                   *tsc(ll,ibio) 
                        sbet(ll+1,ibio) = sd(ll+1,ibio)+(vtmit
     1                   *tsa(ll+1,ibio)
     1                   /tsal(ll,ibio))*sbet(ll,ibio) 
         end do
      end do
        END DO

        DO ibio = 3,nbio
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
                    sn(nwet+ldown,ibio) = sbet(nwet+ldown,ibio)
     1               /tsal(nwet+ldown,ibio) 
             end if
         end if
      end do
c#ifdef MPI
      end do
c#endif
          END DO
        DO ibio = 3,nbio
      do k = ilo,1,-1
!CDIR VECTOR NODEP
         do lh = 1,llay(0,k)
             ll = llay(lh,k)
                sn(ll,ibio) = (sbet(ll,ibio)+vtmit
     1            *tsc(ll,ibio)*sn(ll+1,ibio))
     2            /tsal(ll,ibio) 
         end do
      end do
          END DO

        DO ibio = 3,nbio
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


        IF( XYadv(ll,ibio) .GT.  1. ) XYadv(ll,ibio) = 1E-3
        IF( svert(ll,ibio) .GT.  1. ) svert(ll,ibio) = 1E-3
        IF( XYadv(ll,ibio) .LT. -1. ) XYadv(ll,ibio) = -1E-3
        IF( svert(ll,ibio) .LT. -1. ) svert(ll,ibio) = -1E-3

cc        Tcin(ll,ibio)= XYadv(ll,ibio)+svert(ll,ibio)
        Tcin(ll,ibio) =  
     &     sn(ll,ibio)
     &   + XYadv(ll,ibio)
     &   + svert(ll,ibio)

!      IF( Tcin(ll,ibio) .GT. 10E20 .OR. Tcin(ll,ibio) .LT. -10E20 ) THEN
!        print*, ll, ibio, VNAMOUT( ibio ), sn(ll,ibio),XYadv(ll,ibio),svert(ll,ibio),Tcin(ll,ibio)
!      END IF
                end if
             end do
         end if
c#ifdef MPI
      end do
c#endif
      END DO
      return 
      end subroutine stromTCin
