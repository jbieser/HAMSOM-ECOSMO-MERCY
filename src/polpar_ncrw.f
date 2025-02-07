!----------------------------------------------------------------------
        subroutine ncwriter(pathnc_out, resfilek)     
!----------------------------------------------------------------------
!  write daily data of parameters' 3D arrays to .nc file
! pathnc_out - path to folder where to write
! resfilek - name of file
!	IMPLICIT NONE
!----------------------------------------------------------------------
!====================================================================================
!         include '/usr/local/netcdf-4.1.0/include/netcdf.inc'
       include 'C_model_inc'
       parameter (nx=n-1,ny=m-1,ilop1=ilo+1,khor1=khor+1) 
!  netCDF file id
      integer                  :: ncid=55 !, public
      integer                  :: len
!  dimension ids
      integer                  :: lon_dim,lat_dim,z_dim,dep_dim
      integer                  :: time_dim
      integer, parameter       :: dim1=1,dim3=3 !,dim4=4
      integer                  :: dims(dim3)
cxx      integer                  :: time_len=NF_UNLIMITED
!  variable ids
      integer          :: lon_id,lat_id,z_id,dep_id,time_id !, private
      integer          :: start(4),edges(4)           !, private
      integer          :: tc1_id, tc2_id,tc3_id,tc4_id,tc5_id
      integer          :: tc6_id,
     & tc7_id,tc8_id,tc9_id,tc10_id,tc11_id 

       integer :: ircode, id1, id2, id3, vardim2(2), vardim3(3) 
ccc     &     m,n,ilo !array dimenions 
      integer :: idv1, idv2, idv3, idv4, idv5, idv6, 
     & idv7, idv8, 
     &   idv9, idv10, 
     &   idv11, idv12, idv13, idv14, idv15, idv16, 
     & idv17, idv18, 
     &    idv19, idv20, idv21, idv22,  idv14a

       dimension umit(ndrei),vmit(ndrei),tcmit(ndrei),acmit(ndrei) 
       dimension wcmit(ndrei),szmit(ndrei),dummy1(ndrei)                        
       dimension scmit(ndrei),szahl(ndrei)
       dimension frimit(m,n),hismit(m,n) 
       dimension hisrmit(m,n),uimit(m,n),vimit(m,n),tismit(m,n) 
       dimension sinput(m,n),sinmit(m,n)
       dimension qoi(m,n),qii(m,n),qois(m,n),qiis(m,n) 
       dimension evapmit(khor),precmit(khor)      
 
         dimension dudo3(m,n,ilo),dudo2(m,n) 
         dimension ar_lon(m), ar_lat(m),ar_z(ilo)

        character (len=7) ncout_fil 
cc        character pathnc_out*28, resfilek*12
cc          character pathnc_out*70,resfilek*12
            character pathnc_out*41,resfilek*12

        character       units,long_name
!	character(len=*), optional          :: units,long_name

      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     &            indend(n),isornr(n),isorsr(n),islab(n)
      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)
      common/T_4d/ Tc(ndrei,nchem),Tc_mit(ndrei,nchem)
      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo), 
     * prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),d_(ilo), 
     * qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
      common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,dtrdln(m),dtdlr


! arrays for coordinates
          ar_lon(1)=dln(1) 
          ar_lat(1)=dl
          do ik=2,m
          ar_lon(ik)=ar_lon(ik-1)+dln(ik)
          ar_lat(ik)=dl*real(ik)
         enddo

        do ik=1,ilo
        ar_z(ik)=dz(ik)
        enddo

! CREATING NEW .nc FILE
      write (*,*) pathnc_out,resfilek, ncid
      print*,'resfilek',resfilek
      write (*,*) pathnc_out//resfilek

!       not initialized!
      NOERR = -1 
      iret=nf_create (pathnc_out//resfilek,nf_clobber,ncid)
cc                    call Handle_Err(iret)
 
         print *,'opened netcdf file = ', resfilek,iret,ncid
         print *,iret, NF_NOERR

! DEFINITIONS OF DIMENSIONS !
!  define dimensions

      iret = nf_def_dim(ncid, 'x_dir', n,   id2)
          if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_def_dim(ncid, 'y_dir', m,   id1)
      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_def_dim(ncid, 'dep', ilo, id3)
      if(iret.ne.NF_NOERR) call Handle_Err(iret)


!      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
!      call Handle_Err(iret)

!  define variables

! 1D
c      iret = nf_def_var(ncid,'x_dir',NF_REAL,1,id2,lon_id)
c      if(iret.ne.NF_NOERR) call Handle_Err(iret)
c      iret = nf_def_var(ncid,'y_dir',NF_REAL,1,id1,lat_id)
c      if(iret.ne.NF_NOERR) call Handle_Err(iret)
c      iret = nf_def_var(ncid,'dep',NF_REAL,1,id3,dep_id)
c      if(iret.ne.NF_NOERR) call Handle_Err(iret)
 
!      dims(1) = time_dim
!      iret = nf_def_var(ncid,'time',NF_REAL,1,dims,time_id)
!      call Handle_Err(iret)

!  3(4)D (x,y,z,(t))

        vardim3(1) = id1
        vardim3(2) = id2
        vardim3(3) = id3

         iret = nf_def_var(ncid,'tc1',NF_REAL,3,vardim3,tc1_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
         iret = nf_def_var(ncid,'tc2',NF_REAL,3,vardim3,tc2_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
         iret = nf_def_var(ncid,'tc3',NF_REAL,3,vardim3,tc3_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
         iret = nf_def_var(ncid,'tc4',NF_REAL,3,vardim3,tc4_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
         iret = nf_def_var(ncid,'tc5',NF_REAL,3,vardim3,tc5_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
         iret = nf_def_var(ncid,'tc6',NF_REAL,3,vardim3,tc6_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
         iret = nf_def_var(ncid,'tc7',NF_REAL,3,vardim3,tc7_id)
c         if(iret.ne.NF_NOERR) call Handle_Err(iret)
            iret = nf_def_var(ncid,'tc8',NF_REAL,3,vardim3,tc8_id)
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)
           iret = nf_def_var(ncid,'tc9',NF_REAL,3,vardim3,tc9_id)
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)
          iret = nf_def_var(ncid,'tc10',NF_REAL,3,vardim3,tc10_id)
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)
           iret = nf_def_var(ncid,'tc11',NF_REAL,3,vardim3,tc11_id)
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)

!  assign attributes

!  coordinates
c        iret = nf_put_att_text(ncid,lon_id,'long_name',30,
c     &         'distance from western boundary')
c          iret = nf_put_att_text(ncid,lon_id,'units',2,'km')
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)  
c
c        iret = nf_put_att_text(ncid,lat_id,'long_name',31,
c     &         'distance from northern boundary')
c          iret = nf_put_att_text(ncid,lat_id,'units',2,'km')
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)


c        iret = nf_put_att_text(ncid,dep_id,'long_name',5,'depth')
c          iret = nf_put_att_text(ncid,dep_id,'units',5, 'meter')
c          if(iret.ne.NF_NOERR) call Handle_Err(iret)
!        iret = nf_put_att_text(ncid,time_id,'long_name',4,'time')
!          iret = nf_put_att_text(ncid,time_id,'units',7, 'seconds')
!                    call Handle_Err(iret)
!   variables
      iret = nf_put_att_text(ncid,tc1_id,'long_name',15,
     & 'Hg II dissolved')
      iret = nf_put_att_text(ncid,tc1_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc1_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc1_id,'_FillValue',NF_FLOAT,1,-999.)
cc          if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc2_id,'long_name',19,
     & 'Hg II phytoplankton')
      iret = nf_put_att_text(ncid,tc2_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc2_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc2_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc3_id,'long_name',17,
     & 'Hg II zooplankton')
      iret = nf_put_att_text(ncid,tc3_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc3_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc3_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc4_id,'long_name',17,
     & 'Hg II particulate')
      iret = nf_put_att_text(ncid,tc4_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc4_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc4_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc5_id,'long_name',24,
     & 'Methyl mercury dissolved')
      iret = nf_put_att_text(ncid,tc5_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc5_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc5_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)

      iret = nf_put_att_text(ncid,tc6_id,'long_name',22,
     & 'Methyl mercury phytopl')
      iret = nf_put_att_text(ncid,tc6_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc6_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc6_id,'_FillValue',NF_FLOAT,1,-999.)
c      if(iret.ne.NF_NOERR) call Handle_Err(iret)

      iret = nf_put_att_text(ncid,tc7_id,'long_name',20,
     & 'Methyl mercury zoopl')
      iret = nf_put_att_text(ncid,tc7_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc7_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc7_id,'_FillValue',NF_FLOAT,1,-999.)
c c     if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc8_id,'long_name',26,
     & 'Methyl mercury particulate')
      iret = nf_put_att_text(ncid,tc8_id,'units',4,'ng/l')
      iret=nf_put_att_real(ncid,tc8_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc8_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc9_id,'long_name',3,'hg0')
      iret = nf_put_att_text(ncid,tc9_id,'units',4,'pg/l')
      iret=nf_put_att_real(ncid,tc9_id,'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc9_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc10_id,'long_name',15,
     & 'Hg II sediments')
      iret = nf_put_att_text(ncid,tc10_id,'units',4,'pg/l')
      iret=nf_put_att_real(ncid,tc10_id,
     & 'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc10_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
      iret = nf_put_att_text(ncid,tc11_id,'long_name',24,
     & 'Methyl mercury sediments')
      iret = nf_put_att_text(ncid,tc11_id,'units',4,'pg/l')
        iret=nf_put_att_real(ncid,tc11_id,
     &'missing_value',NF_FLOAT,1,-999.)
      iret=nf_put_att_real(ncid,tc11_id,'_FillValue',NF_FLOAT,1,-999.)
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
!  global attributes
      iret=nf_put_att_text(ncid,NF_GLOBAL,'Title',47,
     & 'ECOSMO:North&Baltic Sea Mercury Schrum & Bieser')
      iret = nf_put_att_text(ncid,NF_GLOBAL,'Created_at',13,
     & 'UiB-GFI & HZG')
cc      if(iret.ne.NF_NOERR) call Handle_Err(iret)
! end of DEFINITIONS
          ircode= nf_enddef(ncid)
cc          if(ircode.ne.NF_NOERR) call Handle_Err(ircode)
! PROVIDE VARIABLES
! 1D
           iret=nf_put_var_real(ncid,lon_id,ar_lon)
cc           if(iret.ne.NF_NOERR) call Handle_Err(iret)
           iret=nf_put_var_real(ncid,lat_id,ar_lat)
cc           if(iret.ne.NF_NOERR) call Handle_Err(iret)
           iret=nf_put_var_real(ncid,dep_id,ar_z)
cc           if(iret.ne.NF_NOERR) call Handle_Err(iret)
! 3D
              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,1)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc1_id,dudo3)
c                   call Handle_Err(ircode)

             do itit=1,ndrei
                 dummy1(itit)= Tc_mit(itit,2)
             enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc2_id,dudo3)
c                    call Handle_Err(ircode)

              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,3)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc3_id,dudo3)
c                   call Handle_Err(ircode)

              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,4)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc4_id,dudo3)
c                   call Handle_Err(ircode)

              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,5)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc5_id,dudo3)
c                   call Handle_Err(ircode)

              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,6)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc6_id,dudo3)
c                   call Handle_Err(ircode)

              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,7)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)   
            ircode=nf_put_var_real(ncid,tc7_id,dudo3)
c                   call Handle_Err(ircode)

              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,8)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)
            ircode=nf_put_var_real(ncid,tc8_id,dudo3)
c                   call Handle_Err(ircode)
              do itit=1,ndrei
                  dummy1(itit)= Tc_mit(itit,9)
              enddo
            call deco1d3d (dudo3,dummy1,ndrei,999.)
            ircode=nf_put_var_real(ncid,tc9_id,dudo3)
c                   call Handle_Err(ircode)

!              do itit=1,ndrei
!                  dummy1(itit)= Tc_mit(itit,10)
!              enddo
!            call deco1d3d (dudo3,dummy1,ndrei,999.)
!            ircode=nf_put_var_real(ncid,tc10_id,dudo3)
c                   call Handle_Err(ircode)

!              do itit=1,ndrei
!                  dummy1(itit)= Tc_mit(itit,11)
!              enddo
!            call deco1d3d (dudo3,dummy1,ndrei,999.)
!            ircode=nf_put_var_real(ncid,tc11_id,dudo3)
c                   call Handle_Err(ircode)

! CLOSING .nc FILE

           ircode=nf_close(ncid)
           print *,'closed netcdf file = ', file_id


      return 
      end 

!======================================================================
        subroutine nc_read2(pathnc_out, resfilek)
!----------------------------------------------------------------------
! read  daily data of parameters' 3D arrays from .nc file
! pathnc_out - path to folder where to write
! resfilek - name of file
!======================================================================
!	IMPLICIT NONE
CC         include '/usr/include/netcdf.inc'
!        include '/usr/local/netcdf-4.1.0/include/netcdf.inc'
        include 'C_model_inc'
      parameter (nx=n-1,ny=m-1,ilop1=ilo+1,khor1=khor+1) 

!  netCDF file id
      integer                  :: ncid !, public
      integer                  :: len
!  dimension ids
      integer                  :: lon_dim,lat_dim,z_dim,dep_dim
      integer                  :: time_dim
      integer, parameter       :: dim1=1,dim3=3 !,dim4=4
      integer                  :: dims(dim3)
cxx      integer                  :: time_len=NF_UNLIMITED
!  variable ids
      integer          :: lon_id,lat_id,z_id,dep_id,time_id !, private
      integer          :: start(4),edges(4)           !, private
      integer          :: tc1_id, tc2_id, tc3_id, tc4_id, tc5_id

       integer :: ircode, id1, id2, id3, vardim2(2), vardim3(3) 
ccc     &     m,n,ilo !dimentions of arrays
      integer :: idv1, idv2, idv3, idv4, idv5, idv6, idv7, idv8, 
     &   idv9, idv10, 
     &   idv11, idv12, idv13, idv14, idv15, idv16, idv17, idv18, 
     &    idv19, idv20, idv21, idv22,  idv14a

       dimension umit(ndrei),vmit(ndrei),tcmit(ndrei),acmit(ndrei) 
       dimension wcmit(ndrei),szmit(ndrei)                        
       dimension scmit(ndrei),szahl(ndrei)
       dimension frimit(m,n),hismit(m,n) 
       dimension hisrmit(m,n),uimit(m,n),vimit(m,n),tismit(m,n) 
       dimension sinput(m,n),sinmit(m,n)
       dimension qoi(m,n),qii(m,n),qois(m,n),qiis(m,n) 
       dimension evapmit(khor),precmit(khor)      
 
         dimension dudo3(m,n,ilo),dudo2(m,n) 
         dimension ar_lon(m), ar_lat(m),ar_z(ilo)

        character (len=7) ncout_fil 
c--------------------Attention pathout change
!i!       character  pathnc_out*28, resfilek*12
cc        character pathnc_out*70,resfilek*12
          character pathnc_out*41,resfilek*12   

        character       units,long_name
!	character(len=*), optional          :: units,long_name

      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     &            indend(n),isornr(n),isorsr(n),islab(n)
      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)
      common/T_4d/ Tc(ndrei,nbio),Tc_mit(ndrei,nbio)
!cxxxxxxxxxxxxxxxxx   Attention, I change pprd into prd !!! (was not important)
      common/num/dz(ilo),av(ilo),ah(ilo),dh(ilo),pd(ilo), 
     * prd(ilo),pr2d(ilo),r2d(ilo),tkw(ilo),tau(ilo),d_(ilo), 
     * qa(ilop1),qbet(ilop1),qn(ilop1),rd(ilop1) 
       common /gitter/ dt,r,g,dl,dlr,dlrh,dln(m),rdln(m),rad,dth,
     * dlvo(m),dlvu(m),gh,rdt,dt2,r4,dtrdln(m),dtdlr


 !     dimension d3(m,n,ilo),d1(mwet) 
 
! arrays for coordinates
          ar_lon(1)=dln(1) 
          ar_lat(1)=dl
          do ik=2,m
          ar_lon(ik)=ar_lon(ik-1)+dln(ik)
          ar_lat(ik)=dl*real(ik)
         enddo

        do ik=1,ilo
        ar_z(ik)=dz(ik)
        enddo

! OPEN .nc FILE
!--------------------1234567890123456789012345678
  !     ircode=nf_create (pathnc_out//resfilek,nf_clobber,ncid)
      print*,'vor dem print'
      print*,  pathnc_out//resfilek
cc      pause 1
        ircode = NF_OPEN(pathnc_out//resfilek,NF_CLOBBER,ncid)
                    call Handle_Err(ircode)
         print *,'opened netcdf INITIAL Tc() file = ', resfilek

! DEFINITIONS OF DIMENSIONS !
!  define dimensions
      iret = nf_inq_dimid(ncid, 'x_dir',    id2)
      call Handle_Err(iret)
      iret = nf_inq_dimid(ncid, 'y_dir',   id1)
      call Handle_Err(iret)
      iret = nf_inq_dimid(ncid, 'dep',  id3)
      call Handle_Err(iret)
!  define variables

         iret = NF_INQ_VARID(ncid,'tc1', tc1_id)
          call Handle_Err(iret)
         iret = NF_INQ_VARID(ncid,'tc2', tc2_id)
          call Handle_Err(iret)
         iret = NF_INQ_VARID(ncid,'tc3', tc3_id)
          call Handle_Err(iret)
         iret = NF_INQ_VARID(ncid,'tc4', tc4_id)
          call Handle_Err(iret)
         iret = NF_INQ_VARID(ncid,'tc5', tc5_id)
          call Handle_Err(iret)
         iret = NF_INQ_VARID(ncid,'tc6', tc6_id)
          call Handle_Err(iret)
         iret = NF_INQ_VARID(ncid,'tc7', tc7_id)
          call Handle_Err(iret)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       ircode=NF_GET_VAR_REAL(ncid,idv8, dudo3)
!                    call Handle_Err(ircode)
!            call comp3d1d (dudo3,scmit,ndrei) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         iret=NF_GET_VAR_REAL(ncid,tc1_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,1)=tcmit(itit) 
            enddo

         iret=NF_GET_VAR_REAL(ncid,tc2_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,2)=tcmit(itit) 
            enddo

         iret=NF_GET_VAR_REAL(ncid,tc3_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,3)=tcmit(itit) 
            enddo

         iret=NF_GET_VAR_REAL(ncid,tc4_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,4)=tcmit(itit) 
            enddo

         iret=NF_GET_VAR_REAL(ncid,tc5_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,5)=tcmit(itit) 
!                write(*,*) Tc(itit,5),tcmit(itit), dudo3(itit,2,2) 
            enddo

         iret=NF_GET_VAR_REAL(ncid,tc6_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,6)=tcmit(itit) 
            enddo

         iret=NF_GET_VAR_REAL(ncid,tc7_id,dudo3)
                 call Handle_Err(ircode)
                 call comp3d1d (dudo3,tcmit,ndrei) 
            do itit=1,ndrei
                Tc(itit,7)=tcmit(itit) 
            enddo

! CLOSING .nc FILE
           ircode=nf_close(ncid)
 !          print *,'closed netcdf file = ', file_id
 !     pause 104

      return 
      end 

!======================================================================
!----------------------------------------------------------------------
        SUBROUTINE Handle_Err( iret )
!c---------------------------------------------------------------------
!c     check errors when writing/reading .NC files
!c---------------------------------------------------------------------

               integer iret
cc        include 'netcdf.inc' 
cc        include '/usr/include/netcdf.inc'
!          include '/usr/local/netcdf-4.1.0/include/netcdf.inc'

             if ( iret.ne.NF_NOERR ) then
             write(*,*) NF_STRERROR( iret )
             stop 'Stopped'
             end if

          return 
        END SUBROUTINE Handle_Err

!======================================================================
!======================================================================
!----------------------------------------------------------------------
        subroutine ncreader(pathnc_in, outfilek,szahl)
!c---------------------------------------------------------------------
! read daily outputs of hydrophys.data (for old Barents Sea version)
!C_--------------------------------------------------------------------

cc         include '/usr/include/netcdf.inc'
!           include '/usr/local/netcdf-4.1.0/include/netcdf.inc'

       include 'C_model_inc'
      parameter (nx=n-1,ny=m-1,ilop1=ilo+1,khor1=khor+1) 
      integer :: ircode, id1, id2, id3, vardim2(2), vardim3(3)  
       integer ::  idv1, idv2, idv3, idv4, idv5, idv6, idv7, idv8, 
     *   idv9, idv10, 
     *   idv11, idv12, idv13, idv14, idv15, idv16, idv17, idv18, 
     *    idv19, idv20, idv21, idv22,  idv14a
       dimension umit(ndrei),vmit(ndrei),tcmit(ndrei),acmit(ndrei) 
       dimension wcmit(ndrei),szmit(ndrei)                        
       dimension scmit(ndrei),szahl(ndrei)
       dimension frimit(m,n),hismit(m,n) 
       dimension hisrmit(m,n),uimit(m,n),vimit(m,n),tismit(m,n) 
       dimension sinput(m,n),sinmit(m,n)
       dimension qoi(m,n),qii(m,n),qois(m,n),qiis(m,n) 
       dimension evapmit(khor),precmit(khor)      
 
         dimension dudo3(m,n,ilo),dudo2(m,n) 
       character (len=7) ncout_fil 
!       character pathnc_in*52, outfilek*12
        character pathnc_in*41,outfilek*12


      common uc(ndrei),vc(ndrei),stc(ndrei),avc(ndrei),z(m,n)
      common zac(khor),wobc(khor),stuvc(khor),fricu(khor),cxc(khor)
      common cyc(khor),pac(khor),txc(khor),tyc(khor)
      common stpc(ndrei),sac(ndrei),tec(ndrei)
      common pres(ilo),wc(ndrei),fricv(khor)
      common /ind/ iwet(khor1),ldep(khor),lazc(khor), 
     *            indend(n),isornr(n),isorsr(n),islab(n)
      common /radis/fqgmit(m,n),fqrmit(khor),fqsmit(khor),fqlmit(khor)


        ircode = NF_OPEN(pathnc_in//outfilek,NF_CLOBBER,ncid)
!	ircode = NF_OPEN('lo1600101.nc',NF_CLOBBER,ncid)
                    call Handle_Err(ircode)
        print *,'opened netcdf file = ', outfilek
 
! DIMENSIONS ! 

        ircode = NF_INQ_DIMID(ncid,'y_dir',id1)
                    call Handle_Err(ircode) 
        ircode = NF_INQ_DIMID(ncid,'x_dir',id2)
                    call Handle_Err(ircode) 
        ircode = NF_INQ_DIMID(ncid,'z_dir',id3)
                    call Handle_Err(ircode) 
        print *,'y_dir = ', id1,'x_dir = ', id2,'z_dir = ', id3

! VARIABLES !
        ircode = NF_INQ_VARID(ncid,'zeta',idv1) 
                    call Handle_Err(ircode)
        ircode = NF_INQ_VARID(ncid,'u-ocean',idv2)
                    call Handle_Err(ircode)
        ircode = NF_INQ_VARID(ncid,'v-ocean',idv3) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'w-ocean',idv4) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'A-vert',idv5) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'s-zahl',idv6) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'temp',idv7) 
                    call Handle_Err(ircode)
          ircode = NF_INQ_VARID(ncid,'salt',idv8) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'ice-c',idv9) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'ice-h',idv10) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'ice-hr',idv11) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'ice-t',idv12) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'ice-u',idv13) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'ice-v',idv14) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'delsal',idv14a) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'fqg',idv15) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'fqr',idv16) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'fqs',idv17) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'fql',idv18) 
                    call Handle_Err(ircode)
        ircode = NF_INQ_VARID(ncid,'qois',idv19) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'qiis',idv20) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'evap',idv21) 
                    call Handle_Err(ircode)
         ircode = NF_INQ_VARID(ncid,'prec',idv22) 
                    call Handle_Err(ircode)
!        print *,' v y_dir = ', idv1, idv2, idv3, idv4, idv5, idv6, idv7, idv8, idv9, idv10,&
!		    idv11, idv12, idv13, idv14, idv15, idv16, idv17, idv18, &
!		    idv19, idv20, idv21, idv22,  idv14a
!VARIABLES VALUES !

       ircode=NF_GET_VAR_REAL(ncid,idv1, dudo2)
                    call Handle_Err(ircode)
                call comp2d1d (zac,dudo2) 

       ircode=NF_GET_VAR_REAL(ncid,idv2, dudo3)
                    call Handle_Err(ircode)
              call comp3d1d (dudo3,uc,ndrei) 

       ircode=NF_GET_VAR_REAL(ncid,idv3, dudo3)
                    call Handle_Err(ircode)
                 call comp3d1d (dudo3,vc,ndrei) 

       ircode=NF_GET_VAR_REAL(ncid,idv4, dudo3)
                    call Handle_Err(ircode)
            call comp3d1d (dudo3,wc,ndrei) 

            wc=wc/1000.

       ircode=NF_GET_VAR_REAL(ncid,idv5, dudo3)
                    call Handle_Err(ircode)
                 call comp3d1d (dudo3,avc,ndrei) 


           avc=avc/10000.

       ircode=NF_GET_VAR_REAL(ncid,idv6, dudo3)
                    call Handle_Err(ircode)
                 call comp3d1d (dudo3,szahl,ndrei) 
  !          goto 2335
!       ircode=NF_GET_VAR_REAL(ncid,idv7, dudo3)
!                    call Handle_Err(ircode)
               call comp3d1d (dudo3,tcmit,ndrei) 

       ircode=NF_GET_VAR_REAL(ncid,idv8, dudo3)
                    call Handle_Err(ircode)
            call comp3d1d (dudo3,scmit,ndrei) 
!----
       ircode=NF_GET_VAR_REAL(ncid,idv9, frimit )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv10,hismit )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv11,hisrmit )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv12,tismit )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv13,uimit)
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv14,vimit)
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv14a,sinmit )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv15,fqgmit)
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv16, dudo2)
                    call Handle_Err(ircode)
             call comp2d1d (fqrmit,dudo2) 

       ircode=NF_GET_VAR_REAL(ncid,idv17, dudo2)
                    call Handle_Err(ircode)
               call comp2d1d (fqsmit,dudo2) 

       ircode=NF_GET_VAR_REAL(ncid,idv18, dudo2)
                    call Handle_Err(ircode)
               call comp2d1d (fqlmit,dudo2) 

       ircode=NF_GET_VAR_REAL(ncid,idv19,qois )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv20,qiis )
                    call Handle_Err(ircode)

       ircode=NF_GET_VAR_REAL(ncid,idv21, dudo2)
                    call Handle_Err(ircode)
                  call comp2d1d (evapmit,dudo2) 

       ircode=NF_GET_VAR_REAL(ncid,idv22, dudo2)
                    call Handle_Err(ircode)
!	  do ii=1,m
!	    do jj=1,n
!		write (*,*) ii,jj,dudo2(ii,jj)
!	    enddo
!	  enddo
             call comp2d1d (precmit,dudo2) 
2335   continue

        ircode=nf_close(ncid)
!         print *,'closed netcdf file = ', file_id
!       pause 1234
!====================================================================================


      return 
      end 
