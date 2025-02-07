#! /bin/csh -f

# Mistral job submission SLURM
#
#SBATCH --partition=compute
##SBATCH --dependency=afterok:34015139

#SBATCH --account=gg0305
#SBATCH --time=08:00:00
#SBATCH --mail-user=johannes.bieser@hereon.de
#SBATCH --output=/home/g/g260102/MECOSMO/model/v3.1/MECOSMO.slurm-%j.out
#SBATCH --error=/home/g/g260102/MECOSMO/model/v3.1/ERR_MECOSMO.slurm-%j.out
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --sockets-per-node=2
#SBATCH --cores-per-socket=12

#SBATCH --job-name=M_krbo15
#####################################################################################
set YEAR = 2013

echo $LD_LIBRARY_PATH
##setenv LD_LIBRARY_PATH "/opt/nco/pgi/openmpi/ib/lib:/opt/netcdf/4.2.1.1/pgi/openmpi/ib/lib:/opt/openmpi/pgi/ib/lib:/opt/pgi/linux86-64/13.2/libso:/storage/M3HOME/lib"
setenv OMP_NUM_THREADS 2
setenv OMP_STACKSIZE 200M
setenv KMP_STACKSIZE 200M
unset OMP_STACKSIZE
unset KMP_STACKSIZE
setenv KMP_BLOCKTIME 200 #60s #200ms
setenv KMP_AFFINITY physical
setenv KMP_LIBRARY turnaround #throughput #turnaround serial
setenv KMP_SETTINGS 1 #0
setenv KMP_VERSION off #on

unlimit stacksize
limit
###################################
#
# Run script for ECOSMO_Hg model
#
###################################

#setenv LD_LIBRARY_PATH /data/M3HOME/lib
setenv PROMPTFLAG N	#IOAPI USER PROMPT

setenv MECOSMO   ${HOME}/MECOSMO/model/v3.1
#setenv DATA      /work/gg0305/g260095/ecosmo/north_b/f_out_1990
setenv DATA	 /work/gg0305/g260102/MECOSMO/data/hourly_new/
setenv PATHTOP   input/geo_input
setenv PATHOUT 	 /work/gg0305/g260102/MECOSMO/res
setenv PATHPHYIN ${DATA}/
setenv PATHBIOIN ${DATA}/ 
#setenv PATHNC_IN ${MECOSMO}/2004
#setenv PATHNC_OUT $PATHOUT/2004res

#set IYEARS_START = 90	#04=2004
#et YEAR = 2002

#set YEAR = 2004
#@ YEAR = $YEAR + $IYEARS_START
set DAY = 001
setenv SDATE ${YEAR}${DAY}
setenv EYEAR $YEAR #SINGLE YEAR RUN

setenv PPP t01
setenv POU m01				#output character mercury daily means
setenv BIOFILE ${PATHBIOIN}/${PPP}	#initial input filename
setenv PHYFILE ${PATHPHYIN}/${PPP}	#initial input filename
setenv FLUFILE ${PATHPHYIN}/${PPP}
setenv GRIDINF ${MECOSMO}/${PPP}grdinfo
setenv CONTROL ${PPP}control
setenv TOPOFILE ${MECOSMO}/topo_nsm_mod		#input bathymetry
setenv OUTFILE  ${POU}000000.nc
setenv CHEMOUT ${PATHOUT}/chemtest.nc

#profile IC
#etenv CHEM_IC ${MECOSMO}/ECOSMO_IC_low

#p = serial
#x = parallel chem+bio
#xf = full paralel 

set restart = 1

#fix 1
#  fixed river loads
#  high kmb
#  modulate rates in e-function only using remin rate
#  all processes active

setenv CHEMOUT2D /work/gg0305/g260102/MECOSMO/res/CHEM2D_
setenv CHEMOUT3D /work/gg0305/g260102/MECOSMO/res/CHEM3D_

#restart files based on Kuss et al., 2017 observations
#setenv CHEMOUT2D ${PATHOUT}/IC_complete_based_on_Kuss_et_al_2D_
#setenv CHEMOUT3D ${PATHOUT}/IC_complete_based_on_Kuss_et_al_3D_


if ( $restart == 1 ) then
  set IC_YEAR = $YEAR
  @ IC_YEAR = $IC_YEAR - 1
###  set IC_YEAR = 2016
  setenv CHEM_IC ${CHEMOUT3D}${IC_YEAR}
  setenv CHEM_IC2D ${CHEMOUT2D}${IC_YEAR}
else
#Fixed basic IC conditions (x5 multiplier for low field in cout.f
  setenv CHEM_IC MECOSMO_IC_newlow
  setenv CHEM_IC2D MECOSMO_IC2D_newlow
  setenv BIO_IC  MECOSMO_IC_newlow

  setenv CHEM_IC IC3D_Py
  setenv CHEM_IC2D IC2D_Py
  setenv BIO_IC IC3D_Py
endif

################

setenv CHEMOUT2D /scratch/g/g260102/CHEM2D_
setenv CHEMOUT3D /scratch/g/g260102/CHEM3D_

setenv NSRIVERS ${MECOSMO}/ns.river.dat.2004
setenv BSRIVERS ${MECOSMO}/bs.river.dat.2011

setenv IPOLMAT ${MECOSMO}/matrix.9
setenv TESTOUT ${MECOSMO}/testout.nc

setenv FBBACKUP lo1backup

setenv INIFILEK re1040000.nc	#inifiles


setenv ITAGA 1
setenv MONA  1
setenv ITAGE 31
setenv MONE 12
setenv IYEARS_START 04	#04=2004


rm ${CHEMOUT2D}${YEAR}
rm ${CHEMOUT3D}${YEAR}

setenv ATMPATH  /work/gg0305/g260102/MECOSMO/data
setenv CONCFILE ${ATMPATH}/Hg_CCLM_GEM/CONC2D.CD72_Hg_CCLM_GEM_yamo_cb05tump_ae6_aq_
setenv DRYDEP   ${ATMPATH}/Hg_CCLM_GEM/DRYDEP.CD72_Hg_CCLM_GEM_yamo_cb05tump_ae6_aq_
setenv WETDEP   ${ATMPATH}/Hg_CCLM_GEM/WETDEP1.CD72_Hg_CCLM_GEM_yamo_cb05tump_ae6_aq_

#2009 meteo
setenv METCRO2D		${ATMPATH}/met/#/storage/M3HOME/data/lmmcip/CD72/	#MCPI 2D meteorological fields
setenv METCRO3D		${ATMPATH}/met/#/storage/M3HOME/data/lmmcip/CD72/	#MCPI 3D meteorological fields

#2006 meteo
setenv METCRO2D		/work/gg0305/g260102/LMMCIP4/output/CD72/
setenv METCRO3D		/work/gg0305/g260102/LMMCIP4/output/CD72/

#production
${MECOSMO}/polpar_main
exit(0)
