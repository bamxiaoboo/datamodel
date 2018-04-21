#!/bin/csh -f

# === Note by Li Ruizhe ===
# Required paramters:
#   NAMELIST_DST_DIR
#   DATA_DST_DIR
#   DATA_SRC_DIR
#
#   RUNTYPE
#   GRID
#   RUN_REFCASE
#   RUN_START_DATE
#   RAMP_CO2_START_YMD
#   DOUT_L_MSNAME
# =========================

cd $NAMELIST_DST_DIR

# --- set grid related variables
if ($GRID =~ 128x60*)     set params = (128  60 26 .true. 1200)
if ($GRID =~ 360x180*)     set params = (360  180 26 .true. 60)

# --- set grid resolution variables
set plon =  $params[1]; set plat = $params[2] ; set plev = $params[3]

#--- set namelist parameters dependent on grid resolution
set flxave = $params[4]


set bndtvo      = Ozone_CMIP5_ACC_SPARC_1850-2099_RCP2.6_T3M_O3.nc
set bndtvs      = amip_bc.gamil.128x060.1951-2010.nc
set absdata     = abs_ems_factors_fastvx.052001.nc
set bndtvaer    = Aerosol1850-2105RCP26gamil.nc
set datinit     = amip_ic.gamil.128x060.0011-01-01-00000.nc
set lnd_pft     = pft-physiology
set lnd_surface = surface-data.128x060.nc


link_data "$DATA_SRC_DIR/ozone/$bndtvo" "$DATA_DST_DIR/"
link_data "$DATA_SRC_DIR/boundary_condition/$bndtvs" "$DATA_DST_DIR/"
link_data "$DATA_SRC_DIR/radiation/$absdata" "$DATA_DST_DIR/"
link_data "$DATA_SRC_DIR/aerosol/$bndtvaer" "$DATA_DST_DIR/"
link_data "$DATA_SRC_DIR/initial_condition/$datinit" "$DATA_DST_DIR/"
link_data "$DATA_SRC_DIR/lnd/$lnd_pft" "$DATA_DST_DIR/"
link_data "$DATA_SRC_DIR/lnd/$lnd_surface" "$DATA_DST_DIR/"

set nsrest = 0
if ($RUN_TYPE == 'restart' || $RUN_TYPE == 'hybrid') then
    set nsrest = 1
cat >! $DATA_DST_DIR/lnd.$CASE_NAME.rpointer << EOF
./$ORIGINAL_CASE_NAME.clm2.r.$RUN_RESTART_DATE-$RUN_RESTART_SECOND
EOF
endif

cat >! ${MODEL_NAME}.stdin << EOF
&atmexp
! ------------------------------------------------------------------------------
! Experiment parameters
  caseid                = '$CASE_NAME'
  nsrest                = $nsrest
! ------------------------------------------------------------------------------
  aqua_planet           = .false.
! ------------------------------------------------------------------------------
  num_x_proc            = $num_x_proc
  num_y_proc            = $num_y_proc
! ------------------------------------------------------------------------------
! Input data parameters
  absems_data           = "$absdata"
  bndtvo                = "$bndtvo"
  ozncyc                = .false.
  bndtvs                = "$bndtvs" 
  sstcyc                = .false.
  ncdata                = "$datinit"
  bndtvaer              = "$bndtvaer"
! ------------------------------------------------------------------------------
! Output parameters
  mss_irt               = 0
  inithist              = "yearly"
  nhtfrq                = 0, -24, -6, -3
  mfilt                 = 1, 30, 120, 240
! ------------------------------------------------------------------------------
! Forcing parameters
  scenario_scon         = "CMIP5"
  scenario_ghg          = "CMIP5"
! ------------------------------------------------------------------------------
! Orbital parameters
  eccen                 = 0.016715
  obliq                 = 23.441
  mvelp                 = 102.7
!------------------------------------------------------------------------------
! daily output
  fincl2                = "TREFHTMN", "TREFHTMX", "TREFHT", "PRECT", "PS", "PSL", "CLDTOT", "PRECC", "PRECSC", "PRECSL", "LHFLX", "SHFLX", "FLDS", "FLUS", "FSDS", "FSUS", "FLUTOA", "T", "Q", "RELHUM", "OMEGA", "U", "V", "Z3", "WSPEED"
! ------------------------------------------------------------------------------
! 6-hour output
  fincl3                = "T:I", "Q:I", "U:I", "V:I", "PS:I", "PSL:I"
! ------------------------------------------------------------------------------
! 3-hour
! for COSP
  fincl4                = "PRECT", "TREFHT:I", "LHFLX", "SHFLX", "FLDS", "FLUS", "FLDSC", "FSDS", "FSUS", "FSUSC", "FSDSC", "PRECC", "PRECSC", "PRECSL", "PS:I", "CLDTOT", "SOLSD"
/

&clmexp
  finidat               = ""
  fpftcon               = "$lnd_pft"
  fsurdat               = "$lnd_surface"
/

EOF
