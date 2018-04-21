#!/bin/bash

Env=$1
Srclist=$3
source $Env

touch $Srclist
cat > $Srclist << EOF
$CODEROOT/libs/shr
$CODEROOT/libs/timing
$CODEROOT/atm/GAMIL2.0_128x60/src/physics/cam1
$CODEROOT/atm/GAMIL2.0_128x60/src/physics/cam1/echam_cu
$CODEROOT/atm/GAMIL2.0_128x60/src/dynamics/eul
$CODEROOT/atm/GAMIL2.0_128x60/src/control
$CODEROOT/atm/GAMIL2.0_128x60/src/couple/c_coupler
$CODEROOT/atm/GAMIL2.0_128x60/src/utils
$CODEROOT/atm/GAMIL2.0_128x60/src/advection/slt
$CODEROOT/atm/GAMIL2.0_128x60/src/ocnsice/dom
$CODEROOT/lnd/CLM2/src/main
$CODEROOT/lnd/CLM2/src/biogeophys
$CODEROOT/lnd/CLM2/src/mksrfdata
$CODEROOT/lnd/CLM2/src/ecosysdyn
$CODEROOT/lnd/CLM2/src/riverroute
$CODEROOT/sice/CSIM4
EOF

exit 0
