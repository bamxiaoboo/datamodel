#!/bin/bash

export env=${1}
export exedir=${2}
export makefile=${3}
export ntasks=${4}
export nthrds=${5}
export grid=${6}

source $env

cd $exedir/obj


touch .tmp
cat > .tmp << EOF; cmp -s .tmp preproc.h || cp -f .tmp preproc.h
#ifndef PREPROC_SET
#define PREPROC_SET
#define COUP_CAM
#define LSMLON  128
#define LSMLAT  60
#endif
 
EOF


touch .tmp
cat > .tmp << EOF; cmp -s .tmp params.h || cp -f .tmp params.h
#ifndef PARAMS_SET
#define PARAMS_SET
#define PCNST   1
#define PNATS   4
#define PLEV    26
#define PLEVR   26
#define PLON    128
#define PCOLS   16
#define PLAT    60
#define PTRM    42
#define PTRN    42
#define PTRK    42
#undef STAGGERED
#undef COUP_SOM
 
 
#endif
EOF

spmd="#define SPMD"

if [ $ntasks -gt 1 ]; then
	spmd="#define SPMD"
fi

cat > .tmp << EOF; cmp -s .tmp misc.h || mv -f .tmp misc.h
#ifndef MISC_SET
#define MISC_SET
#define DATA_SICE
#define DATA_OCN
#define ONLINE_LSM
$spmd
#undef PERGRO
#define NTASK ${ntasks}
#endif
EOF

gmake -j $GMAKE_J -f $makefile || exit 1

exit 0
