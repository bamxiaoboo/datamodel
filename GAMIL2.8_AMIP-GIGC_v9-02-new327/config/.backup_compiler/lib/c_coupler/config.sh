#!/bin/csh -f

set start_date_num = `echo $START_DATE | sed -e 's/-//g'`  # remove "-"
set stop_date_num = `echo $STOP_DATE | sed -e 's/-//g'`  # remove "-"

cat > ${CASEROOT}/CCPL_dir/config/all/env_run.xml << EOF
<?xml version="1.0" ?>
<Time_setting
    case_name="$CASE_NAME"
    model_name="ideal_model_for_CCPL2"
    run_type="$RUN_TYPE"
    leap_year="$LEAP_YEAR"
    start_date="$start_date_num"
    start_second="$START_SECOND"
    reference_date="00010101"
    rest_freq_unit="$REST_FREQ_UNIT"
    rest_freq_count="$REST_FREQ_COUNT"
    stop_option="date"
    stop_date="$stop_date_num"
    stop_second="$STOP_SECOND"
    stop_n="1"
/>
EOF
