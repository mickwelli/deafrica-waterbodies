#!/bin/bash
#PBS -P u46
#PBS -q normal
#PBS -N aus_seasonal
#PBS -l walltime=04:00:00
#PBS -l mem=32GB
#PBS -l jobfs=1GB
#PBS -l ncpus=16
#PBS -l wd
#PBS -l storage=gdata/v10+gdata/r78+gdata/fk4
#PBS -M vanessa.newey@ga.gov.au
#PBS -m abe

module use /g/data/v10/public/modules/modulefiles/
module load dea
module load parallel

parallel --delay 2 --retries 3 --load 100%  --colsep ',' python par_create_wb_seasonal_means_aus.py config_seasonal.ini ::: {1..16},16
wait;
