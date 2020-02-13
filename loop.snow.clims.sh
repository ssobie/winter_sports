#!/bin/bash

varname='swe'

#qsub -N "acc" -v gcm='ACCESS1-0',varname=$varname run.snow.clims.pbs
qsub -N "can" -v gcm='CanESM2',varname=$varname run.snow.clims.pbs
qsub -N "ccs" -v gcm='CCSM4',varname=$varname run.snow.clims.pbs
qsub -N "cnr" -v gcm='CNRM-CM5',varname=$varname run.snow.clims.pbs
qsub -N "csi" -v gcm='CSIRO-Mk3-6-0',varname=$varname run.snow.clims.pbs
qsub -N "gfd" -v gcm='GFDL-ESM2G',varname=$varname run.snow.clims.pbs
qsub -N "hcc" -v gcm='HadGEM2-CC',varname=$varname run.snow.clims.pbs
qsub -N "hes" -v gcm='HadGEM2-ES',varname=$varname run.snow.clims.pbs
qsub -N "inm" -v gcm='inmcm4',varname=$varname run.snow.clims.pbs
qsub -N "mir" -v gcm='MIROC5',varname=$varname run.snow.clims.pbs
qsub -N "mpi" -v gcm='MPI-ESM-LR',varname=$varname run.snow.clims.pbs
qsub -N "mri" -v gcm='MRI-CGCM3',varname=$varname run.snow.clims.pbs



