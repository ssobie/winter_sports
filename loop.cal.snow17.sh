#!/bin/bash
        
qsub -N "shovel" -v site='shovelnose_mountain' run.snow17.calibration.pbs
qsub -N "brook" -v site='brookmere' run.snow17.calibration.pbs
qsub -N "light" -v site='lightning_lake' run.snow17.calibration.pbs
qsub -N "call" -v site='callaghan' run.snow17.calibration.pbs
qsub -N "orchid" -v site='orchid_lake' run.snow17.calibration.pbs
qsub -N "palisade" -v site='palisade_lake' run.snow17.calibration.pbs
qsub -N "grouse" -v site='grouse_mountain' run.snow17.calibration.pbs
qsub -N "dog" -v site='dog_mountain' run.snow17.calibration.pbs
qsub -N "stave" -v site='stave_lake' run.snow17.calibration.pbs
qsub -N "nahat" -v site='nahatlatch' run.snow17.calibration.pbs
qsub -N "wahl" -v site='wahleach' run.snow17.calibration.pbs
qsub -N "kles" -v site='klesilkwa' run.snow17.calibration.pbs
qsub -N "hamil" -v site='hamilton_hill' run.snow17.calibration.pbs
qsub -N "upper" -v site='upper_squamish' run.snow17.calibration.pbs
qsub -N "spuz" -v site='spuzzum_creek' run.snow17.calibration.pbs
qsub -N "chill" -v site='chilliwack_river' run.snow17.calibration.pbs
qsub -N "tenq" -v site='tenquille_lake' run.snow17.calibration.pbs


