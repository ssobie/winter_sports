#!/usr/bin/env python
"""
Save as get-tp.py, then run "python get-tp.py".
 
Input file : None
Output file: tp_20180101.nc
"""
from ecmwfapi import ECMWFDataServer
 
server = ECMWFDataServer()
server.retrieve({
    "class"  : "ei",
    "dataset": "interim",
    "date"   : "2018-01-01",
    "expver" : "1",
    "grid"   : "0.75/0.75",
    "levtype": "sfc",
    "param"  : "228.128",
    "step"   : "12",
    "stream" : "oper",
    "time"   : "00:00:00/12:00:00",
    "type"   : "fc",
    "format" : "netcdf",
    "target" : "/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/pr/tp_20180101.nc",
})
