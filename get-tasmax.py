#!/usr/bin/env python
"""
Save as get-tp.py, then run "python get-tp.py".
 
Input file : None
Output file: tasmax_20180101.nc
"""
from ecmwfapi import ECMWFDataServer
import datetime

starttime = datetime.datetime.strptime('2017-06-01','%Y-%m-%d')
endtime = datetime.datetime.strptime('2018-10-31','%Y-%m-%d')
date_series = [starttime + datetime.timedelta(days=x) for x in range(0,(endtime-starttime).days+1)]
for date in date_series:
    date_string = date.strftime("%Y-%m-%d")
    print(date_string)
    output_file = "/storage/data/projects/rci/data/winter_sports/ERA_INTERIM/tasmax/download/tasmax_12_"+date_string+".nc"
    serv_list = {
        "class"  : "ei",
        "dataset": "interim",
        "date"   : date_string,
        "expver" : "1",
        "grid"   : "0.75/0.75",
        "levtype": "sfc",
        "param"  : "201.128",
        "step"   : "3/6/9/12",
        "stream" : "oper",
        "time"   : "12:00:00",
        "type"   : "fc",
        "format" : "netcdf",
        "target" : output_file,
    }
    server = ECMWFDataServer()
    server.retrieve(serv_list)
