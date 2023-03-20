"""
write a script that takes any time and then go grab the .mat files needed to build a sac of that time and interest
06252019_154340_UTC.mat
local time str in this format year-month-dayThour:minute:second

questions:
do you think it would be usefull to use USGS api to pull data?


date = filename.split('/')[-1].split('.')[0]
month = date[0:2]
day = date[2:4]
year = date[4:8]
hour = date.split('_')[1][0:2]
minute = date.split('_')[1][2:4]
second = date.split('_')[1][4:6]
time = UTCDateTime(int(year), int(month), int(day), int(hour), int(minute), int(second))
use this syntax to try and streamline process of the path to event method
"""

import os
from obspy.core import UTCDateTime
#from libcomcat.search import get_event_by_id, search just try and use their api


def time_converter():
    local_time_str = UTCDateTime(input("What is the start time of the event? "))
    time_after = int(input("What time after the event in seconds are you wanting? "))
    utc_time = local_time_str + time_after
    new_local_time = utc_time.datetime
    new_local_time = new_local_time.strftime('%Y-%m-%dT%H:%M:%S')

    return local_time_str, new_local_time


def Path_to_Event():
    print(os.getcwd())
    os.chdir('../Data_sets/')
    print(os.getcwd())

    requested_year = int(input("What year do you want to look for? "))
    requested_month = int(input("What month do you want to look for? "))
    requested_date = int(input("What date do you want to look for? "))
    time_start_and_time_end = time_converter()

    path = os.path.join(str(requested_year), str(requested_month), str(requested_date))
    print(time_start_and_time_end[0], time_start_and_time_end[1])
    os.chdir(path)
    print(os.getcwd())


def EventGet():
    Path_to_Event()


"""
psuedocode

make sure cwd is the same as where files are stored
make sure that directory exists
get to specific year folder 
go to specific month folder
get to specifc date folder
get range of times before event and after 

could I use utc time to convert from start time to end time? 

maybe write if statement that can change minutes to hours if over certain value and then convert
that added utc time to normal time? 

convert into sac file 





"""
