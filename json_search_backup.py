import json
from dateutil.parser import parse
import os
import scipy
from obspy.core import UTCDateTime
from obspy.io.sac import SACTrace
import numpy as np
import glob

# event_start_file ="/04072018_071603_UTC.mat"
# directory = "/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/2018/4/07"
# os.chdir("/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/2018/4/07")
#
# files = glob.glob(os.path.join(directory, '*.mat'))
#
# try:
#     index = files.index(os.path.join(directory, event_start_file))
# except ValueError:
#     index = -1
#
# after_file = files[index+1:]
#
# print(after_file)


def mat_to_sac(path, event_time):
    event_start_file = event_time
    directory = os.getcwd()
    files = glob.glob(os.path.join(directory, '*.mat'))
    try:
        index = [os.path.normpath(f) for f in files].index(os.path.join(directory, event_start_file))
    except ValueError:
        index = -1
    x = 5
    after_file = files[index+1:index+1+x]

    for file in after_file:
        OSU1_coords = {"W": [36.1344, -97.0819], "N": [36.1342, -97.0813], 'S': [36.1347, -97.0814]}
        stationnames = list(OSU1_coords.keys())

        data_load = scipy.io.loadmat(file)

        data_mat = data_load['SP'][0, 0]['rData']

        OSU1_W = []
        OSU1_N = []
        OSU1_S = []
        for line in data_mat:
            OSU1_W.append(line[0])
            OSU1_S.append(line[1])
            OSU1_N.append(line[2])
        OSU1_data = [OSU1_W, OSU1_N, OSU1_S]

        date = file.split('/')[-1].split('.')[0]
        month = date[0:2]
        day = date[2:4]
        year = date[4:8]
        hour = date.split('_')[1][0:2]
        minute = date.split('_')[1][2:4]
        second = date.split('_')[1][4:6]
        time = UTCDateTime(int(year), int(month), int(day), int(hour), int(minute), int(second))


        for i in range(3):
            data = np.array(OSU1_data[i])

            sampling_rate = 1000
            starttime = time
            long = OSU1_coords[stationnames[i]][1]
            lat = OSU1_coords[stationnames[i]][0]
            channel = stationnames[i]
            stationname = 'OSU1'

            header = {'kstnm': stationname, 'kcmpnm': channel, 'stla': lat, 'stlo': long, 'delta': 1 / sampling_rate,

                      'nzyear': starttime.year, 'nzjday': starttime.julday, 'nzhour': starttime.hour,
                      'nzmin': starttime.minute,
                      'nzsec': starttime.second, 'nzmsec': int(starttime.microsecond / 1000)}

            sac = SACTrace(data=data, **header)
            sac_trace = sac.to_obspy_trace(sac)
            sac_trace.data = data
            sac_trace.stats.sac = header
            print(sac_trace.stats)
            sac_trace.plot()


with open("1000km_range_selected_mag_sorted.json") as f:
    events = json.load(f)

count = 0
limit = 7

for event_id in events:
    if count == limit:
        break
    event_parse = parse(events[event_id]["Time"])
    print(event_parse)
    year = str(event_parse.year)
    month = str(event_parse.month)
    day = str(event_parse.day).zfill(2)
    hour, minute, second = str(event_parse.hour).zfill(2), str(event_parse.minute).zfill(2), \
        str(event_parse.second).zfill(2)
    file_path = os.path.join("/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/", year, month, day)

    if os.path.exists(file_path):
        os.chdir(file_path)
        event_file = f'{month.zfill(2)}{day}{year}_{hour}{minute}{second}_UTC.mat'
        event_file_joined = os.path.join("/Users/stephenyoung/EQ_Stuff/EQ_Files/Data_sets/", year, month, day,
                                         event_file)
        mat_to_sac(file_path, event_file_joined)
        print(event_file_joined)
        print(os.getcwd())
    else:
        print('directory does not exist')
    count += 1


"""
add time as well and suffix on file name... try to use try and except. 
also maybe some sort of limit to certain files...
"""
