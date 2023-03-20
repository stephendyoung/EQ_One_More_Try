from datetime import datetime
from libcomcat.search import count, get_event_by_id, search
import json

time_events = search(starttime=datetime(2017, 1, 1, 0, 0), endtime=datetime(2022, 12, 1, 0, 0), maxradiuskm=1000,
                     latitude=36.1344, longitude=-97.0819, minmagnitude=0, maxmagnitude=9.5)


event_id_dict = {}

for event in time_events:
    current_event = event
    current_event_properties = current_event.properties
    utc_time = current_event["time"]
    adj_utc = utc_time / 1000
    readable_time = str(datetime.fromtimestamp(adj_utc))
    event_id_dict[current_event.id] = {
        "Name": (current_event["place"]),
        "Time": readable_time,
        "Magnitude": (current_event["mag"]),
        "Depth": current_event.depth
    }

sorted_dict = dict(sorted(event_id_dict.items(), key=lambda x: x[1]['Magnitude'], reverse=True))

with open("1000km_range_selected_mag.json", "w") as f:
    json.dump(event_id_dict, f, indent=4)
with open("1000km_range_selected_mag_sorted.json", "w") as f:
    json.dump(sorted_dict, f, indent=4)

"""
important pieces of code 
------------------------------------------------

*****************************************************
Starting piece that just looks at all times events

for i in time_events:
    print(f'{i}')
coyle_ok = time_events[4]  # 2017-4-4 23:49:52 5km NE of coyle OK
properties = coyle_ok.properties
print(f'{coyle_ok_time}: {coyle_ok["time"]}')
*****************************************************
Pulls all properties and also includes the adjusted date time code

# ci3144585 = time_events[4]
# properties = ci3144585.properties
# for prop in properties:
#     print(f'{prop}: {ci3144585[prop]}')

# utc_time = coyle_ok["time"]
# adj_utc = utc_time/1000
#
# print(f'DT:{datetime.fromtimestamp(adj_utc)}')
*****************************************************
First iteration of writing to json but doesn't have nested part 
event_dict = {"Name": [], "Time": [], "Magnitude": [], "Depth": []}

event_dict["Name"].append(coyle_ok["place"])
event_dict["Time"].append(coyle_ok["time"]) this very specific event that I was looking at

rows = []
for i in range(len(event_dict["Name"])):
    row = {
    "ID": event_dict["ID"][i],
    "Time": event_dict["Time"][i],
    "Magnitude": event_dict["Magnitude"][i],
    "Depth": event_dict["Depth"][i]
    }
    rows.append(row)

with open("40km_range_selected_mag.json", "w") as f:
    json.dump(rows, f, indent=4)
*****************************************************
Prints the keys and values for the dictionary... not sure if it works still 

print(" | ".join(event_id_dict.keys()))
for i in range(len(max(event_id_dict.values(), key=len))):
    print(" | ".join(str(event_id_dict.get(key, "")[i]) for key in event_dict.keys()))
    
*****************************************************
Gets all properties from all events: 
  
for prop in current_event_properties:
    print(f'{prop}: {current_event["time"]}')
    
"""