from datetime import datetime
import matplotlib.pyplot as plt
# from usgs-libcomcat.search import count, get_event_by_id, search


box_events = search(starttime=datetime(2017, 1, 1, 0, 0), endtime=datetime(2022, 12, 1, 0, 0), minlatitude=35.7740,
maxlatitude=36.4948, minlongitude=-97.4423, maxlongitude=-96.7215)

