#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from obspy import read
from obspy.io.sac import SACTrace
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime

font = {'weight': 'bold', 'size': 22}
plt.rc('font', **font)

# read sacs into an obspy stream
st_1 = read('../201944222.INF1.STW.TW.sac')

# plot data
st_1.plot()

# In[ ]:


# read specific time around bolide
str_time = UTCDateTime('2019-04-04T22:32:45')
lenth = 60
st_2 = read('../201944222.*.sac', starttime=str_time,
            endtime=str_time + lenth)

# Calibration
st_2[1].data = st_2[1].data * 2.21
st_2[2].data = st_2[2].data * 1.8    # had to change line 32 from 1 to 0 and line 33 from 2 to 0

print(st_2[0].stats)
# plot data
st_2.plot()

# filter data

st_2.filter('bandpass', freqmin=0.3, freqmax=10)

# plot data
st_2.plot()

# spectrogram
for tr in st_2:
    # build figure
    fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
    ax1 = fig.add_axes([0.1, 0.75, 0.7, 0.2])  # [left bottom width height]
    ax1.set_title(tr.stats.channel)
    ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.60], sharex=ax1)
    ax3 = fig.add_axes([0.83, 0.1, 0.03, 0.6])

    # make time vector
    t = np.arange(tr.stats.npts) / tr.stats.sampling_rate

    # title
    ax1.set_title(f'{tr.stats.station}  Start Time = ' + str(tr.stats.starttime))

    # plot time trace
    ax1.plot(t, tr.data, 'k')

    # set y limit for time trace
    # ax1.set_ylim(-1,1)

    # plot spectrogram (bottom subfigure)
    fig = tr.spectrogram(log=True, show=False, dbscale=False, wlen=5, cmap=plt.cm.viridis, axes=ax2)

    # set axis labels
    ax1.set_ylabel('Pressure [pa]')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Frequency [Hz]')

    # set color limits
    mappable = ax2.collections[0]
    ax2.set_ylim(0.1, 20)
    mappable.set_clim(0, 1)

    # colorbar
    plt.colorbar(mappable=mappable, cax=ax3)

    # save
    plt.savefig('./' + str(tr.stats.station) + '-spectro.png')

# save sacs
for tr in st_2:
    tr.write(tr.stats.station + '.trimmed.sac', format='sac')

# In[ ]:


# In[ ]:
