# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
from future.builtins import *  # NOQA
import numpy as np
import ctypes
from obspy import Trace, Stream, UTCDateTime
from obspy.core.compatibility import from_buffer
from ConfigParser import ConfigParser
from io import BytesIO
import os
from struct import unpack

class ChannelGroup(object):

    def __init__(self, ids, freq, channels):
        self.__id = ids
        self.__freq = freq
        self.__channels = channels

    def get_id(self):
        return self.__id

    def set_id(self, ids):
        self.__id = ids

    def get_freq(self):
        return self.__freq

    def set_freq(self, freq):
        self.__freq = freq

    def get_channels(self):
        return self.__channels

    def set_channels(self, channels):
        self.__channels = channels

    ids = property(get_id, set_id)
    freq = property(get_freq, set_freq)
    channels = property(get_channels, set_channels)


class Channel(object):

    def __init__(self, index, name, freq):
        self.__name = name
        self.__freq = freq
        self.__index = index

    def get_name(self):
        return self.__name

    def set_name(self, name):
        self.__name = name

    def get_freq(self):
        return self.__freq

    def set_freq(self, freq):
        self.__freq = freq

    def get_index(self):
        return self.__index

    def set_index(self, index):
        self.__index = index

    name = property(get_name, set_name)
    freq = property(get_freq, set_freq)
    index = property(get_index, set_index)

def _is_sdas(filename):
    try:
        with open(filename, 'rb') as f:
            cfg = ConfigParser()
            cfg.readfp(BytesIO(f.read(98)))
            return cfg.has_option('HEADER', 'OFFSET_TO_DATA')
    except:
        return False

def _read_sdas(filename, **kwargs):
    groups = []
    channel_list = []
    traces = []
    with open(filename, 'rb') as f:
        cfg = ConfigParser()
        cfg.readfp(BytesIO(f.read(98)))
        header_size = cfg.getint('HEADER', 'HEADER_SIZE')
        offset_to_data = cfg.getint('HEADER', 'OFFSET_TO_DATA')
        f.seek(0)
        cfg.readfp(BytesIO(f.read(header_size - 98)))
        header = {}
        header['station'] = cfg.get('SYSTEM', 'NAME')
        group_count = cfg.getint('SYSTEM', 'N_GROUP')
        stream = cfg.get('FILE', 'STREAM')
        fragment_duration = cfg.getint('STREAM' + stream, 'REC_SIZE_SEC')
        duration = cfg.getint('STREAM' + stream, 'FILE_SIZE_SEC')
        channels = cfg.get('STREAM' + stream, 'CH#').split(',')
        for x in range(1, group_count + 1):
            freq = cfg.getint('GROUP' + str(x), 'FREQ')
            group_channels = cfg.get('GROUP' + str(x), 'CH#')
            group = ChannelGroup(x, freq, group_channels.split(','))
            groups.append(group)
        for i, channel in enumerate(channels):
            ch_name = cfg.get('CH' + channel, 'NAME')
            for g in groups:
                if channel in g.channels:
                    freq = g.freq
            channel_list.append(Channel(i + 1, ch_name, freq))
        f.seek(offset_to_data, os.SEEK_SET)
        buff = f.read(256)
        day, month, year, hour, minute, second, ms = unpack('<hhhhhhh', buff[8:22])
        starttime = UTCDateTime(year=year, month=month, day=day, hour=hour, minute=minute, second=second, microsecond=ms)
        header['starttime'] = starttime
        for channel in channel_list:
            fragment_points = channel.freq * fragment_duration
            fragment_size = fragment_points * 2
            block_size = fragment_size * len(channel_list)
            fragment_count = duration / fragment_duration
            f.seek(offset_to_data, os.SEEK_SET)
            channel_data = np.array([], dtype=np.dtype('<i4'))
            for x in range(0, fragment_count):
                f.read(256)
                channel_offset = fragment_size * (channel.index - 1)
                f.seek(channel_offset, os.SEEK_CUR)
                data = f.read(fragment_size)
                dt = np.dtype(np.uint16)
                dt = dt.newbyteorder('<')
                data = np.ascontiguousarray(from_buffer(data, dtype=dt))
                channel_data = np.hstack((channel_data, data))
                f.seek(block_size - fragment_size * channel.index, os.SEEK_CUR)
            header['sampling_rate'] = channel.freq
            header['channel'] = channel.name
            tr = Trace(channel_data, header=header)
            traces.append(tr)
    return Stream(traces)


