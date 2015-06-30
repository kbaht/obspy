"""
GeoInstr Smart24 continuous file bindings to ObsPy core module.

:copyright:
    Yakutsk Division of Geophysic Survey of the Siberian Branch of the RAS, Vladimir Nogovitsyn
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

from struct import unpack
import ctypes
from obspy.core.compatibility import from_buffer
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np
import os


class FrameHeader(ctypes.LittleEndianStructure):
    _fields_ = [("frame_type", ctypes.c_int32),
                ("trailer_offset", ctypes.c_int32),
                ("frame_creator", ctypes.c_char * 8),
                ("frame_dest", ctypes.c_char * 8),
                ("seq_number", ctypes.c_long),
                ("series", ctypes.c_int32),
                ("auth_key", ctypes.c_int32),
                ("channels_number", ctypes.c_int32),
                ("frame_length", ctypes.c_int32),
                ("nominal_time", ctypes.c_char * 20),
                ("string_count", ctypes.c_int32)]


class ChannelDescription(ctypes.LittleEndianStructure):
    _fields_ = [("authentication", ctypes.c_byte),
                ("transformation", ctypes.c_byte),
                ("sensor_type", ctypes.c_byte),
                ("flag", ctypes.c_byte),
                ("site_name", ctypes.c_char*5),
                ("channel_name", ctypes.c_char*3),
                ("lc", ctypes.c_char*2),
                ("data_format", ctypes.c_char*2),
                ("calib", ctypes.c_float),
                ("calper", ctypes.c_float)]


class ChannelStatus(ctypes.LittleEndianStructure):
    _fields_ = [("status", ctypes.c_char * 8),
                ("gps_sync", ctypes.c_char * 20),
                ("clock_diff", ctypes.c_int32),
                ("latitude", ctypes.c_float),
                ("longitude", ctypes.c_float),
                ("altitude", ctypes.c_float),
                ("lsb", ctypes.c_float)]

class ChannelSubframe(ctypes.LittleEndianStructure):
    _fields_ = [("channel_length", ctypes.c_int32),
                ("auth_offset", ctypes.c_int32),
                ("description", ChannelDescription),
                ("timestamp", ctypes.c_char * 20),
                ("time_length", ctypes.c_int32),
                ("npts", ctypes.c_int32),
                ("status_data_size", ctypes.c_int32),
                ("status_data", ChannelStatus),
                ("data_size", ctypes.c_int32)]


def _is_smart24(filename):
    """
    Checks whether a file is Smart24 waveform data or not.

    :type filename: str
    :param filename: Smart24 continuous file to be checked.
    :rtype: bool
    :return: ``True`` if a Smart24 waveform file.
    """
    with open(filename, 'rb') as fh:
        buff = unpack('<ii8s8s', fh.read(24))
    t = buff[0]
    contains = buff[3]
    if t == 5 and contains[:5] == 'FILE_':
        if contains != 'FILE_CON':
            raise ValueError('Read only continuous recording files')
        else:
            return True
    return False


def _read_smart24(filename):
    """
    Reads an Smart24 waveform file and returns a Stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.

    :type filename: str
    :param filename: Smart24 file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream with Traces specified by given file.
    """
    stat = os.stat(filename)
    frame_header = FrameHeader()
    channels_string = []
    traces = []
    with open(filename, 'rb') as fh:
        fh.readinto(frame_header)
        if frame_header.auth_key != 0:
            raise ValueError('Authentication not supported')
        for i in range(frame_header.string_count / 10):
            channels_string.append(unpack('<5s3s2s', fh.read(10)))
        if frame_header.string_count % 4 != 0:
            size = int(frame_header.string_count / 4 + 1) * 4
        fh.read(size - frame_header.string_count)
        starttime = UTCDateTime(frame_header.nominal_time.replace(' ', 'T'))
        while True:
            if fh.tell() == stat.st_size:
                break
            ch_subframe = ChannelSubframe()
            fh.readinto(ch_subframe)
            if ch_subframe.description.authentication != 0:
                raise ValueError('Authentication not supported')
            if ch_subframe.description.transformation != 0:
                raise ValueError('Transformation not supported')
            if ch_subframe.time_length != 60000:
                raise ValueError('Subframe Time Length always equals 60,000')
            if ch_subframe.status_data_size != 48:
                raise ValueError('Channel Status Size is 48 bytes. ')
            data = from_buffer(fh.read(ch_subframe.data_size), dtype='<i4')
            data = np.ascontiguousarray(data)
            fh.read(12)
            tr = Trace(data)
            tr.stats.sampling_rate = ch_subframe.npts / 60
            tr.stats.channel = ch_subframe.description.channel_name
            tr.stats.station = ch_subframe.description.site_name
            tr.stats.location = ch_subframe.description.lc
            tr.stats.starttime = starttime
            tr.stats.calib = ch_subframe.description.calib

            tr.stats.smart24 = {'latitude': ch_subframe.status_data.latitude,
                                'longitude': ch_subframe.status_data.longitude,
                                'altitude': ch_subframe.status_data.altitude,
                                'lsb': ch_subframe.status_data.lsb,
                                'calper': ch_subframe.description.calper}
            traces.append(tr)
            if len(traces) % frame_header.channels_number == 0:
                starttime += 60
    st =  Stream(traces)
    return st.merge()