# -*- coding: utf-8 -*-
"""
Baykal XX bindings to ObsPy core module.

:copyright:
    Vladimir Nogovitsyn
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from future.utils import native

from struct import unpack
import numpy as np
import ctypes
from obspy import Trace, Stream, UTCDateTime
from obspy.core.compatibility import from_buffer


class GeneralHeader53(ctypes.LittleEndianStructure):
    _fields_ = [("number_of_channel", ctypes.c_uint16),
                ("test_type", ctypes.c_uint16),
                ("version", ctypes.c_uint16),
                ("day", ctypes.c_uint16),
                ("month", ctypes.c_uint16),
                ("year", ctypes.c_uint16),
                ("satellites", ctypes.c_uint16),
                ("end_invalid", ctypes.c_uint16),
                ("synchronize", ctypes.c_uint16),
                ("digits", ctypes.c_uint16),
                ("start_invalid", ctypes.c_uint16),
                ("constant", ctypes.c_uint16),
                ("acq_version", ctypes.c_uint16),
                ("system_freq", ctypes.c_uint16),
                ("max_criteria", ctypes.c_uint16),
                ("satellite_instart", ctypes.c_uint16),
                ("station_name", ctypes.c_char * 16),
                ("sample_period", ctypes.c_double),
                ("time_begin", ctypes.c_double),
                ("correction", ctypes.c_double),
                ("latitude", ctypes.c_double),
                ("longitude", ctypes.c_double),
                ("before_synch", ctypes.c_uint64),
                ("synch_point", ctypes.c_uint64),
                ("after_synch", ctypes.c_uint64),
                ("synch_point_start", ctypes.c_uint32),
                ("reserved", ctypes.c_uint32)]


# class GeneralHeader60(ctypes.LittleEndianStructure):
# _fields_ = [("number_of_channel", ctypes.c_uint16),
# ("reserved1", ctypes.c_uint16),
# ("version", ctypes.c_uint16),
#                 ("reserved2", ctypes.c_uint16 * 6),
#                 ("digits", ctypes.c_uint16),
#                 ("reserved3", ctypes.c_uint16),
#                 ("frequency", ctypes.c_uint16),
#                 ("reserved4", ctypes.c_uint16 * 4),
#                 ("station_name", ctypes.c_char * 16),
#                 ("reserved5", ctypes.c_double * 3),
#                 ("latitude", ctypes.c_double),
#                 ("longitude", ctypes.c_double),
#                 ("reserved6", ctypes.c_uint64 * 2),
#                 ("time_begin", ctypes.c_uint64),
#                 ("reserved7", ctypes.c_uint16 * 4)]


class ChannelHeader(ctypes.LittleEndianStructure):
    _fields_ = [("channel_number", ctypes.c_short),
                ("reserved", ctypes.c_char * 6),
                ("channel_name", ctypes.c_char * 24),
                ("channel_mark", ctypes.c_char * 24),
                ("channel_rate", ctypes.c_double),
                ("reserved2", ctypes.c_double)]


def _is_xx(filename):
    """
    Checks whether a file is Baykal XX waveform data or not.

    :type filename: str
    :param filename: XX file to be checked.
    :rtype: bool
    :return: ``True`` if a XX waveform file.
    """
    if _get_xx_version(filename):
        return True
    return False


def _get_xx_version(filename):
    """
    Returns version of Baykal XX  waveform data.

    :type filename: str
    :param filename: XX v53 or v60 file to be checked.
    :rtype: str or False
    :return: version string of XX waveform data or ``False`` if unknown.
    """
    header = GeneralHeader53()
    with open(filename, "rb") as fh:
        fh.readinto(header)
    if header.version == 53 or header.version == 60:
        return header.version
    else:
        return False


def _read_xx(filename, **kwargs):
    """
    Reads an Baykal XX waveform file and returns a Stream object.

    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.

    :type filename: str
    :param filename: XX file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream with Traces specified by given file.
    """
    version = _get_xx_version(filename)
    if version == 53:
        return _read_xx53(filename)
    elif version == 60:
        return None


def _read_xx53(filename):
    """
    Reads an Baykal XX v53 waveform file and returns a Stream object.

    :type filename: str
    :param filename: XX v53 file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream with Traces specified by given file.
    """
    general = GeneralHeader53()
    headers = []
    traces = []
    with open(filename, 'rb') as f:
        f.readinto(general)
        starttime = UTCDateTime(general.year, general.month, general.day)
        header = {}
        header['sampling_rate'] = 1 / general.sample_period
        header['station'] = general.station_name
        header['starttime'] = starttime + general.time_begin
        qty = general.number_of_channel
        for x in range(0, qty):
            ch = ChannelHeader()
            f.readinto(ch)
            headers.append(ch)
        data = f.read(-1)
        data = from_buffer(data, dtype=np.dtype('<i4'))
        nd = np.int32(len(data) / qty)
        data = data.reshape(nd, qty)
        data = np.hsplit(data, qty)
    for x in range(0, qty):
        shape_data = np.ascontiguousarray(data[x].reshape(nd))
        header['channel'] = headers[x].channel_name
        tr = Trace(shape_data, header=header)
        traces.append(tr)
    return Stream(traces)

