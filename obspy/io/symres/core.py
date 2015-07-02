# -*- coding: utf-8 -*-
"""
Symmetric Research *.out bindings to ObsPy core module.
:copyright:
    Ayyyna Nogovitsyna
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from struct import unpack
import numpy as np
import ctypes
from obspy import Trace, Stream, UTCDateTime
from obspy.core.compatibility import from_buffer
import os

def _is_symres(filename):
    """
    Checks whether a file is Symmetric Research *.out waveform data or not.
    :type filename: str
    :param filename: Symmetric Research *.out file to be checked.
    :rtype: bool
    :return: ``True`` if a Symmetric Research *.out waveform file.
    """
    vdaq = os.path.dirname(filename) + "\\vdaq.txt"
    if filename[-4:].capitalize() == ".out" and os.path.isfile(vdaq):
        return True
    return False


def _read_symres(filename, **kwargs):
    """
    Reads an Symmetric Research *.out waveform file and returns a Stream object.
    .. warning::
        This function should NOT be called directly, it registers via the
        ObsPy :func:`~obspy.core.stream.read` function, call this instead.
    :type filename: str
    :param filename: Symmetric Research *.out file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :returns: Stream with Traces specified by given file.
    """
    data = {}
    vdaq = os.path.dirname(filename) + "\\vdaq.txt"
    with open(vdaq, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.split(":")
        data[line[0]] = line[1]
    file = filename[-18:]
    starttime = UTCDateTime(int(file[0:4]), int(file[4:6]), int(file[6:8]), int(file[8:10]),
                            int(file[10:12]), int(file[12:14]))
    header = {}
    header['sampling_rate'] = float(data["SampleRate"])
    header['station'] = data["A-DInfo"].strip()
    header['starttime'] = starttime
    qty = int(data["Channels"])
    traces = []
    with open(filename, 'rb') as f:
        counts = f.read(-1)
        counts = from_buffer(counts, dtype=np.dtype('<i4'))
    nd = np.int32(len(counts) / qty)
    counts = counts.reshape(nd, qty)
    counts = np.hsplit(counts, qty)
    for x in range(0, qty):
        shape_data = np.ascontiguousarray(counts[x].reshape(nd))
        header['channel'] = data["Ch" + str(x) + "ID"].strip()
        tr = Trace(shape_data, header=header)
        traces.append(tr)
    if traces[0].stats.channel == 'TIME':
        sec_hi = traces[0].data[0]
        sec_low = traces[0].data[1]
        msec = traces[0].data[2]
        corr = traces[0].data[3]
        sync = (sec_hi << 16) + sec_low + msec / 1000.0
        if corr > 350:
            sync -= 1
        for t in traces:
            t.stats.starttime = UTCDateTime(sync)
    return Stream(traces)
