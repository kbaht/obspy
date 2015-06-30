# -*- coding: utf-8 -*-
"""
obspy.io.smart24 - GeoInstr Smart24 continuous file read support for ObsPy
=========================================
This module provides read support for Smart24 waveform data.
Authentication and transformation is not supported.

:copyright:
    Yakutsk Division of Geophysic Survey of the Siberian Branch of the RAS, Vladimir Nogovitsyn
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA


if __name__ == '__main__':
    import doctest

    doctest.testmod(exclude_empty=True)