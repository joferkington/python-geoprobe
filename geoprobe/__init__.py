"""
A python module to read and write geoprobe format volumes, horizons,
2d data, and faults

Reads and (mostly) writes seismic data from and to files written by
Landmark's (a subsidiary of Halliburton) Geoprobe software. This
implementation is based on reverse-engineering the file formats, and as
such, is certainly not complete. However, things seem to work.

As a basic example:

>>> from geoprobe import volume
>>> vol = volume('/path/to/geoprobe/volume/file.vol')
>>> print vol.xmin, vol.ymin  # Model coordinate min and max
>>> test_xslice = vol.data[vol.nx/2,:,:] # a memmapped numpy array
"""
__author__ = 'Joe Kington <jkington@geology.wisc.edu>'
__license__ = 'MIT License'
__version__ = '0.3.1'


from horizon import horizon
from volume import volume, Volume
from ezfault import ezfault
from data2d import data2d
import utilities
from colormap import colormap
from swfault import swfault
from tsurf import tsurf

from volume import isValidVolume

__all__ = ['horizon', 'volume', 'Volume', 'ezfault', 'data2d', 'utilities',
           'isValidVolume', 'colormap', 'swfault', 'tsurf']
