import struct
import numpy as np
from common import BinaryFile
from _2dHeader import headerDef as _headerDef
from _2dHeader import headerLength as _headerLength

class data2d(object):
    """Reads geoprobe 2D data files.
    data2d.data is a 2D numpy array containg the seismic data
    data2d.x is a list of the x-coordinates of each trace
    data2d.y is a list of the y-coordinates of each trace"""
    def __init__(self, filename):
        """
        Input:
            filename: The name of the 2D data file
        """
        self._infile = BinaryFile(filename, 'r')
        self._readHeader()
        self._readTraces()

    def _readHeader(self):
        for varname, props in _headerDef.iteritems():
            offset, fmt = props['offset'], props['type']
            self._infile.seek(offset)
            var = self._infile.readBinary(fmt)
            setattr(self, varname, var)

    def _readTraces(self):
        dtype = [('x', '>f4'), ('y', '>f4'), ('traces', '%i>u1'%self._numSamples)]
        self._infile.seek(_headerLength)
        data = np.fromfile(self._infile, dtype=dtype, count=self._numTraces)
        self.x = data['x']
        self.y = data['y']
        self.data = data['traces']

    @property
    def numTraces(self):
        return self.data.shape[0]

    @property
    def numSamples(self):
        return self.data.shape[1]
