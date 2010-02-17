import struct
import numpy as np
from common import BinaryFile
from _2dHeader import headerDef as _headerDef
from _2dHeader import headerLength as _headerLength

class data2d(object):
    """Reads geoprobe 2D data files."""
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
        self._infile.seek(_headerLength)
        self.x, self.y = [], []
        traces = np.zeros((self._numTraces, self._numSamples), dtype=np.uint8)
        for i in xrange(self._numTraces):
            try:
                x,y,traceNum = self._infile.readBinary('>3f')
                self.x.append(x)
                self.y.append(y)
                trace = np.fromfile(self._infile, dtype=np.uint8, count=self._numSamples)
                traces[i:i+self._numSamples] = trace
            except MemoryError, struct.error:
                break
        self.traces = traces
    
