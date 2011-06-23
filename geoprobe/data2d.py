import struct
import numpy as np
from common import BinaryFile
from common import format_headerDef_docs
from _2dHeader import headerDef as _headerDef
from _2dHeader import headerLength as _headerLength

class data2d(object):
    # Not a "normal" docstring so that "useful attributes" is set 
    # at initialization
    __doc__ = """
    Reads geoprobe 2D data files.

    Useful Attributes:
        data: array containg the seismic data as uint8's
        x: a list of the x-coordinates of each trace
        y: a list of the y-coordinates of each trace\n%s
    """ % format_headerDef_docs(_headerDef)

    def __init__(self, filename):
        """
        Input:
            filename: The name of the 2D data file
        """
        self._infile = BinaryFile(filename, 'r')
        self._readHeader()
        self._readTraces()

    def _readHeader(self):
        """
        Read the values stored in the file header and set each one as
        an attribue of the data2d object.
        """
        for varname, props in _headerDef.iteritems():
            offset, fmt = props['offset'], props['type']
            self._infile.seek(offset)
            var = self._infile.readBinary(fmt)
            setattr(self, varname, var)

    @property
    def headerValues(self):
        """
        A dict of all the values stored in the file header
        (Each of these is also an attribute of any data2d object)
        """
        output = {}
        for varname in _headerDef.keys():
            output[varname] = getattr(self, varname)
        return output
        
    def _readTraces(self):
        """
        Read all traces (everything other than the file header) 
        from self._infile
        """
        dtype = [('x', '>f4'), ('y', '>f4'), ('tracenum', '>f4'), 
                ('traces', '%i>u1'%self._numSamples)]
        self._infile.seek(_headerLength)
        data = np.fromfile(self._infile, dtype=dtype, count=self._numTraces)
        self.x = data['x']
        self.y = data['y']
        self.tracenumbers = data['tracenum']
        self.data = data['traces']

    @property
    def numTraces(self):
        """
        The number of traces stored in the file.
        Equivalent to self.data.shape[0]
        """
        return self.data.shape[0]

    @property
    def numSamples(self):
        """
        The number of samples in each trace stored in the file
        Equivalent to self.data.shape[1]
        """
        return self.data.shape[1]
