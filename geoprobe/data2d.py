import numpy as np
from common import BinaryFile
from common import format_headerDef_docs
from _2dHeader import headerDef as _headerDef
from _2dHeader import headerLength as _headerLength

class data2d(object):
    __doc__ = """
    Reads geoprobe 2D data files.

    Useful Attributes:
        data: array containg the seismic data as uint8's
        x: a list of the x-coordinates of each trace
        y: a list of the y-coordinates of each trace
        z: a list of the z-coordinates of each trace\n%s
    """ % format_headerDef_docs(_headerDef)
    # Not a "normal" docstring so that "useful attributes" is set
    # at initialization

    def __init__(self, arg, x=None, y=None, tracenumbers=None, copyFrom=None):
        """
        Input:
            filename: The name of the 2D data file
        """
        if isinstance(arg, basestring):
            filename = arg
            infile = BinaryFile(filename, 'rb')
            self._read_header(infile)
            self._read_traces(infile)
            infile.close()
        else:
            if x is None or y is None:
                raise ValueError('When creating a new data2d object, '\
                        'x and y must both be specified.')
            self._new_file(arg, x, y, tracenumbers, copyFrom)

    def _new_file(self, data, x, y, tracenumbers=None, copyFrom=None,
                  starttime=0.0, endtime=None, samplerate=1.0):
        self.data, self.x, self.y = data, x, y
        if tracenumbers is None:
            tracenumbers = np.arange(1, self.numtraces+1, dtype='>f4')
        self.tracenumbers = tracenumbers
        # Set up the header Dictionary
        if copyFrom is not None:
            # Assume the string is the filename of a 2d data file
            if isinstance(copyFrom, basestring):
                copyFrom = data2d(copyFrom)
            try:
                self.headerValues = copyFrom.headerValues
            except AttributeError:
                raise TypeError('This does not appear to be a valid geoprobe'\
                                ' 2d data object')
        else:
            # Set default attributes
            for varname, info in _headerDef.iteritems():
                setattr(self, varname, info['default'])
            self._numtraces, self._numsamples = data.shape

            if endtime is None:
                endtime = starttime + samplerate * self.numsamples
            self.starttime = starttime
            self.samplerate = samplerate
            self.endtime = endtime

    def write(self, outfile):
        if isinstance(outfile , basestring):
            outfile = BinaryFile(outfile, 'wb')
        self._write_header(outfile)
        self._write_traces(outfile)
        outfile.close()

    def _read_header(self, infile):
        """
        Read the values stored in the file header and set each one as
        an attribue of the data2d object.
        """
        for varname, info in _headerDef.iteritems():
            offset, fmt = info['offset'], info['type']
            infile.seek(offset)
            var = infile.readBinary(fmt)
            setattr(self, varname, var)

    def _write_header(self, outfile):
        """Write the values in self.headerValues to "outfile"."""
        for varname, info in _headerDef.iteritems():
            value = getattr(self, varname, info['default'])
            outfile.seek(info['offset'])
            outfile.writeBinary(info['type'], value)

    def _getHeaderValues(self):
        """
        A dict of all the values stored in the file header
        (Each of these is also an attribute of any data2d object)
        """
        # Return the current instance attributes that are a part of the
        # header definition
        values = {}
        for key in _headerDef.keys():
            # If it's been deleted for some reason, return the default value
            default = _headerDef[key]['default']
            values[key] = getattr(self, key, default)
        return values

    def _setHeaderValues(self, input):
        for key, value in input.iteritems():
            # Only set things in input that are normally in the header
            if key in _headerDef:
                setattr(self, key, value)

    headerValues = property(_getHeaderValues, _setHeaderValues)

    # Unused... Damnit, I need to decide what I'm doing here...
    def _fix_axes(self, data):
        """Reverses the z axis if dz is negative. This ensures that
        self.data[:,0] always corresponds to self.zmin."""
        if self.dz < 0:
            data = data[:, ::-1]
        return data

    def _read_traces(self, infile):
        """
        Read all traces (everything other than the file header)
        from "infile".
        """
        dtype = [('x', '>f4'), ('y', '>f4'), ('tracenum', '>f4'),
                ('traces', '%i>u1'%self._numsamples)]
        infile.seek(_headerLength)
        data = np.fromfile(infile, dtype=dtype, count=self._numtraces)
        self.x = data['x']
        self.y = data['y']
        self.tracenumbers = data['tracenum']
        self.data = data['traces']

    def _write_traces(self, outfile):
        dtype = [('x', '>f4'), ('y', '>f4'), ('tracenum', '>f4'),
                 ('traces', '%i>u1'%self._numsamples)]
        outfile.seek(_headerLength)
        data = np.empty(self.numtraces, dtype=dtype)
        data['x'] = self.x
        data['y'] = self.y
        data['traces'] = self.data
        data['tracenum'] = self.tracenumbers
        data.tofile(outfile, sep='')

    def _get_z(self):
        """Z-values (time/depth) for each trace."""
        try:
            return self._z
        except AttributeError:
            self._z = np.linspace(self.zmin, self.zmax, self.numsamples)
            return self._z
    def _set_z(self, value):
        self._z = value
    z = property(_get_z, _set_z)

    def _bounds(self):
        start = self.z0
        # Fix this!!! "abs" is temporary!!!!
        stop = start + abs(self.dz) * (self.numsamples - 1)
        return start, stop

    def _get_scaled_data(self):
        """Trace array in its original units."""
        return self.data * self.dv + self.v0
    def _set_scaled_data(self, value):
        self.v0 = value.min()
        self.dv = value.ptp() / 256.0
        self.data = (value - self.v0) / self.dv
    scaled_data = property(_get_scaled_data, _set_scaled_data)


    @property
    def zmin(self):
        return min(self._bounds())

    @property
    def zmax(self):
        return max(self._bounds())

    @property
    def numtraces(self):
        """
        The number of traces stored in the file.
        Equivalent to self.data.shape[0]
        """
        return self.data.shape[0]

    @property
    def numsamples(self):
        """
        The number of samples in each trace stored in the file
        Equivalent to self.data.shape[1]
        """
        return self.data.shape[1]
