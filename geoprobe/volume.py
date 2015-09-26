__license__   = "MIT License <http://http://www.opensource.org/licenses/mit-license.php>"
__copyright__ = "2009, Free Software Foundation"
__author__    = "Joe Kington <jkington@wisc.edu>"

import os
import numpy as np

# Dictonary of header values and offsets for a geoprobe volume
from _volHeader import headerDef as _headerDef
from _volHeader import headerLength as _headerLength

# Common methods
from common import BinaryFile
from common import format_headerDef_docs

import utilities

# Factory function for creating new Volume objects...
def volume(input, copyFrom=None, rescale=True, voltype=None):
    """
    Read an existing geoprobe volue or make a new one based on input data
        Input:
            input: Either the path to a geoprobe volume file or data to
                    create a geoprobe object from (either a numpy array or
                    an object convertable into one). The following keyword
                    arguments only apply when creating a new volume from
                    data, not reading from a file.
            copyFrom (default: None): A geoprobe volume object or path to a
                    geoprobe volume file to copy relevant header values
                    from for the new volume object (used only when creating
                    a new volume from a numpy array).
            rescale (default: True): Boolean True/False.  If True, the data
                    will be rescaled so that data.min(), data.max()
                    correspond to 0, 255 in volume.data before being
                    converted to 8-bit values.  If False, the data data
                    will not be rescaled before converting to type uint8.
                    (i.e. 0.9987 --> 0, 257.887 --> 1, etc due to rounding
                    and wrap around)
            voltype (default: None): Explicitly set the type of volume.
                    By default, the volume type will be guessed from the
                    input file. Valid options are: "hdf", "geoprobe_v2"
    """
    typestrings = {'hdf':HDFVolume, 'geoprobe_v2':GeoprobeVolumeV2}
    if voltype is None:
        voltype = formats[0]
    else:
        if voltype in typestrings:
            voltype = typestrings[voltype]

    # What were we given as input?
    if isinstance(input, basestring):
        # Assume strings are filenames of a geoprobe array
        for vol_format in formats:
            if _check_validity(vol_format, input):
                vol = vol_format()
                vol._readVolume(input)
                return vol
        else:
            raise IOError('This does not appear to be a valid geoprobe file!')

    else:
        # If it's not a string, just assume it's a numpy array or
        # convertable into one and try to make a new volume out of it
        vol = voltype()
        vol._newVolume(input, copyFrom=copyFrom, rescale=rescale)
        return vol

def _check_validity(vol_format, filename):
    try:
        volfile = vol_format.format_type(filename, 'rb')
        is_valid = volfile.is_valid()
        volfile.close()
        return is_valid
    except (IOError, EOFError):
        return False

def isValidVolume(filename):
    """Tests whether a filename is a valid geoprobe file. Returns boolean
    True/False."""
    return any([_check_validity(format, filename) for format in formats])

#-- Base Volume Class ---------------------------------------------------------
class Volume(object):
    # Not a "normal" docstring so that "useful attributes" is set at runtime
    __doc__ = """
    Reads and writes geoprobe volumes

    Useful attributes set at initialization:\n%s
    """ % format_headerDef_docs(_headerDef)

    format_type = None

    def __init__(self):
        """See the "volume" function's docstring for constructor information"""
        pass

    def __getitem__(self, value):
        """Expects slices in model coordinates. Returns the equivalent slice
        of volume.data as a numpy array."""
        new_values = self._convert_slice_to_indicies(value)
        dataslice = self.data.__getitem__(new_values)
        return dataslice

    def _convert_slice_to_indicies(self, value):
        """Converts a slice or integer passed into __getitem__ into indicies"""
        value = list(value)
        new_values = []
        for item, axis in zip(value, ('x', 'y', 'z')):
            if isinstance(item, int) or isinstance(item, float):
                new_values.append(int(self.model2index(item, axis=axis)))
            elif isinstance(item, slice):
                newslice = []
                for val in [item.start, item.stop, item.step]:
                    if val is not None:
                        newslice.append(int(self.model2index(val, axis=axis)))
                    else:
                        newslice.append(None)
                new_values.append(slice(*newslice))
            else:
                new_values.append(item)
        if len(new_values) == 1:
            new_values.extend([None, None])
        new_values = tuple(new_values)
        return new_values


    def _readVolume(self,filename):
        """
        Reads the header of a geoprobe volume and sets attributes based on it
        """
        self._filename = filename
        self._infile = self.format_type(filename, 'rb')
        self.headerValues = self._infile.read_header()

    def _newVolume(self, data, copyFrom=None, rescale=True, fix_axes=True):
        """Takes a numpy array and makes a geoprobe volume. This
        volume can then be written to disk using the write() method."""

        data = np.asarray(data)

        #Set up the header Dictionary
        if copyFrom is not None:
            # Assume the string is the filname of a geoprobe volume
            if isinstance(copyFrom, basestring):
                copyFrom = volume(copyFrom)
            try:
                self.headerValues = copyFrom.headerValues
            except AttributeError:
                raise TypeError('This does not appear to be a valid'\
                                ' geoprobe volume object')
        else:
            # Set default attributes
            for varname, info in _headerDef.iteritems():
                setattr(self, varname, info['default'])
            (self.originalnx, self.originalny, self.originalnz) = data.shape

        if rescale:
            self.v0 = data.min()
            self.dv = data.ptp() / 255.0
            data -= self.v0
            data /= self.dv

        if fix_axes:
            data = self._fixAxes(data)

        self._data = data

    def _fixAxes(self, data):
        """
        Reverses the x, y, and z axes if dx, dy, or dz (respectively) are
        negative.  This ensures that self.data[0,0,0] always corresponds
        to self.xmin, self.ymin, self.zmin.
        """
        if self.dz < 0:
            data = data[:,:,::-1]
        if self.dy < 0:
            data = data[:,::-1,:]
        if self.dx < 0:
            data = data[::-1,:,:]
        return data

    def load(self):
        """Reads an entire Geoprobe volume into memory and returns
        a numpy array contining the volume."""
        try:
            dat = self._infile.read_data()
            dat = self._fixAxes(dat)
        except AttributeError:
            # If there's no self._infile attribute, then the volume
            # has been initialized from an array in memory, and we don't
            # need to load the data...
            dat = self.data
        return dat

    def write(self, filename):
        """Writes a geoprobe volume to disk."""
        outfile = self.format_type(filename, 'wb')
        outfile.write_header(self.headerValues)
        outfile.write_data(self.data)

    def crop(self, xmin=None, xmax=None, ymin=None, ymax=None,
            zmin=None, zmax=None, copy_data=True):
        """Crops a volume to the limits specified (in model coordinates)
        by the input arguments.  Returns a new volume instance containing
        data for the cropped region."""

        # Set defaults (verbose, but I want to avoid metaprogramming...)
        if (xmin is None) or (xmin < self.xmin):
            xmin = self.xmin
        if (xmax is None) or (xmax > self.xmax):
            xmax = self.xmax

        if (ymin is None) or (ymin < self.ymin):
            ymin = self.ymin
        if (ymax is None) or (ymax > self.ymax):
            ymax = self.ymax

        if (zmin is None) or (zmin < self.zmin):
            zmin = self.zmin
        if (zmax is None) or (zmax > self.zmax):
            zmax = self.zmax

        # Convert input units to indicies...
        xstart, xstop = self.model2index([xmin, xmax], axis='x')
        ystart, ystop = self.model2index([ymin, ymax], axis='y')
        zstart, zstop = self.model2index([zmin, zmax], axis='z')

        # Crop data
        data = self.data[xstart:xstop, ystart:ystop, zstart:zstop]
        if copy_data:
            data = data.copy()

        # Make a new volume instance and set it's mininum model coords
        vol = type(self)()
        vol._newVolume(data, copyFrom=self, rescale=False, fix_axes=False)
        vol.xmin, vol.ymin, vol.zmin = xmin, ymin, zmin

        return vol

    #-- data property --------------------------------------------------------
    def _getData(self):
        """A 3D numpy array of the volume data.  Contains raw uint8 values.
        If the volume object is based on a file, this is a memory-mapped-file
        array."""
        # Prevents _getData from being called more than once
        try:
            return self._data
        except AttributeError:
            self._data = self._fixAxes(self._infile.memmap_data())
            return self._data

    def _setData(self, newData):
        """Set self.data without making a copy in-memory"""
        newData = np.asanyarray(newData, dtype=np.uint8)
        # Make sure we're dealing with a 3D numpy array
        try:
            self._nx, self._ny, self._nz = newData.shape
        except (ValueError, AttributeError):
            raise TypeError('Data must be a 3D numpy array')

        # We don't update dv and d0 here.  This is to avoid overwriting the
        # "real" dv and d0 when you're manipulating the data in some way. When
        # making a new volume object from an existing numpy array (in
        # _newVolume), dv and d0 are set
        self._data = newData

    data = property(_getData, _setData)
    #--------------------------------------------------------------------------


    #-- headerValues property (avoiding @headerValues.setter to keep
    #-- compatibility w/ python2.5)----
    def _getHeaderValues(self):
        """
        A dict with a key, value pair for each value in the geoprobe volume
        file header.  These are stored as class attributes in each volume
        instance.  Do not set vol.headerValues['something'] = newValue. This
        will silently fail to update vol.something. However, setting
        vol.headerValues = {'something':newValue} will work fine.  Normally,
        header values should be set directly through the volume's attributes,
        e.g. vol.something = newValue.
        """
        # Return the current instance attributes that are a part of the
        # header definition
        values = {}
        for key in _headerDef:
            # If it's been deleted for some reason, return the default value
            default = _headerDef[key]['default']
            values[key] = getattr(self, key, default)
        return values

    def _setHeaderValues(self, input):
        # This needs to raise an error if set via headerValues['x'] = y!
        # Don't know how to do it easily, though...
        for key, value in input.iteritems():
            # Only set things in input that are normally in the header
            if key in _headerDef:
                setattr(self, key, value)

    headerValues = property(_getHeaderValues, _setHeaderValues)
    #--------------------------------------------------------------------------


    #-- xmin, xmax, etc -------------------------------------------------------
    def _getVolumeBound(self, axis, max=True):
        n = [self.nx, self.ny, self.nz][axis]
        d = [self.dx, self.dy, self.dz][axis]
        offset = [self.x0, self.y0, self.z0][axis]
        stop = offset + d * (n-1)
        start = offset
        if ((max is True) & (d>0)) or ((max is False) & (d<0)):
            return stop
        else:
            return start

    def _setVolumeBound(self, value, axis, max=True):
        axisLetter = ['x','y','z'][axis]
        n = [self.nx, self.ny, self.nz][axis]
        d = [self.dx, self.dy, self.dz][axis]
        if ((max is True) & (d>0)) or ((max is False) & (d<0)):
            value = value + (n-1) * abs(d)
        setattr(self, axisLetter+'0', value)

    xmin = property(lambda self: self._getVolumeBound(axis=0, max=False),
            lambda self, value: self._setVolumeBound(value, axis=0, max=False),
            doc="Mininum x model coordinate")

    ymin = property(lambda self: self._getVolumeBound(axis=1, max=False),
            lambda self, value: self._setVolumeBound(value, axis=1, max=False),
            doc="Mininum y model coordinate")

    zmin = property(lambda self: self._getVolumeBound(axis=2, max=False),
            lambda self, value: self._setVolumeBound(value, axis=2, max=False),
            doc="Mininum z model coordinate")

    xmax = property(lambda self: self._getVolumeBound(axis=0, max=True),
            lambda self, value: self._setVolumeBound(value, axis=0, max=True),
            doc="Maximum x model coordinate")

    ymax = property(lambda self: self._getVolumeBound(axis=1, max=True),
            lambda self, value: self._setVolumeBound(value, axis=1, max=True),
            doc="Maximum y model coordinate")

    zmax = property(lambda self: self._getVolumeBound(axis=2, max=True),
            lambda self, value: self._setVolumeBound(value, axis=2, max=True),
            doc="Maximum z model coordinate")
    #--------------------------------------------------------------------------


    #-- nx, ny, nz ------------------------------------------------------------
    @property
    def nx(self):
        """The number of x values in the array (read-only)"""
        return self.data.shape[0]
    @property
    def ny(self):
        """The number of y values in the array (read-only)"""
        return self.data.shape[1]
    @property
    def nz(self):
        """The number of z values in the array (read-only)"""
        return self.data.shape[2]
    #--------------------------------------------------------------------------


    # Need to fix these... should be 3x2 instead of 2x3...
    # fix transform when I do!

    #-- worldCoords property-------------------------------------
    def _getWorldCoords(self):
        """A 2x3 array containing 3 points in world coordinates corresponding
        to the points in modelCoords.
        [[x1, x2, x3],
         [y2, y2, y3]]"""
        return np.asarray(self.georef[:6]).reshape(2,3)
    def _setWorldCoords(self, newValues):
        self.georef[:6] = np.asarray(newValues).flatten()
    worldCoords = property(_getWorldCoords, _setWorldCoords)
    #------------------------------------------------------------

    #-- modelCoords property-------------------------------------
    #   georef format: Xw1,Xw2,Xw3,Yw1,Yw2,Yw3,Ym1,Ym2,Ym3,Xm1,Xm2,Xm3
    #   This means that we need to flipud the modelCoords!!
    #   (looks suspicious, double-check)
    def _getModelCoords(self):
        """A 2x3 array containing 3 points in model coordinates corresponding
        to the points in worldCoords.
        [[x1, x2, x3],
         [y2, y2, y3]]"""
        return np.flipud( np.asarray(self.georef[6:]).reshape(2,3) )
    def _setModelCoords(self, newValues):
        self.georef[6:] = np.asarray(newValues).flatten()
    modelCoords = property(_getModelCoords, _setModelCoords)
    #-----------------------------------------------------------

    @property
    def dxW(self):
        """X-spacing interval in world coordinates"""
        X,Y = self.model2world([self.xmin, self.xmax], [self.ymin, self.ymin])
        dist = np.sqrt(X.ptp()**2 + Y.ptp()**2)
        return dist / self.nx
    @property
    def dyW(self):
        """Y-spacing interval in world coordinates"""
        X,Y = self.model2world([self.xmin, self.xmin], [self.ymin, self.ymax])
        dist = np.sqrt(X.ptp()**2 + Y.ptp()**2)
        return dist / self.ny

    def YSlice(self, Ypos):
        """Takes a slice of the volume at a constant y-value (i.e. the
        slice is in the direction of the x-axis) This is a convience
        function to avoid calling volume.model2index before slicing
        and transposing (for easier plotting) after.
        Input:
            Ypos: Y-Value given in model coordinates
        Output:
            A 2D (NZ x NX) numpy array"""
        Ypos = self.model2index(Ypos, axis='y')
        return self.data[:,Ypos,:].transpose()

    def XSlice(self, Xpos):
        """Takes a slice of the volume at a constant x-value (i.e. the
        slice is in the direction of the y-axis) This is a convience
        function to avoid calling volume.model2index before slicing
        and transposing (for easier plotting) after.
        Input:
            Xpos: X-Value given in model coordinates
        Output:
            A 2D (NZ x NY) numpy array"""
        Xpos = self.model2index(Xpos)
        return self.data[Xpos,:,:].transpose()

    def ZSlice(self, Zpos):
        """Takes a slice of the volume at a constant z-value (i.e.
        a depth slice) This is a convience function to avoid calling
        volume.model2index before slicing and transposing (for
        easier plotting) after.
        Input:
            Zpos: Z-Value given in model coordinates
        Output:
            A 2D (NY x NX) numpy array"""
        Zpos = self.model2index(Zpos, axis='z')
        return self.data[:,:,Zpos].transpose()

    def extract_section(self, x, y, zmin=None, zmax=None, coords='model'):
        """Extracts an "arbitrary" section defined by vertices in *x* and *y*.
        Input:
            *x*: A sequence of x coords in the coordinate system specified by
                *coords*. (*coords* defaults to "model")
            *y*: A sequence of y coords in the coordinate system specified by
                *coords*
            *zmin*: The minimum "z" value for the returned section
                (in depth/time)
            *zmax*: The maximum "z" value for the returned section
                (in depth/time)
            *coords*: Either "model" or "world", specifying whether *x* and *y*
                are given in model or world coordinates."""
        if coords == 'world':
            x, y = self.world2model(x, y)
        elif coords != 'model':
            raise ValueError('"coords" must be either "world" or "model".')
        if zmin is None:
            zmin = self.zmin
        if zmax is None:
            zmax = self.zmax
        zmax = self.model2index(zmax, axis='z') + 1
        zmin = self.model2index(zmin, axis='z')

        # If zmin and zmax are out of bounds, things will work fine, but the
        # returned array won't represent data between zmin and zmin, only the
        # data between self.zmin and self.zmax. Best to be noisy!
        if (zmin < 0) or (zmax > self.nz):
            msg = '"zmin" and "zmax" must be between %.1f and %.1f'
            raise ValueError(msg % (self.zmin, self.zmax))

        x, y = self.model2index(x, y)
        section, xi, yi = utilities.extract_section(self.data, x, y, zmin, zmax)
        xi, yi = self.index2model(xi, yi)
        if coords == 'world':
            xi, yi = self.model2world(xi, yi)
        return section, xi, yi

    def model2index(self, *coords, **kwargs):
        """Converts model coordinates into array indicies
        Input:
            *coords: x,y,z model coordinates to be converted. Will accept numpy
                arrays.  If only one coordinate is specified, it is assumed to
                be x, and if two are specified, they are assumed to be x,y.
                To override this choice (i.e. convert just y or z) use the
                axis keyword argument
            axis (default None): The axis that the given model coordinate
                is from.  Use 0,1,2 or 'x','y','z'.
        Output:
            A tuple of converted coordinates (or just the coord if only one
            was input)
        Examples:
            XI = volume.model2index(X)
            XI,YI = volume.model2index(X,Y)
            XI,YI,ZI = volume.model2index(X,Y,Z)
            ZI = volume.model2index(Z, axis=2) # (or axis='z')
            YI,ZI = volume.model2index(Y,Z,axis='y')
        """
        return self._modelIndex(*coords, **kwargs)

    def index2model(self, *coords, **kwargs):
        """Converts array indicies to model coordinates
        Input:
            *coords: x,y,z array indicies to be converted. Will accept numpy
                arrays.  If only one index is specified, it is assumed to
                be x, and if two are specified, they are assumed to be x,y.
                To override this choice (i.e. convert just y or z) use the
                axis keyword argument
            axis (default None): The axis that the given model coordinate
                is from.  Use 0,1,2 or 'x','y','z'.
        Output:
            A tuple of converted indicies (or just the index if only one was
            input)
        Examples:
            X = volume.index2model(XI)
            X,Y = volume.index2model(XI,YI)
            X,Y,Z = volume.index2model(XI,YI,ZI)
            Z = volume.index2model(ZI, axis=2) # (or axis='z')
        """
        kwargs['inverse'] = True
        return self._modelIndex(*coords, **kwargs)

    def _modelIndex(self, *coords, **kwargs):
        """
        Consolidates model2index and index2model. See docstrings of
        index2model and model2index for more detail.
        Input:
            value: The model coordinate or volume index to be converted
            axis (default 0): Which axis the coordinate represents.
                Can be either 0,1,2 or 'X','Y','Z'
            inverse (defalut False): True to convert index to model coords
        Output:
            converted coordinate(s)
        """
        #-- Again, this looks needlessly complex... ------------------
        # Unfortunately, I can't figure out a more simple method while
        # still retaining the flexibility I want this function to have
        # Also, I'd really rather not use **kwagrs here, but I have to
        # use *args first, and I can't have *args and the axis=None.

        def convert(value,axis,inverse):
            """Actually convert the coordinates"""
            value = np.asarray(value)

            # Select the proper starting point and step value
            #   The "axis % 3" is to allow some odd but useful stuff...
            #   E.g. converting Y,Z pairs
            axis = axis % 3
            mins = [self.xmin, self.ymin, self.zmin]
            Ds = [abs(self.dx), abs(self.dy), abs(self.dz)]
            minimum, d = mins[axis], Ds[axis]

            # Convert the coordinates
            if inverse:  # index2model
                return value * d + minimum
            else: # model2index
                idx = (value - minimum) / d
                if int_conversion:
                    return idx.astype(int)
                else:
                    return idx

        #-- Handle user input -------------------------------------------------

        # Set the default values of axis and inverse
        axis = kwargs.get('axis', 0)
        inverse = kwargs.get('inverse', False)
        int_conversion = kwargs.get('int_conversion', True)

        # Make sure we have a valid axis
        if axis not in [0,1,2,'x','y','z','X','Y','Z']:
            raise ValueError('"axis" must be one of 0,1,2 or "x","y","z"')

        # Allow both 0,1,2 and 'x','y','z' (or 'X','Y','Z') for axis
        if isinstance(axis, basestring):
            axis = axis.upper()
            axis = {'X':0, 'Y':1, 'Z':2}[axis]

        # Handle calling f(x), f(x,y), f(x,y,z), f(z,axis=2), etc
        converted = [convert(x, i+axis, inverse) for i,x in enumerate(coords)]

        # If there's just one value, return it, otherwise return a sequence
        if len(converted) == 1:
            return converted[0]
        else:
            return converted

    def model2world(self,crossline,inline=None):
        """
        Converts model coordinates to world coordinates.  Accepts either a Nx2
        list/numpy.array or 2 seperate lists/numpy.array's with crossline,
        inline coordinates.  Returns 2 numpy arrays X and Y with world
        coordinates. If a Nx2 or 2xN array was given, returns a single Nx2 or
        2xN array instead of seperate x & y 1D arrays.
        """
        return self._transformCoords(crossline, inline, self.transform)

    def world2model(self, x, y=None):
        """
        Converts world coordinates to model coordinates.  Accepts either a Nx2
        list/numpy.array or 2 seperate lists/numpy.array's with x, y
        coordinates.  Returns 2 numpy arrays X and Y with model coordinates.
        If a Nx2 or 2xN array was given, returns a single Nx2 or 2xN array
        instead of seperate x & y 1D arrays.
        """
        return self._transformCoords(x, y, self.invtransform)

    def _transformCoords(self,x,y,transform):
        """
        Consolidates model2world and world2model.  Takes x,y and a transform
        matrix and ouputs transformed coords X and Y (both are Nx1).  If y is
        None, assumes x is a Nx2 matrix where y is the 2nd column
        """
        #-- Process input ------------------
        return_array, transpose = False, False
        x = np.squeeze(np.asarray(x))

        # If only one array is given...
        if y is None:
            return_array = True
            # Assume x is 2xN or Nx2 and slice appropriately
            if 2 in x.shape:
                if x.shape[0] == 2:
                    x,y = x[0,:], x[1,:]
                elif x.shape[1] == 2:
                    transpose = True
                    x,y = x[:,0], x[:,1]
                x = np.squeeze(x)
            else:
                raise ValueError('If Y is not given, X must be an'
                                 ' Nx2 or 2xN array!')
        y = np.squeeze(np.asarray(y))
        shape = x.shape
        x, y = x.flatten(), y.flatten()

        #-- Convert Coordinates -------------------------------------
        try:
            # Make a G-matrix from the input coords
            dataIn = np.vstack((x,y,np.ones(x.size)))
        except ValueError:
            raise ValueError('X and Y inputs must be the same size!!')
        # Output world coords
        dataOut = np.dot(transform,dataIn)

        #-- Output Data with Same Shape as Input --------------------
        if return_array:
            if transpose:
                dataOut = dataOut.T
            return dataOut
        else:
            # Return x, y
            X,Y = dataOut[0,:], dataOut[1,:]
            if X.size == 1:
                return X[0], Y[0]
            else:
                return X.reshape(shape), Y.reshape(shape)

    @property
    def transform(self):
        """A 2x3 numpy array describing an affine transformation between
        model and world coordinates. (Read only property. Calculated from
        volume.modelCoords and volume.worldCoords.)"""
        # Detailed explanation of inversion...
        # Ok:
        #  Ax_m + By_m + C = x_w
        #  Ex_m + Fy_m + G = y_w
        # where x_m, y_m are the model coords
        #       x_w, y_w are the world coords
        #       A,B,C,E,F,G are the transformation between them (m)
        #
        # This can be expressed in matrix form as:
        # G = |x_m1, y_m1, 1|  m = |A, E|  d = |x_w1, y_w1|
        #     |x_m2, y_m2, 1|      |B, F|      |x_w2, y_w2|
        #     | .     .    .|      |C, G|      | .     .  |
        #     |x_mn, y_mn, 1|                  |x_wn, y_wn|
        #  Where:
        #     G*m = d
        #
        # Then we solve for m, as we know G and d
        #   m = d*inv(G)
        G = np.vstack((self.modelCoords, [1,1,1]))
        d = self.worldCoords
        m = np.dot(d,np.linalg.inv(G))
        return m

    @property
    def invtransform(self):
        """A 2x3 numpy array to transform between world and
        model coordinates. (Read only property. Calculated from
        volume.modelCoords and volume.worldCoords.)"""
        # See explanation in self.transform
        G = np.vstack((self.worldCoords, [1,1,1]))
        d = self.modelCoords
        m = np.dot(d,np.linalg.inv(G))
        return m

class GeoprobeVolumeFileV2(object):
    """Low-level operations for reading and writing to Geoprobe Volume format
    version 2.0 seismic data files."""
    def __init__(self, filename, mode):
        self.filename = filename
        self.mode = mode
        self._file = BinaryFile(self.filename, self.mode)

    def is_valid(self):
        """Returns a boolean indicating whether this is a valid file."""
        header = self.read_header()
        nx, ny, nz = [header[item] for item in ('_nx', '_ny', '_nz')]
        volSize = os.stat(self.filename).st_size
        predSize = nx*ny*nz + _headerLength

        if (header['magicNum'] != 43970) or (volSize != predSize):
            return False
        else:
            return True


    def read_header(self):
        """Reads and returns the header of a geoprobe volume."""
        header = dict()
        for varname, info in _headerDef.iteritems():
            self._file.seek(info['offset'])
            value = self._file.readBinary(info['type'])
            header[varname] = value
        return header

    def read_data(self):
        """Reads an entire Geoprobe volume into memory and returns
        a numpy array contining the volume."""
        header = self.read_header()
        nx, ny, nz = [header[item] for item in ('_nx', '_ny', '_nz')]
        dat = np.fromfile(self.filename, dtype=np.uint8)
        dat = dat[_headerLength:]
        dat = dat.reshape((nz, ny, nx)).T
        return dat

    def memmap_data(self):
        """Return an object similar to a memory-mapped numpy array."""
        header = self.read_header()
        nx, ny, nz = [header[item] for item in ('_nx', '_ny', '_nz')]
        dat = np.memmap(filename=self.filename, mode='r',
            offset=_headerLength, order='F',
            shape=(nx, ny, nz)
            )
        return dat

    def write_header(self, header):
        """Write the values in the dict "header" to the file."""
        for varname, info in _headerDef.iteritems():
            value = header.get(varname, info['default'])
            self._file.seek(info['offset'])
            self._file.writeBinary(info['type'], value)

    def write_data(self, data):
        """Writes a geoprobe volume to disk."""
        self._file.seek(_headerLength)
        # Write to file in Fortran order...
        #  Not using data.ravel('F'), as it appears to make a copy if the
        #  array is C-ordered (doubling memory usage). Instead, we just
        #  write the transpose with the tofile method. (tofile writes in
        #  C-order regardless of the order of the input array, thus
        #  requring the transpose for both F and C ordered arrays)
        data.T.tofile(self._file)

    def close(self):
        return self._file.close()

class HDFVolumeFile(object):
    """Low level operations for reading and writing to hdf5 files."""
    dataset_name = '/volume'
    def __init__(self, filename, mode):
        import h5py
        self.filename = filename
        self.mode = mode
        if 'b' in self.mode:
            self.mode = self.mode.replace('b', '')
        self._file = h5py.File(self.filename, self.mode)
        if 'w' not in self.mode:
            self.dataset = self._file[self.dataset_name]

    def is_valid(self):
        return self.dataset_name in self._file

    def read_header(self):
        return self._file.attrs

    def read_data(self):
        out = np.empty(self.dataset.shape, self.dataset.dtype)
        self.dataset.read_direct(out)
        return out

    def memmap_data(self):
        return self.dataset

    def write_header(self, header):
        for key, value in header.iteritems():
            self._file.attrs[key] = value

    def write_data(self, data):
        header = self.read_header()
        dx, dy, dz = [header[item] for item in ['dx', 'dy', 'dz']]
        x0, y0, z0 = [header[item] for item in ['x0', 'y0', 'z0']]

        # Because h5py doesn't support negative steps when indexing, we have to
        # reverse the x0, y0, z0 before writing and make all d's positive!
        if dz < 0:
            z0 = z0 + (data.shape[2] - 1) * dz
            data = data[:,:,::-1]
        if dy < 0:
            y0 = y0 + (data.shape[1] - 1) * dy
            data = data[:,::-1,:]
        if dx < 0:
            x0 = x0 + (data.shape[0] - 1) * dx
            data = data[::-1,:,:]

        self.write_header(dict(dx=abs(dx), dy=abs(dy), dz=abs(dz),
                          x0=x0, y0=y0, z0=z0))
        self.dataset = self._file.create_dataset(self.dataset_name, data=data)

    def close(self):
        return self._file.close()

#-- Valid file formats --------------------------------------------------------
class GeoprobeVolumeV2(Volume):
    format_type = GeoprobeVolumeFileV2

class HDFVolume(Volume):
    format_type = HDFVolumeFile

formats = [GeoprobeVolumeV2]

# If h5py is available, add HDFVolume's to the list of available formats
try:
    import h5py
    formats.append(HDFVolume)
except ImportError:
    pass


