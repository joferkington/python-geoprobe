
__license__   = "MIT License <http://http://www.opensource.org/licenses/mit-license.php>"
__copyright__ = "2009, Free Software Foundation"
__author__    = "Joe Kington <jkington@wisc.edu>"

import struct, os, array 
import numpy as np

#Dictonary of header values and offsets for a geoprobe volume
from _header import headerDef as _headerDef
from _header import headerLength as _headerLength
    

def isGeoprobeVolume(filename):
    try: testvol = vol(filename)
    except: return False
    
    volID = testvol.magic
    volSize = os.stat(filename).st_size
    predSize = testvol.nx*testvol.ny*testvol.nz + _headerLength

    if (volID!=43970) or (volSize!=predSize): return False
    else: return True


class volume(object):
    """Reads and writes geoprobe volumes"""

    # TODO: Implement a clip method

    def __init__(self, input, copyFrom=None, rescale=True):
        """Read an existing geoprobe volue or make a new one based on input data
            Input: 
                input: Either the path to a geoprobe volume file or data to create
                        a geoprobe object from (either a numpy array or an object
                        convertable into one). The following keyword arguments only
                        apply when creating a new volume from data, not reading from
                        a file.
                copyFrom (default: None): A geoprobe volume object or path to a geoprobe
                        volume file to copy relevant header values from for the new 
                        volume object (used only when creating a new volume from a
                        numpy array).
                rescale (default: True): Boolean True/False.  If True, the data will be
                        rescaled so that data.min(), data.max() correspond to 0, 255 in
                        volume.data before being converted to 8-bit values.  If False, 
                        the data data will not be rescaled before converting to type 
                        uint8. (i.e. 0.9987 --> 0, 257.887 --> 1, etc due to rounding
                        and wrap around)""" 

        # What were we given as input?  
        if isinstance(input, str):  
            # Assume strings are filenames of a geoprobe array
            self._readVolume(input)
        else:
            # If it's not a string, just assume it's a numpy array or convertable
            # into one and try to make a new volume out of it
            self._newVolume(input, copyFrom, rescale)
 
    def _readVolume(self,filename):
        """Reads the header of a geoprobe volume and sets attributes based on it"""
    
        #Probably best to let the calling program handle IO errors
        self._filename = filename
        self._infile = _volumeFile(filename,'rb')

        #--Set up proper attributes (setting headerValues handles this)
        self.headerValues = self._infile.readHeader()
       
    def _newVolume(self,data,copyFrom=None,rescale=True):
        """Takes a numpy array and makes a geoprobe volume. This
        volume can then be written to disk using the write() method."""

        data = np.asarray(data)

        #Set up the header Dictionary
        if copyFrom:
            if isinstance(copyFrom, str): # Assume the string is the filname of a geoprobe volume
                copyFrom = volume(copyFrom)
            try:
                self.headerValues = copyFrom.headerValues
            except AttributeError:
                raise TypeError('This does not appear to be a valid geoprobe volume object')
        else:
            # Set default attributes
            for id, values in _headerDef.iteritems():
                setattr(self, id, values['default'])
            (self.originalNx, self.originalNy, self.originalNz) = data.shape

        if rescale:
            self.v0 = data.min()
            self.dv = data.ptp() / 255.0
            data -= self.v0
            data /= self.dv

        self.data = self._fixAxes(data)

    def _fixAxes(self,data):
        """Transposes the axes of geoprobe volume numpy
        array so that axis1=x, axis2=y, and axis3=z. Also
        flips the z-axis, if needed."""
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
        dat = np.fromfile(self._filename, dtype=np.uint8)
        dat = dat[_headerLength:]
        dat = dat.reshape((self.nz, self.ny, self.nx)).T
        return dat

    def write(self, filename):
        """Writes a geoprobe volume to disk using memmapped arrays"""
        self._infile = _volumeFile(filename, 'w')
        self._infile.writeHeader(self)
        self._infile.close()

        self._filemap = np.memmap(filename, offset=_headerLength, 
                dtype=np.uint8, mode='r+', shape=(self.nx,self.ny,self.nz), 
                order='Fortran')

        # Reverse the various axes, if necessary
        #data = self._fixAxes(self.data)

        self._filemap[:,:,:] = self.data[:,:,:]
        self._filemap.flush()

    #-- data property ------------------------------------------------
    def _getData(self):
        """A 3D numpy array of the volume data.  Contains raw uint8 values.
        If the volume object is based on a file, this is a memory-mapped-file
        array."""
        try:
            # Have we already done this? 
            # (Isn't done on initilization to avoid large virtual mem usage)
            dat = self._data
        except AttributeError:
            #If not, we're reading from a file, so make a memory-mapped-file numpy array
            dat = np.memmap(filename=self._filename, mode='r',
                offset=_headerLength, order='F', 
                shape=(self._nx, self._ny, self._nz) 
                )
            dat = self._fixAxes(dat)
            self._data = dat
        return dat
    def _setData(self, newData):
        newData = np.asarray(newData).astype(np.uint8)
        try:
            self._nx, self._ny, self._nz = newData.shape
        except ValueError, AttributeError:
            raise TypeError('Data must be a 3D numpy array')
        # We don't update dv and d0 here.  This is to avoid overwriting the "real"
        # dv and d0 when you're manipulating the data in some way. When making a 
        # new volume object from an existing numpy array (in _newVolume), dv and d0 are set
        self._data = newData
    data = property(_getData, _setData)
    #-------------------------------------------------------------------


    #-- headerValues property (avoiding @headerValues.setter to keep compatibility w/ python2.5)----
    def _getHeaderValues(self):
        """A dict with a key, value pair for each value in the geoprobe volume file 
        header.  These are stored as class attributes in each volume instance.
        Do not set vol.headerValues['something'] = newValue. This will silently fail to 
        update vol.something. However, setting vol.headerValues = {'something':newValue} 
        will work fine.  Normally, header values should be set directly through the volume's
        attributes, e.g. vol.something = newValue."""
        # Return the current instance attributes that are a part of the header definition  
        values = {}
        for key in _headerDef:
            # If it's been deleted for some reason, return the default value
            default = _headerDef[key]['default']  
            values[key] = getattr(self, key, default)
        return values
    def _setHeaderValues(self, input):
        # This needs to raise an error if set via headerValues['x'] = y! Don't know how to do it easily, though...
        for key, value in input.iteritems():
            # Only set things in input that are normally in the header
            if key in _headerDef:
                setattr(self, key, input[key])
    headerValues = property(_getHeaderValues, _setHeaderValues)
    #------------------------------------------------------------------------------------------


    #-- xmin, xmax, etc -----------------------------------------------------------------------
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
        offset = [self.x0, self.y0, self.z0][axis]
        if ((max is True) & (d>0)) or ((max is False) & (d<0)):
            value = value + (n-1) * d
        setattr(self, axisLetter+'0', value)

    xmin = property(lambda self: self._getVolumeBound(axis=0, max=False), 
                    lambda self, value: self._setVolumeBound(value, 0, max=False),
                    doc="Mininum x model coordinate")
    ymin = property(lambda self: self._getVolumeBound(axis=1, max=False), 
                    lambda self, value: self._setVolumeBound(value, axis=1, max=False),
                    doc="Mininum y model coordinate")
    zmin = property(lambda self: self._getVolumeBound(axis=2, max=False), 
                    lambda self, value: self._setVolumeBound(value, axis=2, max=False),
                    doc="Mininum z model coordinate")
    xmax = property(lambda self: self._getVolumeBound(axis=0, max=True), 
                    lambda self, value: self._setVolumeBound(value, 0, max=True),
                    doc="Maximum x model coordinate")
    ymax = property(lambda self: self._getVolumeBound(axis=1, max=True), 
                    lambda self, value: self._setVolumeBound(value, 1, max=True),
                    doc="Maximum y model coordinate")
    zmax = property(lambda self: self._getVolumeBound(axis=2, max=True), 
                    lambda self, value: self._setVolumeBound(value, 2, max=True),
                    doc="Maximum z model coordinate")
    #-----------------------------------------------------------------------------------


    #-- nx, ny, nz ---------------------------------------------------------------------
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
    #-----------------------------------------------------------------------------------


    # Need to fix these... should be 3x2 instead of 2x3... fix transform when I do.

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
    #   This means that we need to flipud the modelCoords!! (looks suspicious, double-check)
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
        Ypos = self.model2index(Ypos)
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
        Zpos = self.model2index(Zpos)
        return self.data[:,:,Zpos].transpose()

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
            A tuple of converted coordinates (or just the coord if only one was input)
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
            A tuple of converted indicies (or just the index if only one was input)
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
            axis = axis % 3 # This is to allow some odd but useful stuff... E.g. converting Y,Z pairs
            mins = [self.xmin, self.ymin, self.zmin]
            Ds = [abs(self.dx), abs(self.dy), abs(self.dz)]
            min, d = mins[axis], Ds[axis]

            # Convert the coordinates
            if inverse:  # index2model
                return value * d + min
            else: # model2index
                idx = (value - min) / d
                return idx.astype(np.int)

        #-- Handle user input ------------------------------------------------------ 

        # Set the default values of axis and inverse (unfortunately, have to use **kwargs)
        axis = kwargs.get('axis', 0)
        inverse = kwargs.get('inverse', False)

        # Make sure we have a valid axis
        if axis not in [0,1,2,'x','y','z','X','Y','Z']:  
            raise ValueError('"axis" must be one of 0,1,2 or "x","y","z"')

        # Allow both 0,1,2 and 'x','y','z' (or 'X','Y','Z') for axis
        if isinstance(axis, str): 
            axis = axis.upper()
            axis = {'X':0, 'Y':1, 'Z':2}[axis]

        # Handle calling f(x), f(x,y), f(x,y,z), f(z,axis=2), etc 
        converted = [convert(x, i+axis, inverse) for i,x in enumerate(coords)]

        # If there's just one value, return it, otherwise return a tuple
        if len(converted) == 1: return converted[0]
        else: return converted

    def model2world(self,crossline,inline=None):
        """Converts model coordinates to world coordinates.  
        Accepts either a Nx2 list/numpy.array or 2 seperate
        lists/numpy.array's with crossline, inline coordinates.
        Returns 2 numpy arrays X and Y with world coordinates."""
        return self._transformCoords(crossline, inline, self.transform)

    def world2model(self, x, y=None):
        """Converts world coordinates to model coordinates.  
        Accepts either a Nx2 list/numpy.array or 2 seperate
        lists/numpy.array's with x, y coordinates.
        Returns 2 numpy arrays X and Y with model coordinates."""
        return self._transformCoords(x, y, self.invtransform)

    def _transformCoords(self,x,y,transform):
        """Consolidates model2world and world2model.
        Takes x,y and a transform matrix and ouputs
        transformed coords X and Y (both are Nx1).  
        If y is None, assumes x is a Nx2 matrix where
        y is the 2nd column"""

        #-- Process input ------------------
        x = np.squeeze(np.asarray(x))
        # If only one array is given...
        if y is None:
            # Assume x is 2xN or Nx2 and slice appropriately
            if 2 in x.shape:
                if x.shape[0] == 2:    x,y = x[0,:], x[1,:]
                elif x.shape[1] == 2:  x,y = x[:,0], x[:,0]
                x = np.squeeze(x)
            else:  
                raise ValueError('If Y is not given, X must be an Nx2 or 2xN array!')
        y = np.squeeze(np.asarray(y))
        # Check size of input
        ndata = x.size
        if x.size != y.size:
            raise ValueError('X and Y inputs must be the same size!!')

        #-- Convert Coordinates ------------------
        dataIn = np.vstack((x,y,np.ones(ndata))) # Input model coords in the form of a G-matrix
        dataOut = np.dot(transform,dataIn)  # Output world Coords
        
        return dataOut[0,:], dataOut[1,:]  # X,Y
    
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
   




#-- Raw reading and writing -------------------------------------------
class _volumeFile(file):
    """Raw reading and writing of values in the volume header"""
    def __init__(self, *args, **kwargs):
        """Initialization is the same as a normal file object
        %s""" % file.__doc__
        file.__init__(self, *args, **kwargs)

    def readBinary(self,fmt):
        """Read and unpack a binary value from the file based
        on string fmt (see the struct module for details).
        """
        size = struct.calcsize(fmt)
        data = self.read(size)
        data = struct.unpack(fmt, data)

        # Unpack the tuple if it only has one value
        if len(data) == 1: data = data[0]

        # Join strings and strip trailing zeros
        if 's' in fmt:
            data = ''.join(data)
            data = data.strip('\x00')

        return data

    def writeBinary(self, fmt, dat):
        """Pack and write data to the file according to string fmt."""
        # If it's a string, pad with \x00 to the appropriate length
        if 's' in fmt:
            length = int(fmt.replace('s',''))
            dat.ljust(length, '\x00')
            dat = struct.pack(fmt, dat)
        # Hackish, but it works...
        else:
            try: dat = struct.pack(fmt, *dat) 
            except TypeError: dat = struct.pack(fmt, dat) # If it's not a sequence, don't expand
        self.write(dat)

    def readHeader(self):
        """Read header values from a geoprobe volume file.  Reads
        values based on the dict defined in _header.py. Returns
        a dict of name:value pairs."""
        headerValues = {}
        for value, props in _headerDef.iteritems():
            offset = props['offset']
            fmt = props['type']
            self.seek(offset)
            headerValues[value] = self.readBinary(fmt)
        return headerValues

    def writeHeader(self, volInstance):
        """Write the header values contained in a volume instance
        "volInstance" at the offsets defined in _header.py.""" 
        for key, props in _headerDef.iteritems():
            default = props['default']
            fmt = props['type']
            offset = props['offset']
            currentValue = getattr(volInstance, key, default)
            self.seek(offset)
            self.writeBinary(fmt, currentValue)

    
