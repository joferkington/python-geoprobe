"""Provides a dictonary containing the name, offset, type, and default value 
for each known value in a GeoProbe volume header"""

__license__   = "MIT License <http://http://www.opensource.org/licenses/mit-license.php>"
__copyright__ = "2009, Free Software Foundation"
__author__    = "Joe Kington <jkington@wisc.edu>"

"""Reversed engineered sometime in early 2009 by Joe Kington to make my dissertation a bit easier.
Geoprobe volumes consist of a 3072 byte header, with various essential metadata.  

Voxels are stored as 1 byte values (i.e. 0-255 or -127 to +127) in fortran order (x cycles fastest, 
then y, then z. In other words, the volume is stored as a collection of z slices).  The bottom 
slice is typically written first (thus the default dz of -1).

There's definitely some stuff I don't understand (thus the "_unknown*" variables).  However, this
seems to work for both reading and writing.  No guarentees that it will handle every geoprobe 
volume you'll ever come across, though!"""


#--Header Definition--------------------------------------------------------------
#    default:None indicates that the vaule must be obtained from the input data
headerLength = 3072 #In bytes
headerDef = {
        'magicNum':         {'offset':0,    'type':'>i',    'default':43970 },      # Number indicating a geoprobe volume (default --> GeoProve Verion 2)
        '_unknown1':        {'offset':4,    'type':'>4i',   'default':[1,0,0,0] },  # Unknown, important?
        'path':             {'offset':20,   'type':'300s',  'default':300*" " },    # Path of geoprobe volume (why is this in here?)
        '_unknown2':        {'offset':320,  'type':'>4i',   'default':[8,0,0,8] },  # Unknown, important? (8bit? But why two of them?)
        '_nx':              {'offset':336,  'type':'>I',    'default':None },       # Number of x-values. This is _nx to allow the property volume.nx, which returns vol.data.shape[0]
        '_ny':              {'offset':340,  'type':'>I',    'default':None },       # Number of y-values. (see above)
        '_nz':              {'offset':344,  'type':'>I',    'default':None },       # Number of z-values. (see above)
        '_unknown3':        {'offset':348,  'type':'>i',    'default':0    },       # Unknown, padding??
        'v0':               {'offset':352,  'type':'>f',    'default':0    },       # Voxel value calibration factor
        'x0':               {'offset':356,  'type':'>f',    'default':0    },       # X-axis calibration factor (e.g. x = i*dx + x0, where i is the index value)
        'y0':               {'offset':360,  'type':'>f',    'default':0    },       # Y-axis calibration factor
        'z0':               {'offset':364,  'type':'>f',    'default':0    },       # Z-axis calibration factor
        'dv':               {'offset':368,  'type':'>f',    'default':1    },       # Voxel value scaling factor
        'dx':               {'offset':372,  'type':'>f',    'default':1    },       # X-axis scaling factor
        'dy':               {'offset':376,  'type':'>f',    'default':1    },       # Y-axis scaling factor
        'dz':               {'offset':380,  'type':'>f',    'default':1    },       # Z-axis scaling factor
        'vUnit':            {'offset':384,  'type':'16s',   'default':'Unknown' },  # Physical units for voxels
        'xUnit':            {'offset':400,  'type':'16s',   'default':'Unknown' },  # Physical units for the X-axis
        'yUnit':            {'offset':416,  'type':'16s',   'default':'Unknown' },  # Physical units for the Y-axis
        'zUnit':            {'offset':432,  'type':'16s',   'default':'Unknown' },  # Physical units for the Z-axis
        'vDescrip':         {'offset':448,  'type':'16s',   'default':'Unknown' },  # Voxel description
        'xDescrip':         {'offset':464,  'type':'16s',   'default':'Unknown' },  # X-axis description
        'yDescrip':         {'offset':480,  'type':'16s',   'default':'Unknown' },  # Y-axis description
        'zDescrip':         {'offset':496,  'type':'16s',   'default':'Unknown' },  # Z-axis description
        '_unknown4':        {'offset':512,  'type':'>536i', 'default':536*[0] },    # Padding??
        '_unknown5':        {'offset':2584, 'type':'8s',    'default':'\xaa\xff\xff#\xaa\x00\x00#' }, #Probably important!  No idea what it is, though...
        'georef':           {'offset':2592, 'type':'>12d',  'default':[0,0,1,0,1,1,0,1,1,0,0,1] },    # 3 sets of points for georeferencing. Order: worldX1, worldX2, worldX3, worldY1, worldY2, worldY3, modelY1, modelY2, modelY3, modelX1, modelX2, modelX3
        '_unknown6':        {'offset':2688, 'type':'>14i',  'default':14*[0] },     # Padding??
        'originalNx':       {'offset':2744, 'type':'>I',    'default':None },       # Original dimensions of the x-axis. No idea why these are here
        'originalNy':       {'offset':2748, 'type':'>I',    'default':None },       # Original dimensions of the y-axis
        'originalNz':       {'offset':2752, 'type':'>I',    'default':None },       # Original dimensions of the z-axis
        '_unknown7':        {'offset':2756, 'type':'>6i',   'default':6*[0] },      # Padding??
        'segmentName':      {'offset':2762, 'type':'50s',   'default':50*" " },     # Seems to be the pathname relative to the base geoprobe project directory
        '_unknown8':        {'offset':2812, 'type':'>i',    'default':0      },     # Padding??
        'seisworksProject': {'offset':2816, 'type':'256s',  'default':None }        # Name of the associated seisworks project (?)
}
