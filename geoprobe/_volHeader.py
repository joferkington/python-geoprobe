"""Provides a dictonary containing the name, offset, type, and default value
for each known value in a GeoProbe volume header"""

__license__   = "MIT License <http://http://www.opensource.org/licenses/mit-license.php>"
__copyright__ = "2009, Free Software Foundation"
__author__    = "Joe Kington <jkington@wisc.edu>"

"""
Reversed engineered sometime in early 2009 by Joe Kington to make my
dissertation a bit easier.  Geoprobe volumes consist of a 3072 byte header,
with various essential metadata.

Voxels are stored as 1 byte values (i.e. 0-255 or -127 to +127) in fortran
order (x cycles fastest, then y, then z. In other words, the volume is stored
as a collection of z slices).  The bottom slice is typically written first
(thus the default dz of -1).

There's definitely some stuff I don't understand (thus the "_unknown*"
variables).  However, this seems to work for both reading and writing.  No
guarentees that it will handle every geoprobe volume you'll ever come across,
though!
"""


#--Header Definition--------------------------------------------------------------
#    default:None indicates that the vaule must be obtained from the input data
headerLength = 3072 #In bytes
headerDef = {
    'magicNum':  {'offset':0, 'type':'>i', 'default':43970,
                  'doc':'Number indicating a geoprobe volume'},

    '_unknown1': {'offset':4, 'type':'>4i', 'default':[1,0,0,0],
                  'doc':'Unknown, important?'},

    'path':      {'offset':20, 'type':'300s', 'default':300*" ",
                  'doc':'Path of original geoprobe volume'},

    '_unknown2': {'offset':320, 'type':'>4i', 'default':[8,0,0,8],
                  'doc':'Unknown, important? (8bit? but why two of them?)'},

    '_nx':       {'offset':336, 'type':'>i', 'default':None,
                  'doc':'Number of x-values. this is _nx to allow the' \
                        ' property volume.nx, which returns vol.data.shape[0]'},

    '_ny':       {'offset':340, 'type':'>i', 'default':None,
                  'doc':'Number of y-values. (see above)'},

    '_nz':       {'offset':344, 'type':'>i', 'default':None,
                  'doc':'Number of z-values. (see above)'},

    'v0':        {'offset':352, 'type':'>f', 'default':0,
                  'doc':'Voxel value calibration factor (Physical voxel' \
                        ' values = V*dv + v0, where V is the raw uint8' \
                        ' value in volume.data)'},

    'x0':        {'offset':356, 'type':'>f', 'default':0,
                  'doc':'X-axis calibration factor (e.g. x = i*dx + x0)'},

    'y0':        {'offset':360, 'type':'>f', 'default':0,
                  'doc':'Y-axis calibration factor (e.g. y = j*dy + y0)'},

    'z0':        {'offset':364, 'type':'>f', 'default':0,
                  'doc':'Z-axis calibration factor (e.g. z = k*dz + z0)'},

    'dv':        {'offset':368, 'type':'>f', 'default':1,
                  'doc':'Voxel value scaling factor (Physical voxel' \
                        ' values = V*dv + v0, where V is the raw uint8' \
                        ' value in volume.data)'},

    'dx':        {'offset':372, 'type':'>f', 'default':1,
                  'doc':'X-axis scaling factor (e.g. x = i*dx + x0)'},

    'dy':        {'offset':376, 'type':'>f', 'default':1,
                  'doc':'Y-axis scaling factor (e.g. y = j*dy + y0)'},

    'dz':        {'offset':380, 'type':'>f', 'default':1,
                  'doc':'Z-axis scaling factor (e.g. z = k*dz + z0)'},

    'vunit':     {'offset':384, 'type':'16s', 'default':'unknown',
                  'doc':'Physical units for voxels'},

    'xunit':     {'offset':400, 'type':'16s', 'default':'unknown',
                  'doc':'Physical units for the x-axis'},

    'yunit':     {'offset':416, 'type':'16s', 'default':'unknown',
                  'doc':'Physical units for the y-axis'},

    'zunit':     {'offset':432, 'type':'16s', 'default':'unknown',
                  'doc':'Physical units for the z-axis'},

    'vdescrip':  {'offset':448, 'type':'16s', 'default':'unknown',
                  'doc':'Voxel description'},

    'xdescrip':  {'offset':464, 'type':'16s', 'default':'unknown',
                  'doc':'X-axis description'},

    'ydescrip':  {'offset':480, 'type':'16s', 'default':'unknown',
                  'doc':'Y-axis description'},

    'zdescrip':  {'offset':496, 'type':'16s', 'default':'unknown',
                  'doc':'Z-axis description'},

    '_unknown3': {'offset':2584, 'type':'8s', 'default':'\xaa\xff\xff#\xaa\x00\x00#',
                  'doc':'Probably important!  no idea what it is, though...'},

    'georef':    {'offset':2592, 'type':'>12d', 'default':[0,0,1,0,1,1,0,1,1,0,0,1],
                  'doc':'3 sets of points for georeferencing. order:' \
                        ' worldx1, worldx2, worldx3, worldy1, worldy2,' \
                        ' worldy3, modely1, modely2, modely3, modelx1,' \
                        ' modelx2, modelx3'},

    'originalnx': {'offset':2744, 'type':'>i', 'default':None,
                   'doc':'Original dimensions of the x-axis'},

    'originalny': {'offset':2748, 'type':'>i', 'default':None,
                   'doc':'Original dimensions of the y-axis'},

    'originalnz': {'offset':2752, 'type':'>i', 'default':None,
                   'doc':'Original dimensions of the z-axis'},

    'segmentname':{'offset':2762, 'type':'50s',   'default':50*" ",
                   'doc':'Pathname relative to the base geoprobe project directory'},

    'seisworksproject':{'offset':2816, 'type':'256s',  'default':'',
                        'doc':'Name of the associated seisworks project'}
}


# The following are almost definitely padding. I'm preserving them here in case the locations are ever needed...
        # '_padding1':        {'offset':348,  'type':'>i',    'default':0    },       # Unknown, padding??
        # '_padding2':        {'offset':512,  'type':'>536i', 'default':536*[0] },    # Padding??
        # '_padding3':        {'offset':2688, 'type':'>14i',  'default':14*[0] },     # Padding??
        # '_padding4':        {'offset':2812, 'type':'>i',    'default':0      },     # Padding??
        # '_padding5':        {'offset':2756, 'type':'>6i',   'default':6*[0] },      # Padding??
