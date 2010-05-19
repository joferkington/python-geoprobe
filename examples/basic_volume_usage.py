"""
A quick example of viewing data stored in a geoprobe volume file.
"""
import os

import numpy as np
import matplotlib.pyplot as plt

import geoprobe

def main():
    # Path to the example data dir relative to the location of this script.
    datadir = os.path.dirname(__file__) + '/data/'

    # Read an existing geoprobe volume
    vol = geoprobe.volume(datadir + 'Volumes/example.vol')

    # Print some info
    print_info(vol)

    # Example plots
    plot(vol)

def plot(vol):
    """Plot the first inline and first crossline in "vol", a geoprobe.volume
    instance."""
    # Plot the first inline in the volume
    plt.figure()
    plt.imshow(vol.XSlice(vol.xmin))
    # Note: instead of vol.XSlice, we could have used vol.data[0,:,:].T
    plt.title('Inline %i' % vol.xmin)
    
    # Plot the first crossline in the volume
    plt.figure()
    plt.imshow(vol.YSlice(vol.ymin))
    # Note: instead of vol.XSlice, we could have used vol.data[:,0,:].T
    plt.title('Crossline %i' % vol.xmin)

    plt.show()

def print_info(vol):
    """Print some basic information about "vol", a geoprobe.volume instance."""
    # Print out some basic information 
    print 'The volume has dimensions of (nx, ny, nz):', vol.nx, vol.ny, vol.nz
    print 'The inline coordinates range from', vol.xmin, 'to', vol.xmax
    print 'The inline spacing is:', vol.dxW, 'world units'
    print 'The crossline coordinates range from', vol.ymin, 'to', vol.ymax
    print 'The crossline spacing is:', vol.dyW, 'world units'
    print 'The depth/time coordinates range from', vol.zmin, 'to', vol.zmax

    # Determine the locations of the corners
    print 'The world coordinates of the corners of the volume are:'
    print '    Lower-left:', vol.model2world(vol.xmin, vol.ymin)
    print '    Upper-left:', vol.model2world(vol.xmin, vol.ymax)
    print '    Upper-right:', vol.model2world(vol.xmax, vol.ymax)
    print '    Lower-right:', vol.model2world(vol.xmax, vol.ymin)

if __name__ == '__main__':
    main()
