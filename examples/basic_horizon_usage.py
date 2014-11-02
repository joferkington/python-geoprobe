"""
A quick example of viewing data stored in a geoprobe horizon file
"""

import os

import matplotlib.pyplot as plt

import geoprobe

def main():
    # Path to the example data dir relative to the location of this script.
    # This is just so that the script can be called from a different directory
    datadir = os.path.dirname(__file__) + '/data/'

    # Read an existing geoprobe horizon
    hor = geoprobe.horizon(datadir + 'Horizons/channels.hzn')

    print_info(hor)
    plot(hor)

def print_info(hor):
    """Print some basic information about "hor", a geoprobe.horizon instance"""
    print 'The horizon has a total of %i points, %i of which are'\
          ' auto-tracked' % (hor.data.size, hor.surface.size)
    print 'The horizon has %i manually picked lines' % len(hor.lines)
    print 'The inline coordinates range from', hor.xmin, 'to', hor.xmax
    print 'The crossline coordinates range from', hor.ymin, 'to', hor.ymax
    print 'The depth/time coordinates range from', hor.zmin, 'to', hor.zmax

def plot(hor):
    """Plot the "filled" z-values and "manual picks" in the geoprobe.horizon
    instance"""
    #-- Plot the "filled" values ----------------------------------------------
    plt.figure()

    # hor.grid is a 2D array of the Z-values stored in the horizon
    plt.imshow(hor.grid, cmap=plt.cm.jet_r,
            extent=(hor.xmin, hor.xmax, hor.ymax, hor.ymin))

    #-- Plot the "manual picks" -----------------------------------------------
    plt.hold(True)

    # Here, "line" is a numpy structured array with fields 'x', 'y', 'z', etc.
    # "line_info" is a 4-tuple of (xdir, ydir, zdir, ID) (and is unused here)
    for line_info, line in hor.lines:
        plt.plot(line['x'], line['y'], 'g-')

    #-- Labels, titles, etc ---------------------------------------------------
    cb = plt.colorbar(orientation='horizontal')
    cb.set_label('Depth in meters below sea level')

    plt.axis('image')
    plt.title('An example horizon file')
    plt.xlabel('Inline')
    plt.ylabel('Crossline')

    plt.show()

if __name__ == '__main__':
    main()
