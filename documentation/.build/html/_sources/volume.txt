Welcome to my project's documentation!
======================================

:mod:`volume` -- Reference
=============================================
Geoprobe volume files contain a 3D array of values discretized and stored as 
unsigned 8-bit integers. Reading from and writing to Geoprobe volumes is supported
through the volume class.

.. autoclass:: geoprobe.volume
.. automethod:: geoprobe.volume.__init__

The raw uint8 array values can be accessed through the volume.data attribute. This 
is a 3D numpy array (or memory-mapped-file array) with shape (nx,ny,nz). The 
attributes nx, ny, and nz are readonly.  To change them, reshape or otherwise 
operate on volume.data.

.. autoattribute:: geoprobe.volume.nx
.. autoattribute:: geoprobe.volume.ny
.. autoattribute:: geoprobe.volume.nz
.. autoattribute:: geoprobe.volume.data

If the volume object is created from an existing volume file, volume.data will be
a memory-mapped-file array.  Otherwise, this will contain whatever array you set
when creating the object (or assign to volume.data afterwards).  When setting this
attribute, the new array will be recast as a uint8 array.  (i.e. the values will
always be between 0-255, and overflow/wrap-around will occur if the new array has
values outside of this range)

To load volume.data from a memory-mapped file entirely into memory (for faster
access) use volume.load

.. automethod:: geoprobe.volume.load

The original values stored in the volume can be recovered with the volume.d0 and
volume.dv properties.  For example:
   >>>rescaled = volume.data * volume.dv + volume + d0

.. attribute:: geoprobe.volume.dv

   Voxel value scaling factor 

.. attribute:: geoprobe.volume.d0

   Voxel value calibration factor

There are three coordinate systems stored in a geoprobe volume file. 
   1) The array indicies. These range from 0 to nx-1, 0 to ny-1, 0 to nz-1.
   2) The model coordinates (e.g. inline and crossline). These range from volume.xmin
      to volume.xmax, volume.ymin to volume.ymax, and volume.zmin to volume.zmax.
   3) The world coordinates.  These are typically a map projection of some sort (e.g.
      utm). However, the only restriction is that they must be related to the model
      coordinates via an affine transformation.  

.. automethod:: geoprobe.volume.index2model
.. automethod:: geoprobe.volume.model2index
.. automethod:: geoprobe.volume.model2world
.. automethod:: geoprobe.volume.world2model
.. autoattribute:: geoprobe.volume.worldcoords
.. autoattribute:: geoprobe.volume.modelcoords
.. autoattribute:: geoprobe.volume.transform
.. autoattribute:: geoprobe.volume.invtransform
.. autoattribute:: geoprobe.volume.modelCoords
.. autoattribute:: geoprobe.volume.worldCoords

The x, y, and z spacing and starting values are defined in the following attributes.

.. attribute:: geoprobe.volume.x0

   x axis calibration factor.  (model x-coordinate = dx * i + x0, where i is the x-axis index)
   Defaults to 0.0

.. attribute:: geoprobe.volume.y0

   y axis calibration factor. (model y-coordinate = dy * i + y0, where i is the y-axis index)
   Defaults to 0.0

.. attribute:: geoprobe.volume.z0

   z axis calibration factor. (model z-coordinate = dz * i + z0, where i is the z-axis index)
   Defaults to 0.0

.. attribute:: geoprobe.volume.dx
   
   x axis scaling factor. (model x-coordinate = dx * i + x0, where i is the x-axis index)
   :attribute: `geoprobe.volume.dx` may be negative, but cannot be 0. Defaults to 1.0

.. attribute:: geoprobe.volume.dy

   y axis scaling factor. (model y-coordinate = dy * i + y0, where i is the y-axis index)
   :attribute: `geoprobe.volume.dy` may be negative, but cannot be 0. Defaults to 1.0

.. attribute:: geoprobe.volume.dz

   z axis scaling factor. (model z-coordinate = dz * i + z0, where i is the z-axis index)
   :attribute: `geoprobe.volume.dz` may be negative, but cannot be 0. Defaults to 1.0


The xmin, xmax, ymin, ymax, zmin, and zmax attributes define the model coords. Setting
these will change x0, y0, or z0 (respectively).

.. autoattribute:: geoprobe.volume.xmin
.. autoattribute:: geoprobe.volume.xmax
.. autoattribute:: geoprobe.volume.ymin
.. autoattribute:: geoprobe.volume.ymax
.. autoattribute:: geoprobe.volume.zmin
.. autoattribute:: geoprobe.volume.zmax

Changing the sign of dx, dy, or dz changes the order in which the array is stored on
disk and will change volume.<x,y,z>min and volume.<x,y,z>max.  However, volume.data 
is always accessed such that volume.data[0,0,0] corresponds to xmin, ymin, zmin and 
volume.data[nx,ny,nz] corresponds to xmax, ymax, zmax.  

The world x and y grid spacing can be obtained from the volume.dxW and volume.dyW attributes.

.. autoattribute:: geoprobe.volume.dxW
.. autoattribute:: geoprobe.volume.dyW

There are also convience methods to quickly extract an axis-aligned slice at a given model
coordinate. 

.. automethod:: geoprobe.volume.XSlice
.. automethod:: geoprobe.volume.YSlice
.. automethod:: geoprobe.volume.ZSlice
