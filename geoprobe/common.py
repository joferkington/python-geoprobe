#! /usr/bin/python
import sys, os
import struct
import weakref
import numpy as np

#-- Raw reading and writing -------------------------------------------
class BinaryFile(file):
    """
    Automatically packs or unpacks binary data according to a format
    when reading or writing.
    """
    def __init__(self, *args, **kwargs):
        """
        Initialization is the same as a normal file object
        %s""" % file.__doc__
        file.__init__(self, *args, **kwargs)

    def readBinary(self,fmt):
        """
        Read and unpack a binary value from the file based
        on string fmt (see the struct module for details).
        """
        size = struct.calcsize(fmt)
        data = self.read(size)
        # Reading beyond the end of the file just returns ''
        if len(data) != size:
            raise EOFError('End of file reached')
        data = struct.unpack(fmt, data)

        for item in data:
            # Strip trailing zeros in strings 
            if isinstance(item, str):
                item = item.strip('\x00')

        # Unpack the tuple if it only has one value
        if len(data) == 1: data = data[0]

        return data

    def writeBinary(self, fmt, dat):
        """Pack and write data to the file according to string fmt."""
        # Try expanding input arguments (struct.pack won't take a tuple)
        try: 
            dat = struct.pack(fmt, *dat) 
        except (TypeError, struct.error): 
            # If it's not a sequence (TypeError), or if it's a 
            # string (struct.error), don't expand.
            dat = struct.pack(fmt, dat) 
        self.write(dat)

class StaticCache(object):
    """
    A decorator for very simple cacheing of values. The point isn't 
    to memonize, just to make RAII a bit easier and avoid calling 
    particularly expensive getters in properties more than once.
    """
    def __init__(self, function):
        self.function = function
    def __call__(self, *args):
        try:
            value =  self.cached_value()
            if value is None: # Due to weakref
                raise AttributeError
        except AttributeError:
            retval = self.function(*args)
            self.cached_value = weakref.ref(retval)
            value = self.cached_value()
        return value


