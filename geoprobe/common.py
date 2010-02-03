#! /usr/bin/python
import sys, os
import struct
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

