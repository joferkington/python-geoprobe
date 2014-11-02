"""Recipies and various code utility functions."""

import struct
import textwrap

#-- Miscellaneous -------------------------------------------------------------
def format_headerDef_docs(headerDef, initial_indent=8, subsequent_indent=12):
    """
    Format the attributes contained in a headerDef for pretty printing
    (i.e. for use in docstrings)
    """
    attribute_docs = ''
    initial_indent *= ' '
    subsequent_indent *= ' '

    for key in sorted(headerDef.keys()):
        value = headerDef[key]
        default = value['default']
        if isinstance(default, basestring):
            default = default.strip()

        doc = '%s: %s (default=%s)' % (key, value['doc'], repr(default))
        doc = textwrap.fill(doc, initial_indent=initial_indent,
                subsequent_indent=subsequent_indent)

        if not key.startswith('_'):
            attribute_docs += doc + '\n'

    return attribute_docs

class cached_property(object):
    """
    A decorator class that ensures that properties are only evaluated once.
    From: <http://code.activestate.com/recipes/363602-lazy-property-evaluation/>
    """
    def __init__(self, calculate_function):
        self._calculate = calculate_function

    def __get__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.func_name, value)
        return value

#-- Raw reading and writing ---------------------------------------------------
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
            if isinstance(item, basestring):
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


