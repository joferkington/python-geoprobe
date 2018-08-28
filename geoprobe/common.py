"""Recipies and various code utility functions."""

from six import string_types

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
        if isinstance(default, string_types):
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
        setattr(obj, self._calculate.__name__, value)
        return value

#-- Raw reading and writing ---------------------------------------------------
def read_binary(infile, fmt):
    """
    Read and unpack a binary value from the file based on string fmt (see the
    struct module for details).
        Input:
            infile: A file-like object to read from.
            fmt: A ``struct`` format string.
        Output:
            A tuple of unpacked data (or a single item if only one item is
            returned from ``struct.unpack``).
    """
    size = struct.calcsize(fmt)
    data = infile.read(size)
    # Reading beyond the end of the file just returns ''
    if len(data) != size:
        raise EOFError('End of file reached')
    data = struct.unpack(fmt, data)

    for item in data:
        # Strip trailing zeros in strings
        if isinstance(item, string_types):
            item = item.strip('\x00')

    # Unpack the tuple if it only has one value
    if len(data) == 1: data = data[0]

    return data

def write_binary(outfile, fmt, dat):
    """
    Pack and write data to the file according to string fmt.
        Input:
            outfile: An open file-like object to write to.
            fmt: A ``struct`` format string.
            dat: Data to pack into binary form.
    """
    # Try expanding input arguments (struct.pack won't take a tuple)
    try:
        dat = struct.pack(fmt, *dat)
    except (TypeError, struct.error):
        # If it's not a sequence (TypeError), or if it's a
        # string (struct.error), don't expand.
        dat = struct.pack(fmt, dat)
    outfile.write(dat)
