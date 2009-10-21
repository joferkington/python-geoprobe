from distutils.core import setup

setup(
    name = 'geoprobe',
    version = 0.1,
    description = "Reads and (partially) writes seismic data in Landmark's Geoprobe format",
    author = 'Joe Kington',
    author_email = 'joferkington@gmail.com',
    url = 'http://code.google.com/p/python-geoprobe/',
    packages = ['geoprobe'],
    package_data = {'geoprobe' : ['__init__.py', 'volume.py', 'horizon.py', 
                                  '_header.py', 'ezfault.py', 'utilities.py']}
)
