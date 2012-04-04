from setuptools import setup

setup(
    name = 'geoprobe',
    version = '0.3',
    description = "Reads and (partially) writes seismic data in Landmark's Geoprobe format",
    author = 'Joe Kington',
    author_email = 'joferkington@gmail.com',
    license = 'MIT',
    url = 'http://code.google.com/p/python-geoprobe/',
    packages = ['geoprobe'],
    package_data = {'geoprobe' : ['__init__.py', 'volume.py', 'horizon.py', 
                                  '_volHeader.py', 'ezfault.py', 'utilities.py',
                                  '_2dHeader.py', 'data2d.py', 'common.py',
                                  'swfault.py', 'tsurf.py']}
)
