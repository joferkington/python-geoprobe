from setuptools import setup

setup(
    name = 'geoprobe',
    version = '0.3.1',
    description = "Reads and (partially) writes seismic data in Landmark's Geoprobe format",
    author = 'Joe Kington',
    author_email = 'joferkington@gmail.com',
    license = 'MIT',
    url = 'https://github.com/joferkington/python-geoprobe',
    packages = ['geoprobe'],
    package_data = {'geoprobe' : ['__init__.py', 'volume.py', 'horizon.py',
                                  '_volHeader.py', 'ezfault.py', 'utilities.py',
                                  '_2dHeader.py', 'data2d.py', 'common.py',
                                  'swfault.py', 'tsurf.py']}
)
