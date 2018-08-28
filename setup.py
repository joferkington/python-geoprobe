from setuptools import setup, find_packages

setup(
    name = 'geoprobe',
    version = '0.4.0',
    description = "Reads and (partially) writes seismic data in Landmark's Geoprobe format",
    author = 'Joe Kington',
    author_email = 'joferkington@gmail.com',
    license = 'MIT',
    url = 'https://github.com/joferkington/python-geoprobe',
    packages = find_packages(),
    install_requires = ['numpy', 'six'],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        ],
)
