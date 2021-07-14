# -*- coding: utf-8 -*-
"""


    Copyright (C) 2013-2019 ITB - CNR

    This file is part of isoSegmenter.

    isoSegmenter is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    isoSegmenter is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with isoSegmenter.  If not, see <http://www.gnu.org/licenses/>.


If you use isoSegmenter in your work, please cite this manuscript:

    Cozzi P, Milanesi L, Bernardi G. Segmenting the Human Genome into
    Isochores. Evolutionary Bioinformatics. 2015;11:253-261.
    doi:10.4137/EBO.S27693


Created on Thu Nov 26 12:23:08 2015

@author: Paolo Cozzi <paolo.cozzi@ibba.cnr.it>

A module to install isoSegmenter through pip

"""

import io
import os
import re

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))


# Get the long description from the README file
# The long_description field is used by PyPI when you are registering a
# package, to build its home page.
# TODO: convert README.md in README.rst
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


# define function to parse versions
# https://packaging.python.org/en/latest/single_source_version.html
def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='GClib',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=find_version('GClib', "__init__.py"),

    description='A program for segmenting genomes into isochores',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/bunop/isoSegmenter',

    # Author details
    author='Paolo Cozzi',
    author_email='paolo.cozzi@ptp.it',

    # Choose your license
    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        # 'Programming Language :: Python :: 3.4',
    ],

    # What does your project relate to?
    keywords='isochore segmentation',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['tests']),

    # Scripts are files containing Python source code, intended to be started
    # from the command line. The scripts option simply is a list of files to
    # be handled in this way
    scripts=[
        'scripts/isoSegmenter.py',
        'scripts/tileImages.py',
        'scripts/isoFamily.py'],

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[
        "numpy",
        "gdmodule>0.58",
        "Pillow",
        "matplotlib",
        "biopython<1.77"],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        'dev': ['check-manifest'],
        'test': ['coveralls', 'nose'],
    },

    # testing modules
    test_suite="test",

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'GClib': ['FreeSerifBold.ttf'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
#    entry_points={
#        'console_scripts': [
#            'sample=sample:main',
#        ],
#    },
)
