"""
Setup module for BIBIMBpy.

Pau Ramos 2022

Based on:
https://github.com/brugalada/BIBIMBpy.git
"""

import re
from codecs import open
from os import path

from setuptools import setup, find_packages

this_folder = path.abspath(path.dirname(__file__))

# Get the long description from the README and HISTORY files
with open(path.join(this_folder, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
with open(path.join(this_folder, 'HISTORY.md'), encoding='utf-8') as f:
    history = f.read()

# Get code version from __init__.py (see https://github.com/dfm/emcee/blob/master/setup.py)
vre = re.compile("__version__ = \"(.*?)\"")
m = open(path.join(this_folder, "bibimbpy", "__init__.py")).read()
code_version = vre.findall(m)[0]

setup(
    name='BIBIMBpy',

    version=code_version,

    description='Backward Integration Basic Interface Module for bars with Python',
    long_description=long_description + "\n\n"
                     + '# Changelog\n'
                     + '\n\n'
                     + history,
    long_description_content_type="text/markdown",

    url="https://github.com/brugalada/BIBIMBpy",

    author='Pau Ramos',
    author_email='pramos@fqa.ub.edu',

    license='LGPLv2.1',

    #packages=['orbits', 'initialize', 'utils'],
    packages=find_packages(exclude=['examples']),
    package_data={'': ['LICENSE', 'HISTORY.md', 'INSTALL.md']},
    include_package_data=True,
    #install_requires=['numpy', 'scipy', 'matplotlib', 'agama'],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2.1 (LGPLv2.1)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],

    entry_points={},
)
