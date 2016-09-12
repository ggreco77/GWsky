#!/usr/bin/env python
#-*- coding: iso-8859-15 -*-

from setuptools import setup, find_packages
from glob import glob
import re

package_name = 'GWSKY'
version_str = re.search(r'^__version__\s+=\s+[\'"](.+)[\'"]',
        open('%s/version.py' % (package_name, )).read(),
        re.M).group(1)

setup(name=package_name,
        version=version_str,
        description='GWsky: tiling the skymap in Fields of View',
        author='Giuseppe Greco',
        author_email='giuseppe.greco@uniurb.it',
        license='BSD',
        packages=find_packages(),
        long_description=open('README.md').read(),
        install_requires=['requests',
            ]
)
