#!/usr/bin/env python
#-*- coding: iso-8859-15 -*-

from distutils.core import setup


exec(open('GWsky/version.py').read())


setup(name='GWsky',
      version=__version__,
      description='GWsky: tiling the skymap in Fields of View',
      author='Giuseppe Greco',
      author_email='giuseppe.greco@uniurb.it',
      license='BSD',
      url='https://github.com/ggreco77/GWsky',
      requires=['astropy', 'healpy','matplotlib', 'pandas', 'numpy'],
      long_description="GWsky is an interactive Python script to generate a sequence of pointings given a specific Field of View (FoV). The script aims to split the large GW sky localization into several independent areas. It defines a sequence of FoVs from a fixed position over the sky, e.g., starting from the highest probability pixel. The results are displayed in Aladin Sky Atlas (http://aladin.u-strasbg.fr/) using the SAMPIntegratedClient class. ",
      classifiers=[
          'Development Status :: v2 - Beta',
          'Programming Language :: Python',
          'License :: BSD',
          'Topic :: Scientific/Engineering :: Astronomy',
      ],
)
