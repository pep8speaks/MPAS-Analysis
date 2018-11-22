#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

from setuptools import setup, find_packages
import warnings

isrelease = True

version = '1.1'

if not isrelease:
    import subprocess
    try:
        pipe = subprocess.Popen(
            ["git", "describe", "--always", "--match", "v[0-9]*"],
            stdout=subprocess.PIPE)
        (version, stderr) = pipe.communicate()
    except:
        warnings.warn("WARNING: Couldn't get git revision, using generic "
                      "version string")

setup(name='mpas_analysis',
      version=version,
      description='Analysis for Model for Prediction Across Scales (MPAS) '
                  'simulations.',
      url='https://github.com/MPAS-Dev/MPAS-Analysis',
      author='MPAS-Analysis Developers',
      author_email='mpas-developers@googlegroups.com',
      license='BSD',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'mpas_analysis': ['config.default'],
                    'mpas_analysis.shared.html': ['templates/*'],
                    'mpas_analysis.test': ['test*/*', 'test*/*/*'],
                    'mpas_analysis.shared.plot':
                        ['ColourMapSuite3/*/*.xml',
                         'SciVisColorColormaps/*.xml'],
                    'mpas_analysis.obs': ['analysis_input_files',
                                          'observational_datasets.xml']},
      install_requires=['numpy', 'scipy', 'matplotlib', 'netCDF4', 'xarray',
                        'dask', 'bottleneck', 'basemap', 'lxml',
                        'pyproj', 'pillow', 'cmocean', 'progressbar2',
                        'requests', 'shapely', 'cartopy'],
      entry_points={'console_scripts':
                    ['mpas_analysis = mpas_analysis.__main__:main',
                     'download_analysis_data = '
                     'mpas_analysis.__main__:download_analysis_data']})
