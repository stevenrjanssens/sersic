#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
	name='sersic',
	version='0.0.0dev',
	description='Sersic profile utils',
	author='Steven Janssens',
	url='http://github.com/stevenrjanssens/sersic',
	packages=find_packages('src'),
	package_dir={'': 'src'},
    install_requires=['numpy', 'scipy', 'pytest']
)
