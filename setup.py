# -*- coding: utf-8 -*-
""" Install PlugFlowReactor package. """
from setuptools import setup

install_requires = [
    'cantera>=2.4.0',
    'numpy>=1.11.0',
    'matplotlib>=2.2.0',
    'pandas>=0.23.4'
]

setup(
    name = 'PlugFlowReactor',
    packages = ['PlugFlowReactor'],
    include_package_data = True,
    install_requires = install_requires,
    author = 'Walter Dal\'Maz Silva',
    author_email = 'waltermateriais@gmail.com',
    description = 'Plug-flow reactor implemented with Cantera',
    keywords = 'kinetics reactor chemistry cantera transport',
    url = 'https://github.com/waltermateriais/PlugFlowReactor/',
    license = 'UNLICENSE',
    version = '0.1.0',
    zip_safe = False
)
