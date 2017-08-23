#!/usr/bin/env python
from setuptools import setup

setup(
    name = "alhambra",
    version = "0.4.0",
    packages = ['alhambra'],

    install_requires = ['numpy','stickydesign >= 0.5.0','svgwrite','lxml','shutilwhich','peppercompiler', 'ruamel.yaml'],

    include_package_data=True,

    entry_points={
        'console_scripts': [
            'alhambra = alhambra.scripts:alhambra'
            ]
        },
    author = "Constantine Evans",
    author_email = "cgevans@evans.foundation",
    description = "DX Tile Set Designer",
)
