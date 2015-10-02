#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name = "tilesetdesigner",
    version = "0.0.1",
    packages = ['tilesetdesigner'],

    install_requires = ['numpy','stickydesign','svgwrite','lxml'],

    include_package_data=True,
    #cmdclass={'build': build_with_spurious, 'develop': develop_with_spurious},
    entry_points={
        'console_scripts': [
            'tilesetdesigner = tilesetdesigner.scripts:tilesetdesigner'
            ]
        },
    author = "Constantine Evans",
    author_email = "cge@dna.caltech.edu",
    description = "DX Tile Set Designer",
)
