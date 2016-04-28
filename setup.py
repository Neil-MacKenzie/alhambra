#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name = "tilesetdesigner",
    version = "0.2.0",
    packages = ['tilesetdesigner'],

    install_requires = ['numpy','stickydesign >= 0.4.2','svgwrite','lxml','shutilwhich'],

    include_package_data=True,
    #cmdclass={'build': build_with_spurious, 'develop': develop_with_spurious},
    entry_points={
        'console_scripts': [
            'tilesetdesigner = tilesetdesigner.scripts:tilesetdesigner'
            ]
        },
    author = "Constantine Evans",
    author_email = "cgevans@evans.foundation",
    description = "DX Tile Set Designer",
)
