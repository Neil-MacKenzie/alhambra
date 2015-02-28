#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name = "tilesetdesigner",
    version = "0.0.1dev",
    packages = ['tilesetdesigner'],

    install_requires = ['numpy','stickydesign','svgwrite','lxml'],

    include_package_data=True,

    author = "Constantine Evans",
    author_email = "cge@dna.caltech.edu",
    description = "DX Tile Set Designer",
)
