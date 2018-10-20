#!/usr/bin/env python
from setuptools import setup

setup(
    name="alhambra",
    version="1.0.0.dev2",
    packages=['alhambra'],
    install_requires=[
        'numpy', 'stickydesign >= 0.8.0', 'svgwrite', 'lxml', 'shutilwhich',
        'peppercompiler >= 0.1.2', 'ruamel.yaml', 'cssutils'
    ],
    include_package_data=True,
    entry_points={'console_scripts': ['alhambra = alhambra.scripts:alhambra']},
    author="Constantine Evans",
    author_email="cgevans@evans.foundation",
    description="DX Tile Set Designer", )
