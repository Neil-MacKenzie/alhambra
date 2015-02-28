#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.command.build import build
from setuptools.command.develop import develop

class build_with_spurious(build):
    def run(self):
        import os
        print "Compiling spuriousSSM"
        os.system("cc -Wall -O3 SpuriousDesign/spuriousSSM.c -o tilesetdesigner/spuriousSSM -lm")
        
        build.run(self)

class develop_with_spurious(develop):
    def run(self):
        import os
        print "Compiling spuriousSSM"
        os.system("cc -Wall -O3 SpuriousDesign/spuriousSSM.c -o tilesetdesigner/spuriousSSM -lm")
        
        develop.run(self)

setup(
    name = "tilesetdesigner",
    version = "0.0.1dev",
    packages = ['tilesetdesigner'],

    install_requires = ['numpy','stickydesign','svgwrite','lxml'],

    include_package_data=True,
    cmdclass={'build': build_with_spurious, 'develop': develop_with_spurious},
    entry_points={
        'console_scripts': [
            'tilesetdesigner = tilesetdesigner.scripts:tilesetdesigner'
            ]
        },
    author = "Constantine Evans",
    author_email = "cge@dna.caltech.edu",
    description = "DX Tile Set Designer",
)
