# Introduction

Alhambra, formerly tilesetdesigner, is a package for designing DNA tile systems,
currently focused on DX tiles.  It uses stickydesign to create sticky end
sequences, peppercompiler with spuriousSSM to create core sequences, and xgrow
to simulate systems.  It uses an extensible system for tileset design, and a
flexible YAML format for describing the tilesets.

# Installation 

Alhambra is designed to be installed as a Python package.  To install the
current (semi-)stable version from the Python Package Index, you can simply use

    pip install alhambra
	
(If you are reading this after the DNA 24 conference, I would suggest waiting
until 2018-10-15 to do so, as I will be making a few changes, adding documentation,
and putting out a new release for pip that won't be very outdated.  The current pip
version does not include the reduction code.)
	
Alhambra is designed to work with Python 3, but should also work with Python
2.7, though some exceptions may fail, as it makes heavy use of exception
chaining.

To install development versions, you can check out this github repository, and
use `pip -e` or some other method for installing Python packages.  Multiple
versions can be handled through `virtualenv`.

All Alhambra requirements should be handled through setuptools dependencies, but
this is not currently the case for xgrow and xgrowutils.

# Usage

See the documentation in the docs folder.  (Note: If you are reading this after the
DNA 24 conference, I will be putting the documentation online by 2018-10-15, which
should be more convenient to access.)

# Questions

Please send any questions to Constantine Evans, at cevans@evanslabs.org or cge@dna.caltech.edu.
