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
	
To install a prerelease version with reduction (from DNA 24) and other significant
improvements, use

	pip install alhambra=1.0.0a1
	
Once more documentation is present, and more tests are done, this new version will
become the default.
	
Alhambra is designed to work with Python 3, and may fail with Python 2.

To install development versions, you can check out this github repository, and
use `pip -e` or some other method for installing Python packages.  Multiple
versions can be handled through `virtualenv`.

All Alhambra requirements should be handled through setuptools dependencies, but
this is not currently the case for xgrow and xgrowutils.

# Usage

[Documentation is available online on readthedocs.io](https://alhambra.readthedocs.io/en/latest/).  In particular, see 
[the tutorial](https://alhambra.readthedocs.io/en/latest/tutorial.html).  It is also available in the docs/ folder.

Most user-facing functions are on the TileSet class.

# Questions

Please send any questions to Constantine Evans, at cevans@evanslabs.org or cge@dna.caltech.edu.
