vim:textwidth=60:ft=pandoc

# Tasks

* Make sure that functions don't modify their inputs. They
  should deepcopy tilesets. (I think this is already done,
  but I'm not sure. - CE)
* Write scripts to take file inputs and do stuff with them,
  without going into a python interpreter at all. These
  should be integrated into the package/setup.py, as there
  is some way to do this.
* Implement code for generating adapter tiles and adapter
  tile diagrams, with flexibility as far as what adapters to
  make and what scaffold/staples are being used. Note that
  the current way of listing seeds isn't necessarily ideal
  for this, as we'd like to be able to generate sequences
  for a bunch of adapter tiles that give different starting
  values, and then only use a subset of them at once. Think
  about this.
* Figure out what to do about SpuriousSSM. Everything else
  should now be simple to install and configure, but this
  remains a challenge. How can we include it such that
  setup.py will compile it? I've (CE) asked about this on
  stackoverflow, and we'll see what answers we can get. An
  alternative would be to have an easy installation method
  for spuriousSSM separately, by, for example, having a
  configure script and so on, but this would still make
  things more complicated.
* Implement code for quickly and easily creating desired
  guard strands. This should be pretty trivial.
* Implement more tile types; some tile types used in the
  binary counter, for example, are not currently supported.
  This requires three things: a class (eg, in tiletypes.py),
  a comp file (in peppercomps), and an svg base file (in
  seqdiagrambases).
* Clean up old code. I don't think endreorder.py is used at
  all. tileutils.py is probably mostly or entirely replaced
  by tiletypes.py. Other bits of the included code are
  probably pretty bad.
* Add documentation.
* Add both a test suite for testing changes, and checks
  throughout to test consistency and quality of sets during
  design.
