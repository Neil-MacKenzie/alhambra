vim:textwidth=60:ft=pandoc

# Tasks

* Clean up old code. I don't think endreorder.py is used at
  all. tileutils.py is probably mostly or entirely replaced
  by tiletypes.py. Other bits of the included code are
  probably pretty bad.
* Add documentation.
* Add both a test suite for testing changes, and checks
  throughout to test consistency and quality of sets during
  design.
* Add support for controlling adapter tiles and tile
  presence in xgrow.
* Possibly add support for running xgrow directly from the
  file, through some script? Or update xgrow_wrap to do
  this.
* Figure out how to do something about the xgrow import file
  problem.
* Allow leaving out use_adapters, in which case adapters are
  just used in order listed.
* Figure out what to do about the file format. Create
  functions for dealing with it, and possibly for
  combining/updating easily.
* Add a way to handle labels.
