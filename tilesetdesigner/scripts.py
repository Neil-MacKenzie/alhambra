from tilesetdesigner import *
import sys
import yaml

def tilesetdesigner():
    import logging
    import sys
    logging.getLogger().setLevel(logging.INFO)
    sysname = sys.argv[1]
    fname = sys.argv[2]
    out = sys.argv[3]
    sys = design_set(fname, sysname)
    yaml.dump(sys, open(out,'w'))
