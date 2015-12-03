from .designer import *
import yaml
import argparse
import logging
import sys
import os.path

def tilesetdesigner():
    import sys
    
    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", help="tileset input YAML file")
    parser.add_argument("-n", "--name", help="base name for output files")
    parser.add_argument("-o", "--out", help="output file (if not just NAME-out.yaml)")
    parser.add_argument("-r", "--reorderargs", help="extra arguments (as dict) to reorder function", default="{}")
    parser.add_argument("-s", "--spuriousargs", help="extra params (use \"s) for spuriousSSM", default="")
    parser.add_argument("-v", "--verbose", help="increase verbosity", action="store_true", default=False)
    parser.add_argument("-f", "--force", help="overwrite output file", default=False, action="store_true")
    parser.add_argument("-d", "--diagrams", help="make diagrams as well", default=False, action="store_true")
    parser.add_argument("-x", "--xgrow", help="make xgrow files as well", default=False, action="store_true")
    parser.add_argument("--test", help="run quickly for testing, making bad sequences", default=False, action="store_true")
    args = parser.parse_args()

    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)
    
    base = os.path.splitext(os.path.basename(args.inputfile))[0]
    
    if not args.name:
        args.name = base
    else:
        base = args.name
    if not args.out:
        args.out = base+'-out.yaml'

    args.reorderargs=eval(args.reorderargs) # FIXME: highly dangerous
    if args.test:
        args.reorderargs['steps']=5
        args.spuriousargs+=" bored=5"

    # Check that we're not clobbering something here for the output. FIXME: do
    # this for temp files too.
    if os.path.exists(args.out) and not args.force:
        logging.error("Output file already exists, and --force is not enabled.")
        sys.exit(1)

    for ext in ['.sys','.save','.mfe','.fix','.pil','.seqs']:
        if os.path.exists(args.name+ext) and not args.force:
            logging.error("Temporary file {} already exists, and --force is not enabled.".format(args.name+ext))
            sys.exit(1)

    if not os.path.exists(args.inputfile):
        logging.error("Input file does not exist.")
        sys.exit(1)

    fix_paths()

    sys = design_set( \
            args.inputfile, 
            args.name,
            reorderopts=args.reorderargs,
            coreopts={'spurious_pars': args.spuriousargs})

    yaml.dump(sys, open(args.out,'w'))

    if args.diagrams:
        create_abstract_diagrams( sys , base+'-abstract.svg' )
        create_sequence_diagrams( sys, base+'-sequences.svg' )
        if 'createseqs' in sys['seed'].keys():
            create_adapter_sequence_diagrams( sys, base+'-adapters.svg' )

    if args.xgrow:
        from . import stxg
        yaml.dump( stxg.from_yaml_endadj( sys, perfect=True ), 
                   open(base+'-perfectends.stxg', 'w') )
        yaml.dump( stxg.from_yaml_endadj( sys, perfect=True, rotate=True ), 
                   open(base+'-perfectends-rotated.stxg', 'w') )
        yaml.dump( stxg.from_yaml_endadj( sys, perfect=False ), 
                   open(base+'-withends.stxg', 'w') )
        yaml.dump( stxg.from_yaml_endadj( sys, perfect=False, rotate=True ), 
                   open(base+'-withends-rotated.stxg', 'w') )
