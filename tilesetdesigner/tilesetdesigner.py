# vim: set sw=4 ts=4
from random import shuffle
import copy
import tileutils
import stickydesign as sd
import logging


# Color dictionary for xgrow colors...
import pkg_resources
import os.path
rgbv = pkg_resources.resource_stream(__name__, os.path.join('data','rgb.txt'))
xcolors={ " ".join(y[3:]): "rgb({},{},{})".format(y[0],y[1],y[2]) for y in [x.split() for x in rgbv] }

import tiletypes

tfactory = tiletypes.TileFactory()

def design_set(tileset, name, includes=[pkg_resources.resource_filename(__name__,'peppercomps')]):
    """DO EVERYTHING

    :tileset: TODO
    :name: TODO
    :returns: TODO

    """
    import yaml 
    # If tileset is a filename, open and load it. If it's a file, open it.
    # If it's a dict / something else, copy and use it.

    if hasattr(tileset,'read'):
        tileset = yaml.load(tileset)
    elif isinstance(tileset, basestring):
        with open(tileset,'r') as f:
            tileset = yaml.load(f)
    else:
        tileset = copy.deepcopy(tileset)

    # Make some sequences for the sticky ends. But we probably don't need
    # to save the non-reordered version!
    tileset_with_ends_randomorder = create_sticky_end_sequences( tileset )
    
    # Now reorder them.
    tileset_with_ends_ordered = reorder_sticky_ends( tileset_with_ends_randomorder )

    # Now create the strands.
    tileset_with_strands = create_strand_sequences( tileset_with_ends_ordered, name, includes = includes )

    # Now do some output.
    return tileset_with_strands


def create_sticky_end_sequences( tileset, inputs='complements', *options ):
    """Given a tileset dictionary in the usual format, without sticky ends, create sticky
    end sequences via stickydesign, and add them *in random order* to the tileset."""

    # Create a copy of the tileset, so that the original is unchanged.
    tset = copy.deepcopy(tileset)

    # To start, we'll create two sets to hold ends of each type. For each tile,
    # we'll then add the ends to the sets.
    endsets = { 'DT': set(), 'TD': set() }

    # Go through each tile, adding the ends to the sets, and making sure that
    # the ends were not previously added to the other set.
    for tile in tset['tiles']:
        ends = tfactory(tile).tile_ends
        
        for end in ends:
            if end[1] == 'TD':
                endsets['TD'].add( end[0] )
                if end[0] in endsets['DT'] or (tileutils.compname(end[0]) in
                    endsets['DT']):
                    raise ValueError("end {} in {} is TD, but was previously \
                    DT.".format(end[0],tile['name']))
            elif end[1] == 'DT':
                endsets['DT'].add( end[0] )
                if (end[0] in endsets['TD']) or (tileutils.compname(end[0]) in
                    endsets['TD']):
                    raise ValueError("end {} in {} is DT, but was previously \
                    TD.".format(end[0],tile['name']))
            else:
                pass
                # For now, we don't do anything for other end types, which will
                # primarily be fake/hairpins

    # Check to ensure that every end is also used as a complement:
    onetd = set( x[:-1] for x in endsets['TD'] if x[-1]=='/' ) ^ set( x for x in endsets['TD'] if x[-1]!='/' )
    onedt = set( x[:-1] for x in endsets['DT'] if x[-1]=='/' ) ^ set( x for x in endsets['DT'] if x[-1]!='/' )
    if onetd or onedt:
        warnings.warn("FIXME: IMPLEMENT THIS WARNING MESSAGE (END/COMPLEMENT \
        MISSING")

    # Make all ends into non-complements, adding the end name if only the
    # complement was used.
    endsets['TD'] = set( x[:-1] for x in endsets['TD'] if x[-1]=='/' ) | set( x for x in endsets['TD'] if x[-1]!='/' )
    endsets['DT'] = set( x[:-1] for x in endsets['DT'] if x[-1]=='/' ) | set( x for x in endsets['DT'] if x[-1]!='/' )

    # The number of ends of each type to design. FIXME: should deal with
    # preexisting ends in the file.
    numTD = len(endsets['TD'])
    numDT = len(endsets['DT'])

    # Create the ends with stickydesign.
    ends = {'TD': [], 'DT': []}
    ends['TD'] = sd.easyends('TD',5,number=numTD).tolist()
    ends['DT'] = sd.easyends('DT',5,number=numDT).tolist()

    # Check that we found a sufficient number of ends. FIXME: is this necessary,
    # or will stickydesign raise an error?
    assert numTD == len(ends['TD'])
    assert numDT == len(ends['DT'])

    # Shuffle the lists of end sequences, to ensure that they're random order, and that ends
    # used earlier in the set are not always better than those used later.
    shuffle(ends['TD'])
    shuffle(ends['DT'])

    if 'ends' not in tset.keys():
        tset['ends'] = []

    # Create the list of ends and their sequences. FIXME: deal with preexisting
    # ends.
    for endtype in ['TD','DT']:
        for end in endsets[endtype]:
            tset['ends'].append( { 'name': end, 'type': endtype, 'fseq': \
                ends[endtype].pop() } )
   
    # If all inputs to tiles are complements, then make the input array for each
    # tile.
    if inputs == 'complements':
        for tile in tset['tiles']:
            tile['input'] = [ int(x[-1]=='/') for x in tile['ends'] ]

    return tset

def reorder_sticky_ends( tileset, hightemp=0.1, lowtemp=1e-7, steps=15000, update=1000, *options ):
    import endreorder2
    import anneal

    tset = copy.deepcopy(tileset)

    reordersys = endreorder2.EndSystemFseq( tset )

    # FIXME: better parameter control here.
    annealer = anneal.Annealer( reordersys.score, reordersys.mutate )

    newstate = annealer.anneal( reordersys.initstate, hightemp, lowtemp, steps, update )

    # Now take that new state, and apply it to the new tileset.
    for end in tset['ends']:
        if end['type'] in ['DT','TD']:
            eloc = reordersys.enlocs[end['name']]
            end['fseq'] = newstate[0].seqs[eloc[1]].tolist()[eloc[0]]

    return tset

def create_strand_sequences( tileset, basename, includes=None, *options ):
    import DNACircuitCompiler.compiler as compiler
    import DNACircuitCompiler.design.spurious_design as spurious_design
    import DNACircuitCompiler.finish as finish

    tileset = copy.deepcopy(tileset)

    create_pepper_input_files( tileset, basename )

    compiler.compiler( basename, [], basename+'.pil', basename+'.save', \
            fixed_file=basename+'.fix', includes=includes, synth=True )

    spurious_design.design( basename, infilename=basename+'.pil', outfilename=basename+'.mfe',
            verbose=True, struct_orient=True, tempname=basename+'-temp',
            findmfe=False) 
    # FIXME: this is a debugging test to make running faster. Fix.

    finish.finish( basename+'.save', designname=basename+'.mfe', seqsname=basename+'.seqs',
            strandsname=None, run_kin=False, cleanup=False, trials=0, time=0, temp=27, conc=1,
            spurious=False, spurious_time=0 ) #FIXME: shouldn't need so many options.

    tileset_with_strands = load_pepper_output_files( tileset, basename )

    return tileset_with_strands
    

def create_pepper_input_files( tileset, basename ):
    # We first need to create a fixed sequence list/file for pepper.
    with open(basename+'.fix','w') as fixedfile:
        for end in tileset['ends']:
            seq = end['fseq'][1:-1]
            if end['type'] == 'TD': 
                adj = end['fseq'][-1]
                cadj = end['fseq'][1]
            elif end['type'] == 'DT':        
                adj = end['fseq'][-1]
                cadj = end['fseq'][1]
            else:
                print "warning! end {} not recognized".format(end['name'])
            fixedfile.write( "signal e_{0} = {1}\n".format( end['name'], seq.upper() ) ) 
            fixedfile.write( "signal a_{0} = {1}\n".format( end['name'], adj.upper() ) ) 
            fixedfile.write( "signal c_{0} = {1}\n".format( end['name'], cadj.upper() ) )

  # Now we'll create the system file in parts.
    importlist = set()
    compstring = ""
    
    for tile in tileset['tiles']:
        e = [[],[]]
        for end in tile['ends']:
            if (end == 'hp'):
                continue # skip hairpins, etc that aren't designed by stickydesign
            e[0].append( 'e_'+end.replace('/','*') )
            if end[-1]=='/':
                a = 'c_'+end[:-1]+'*'
            else:
                a = 'a_'+end
            e[1].append( a )
        s1 = " + ".join(e[0])
        s2 = " + ".join(e[1])
        tiletype = tile['type']
        if 'extra' in tile.keys():
            tiletype+='_'+tile['extra']
        compstring += "component {} = {}: {} -> {}\n".format(tile['name'],tiletype,s1,s2)
        importlist.add( tiletype )
    
    with open(basename+'.sys','w') as sysfile:
        sysfile.write("declare system {}: ->\n\n".format(basename))
        sysfile.write("import "+", ".join(importlist)+"\n\n")
        sysfile.write(compstring)


def load_pepper_output_files( tileset, basename ):
    import re

    tset = copy.deepcopy(tileset)

    seqsstring = open(basename+'.seqs').read()

    # FIXME: we should do more than just get these sequences. We should also check that our
    # ends, complements, adjacents, etc are still correct. But this is a pretty low priority.

    for tile in tset['tiles']:
        pepperstrands = re.compile('strand '+tile['name']+'-([^ ]+) = ([^\n]+)').findall(seqsstring)
        tile['fullseqs'] = tileutils.order_pepper_strands(pepperstrands)
    
    return tset

def create_sequence_diagrams( tileset, filename, *options ):
    from lxml import etree
    import pkg_resources
    import os.path

    base = etree.parse( \
            pkg_resources.resource_stream(__name__,os.path.join('seqdiagrambases','blank.svg')) \
            )
    baseroot = base.getroot()
    pos = 150
    for tiledef in tileset['tiles']:

        tile = tfactory.parse(tiledef)
        group = tile.sequence_diagram()

        group.attrib['transform'] = 'translate(0,{})'.format(pos)
        pos+=150
        baseroot.append(group)

    base.write( filename )


def create_abstract_diagrams( tileset, filename, *options ):
    import svgwrite
    
    drawing = svgwrite.Drawing( filename )

    pos = 0
    lim = 10000
    
    for tiledef in tileset['tiles']:
        tile = tfactory.parse(tiledef)
        
        group = tile.abstract_diagram(drawing)

        group['transform'] = "translate({},{})".format((pos%lim)*150,pos/lim*150)
        pos+=1
                
        drawing.add( group )
        pos += 1        
    
    drawing.save()

def create_layout_diagrams( tileset, xgrowarray, filename, scale=1, *options ):
    import svgwrite
    import stxg

    b = svgwrite.Drawing( filename )

    c = b.g()

    svgtiles = {}
    
    for tiledef in tileset['tiles']:
        tile = tfactory.parse(tiledef)

        group = tile.abstract_diagram(b)
        svgtiles[tile['name']] = group
        
    
    tilelist=stxg.from_yaml_endadj(tileset)['tiles']
    tilen = [None]+[x['name'] for x in tilelist]
    firstxi = 10000
    firstyi = 10000

    for yi in range(0,xgrowarray.shape[0]):
        for xi in range(0,xgrowarray.shape[1]):
            tn = tilen[xgrowarray[yi,xi]]
            if tn and tn[-5:]=='_left':
                tn=tn[:-5]
            if tn and tn[-4:]=='_top':
                tn=tn[:-4]
            if not (tn in svgtiles.keys()):
                continue
            if xi < firstxi: 
                firstxi = xi
            if yi < firstyi: 
                firstyi = yi
            st = svgtiles[tn].copy()
            st['transform']='translate({},{})'.format(xi*100,yi*100)
            c.add(st)

    c['transform'] = "translate({},{}) scale({})".format(-firstxi*100*scale,-firstyi*100*scale,scale)

    b.add(c)
    b.save()

