# vim: set sw=4 ts=4
from random import shuffle
import copy
import stickydesign as sd
import logging
import warnings
import numpy as np
import re

from stickydesign.energetics_daoe import energetics_daoe
default_energetics = energetics_daoe(temperature=33, mismatchtype='combined', coaxparams=True) 

import pkg_resources
import os

from . import tiletypes
from . import seeds
from . import util
from . import seq
from . import tilesets
#from tilesets import *

seedtypes = seeds.seedtypes

from peppercompiler import compiler as compiler
from peppercompiler.design import spurious_design as spurious_design
from peppercompiler import finish as finish
from peppercompiler.DNA_classes import wc
import ruamel.yaml as yaml
from . import tilesets
        
def design_set(tileset, name='tsd_temp', includes=[pkg_resources.resource_filename(__name__,'peppercomps-j1')], energetics=None, stickyopts={}, reorderopts={}, coreopts={}, keeptemp=False):
    """
    Helper function to design sets from scratch, calling the numerous parts of tilesetdeisgner. You may want to use the tilesetdesigner
    shell script instead.

    As with other functions in tilesetdesigner, this should not clobber inputs.

    :tileset: tileset definition dictionary, or an IO object with a read attribute, or a filename.
    :name: base name for temporary files (default tsd_temp)
    :returns: tileset definition dictionary, with added sequences

    """
    if energetics == None:
        energetics = default_energetics
    
    if hasattr(tileset,'read'):
        tileset = tilesets.load_tileset_dict(tileset)
    else:
        tileset = tilesets.tileset_dict(tileset)

    tileset.check_consistent()
    tileset_with_ends_randomorder, new_ends = create_sticky_end_sequences(
        tileset, energetics=energetics, **stickyopts )
    tileset_with_ends_ordered = reorder_sticky_ends( tileset_with_ends_randomorder,
                                                     newends=new_ends,
                                                     energetics=energetics,
                                                     **reorderopts )
    tileset_with_strands = create_strand_sequences( tileset_with_ends_ordered, name,
                                                    includes = includes, **coreopts )

    if 'guards' in tileset_with_strands.keys():
        tileset_with_strands = create_guard_strand_sequences( tileset_with_strands )

    # FIXME: this is temporary, until we have a better way of deciding.
    if 'createseqs' in tileset_with_strands['seed'].keys():
        tileset_with_strands = create_adapter_sequences( tileset_with_strands )

    if not keeptemp:
        os.remove(name+'.fix')
        os.remove(name+'.mfe')
        os.remove(name+'.pil')
        os.remove(name+'.save')
        os.remove(name+'.seqs')
        os.remove(name+'.sys')

    tileset_with_strands.check_consistent()       
    return tileset_with_strands

def create_sticky_end_sequences( tileset, energetics=None, sdopts={}, _sdopts_2_function=None ):
    """\
Create sticky end sequences for a tileset, using stickydesign.  This new
version should be more flexible, and should be able to keep old sticky ends,
accounting for them in creating new ones.

Parameters
----------

tileset: the tileset to create sticky ends sequences for.  This will be copied
     and returned, not modified.

energetics: the energetics instance to use for the design.  If None (default), 
will use alhambra.designer.default_energetics.

Outputs (tileset, new_ends) where new_ends is a list of new end names that
were designed.
"""
    if not energetics:
        energetics = default_energetics
    # Steps for doing this:

    # Create a copy of the tileset.
    newtileset = copy.deepcopy(tileset)

    # Build a list of ends from the endlist in the tileset.  Do this
    # by creating a named_list, then merging them into it.
    ends = util.named_list() 

    if 'ends' in newtileset.keys():
        ends = util.merge_endlists( ends, newtileset['ends'],
                                    fail_immediate=False, in_place=True )

    # This is the endlist from the tiles themselves.
    if 'tiles' in newtileset.keys(): # maybe you just want ends?
        # this checks for end/complement usage, and whether any
        # previously-describedends are unused
        # FIXME: implement
        #tiletypes.check_end_usage(newtileset['tiles'], ends)
        
        endlist_from_tiles = tiletypes.endlist_from_tilelist(
                                  newtileset['tiles'] )

    ends = util.merge_endlists( ends, endlist_from_tiles, in_place=True )

    # Ensure that if there are any resulting completely-undefined ends, they
    # have their sequences removed.
    for end in ends:
        if 'fseq' in end.keys() and end['fseq']=='nnnnnnn':
            del(end['fseq'])
    
    # Build inputs suitable for stickydesign: lists of old sequences for TD/DT,
    # and numbers of new sequences needed.
    oldDTseqs = [ end['fseq'] for end in ends \
              if 'fseq' in end.keys() and end['type']=='DT' ]
    oldTDseqs = [ end['fseq'] for end in ends \
              if 'fseq' in end.keys() and end['type']=='TD' ]

    newTDnames = [ end['name'] for end in ends \
              if 'fseq' not in end.keys() and end['type']=='TD' ]
    newDTnames = [ end['name'] for end in ends \
              if 'fseq' not in end.keys() and end['type']=='DT' ]

    
    # Deal with energetics, considering potential old sequences.
    # FIXME: EXPLAIN WHAT THIS ABSTRUSE CODE DOES...
    # TODO: tests needs to test this
    targets = []
    if len(oldDTseqs)==0 and len(oldTDseqs)==0:
        targets.append( sd.enhist( 'DT', 5, energetics=energetics)[2]['emedian'] )
        targets.append( sd.enhist( 'TD', 5, energetics=energetics)[2]['emedian'] )
    if len(oldDTseqs)>0:
        targets.append( energetics.matching_uniform(sd.endarray(oldDTseqs,'DT')) )
    if len(oldTDseqs)>0:
        targets.append( energetics.matching_uniform(sd.endarray(oldTDseqs,'TD')) )
    targetint = np.average(targets)
    
    # Create new sequences.
    newTDseqs = sd.easyends( 'TD', 5, number=len(newTDnames),
                             energetics=energetics,
                             interaction=targetint, **sdopts).tolist()

    if _sdopts_2_function:
        sdopts = _sdopts_2_function(sdopts, newTDseqs, targetint, energetics, 'TD')
    
    newDTseqs = sd.easyends( 'DT', 5, number=len(newDTnames),
                             energetics=energetics,
                             interaction=targetint, **sdopts).tolist()

    # FIXME: move to stickydesign
    assert len(newTDseqs) == len(newTDnames)
    assert len(newDTseqs) == len(newDTnames)

    # Shuffle the lists of end sequences, to ensure that they're random order, and that ends
    # used earlier in the set are not always better than those used later.
    shuffle(newTDseqs)
    shuffle(newDTseqs)

    for name,seq in zip(newDTnames,newDTseqs):
        ends[name]['fseq'] = seq
    for name,seq in zip(newTDnames,newTDseqs):
        ends[name]['fseq'] = seq

    ends.check_consistent()

    # Ensure that the old and new sets have consistent end definitions,
    # and that the tile definitions still fit.
    tiletypes.merge_endlists( tileset['ends'], ends )
    tiletypes.merge_endlists(
        tiletypes.endlist_from_tilelist(newtileset['tiles']),
        ends )
    
    # Apply new sequences to tile system.
    newtileset['ends'] = ends

    return (newtileset, newTDnames+newDTnames)    

def reorder_sticky_ends( tileset, newends=[], hightemp=0.1, lowtemp=1e-7, steps=45000, update=1000, energetics=None ):
    """Given a tileset dictionary that includes sticky end sequences, reorder these to
    try to optimize error rates."""
    from . import endreorder
    from . import endreorder_fast
    from . import anneal

    if energetics == None:
        energetics = default_energetics
    
    tset = copy.deepcopy(tileset)

    if 'info' not in tset.keys():
        tset['info'] = {}
    
    reordersys = endreorder_fast.EndSystemFseq( tset, newends, energetics=energetics )
    
    #reordersys_old = endreorder.EndSystemFseq( tset, energetics=energetics )

    # FIXME: better parameter control here.
    annealer = anneal.Annealer( reordersys.score, reordersys.mutate )

    newstate = annealer.anneal( reordersys.initstate, hightemp, lowtemp, steps, update )

    # Now take that new state, and apply it to the new tileset.
    seqs = reordersys.slowseqs( newstate[0] )
    for end in tset['ends']:
        if end['type'] in ['DT','TD']:
            eloc = reordersys.enlocs[end['name']]
            end['fseq'] = seqs[eloc[1]].tolist()[eloc[0]]

    ri = {}

    ri['score'] = reordersys.score(newstate[0])
    
    tset['info']['reorder'] = ri
            
    # Ensure that only ends in newends moved: that all others remain mergeable:
    if newends:
        old_ends_from_new_set = util.named_list( end for end in tset['ends']
                                                 if end['name'] not in newends )
        util.merge_endlists(tileset['ends'],old_ends_from_new_set)
        
    # Ensure system consistency
    tset.check_consistent()
    return tset

def create_strand_sequences( tileset, basename, includes=None, spurious_pars="verboten_weak=1.5", *options ):
    """Given a tileset dictionary with sticky ends sequences, create core sequences for tiles."""

    newtileset = copy.deepcopy(tileset)

    create_pepper_input_files( newtileset, basename )

    compiler.compiler( basename, [], basename+'.pil', basename+'.save', \
            fixed_file=basename+'.fix', includes=includes, synth=True )

    spurious_design.design( basename, infilename=basename+'.pil', outfilename=basename+'.mfe',
                            verbose=True, struct_orient=True, tempname=basename+'-temp',
                            extra_pars=spurious_pars, findmfe=False, cleanup=False ) 

    if 'info' not in newtileset.keys():
        newtileset['info'] = {}

    with open(basename+'-temp.sp') as f:
        a = f.read()
        cdi = {}
        cdi['basename'] = basename
        cdi['score_verboten'] = float(
            re.findall(r'score_verboten\s+score\s+=\s+([+-]?[\d.,]+)',a)[1])
        cdi['score_spurious'] = float(
            re.findall(r'score_spurious\s+score\s+=\s+([+-]?[\d.,]+)',a)[1])
        cdi['score_bonds'] = float(
            re.findall(r'score_bonds\s+score\s+=\s+([+-]?[\d.,]+)',a)[1])
        cdi['score'] = float(
            re.findall(r'weighted score\s+=\s+([+-]?[\d.,]+)',a)[1])
        cdi['spurious_output'] = re.search(
            r"(?<=FINAL\n\n)[\w\W]+weighted.*",a,re.MULTILINE).group(0)

    newtileset['info']['core']=cdi
        
    finish.finish( basename+'.save', designname=basename+'.mfe', seqsname=basename+'.seqs',
            strandsname=None, run_kin=False, cleanup=False, trials=0, time=0, temp=27, conc=1,
            spurious=False, spurious_time=0 ) #FIXME: shouldn't need so many options.

    tileset_with_strands = load_pepper_output_files( newtileset, basename )

    # Ensure:
    util.merge_endlists( tileset['ends'],
                         tiletypes.endlist_from_tilelist(tileset_with_strands['tiles']) )  # Ends still fit
    for tile in tileset_with_strands['tiles']:
        oldtile = tiletypes.gettile( tileset, tile['name'] )
        if 'fullseqs' in oldtile.keys():
            for old,new in zip(oldtile['fullseqs'],tile['fullseqs']):
                seq.merge(old,new)  # old tile sequences remain
        assert oldtile['ends'] == tile['ends']
    
    # Check that old end sequences remain
    util.merge_endlists( tileset['ends'], tileset_with_strands['ends'] ) 
    
    return tileset_with_strands
    
def create_pepper_input_files( tileset, basename ):
    # Are we creating adapters in Pepper?
    if seeds.seedtypes[tileset['seed']['type']].needspepper:
        seedclass = seeds.seedtypes[tileset['seed']['type']]
        createadapts = True
    else:
        createadapts = False

    fixedfile = open(basename+".fix",'w')
    # We first need to create a fixed sequence list/file for pepper.
        # Add fixed sticky end and adjacent tile sequences.
    for end in tileset['ends']:
        if 'fseq' not in end.keys():
            continue
        seq = end['fseq'][1:-1]
        if end['type'] == 'TD': 
            adj = end['fseq'][-1]
            cadj = end['fseq'][0] # FIXME: WAS [1], OFF BY ONE!
        elif end['type'] == 'DT':        
            adj = end['fseq'][0]
            cadj = end['fseq'][-1] # FIXME: WAS [1], OFF BY ONE!
        else:
            print("warning! end {} not recognized".format(end['name']))
        fixedfile.write( "signal e_{0} = {1}\n".format( end['name'], seq.upper() ) ) 
        fixedfile.write( "signal a_{0} = {1}\n".format( end['name'], adj.upper() ) ) 
        fixedfile.write( "signal c_{0} = {1}\n".format( end['name'], cadj.upper() ) )
        # If we are creating adapter tiles in Pepper, add origami-determined sequences
    if createadapts:
        for i,core in enumerate(seedclass.cores, 1):
            fixedfile.write( "signal origamicore_{0} = {1}\n".format(i,core) )

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
        if 'fullseqs' in tile.keys():
            fixedfile.write( "structure {}-tile = ".format(tile['name'])+"+".join(
                [seq.upper() for seq in tile['fullseqs']]) + "\n" )

   
    if createadapts:
        importlist, compstring = seedclass.create_pepper_input_files( tileset['seed'], importlist, compstring )

    with open(basename+'.sys','w') as sysfile:
        sysfile.write("declare system {}: ->\n\n".format(basename))
        sysfile.write("import "+", ".join(importlist)+"\n\n")
        sysfile.write(compstring)


def load_pepper_output_files( tileset, basename ):
    import re
    
    # Are we creating adapters in Pepper?
    if seeds.seedtypes[tileset['seed']['type']].needspepper:
        seedclass = seeds.seedtypes[tileset['seed']['type']]
        createadapts = True

    tset = copy.deepcopy(tileset)

    seqsstring = open(basename+'.seqs').read()

    # FIXME: we should do more than just get these sequences. We should also check that our
    # ends, complements, adjacents, etc are still correct. But this is a pretty low priority.
    # UPDATE for 2017: WOW, LOOK AT ME BEING AN IDIOT IN THE COMMENTS - CGE

    for tile in tset['tiles']:
        pepperstrands = re.compile('strand '+tile['name']+'-([^ ]+) = ([^\n]+)').findall(seqsstring)
        tile['fullseqs'] = tiletypes.order_pepper_strands(pepperstrands)
    
    for adapter in tset['seed']['adapters']:
        pepperstrands = re.compile('strand '+adapter['name']+'-([^ ]+) = ([^\n]+)').findall(seqsstring)
        adapter['fullseqs'] = tiletypes.order_pepper_strands(pepperstrands)

    return tset

def create_guard_strand_sequences( tileset ):
    tset = copy.deepcopy(tileset)


    for guard in tset['guards']:
        tile = tiletypes.gettile( tset, guard[0] )
        guard.append(  wc(tile['fullseqs'][guard[1]-1]) )
    
    return tset

def create_adapter_sequence_diagrams( tileset, filename, *options ):
    from lxml import etree
    import pkg_resources
    import os.path

    base = etree.parse( \
            pkg_resources.resource_stream(__name__,os.path.join('seqdiagrambases','blank.svg')) \
            )
    baseroot = base.getroot()
    pos = 200
    for adapterdef in tileset['seed']['adapters']:

        seedclass = seedtypes[ tileset['seed']['type'] ]
        group = seedclass.create_adapter_sequence_diagram( adapterdef )

        group.attrib['transform'] = 'translate(0,{})'.format(pos)
        pos+=200
        baseroot.append(group)

    base.write( filename )

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

        tile = tiletypes.tfactory.parse(tiledef)
        group = tile.sequence_diagram()

        group.attrib['transform'] = 'translate(0,{})'.format(pos)
        pos+=150
        baseroot.append(group)

    base.write( filename )

def create_adapter_sequences( tileset ):
    seedclass = seedtypes[ tileset['seed']['type'] ]
    if seedclass.needspepper:
        warnings.warn("This set must have adapter sequences created during regular sequence design. You can ignore this if you just created sequences.")
        return tileset
    return seedclass.create_adapter_sequences( tileset )
    

def create_abstract_diagrams( tileset, filename, *options ):
    import svgwrite
    
    drawing = svgwrite.Drawing( filename )

    pos = 0
    lim = 10000
    
    for tiledef in tileset['tiles']:
        tile = tiletypes.tfactory.parse(tiledef)
        
        group = tile.abstract_diagram(drawing)

        group['transform'] = "translate({},{})".format((pos%lim)*150,pos/lim*150)
        pos+=1
                
        drawing.add( group )
        pos += 1        
    
    drawing.save()

def create_layout_diagrams( tileset, xgrowarray, filename, scale=1, *options ):
    import svgwrite
    from . import stxg

    b = svgwrite.Drawing( filename )

    c = b.g()

    svgtiles = {}
    
    for tiledef in tileset['tiles']:
        tile = tiletypes.tfactory.parse(tiledef)

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

