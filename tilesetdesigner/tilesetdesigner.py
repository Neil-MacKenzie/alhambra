from random import shuffle
import copy
import tileutils
import stickydesign as sd


# Color dictionary for xgrow colors...
import pkg_resources
import os.path
rgbv = pkg_resources.resource_stream(__name__, os.path.join('data','rgb.txt'))
xcolors={ " ".join(y[3:]): "rgb({},{},{})".format(y[0],y[1],y[2]) for y in [x.split() for x in rgbv] }

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
        ends = tileutils.get_tile_ends(tile)
        
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

def reorder_sticky_ends( tileset, *options ):
    import endreorder2
    import anneal

    tset = copy.deepcopy(tileset)

    reordersys = endreorder2.EndSystemFseq( tset )

    # FIXME: better parameter control here.
    annealer = anneal.Annealer( reordersys.score, reordersys.mutate )

    newstate = annealer.anneal( reordersys.initstate, 0.1, 1e-8, 4000, 200 )

    # Now take that new state, and apply it to the new tileset.
    for end in tset['ends']:
        if end['type'] in ['DT','TD']:
            eloc = b.enlocs[end['name']]
            end['fseq'] = newstate[0].seqs[eloc[1]].tolist()[eloc[0]]

    return tset

def create_core_sequences( tileset, *options ):
    pass

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
            if (end['type'] != 'DT') and (end['type'] != 'TD'):
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

def create_abstract_diagrams( tileset, filename, *options ):
    import svgwrite
    
    

    b = svgwrite.Drawing( filename )

    pos = 0
    lim = 10000
    
    for tile in tileset['tiles']:
        a = b.g(transform="translate({},{})".format((pos%lim)*150,pos/lim*150))
        
        if tile['type'] in ['tile_daoe_3up','tile_daoe_5up','tile_single']:
            # Tile Box
            a.add( b.rect((100,100),(100,100),stroke="black",fill=xcolors[tile['color']]) )

            # End Names
            a.add( b.text( tile['ends'][0], insert=(150,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][1], insert=(190,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,150)", font_size="14pt") )
            a.add( b.text( tile['ends'][2], insert=(150,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][3], insert=(110,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,150)", font_size="14pt") ) 

            # Tile Name
            a.add( b.text( tile['name'], insert=(150,150), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            
            orient = None
            if tile['type'] == 'tile_daoe_3up':
                orient = ('3','5')
            elif tile['type'] == 'tile_daoe_5up':
                orient = ('5','3')
            
            if orient:
                # Orientation
                a.add( b.text( orient[0], insert=(192,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(108,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))

        if tile['type'] in ['tile_daoe_doublehoriz_35up', 'tile_doublehoriz']:
            # Tile Box
            a.add( b.rect((100,100),(200,100),stroke="black",fill=xcolors[tile['color']]) )

            # End Names
            a.add( b.text( tile['ends'][0], insert=(150,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][1], insert=(250,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][2], insert=(290,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,290,150)", font_size="14pt") )
            a.add( b.text( tile['ends'][3], insert=(250,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][4], insert=(150,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][5], insert=(110,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,150)", font_size="14pt") ) 

            # Tile Name
            a.add( b.text( tile['name'], insert=(200,150), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            
            orient = None
            if tile['type'] == 'tile_daoe_doublehoriz_35up':
                orient = ('3','5')
            elif tile['type'] == 'tile_daoe_5up':
                orient = ('5','3')
            
            if orient:
                # Orientation
                a.add( b.text( orient[0], insert=(192,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(292,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(108,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[0], insert=(208,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            pos += 1
            
        if tile['type'] in ['tile_daoe_doublevert_35up', 'tile_doublevert']:
            # Tile Box
            a.add( b.rect((100,100),(100,200),stroke="black",fill=xcolors[tile['color']]) )

            # End Names
            a.add( b.text( tile['ends'][0], insert=(150,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][1], insert=(190,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,150)", font_size="14pt") )
            a.add( b.text( tile['ends'][2], insert=(190,250), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,250)", font_size="14pt") )
            a.add( b.text( tile['ends'][3], insert=(150,290), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][4], insert=(110,250), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,250)", font_size="14pt") ) 
            a.add( b.text( tile['ends'][5], insert=(110,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,150)", font_size="14pt") ) 

            
            # Tile Name
            a.add( b.text( tile['name'], insert=(150,200), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            
            orient = None
            if tile['type'] == 'tile_daoe_doublevert_35up':
                orient = ('3','5')
            elif tile['type'] == 'tile_daoe_5up':
                orient = ('5','3')
            
            if orient:
                # Orientation
                a.add( b.text( orient[0], insert=(192,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[0], insert=(108,292), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(192,208), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(108,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            pos += 1
                
        b.add(a)
        pos += 1        
    b.save()

def create_layout_diagrams( tileset, xgrowarray, filename, scale=1, *options ):
    import svgwrite
    import stxg

    b = svgwrite.Drawing( filename )

    c = b.g()

    svgtiles = {}

    for tile in tileset['tiles']:
        a = b.g()
        
        if tile['type'] in ['tile_daoe_3up','tile_daoe_5up','tile_single']:
            # Tile Box
            a.add( b.rect((100,100),(100,100),stroke="black",fill=xcolors[tile['color']]) )

            # End Names
            a.add( b.text( tile['ends'][0], insert=(150,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][1], insert=(190,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,150)", font_size="14pt") )
            a.add( b.text( tile['ends'][2], insert=(150,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][3], insert=(110,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,150)", font_size="14pt") ) 

            # Tile Name
            a.add( b.text( tile['name'], insert=(150,150), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            
            orient = None
            if tile['type'] == 'tile_daoe_3up':
                orient = ('3','5')
            elif tile['type'] == 'tile_daoe_5up':
                orient = ('5','3')
            
            if orient:
                # Orientation
                a.add( b.text( orient[0], insert=(192,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(108,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))

        if tile['type'] in ['tile_daoe_doublehoriz_35up', 'tile_doublehoriz']:
            # Tile Box
            a.add( b.rect((100,100),(200,100),stroke="black",fill=xcolors[tile['color']]) )

            # End Names
            a.add( b.text( tile['ends'][0], insert=(150,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][1], insert=(250,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][2], insert=(290,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,290,150)", font_size="14pt") )
            a.add( b.text( tile['ends'][3], insert=(250,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][4], insert=(150,190), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][5], insert=(110,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,150)", font_size="14pt") ) 

            # Tile Name
            a.add( b.text( tile['name'], insert=(200,150), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            
            orient = None
            if tile['type'] == 'tile_daoe_doublehoriz_35up':
                orient = ('3','5')
            elif tile['type'] == 'tile_daoe_5up':
                orient = ('5','3')
            
            if orient:
                # Orientation
                a.add( b.text( orient[0], insert=(192,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(292,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(108,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[0], insert=(208,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
            
        if tile['type'] in ['tile_daoe_doublevert_35up', 'tile_doublevert']:
            # Tile Box
            a.add( b.rect((100,100),(100,200),stroke="black",fill=xcolors[tile['color']]) )

            # End Names
            a.add( b.text( tile['ends'][0], insert=(150,110), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            a.add( b.text( tile['ends'][1], insert=(190,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,150)", font_size="14pt") )
            a.add( b.text( tile['ends'][2], insert=(190,250), text_anchor='middle', alignment_baseline="middle", transform="rotate(90,190,250)", font_size="14pt") )
            a.add( b.text( tile['ends'][3], insert=(150,290), text_anchor='middle', alignment_baseline="middle", font_size="14pt") )
            a.add( b.text( tile['ends'][4], insert=(110,250), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,250)", font_size="14pt") ) 
            a.add( b.text( tile['ends'][5], insert=(110,150), text_anchor='middle', alignment_baseline="middle", transform="rotate(-90,110,150)", font_size="14pt") ) 

            
            # Tile Name
            a.add( b.text( tile['name'], insert=(150,200), text_anchor='middle', alignment_baseline="middle", font_size="14pt" ) )
            
            orient = None
            if tile['type'] == 'tile_daoe_doublevert_35up':
                orient = ('3','5')
            elif tile['type'] == 'tile_daoe_5up':
                orient = ('5','3')
            
            if orient:
                # Orientation
                a.add( b.text( orient[0], insert=(192,108), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[0], insert=(108,292), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(192,208), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))
                a.add( b.text( orient[1], insert=(108,192), text_anchor='middle', alignment_baseline="middle", font_size="9pt"))

        svgtiles[tile['name']] = a
        
    
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


def create_layout_diagram( tileset, array, tiledict, *options ):
    pass


