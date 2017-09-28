# vim: set sw=4 ts=4
import copy
import warnings
import re

import pkg_resources
import os

from . import tiletypes
from . import seeds
from . import util
from . import seq
from .tilesets import TileSet
from . import stickyends

from peppercompiler import compiler as compiler
from peppercompiler.design import spurious_design as spurious_design
from peppercompiler import finish as finish
from peppercompiler.DNA_classes import wc


def design_set(
        tileset,
        name='tsd_temp',
        includes=[pkg_resources.resource_filename(__name__, 'peppercomps-j1')],
        energetics=None,
        stickyopts={},
        reorderopts={},
        coreopts={},
        keeptemp=False):
    """Helper function to design sets from scratch, calling the numerous parts of
    tilesetdesigner. You may want to use the tilesetdesigner shell script
    instead.

    As with other functions in tilesetdesigner, this should not clobber inputs.

    :tileset: tileset definition dictionary, or an IO object with a read
    attribute, or a filename.

    :name: base name for temporary files (default tsd_temp)

    :returns: tileset definition dictionary, with added sequences

    """
    if energetics is None:
        energetics = stickyends.DEFAULT_ENERGETICS

    if hasattr(tileset, 'read'):
        tileset = TileSet.load(tileset)
    else:
        tileset = TileSet(tileset)

    tileset.check_consistent()
    tileset_with_ends_randomorder, new_ends = stickyends.create_sequences(
        tileset, energetics=energetics, **stickyopts)
    tileset_with_ends_ordered = stickyends.reorder(
        tileset_with_ends_randomorder,
        newends=new_ends,
        energetics=energetics,
        **reorderopts)
    tileset_with_strands = create_strand_sequences(
        tileset_with_ends_ordered, name, includes=includes, **coreopts)

    if 'guards' in tileset_with_strands.keys():
        tileset_with_strands = create_guard_strand_sequences(
            tileset_with_strands)

    # FIXME: this is temporary, until we have a better way of deciding.
    if 'createseqs' in tileset_with_strands['seed'].keys():
        tileset_with_strands = create_adapter_sequences(tileset_with_strands)

    if not keeptemp:
        os.remove(name + '.fix')
        os.remove(name + '.mfe')
        os.remove(name + '.pil')
        os.remove(name + '.save')
        os.remove(name + '.seqs')
        os.remove(name + '.sys')

    tileset_with_strands.check_consistent()
    return tileset_with_strands


def create_strand_sequences(tileset,
                            basename,
                            includes=None,
                            spurious_pars="verboten_weak=1.5",
                            *options):
    """Given a tileset dictionary with sticky ends sequences, create core sequences
for tiles.
    """

    newtileset = copy.deepcopy(tileset)

    create_pepper_input_files(newtileset, basename)

    compiler.compiler(
        basename, [],
        basename + '.pil',
        basename + '.save',
        fixed_file=basename + '.fix',
        includes=includes,
        synth=True)

    spurious_design.design(
        basename,
        infilename=basename + '.pil',
        outfilename=basename + '.mfe',
        verbose=True,
        struct_orient=True,
        tempname=basename + '-temp',
        extra_pars=spurious_pars,
        findmfe=False,
        cleanup=False)

    if 'info' not in newtileset.keys():
        newtileset['info'] = {}

    with open(basename + '-temp.sp') as f:
        a = f.read()
        cdi = {}
        cdi['basename'] = basename
        cdi['score_verboten'] = float(
            re.findall(r'score_verboten\s+score\s+=\s+([+-]?[\d.,]+)', a)[1])
        cdi['score_spurious'] = float(
            re.findall(r'score_spurious\s+score\s+=\s+([+-]?[\d.,]+)', a)[1])
        cdi['score_bonds'] = float(
            re.findall(r'score_bonds\s+score\s+=\s+([+-]?[\d.,]+)', a)[1])
        cdi['score'] = float(
            re.findall(r'weighted score\s+=\s+([+-]?[\d.,]+)', a)[1])
        cdi['spurious_output'] = re.search(r"(?<=FINAL\n\n)[\w\W]+weighted.*",
                                           a, re.MULTILINE).group(0)

    newtileset['info']['core'] = cdi

    finish.finish(
        basename + '.save',
        designname=basename + '.mfe',
        seqsname=basename + '.seqs',
        strandsname=None,
        run_kin=False,
        cleanup=False,
        trials=0,
        time=0,
        temp=27,
        conc=1,
        spurious=False,
        spurious_time=0)  # FIXME: shouldn't need so many options.

    tileset_with_strands = load_pepper_output_files(newtileset, basename)

    # Ensure:
    util.merge_endlists(tileset['ends'],
                        tiletypes.endlist_from_tilelist(
                            tileset_with_strands['tiles']))  # Ends still fit
    for tile in tileset_with_strands['tiles']:
        oldtile = tiletypes.gettile(tileset, tile['name'])
        if 'fullseqs' in oldtile.keys():
            for old, new in zip(oldtile['fullseqs'], tile['fullseqs']):
                seq.merge(old, new)  # old tile sequences remain
        assert oldtile['ends'] == tile['ends']

    # Check that old end sequences remain
    util.merge_endlists(tileset['ends'], tileset_with_strands['ends'])

    return tileset_with_strands


def create_pepper_input_files(tileset, basename):
    # Are we creating adapters in Pepper?
    if seeds.seedtypes[tileset['seed']['type']].needspepper:
        seedclass = seeds.seedtypes[tileset['seed']['type']]
        createadapts = True
    else:
        createadapts = False

    fixedfile = open(basename + ".fix", 'w')
    # We first need to create a fixed sequence list/file for pepper.
    # Add fixed sticky end and adjacent tile sequences.
    for end in tileset['ends']:
        if 'fseq' not in end.keys():
            continue
        seq = end['fseq'][1:-1]
        if end['type'] == 'TD':
            adj = end['fseq'][-1]
            cadj = end['fseq'][0]  # FIXME: WAS [1], OFF BY ONE!
        elif end['type'] == 'DT':
            adj = end['fseq'][0]
            cadj = end['fseq'][-1]  # FIXME: WAS [1], OFF BY ONE!
        else:
            print("warning! end {} not recognized".format(end['name']))
        fixedfile.write(
            "signal e_{0} = {1}\n".format(end['name'], seq.upper()))
        fixedfile.write(
            "signal a_{0} = {1}\n".format(end['name'], adj.upper()))
        fixedfile.write(
            "signal c_{0} = {1}\n".format(end['name'], cadj.upper()))
        # If we are creating adapter tiles in Pepper, add origami-determined
        # sequences
    if createadapts:
        for i, core in enumerate(seedclass.cores, 1):
            fixedfile.write("signal origamicore_{0} = {1}\n".format(i, core))

    # Now we'll create the system file in parts.
    importlist = set()
    compstring = ""

    for tile in tileset['tiles']:
        e = [[], []]
        for end in tile['ends']:
            if (end == 'hp'):
                continue
                # skip hairpins, etc that aren't designed by stickydesign
            e[0].append('e_' + end.replace('/', '*'))
            if end[-1] == '/':
                a = 'c_' + end[:-1] + '*'
            else:
                a = 'a_' + end
            e[1].append(a)
        s1 = " + ".join(e[0])
        s2 = " + ".join(e[1])
        tiletype = tile['type']
        if 'extra' in tile.keys():
            tiletype += '_' + tile['extra']
        compstring += "component {} = {}: {} -> {}\n".format(
            tile['name'], tiletype, s1, s2)
        importlist.add(tiletype)
        if 'fullseqs' in tile.keys():
            fixedfile.write("structure {}-tile = ".format(tile['name']) +
                            "+".join([seq.upper()
                                      for seq in tile['fullseqs']]) + "\n")

    if createadapts:
        importlist, compstring = seedclass.create_pepper_input_files(
            tileset['seed'], importlist, compstring)

    with open(basename + '.sys', 'w') as sysfile:
        sysfile.write("declare system {}: ->\n\n".format(basename))
        sysfile.write("import " + ", ".join(importlist) + "\n\n")
        sysfile.write(compstring)


def load_pepper_output_files(tileset, basename):
    import re

    # Are we creating adapters in Pepper?
    # if seeds.seedtypes[tileset['seed']['type']].needspepper:
    #     seedclass = seeds.seedtypes[tileset['seed']['type']]
    #     createadapts = True

    tset = copy.deepcopy(tileset)

    seqsstring = open(basename + '.seqs').read()

    # FIXME: we should do more than just get these sequences. We should also
    # check that our ends, complements, adjacents, etc are still correct. But
    # this is a pretty low priority.  UPDATE for 2017: WOW, LOOK AT ME BEING AN
    # IDIOT IN THE COMMENTS - CGE

    for tile in tset['tiles']:
        pepperstrands = re.compile('strand ' + tile['name'] +
                                   '-([^ ]+) = ([^\n]+)').findall(seqsstring)
        tile['fullseqs'] = tiletypes.order_pepper_strands(pepperstrands)

    for adapter in tset['seed']['adapters']:
        pepperstrands = re.compile('strand ' + adapter['name'] +
                                   '-([^ ]+) = ([^\n]+)').findall(seqsstring)
        adapter['fullseqs'] = tiletypes.order_pepper_strands(pepperstrands)

    return tset


def create_guard_strand_sequences(tileset):
    tset = copy.deepcopy(tileset)

    for guard in tset['guards']:
        tile = tiletypes.gettile(tset, guard[0])
        guard.append(wc(tile['fullseqs'][guard[1] - 1]))

    return tset


def create_adapter_sequence_diagrams(tileset, filename, *options):
    from lxml import etree
    import pkg_resources
    import os.path

    base = etree.parse(
        pkg_resources.resource_stream(
            __name__, os.path.join('seqdiagrambases', 'blank.svg')))
    baseroot = base.getroot()
    pos = 200
    for adapterdef in tileset['seed']['adapters']:

        seedclass = seeds.seedtypes[tileset['seed']['type']]
        group = seedclass.create_adapter_sequence_diagram(adapterdef)

        group.attrib['transform'] = 'translate(0,{})'.format(pos)
        pos += 200
        baseroot.append(group)

    base.write(filename)


def create_sequence_diagrams(tileset, filename, *options):
    from lxml import etree
    import pkg_resources
    import os.path

    base = etree.parse(
        pkg_resources.resource_stream(
            __name__, os.path.join('seqdiagrambases', 'blank.svg')))
    baseroot = base.getroot()
    pos = 150
    for tiledef in tileset['tiles']:

        tile = tiletypes.tfactory.parse(tiledef)
        group = tile.sequence_diagram()

        group.attrib['transform'] = 'translate(0,{})'.format(pos)
        pos += 150
        baseroot.append(group)

    base.write(filename)


def create_adapter_sequences(tileset):
    seedclass = seeds.seedtypes[tileset['seed']['type']]
    if seedclass.needspepper:
        warnings.warn(
            "This set must have adapter sequences created during\
 regular sequence design. You can ignore this if you just created sequences."
        )
        return tileset
    return seedclass.create_adapter_sequences(tileset)


def create_abstract_diagrams(tileset, filename, *options):
    import svgwrite

    drawing = svgwrite.Drawing(filename)

    pos = 0
    lim = 10000

    for tiledef in tileset['tiles']:
        tile = tiletypes.tfactory.parse(tiledef)

        group = tile.abstract_diagram(drawing)

        group['transform'] = "translate({},{})".format((pos % lim) * 150,
                                                       pos / lim * 150)
        pos += 1

        drawing.add(group)
        pos += 1

    drawing.save()


def create_strand_order_list(tileset):
    tf = [tiletypes.tfactory.parse(t) for t in tileset['tiles']]
    return sum((x.orderableseqs for x in tf), [])


def create_layout_diagrams(tileset, xgrowarray, filename, scale=1, *options):
    import svgwrite
    from . import stxg

    b = svgwrite.Drawing(filename)

    c = b.g()

    svgtiles = {}

    for tiledef in tileset['tiles']:
        tile = tiletypes.tfactory.parse(tiledef)

        group = tile.abstract_diagram(b)
        svgtiles[tile['name']] = group

    tilelist = stxg.from_yaml_endadj(tileset)['tiles']
    tilen = [None] + [x['name'] for x in tilelist]
    firstxi = 10000
    firstyi = 10000

    for yi in range(0, xgrowarray.shape[0]):
        for xi in range(0, xgrowarray.shape[1]):
            tn = tilen[xgrowarray[yi, xi]]
            if tn and tn[-5:] == '_left':
                tn = tn[:-5]
            if tn and tn[-4:] == '_top':
                tn = tn[:-4]
            if not (tn in svgtiles.keys()):
                continue
            if xi < firstxi:
                firstxi = xi
            if yi < firstyi:
                firstyi = yi
            st = svgtiles[tn].copy()
            st['transform'] = 'translate({},{})'.format(xi * 100, yi * 100)
            c.add(st)

    c['transform'] = "translate({},{}) scale({})".format(
        -firstxi * 100 * scale, -firstyi * 100 * scale, scale)

    b.add(c)
    b.save()
