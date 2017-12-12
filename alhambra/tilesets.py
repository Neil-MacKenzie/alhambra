import ruamel.yaml as yaml
from ruamel.yaml.comments import CommentedMap
from ruamel.yaml.representer import RoundTripRepresenter
import warnings

import copy
import re

import pkg_resources
import os

from .tiles import TileList
from .ends import EndList

from . import tilestructures
from . import seeds
from . import util
from . import seq

from peppercompiler import compiler as compiler
from peppercompiler.design import spurious_design as spurious_design
from peppercompiler import finish as finish
from peppercompiler.DNA_classes import wc

import numpy as np
import stickydesign as sd

import stickydesign.plots as sdplots
from collections import Counter
from stickydesign import EnergeticsDAOE
from matplotlib import pylab


import collections
from random import shuffle
from datetime import datetime, timezone

import logging

DEFAULT_ENERGETICS = sd.EnergeticsDAOE(
    temperature=33, mismatchtype='combined', coaxparams=True)

DEFAULT_MULTIMODEL_ENERGETICS = [
    sd.EnergeticsDAOE(
        temperature=33, mismatchtype='combined', coaxparams='protozanova'),
    sd.EnergeticsDAOE(
        temperature=33, mismatchtype='combined', coaxparams='pyshni'),
    sd.EnergeticsDAOE(
        temperature=33, mismatchtype='combined', coaxparams='peyret'),
    sd.EnergeticsDAOE(
        temperature=33, mismatchtype='combined', coaxparams=False)]

DEFAULT_MM_ENERGETICS_NAMES = ['Prot', 'Pysh', 'Peyr', 'None']

SELOGGER = logging.getLogger(__name__)

DEFAULT_REGION_ENERGETICS = EnergeticsDAOE(temperature=33,
                                           coaxparams=False,
                                           danglecorr=False)


class TileSet(CommentedMap):
    def __init__(self, val={}):
        CommentedMap.__init__(self, val)

        self['ends'] = EndList(self.get('ends', []))
        self['tiles'] = TileList(self.get('tiles', []))

    @classmethod
    def from_file(cls, name_or_stream, *args, **kwargs):
        # Assume a stream:
        if getattr(name_or_stream, 'read', None) is None:
            return cls(yaml.round_trip_load(open(name_or_stream, 'r'),
                                            *args, **kwargs))
        else:
            return cls(yaml.round_trip_load(name_or_stream, *args, **kwargs))

    def to_file(self, name_or_stream):
        if getattr(name_or_stream, 'read', None) is None:
            return yaml.round_trip_dump(self, open(name_or_stream, 'w'))
        else:
            return yaml.round_trip_dump(self, name_or_stream)

    def dump(self, *args, **kwargs):
        return yaml.round_trip_dump(self, *args, **kwargs)
        
    @property
    def tiles(self):
        return self['tiles']

    def ends():
        doc = """Doc string"""
        def fget(self):
            return self['ends']
    
        def fset(self, value):
            self['ends'] = value
    
        def fdel(self):
            del self['ends']
        return locals()
    ends = property(**ends())
    
    @property
    def seed(self):
        return self.get('seed', None)

    def create_abstract_diagrams(self, filename, *options):
        from lxml import etree
        import os
        import pkg_resources
        base = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join('seqdiagrambases', 'blank.svg')))
        baseroot = base.getroot()

        pos = 0
        for tile in self.tiles:
            group, n = tile.abstract_diagram(self)

            group.attrib['transform'] = "translate({},{})".format(
                (pos % 12) * 22, (pos // 12) * 22)
            pos += n

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
        for tile in tileset.tiles:
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

    def create_layout_diagrams(tileset, xgrowarray, filename, scale=1, *options):
        from lxml import etree
        base = etree.parse(
            pkg_resources.resource_stream(
                __name__, os.path.join('seqdiagrambases', 'blank.svg')))
        baseroot = base.getroot()

        svgtiles = {}

        for tile in tileset.tiles:
            group, n = tile.abstract_diagram(tileset)
            svgtiles[tile['name']] = group

        from . import xgrow
        tilelist = xgrow.generate_xgrow_dict(tileset, perfect=True)['tiles']
        tilen = [None] + [x['name'] for x in tilelist]
        firstxi = 10000
        firstyi = 10000
        import copy
        for yi in range(0, xgrowarray.shape[0]):
            for xi in range(0, xgrowarray.shape[1]):
                tn = tilen[xgrowarray[yi, xi]]
                if tn and tn[-5:] == '_left':
                    tn = tn[:-5]
                if tn and tn[-7:] == '_bottom':
                    tn = tn[:-7]
                if not (tn in svgtiles.keys()):
                    continue
                if xi < firstxi:
                    firstxi = xi
                if yi < firstyi:
                    firstyi = yi
                st = copy.deepcopy(svgtiles[tn])
                st.attrib['transform'] = 'translate({},{})'.format(
                    xi * 10, yi * 10)
                baseroot.append(st)

        base.write(filename)

        
    @property
    def strand_order_list(self):
        return [y for x in self.tiles for y in x.orderableseqs]
    
    def check_consistent(self):
        # * END LIST The end list itself must be consistent.
        # ** Each end must be of understood type
        # ** Each end must have a valid sequence or no sequence
        # ** There must be no more than one instance of each name
        # ** WARN if there are ends with no namecounts
        # * TILE LIST
        # ** each tile must be of understood type (must parse)
        # ** ends in the tile list must be consistent (must merge)
        # ** there must be no more than one tile with each name
        self.tiles.check_consistent()
        endsfromtiles = self.tiles.endlist()

        # ** WARN if any end that appears does not have a complement used or vice versa
        # ** WARN if there are tiles with no name
        # * TILE + END
        # ** The tile and end lists must merge validly
        # (checks sequences, adjacents, types, complements)
        self.ends.merge(endsfromtiles)
        
        # ** WARN if tilelist has end references not in ends
        # ** WARN if merge is not equal to the endlist
        # ** WARN if endlist has ends not used in tilelist
        # * ADAPTERS / SEEDS
        if 'seed' in self.keys():
            # ** seeds must be of understood type
            assert self['seed']['type'] in seeds.seedtypes.keys()
            # ** adapter locations must be valid
            sclass = seeds.seedtypes[self['seed']['type']]

            sclass.check_consistent(self)

            sclass.check_sequence(self)
            # ** each adapter must have no sequence or a consistent sequence
            # *** the RH strand must match the associated tile
            # *** the ends in the sequence must match the ends in the endlist
            # *** the LH sequence must be validly binding to both RH and
            #     origami
            # ** each adapter must have valid definition, which means for us:
            # *** if both tile mimic and ends are specified, they must match
    
    def summary(self):
        self.check_consistent()
        info = {'ntiles': len(self.tiles),
                'nends':  len(self.ends),
                'ntends': len(tilestructures.endlist_from_tilelist(self.tiles)),
                'tns':    " ".join(x['name'] for x in self.tiles if 'name' in x.keys()),
                'ens':    " ".join(x['name'] for x in self.ends if 'name' in x.keys()),
                'name':   " {}".format(self['info']['name']) if \
                                       ('info' in self.keys() and \
                                        'name' in self['info'].keys()) else ""}
        tun = sum( 1 for x in self.tiles if 'name' not in x.keys() )
        if tun > 0:
            info['tns'] += " ({} unnamed)".format(tun)
        eun = sum( 1 for x in self.ends if 'name' not in x.keys() )
        if eun > 0:
            info['ens'] += " ({} unnamed)".format(eun)
        return "TileSet{name}: {ntiles} tiles, {nends} ends, {ntends} ends in tiles.\nTiles: {tns}\nEnds:  {ens}".format(**info)

    def __str__( self ):
        return self.summary()

    def copy(self):
        return copy.deepcopy(self)

    def __deepcopy__(self, memo):
        # type: (Any) -> Any
        res = self.__class__()
        memo[id(self)] = res
        for k in self:
            res[k] = copy.deepcopy(self[k])
            self.copy_attributes(res, deep=True)
        return res

    def dump(self, stream):
        return yaml.round_trip_dump(self, stream)

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
            energetics = DEFAULT_ENERGETICS

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

    def create_end_sequences(tileset, method='default', energetics=None,
                             trials=100, sdopts={}, ecpars={}):
        """Create sticky end sequences for a tileset, using stickydesign.
    This new version should be more flexible, and should be able to
    keep old sticky ends, accounting for them in creating new ones.

    Parameters
    ----------

    tileset: the tileset to create sticky ends sequences for.  This
         will be copied and returned, not modified.

    method: [default 'default'] if 'default', use the default,
    single-model sequence design.  If 'multimodel', use multimodel end
    choice.

    energetics: the energetics instance to use for the design, or list
    of energetics for method='multimodel', in which case the first
    will be the primary.  If None (default), will use
    alhambra.DEFAULT_ENERGETICS, or, if method='multimodel', will use
    alhambra.DEFAULT_MM_ENERGETICS.


    Outputs (tileset, new_ends) where new_ends is a list of new end
    names that were designed.

        """
        info = {}
        info['method'] = method
        info['time'] = datetime.now(tz=timezone.utc).isoformat()
        info['sd_version'] = sd.version.__version__

        if not energetics:
            if method == 'multimodel':
                energetics = DEFAULT_MULTIMODEL_ENERGETICS
            else:
                energetics = DEFAULT_ENERGETICS
        if method == 'multimodel' and not isinstance(energetics,
                                                     collections.Iterable):
            raise ValueError("Energetics must be an iterable for multimodel.")
        elif method == 'multimodel':
            all_energetics = energetics
            energetics = all_energetics[0]
            info['energetics'] = [str(e) for e in all_energetics]
            info['trails'] = trials
        elif method == 'default':
            info['energetics'] = str(energetics)

        # Steps for doing this:

        # Create a copy of the tileset.
        newtileset = tileset.copy()

        # Build a list of ends from the endlist in the tileset.  Do this
        # by creating a NamedList, then merging them into it.
        ends = EndList()

        if newtileset.ends:
            ends.merge(newtileset.ends,
                       fail_immediate=False, in_place=True)

        # This is the endlist from the tiles themselves.
        if newtileset.tiles:  # maybe you just want ends?
            # this checks for end/complement usage, and whether any
            # previously-describedends are unused
            # FIXME: implement
            # tilestructures.check_end_usage(newtileset.tiles, ends)

            endlist_from_tiles = newtileset.tiles.endlist()

        ends.merge(endlist_from_tiles, in_place=True)

        # Ensure that if there are any resulting completely-undefined ends, they
        # have their sequences removed.
        #for end in ends:
        #    if end.fseq and set(end.fseq) == {'n'}:
        #        del(end.fseq)

        # Build inputs suitable for stickydesign: lists of old sequences for TD/DT,
        # and numbers of new sequences needed.
        oldDTseqs = [
            end.fseq for end in ends
            if end.etype == 'DT' and seq.is_definite(end.fseq)
        ]
        oldDTarray = sd.endarray(oldDTseqs, 'DT')
        oldTDseqs = [
            end.fseq for end in ends
            if end.etype == 'TD' and seq.is_definite(end.fseq)
        ]
        oldTDarray = sd.endarray(oldTDseqs, 'TD')

        newTD = [
            end for end in ends
            if end.etype == 'TD' and not seq.is_definite(end.fseq)
        ]
        newDT = [
            end for end in ends
            if end.etype == 'DT' and not seq.is_definite(end.fseq)
        ]

        # Deal with energetics, considering potential old sequences.
        # FIXME: EXPLAIN WHAT THIS ABSTRUSE CODE DOES...
        # TODO: tests needs to test this
        targets = []
        if len(oldDTseqs) == 0 and len(oldTDseqs) == 0:
            targets.append(sd.enhist('DT', 5, energetics=energetics)[2]['emedian'])
            targets.append(sd.enhist('TD', 5, energetics=energetics)[2]['emedian'])
        if len(oldDTseqs) > 0:
            targets.append(
                energetics.matching_uniform(oldDTarray))
        if len(oldTDseqs) > 0:
            targets.append(
                energetics.matching_uniform(oldTDarray))
        targetint = np.average(targets)

        if any(not seq.is_null(end.fseq) for end in newTD):
            TDtemplates = [end.fseq for end in newTD]
        else:
            TDtemplates = None
        if any(not seq.is_null(end.fseq) for end in newDT):
            DTtemplates = [end.fseq for end in newDT]
        else:
            DTtemplates = None

        if method == 'default':
            if TDtemplates or DTtemplates:
                raise NotImplementedError
            # Create new sequences.
            newTDseqs = sd.easyends(
                'TD',
                5,
                number=len(newTD),
                energetics=energetics,
                interaction=targetint,
                **sdopts).tolist()

            newDTseqs = sd.easyends(
                'DT',
                5,
                number=len(newDT),
                energetics=energetics,
                interaction=targetint,
                **sdopts).tolist()

        elif method == 'multimodel':
            SELOGGER.info(
                "starting multimodel sticky end generation " +
                "of TD ends for {} DT and {} TD ends, {} trials.".format(
                    len(newDT), len(newTD), trials))

            newTDseqs = []
            pl = util.ProgressLogger(SELOGGER, trials*2)
            for i in range(0, trials):
                endchooserTD = sd.multimodel.endchooser(all_energetics,
                                                        templates=TDtemplates,
                                                        **ecpars)

                newTDseqs.append(
                    sd.easyends(
                        'TD',
                        5,
                        number=len(newTD),
                        oldends=oldTDseqs,
                        energetics=energetics,
                        interaction=targetint,
                        echoose=endchooserTD,
                        **sdopts))
                pl.update(i)

            if oldTDseqs:
                tvals = [[e.matching_uniform(oldTDarray[0:1])
                         for e in all_energetics] *
                         len(newTDseqs)] * len(newTDseqs)
            else:
                tvals = [[e.matching_uniform(x[0:1])
                          for e in all_energetics] for x in newTDseqs]

            endchoosersDT = [sd.multimodel.endchooser(all_energetics,
                                                      target_vals=tval,
                                                      templates=DTtemplates,
                                                      **ecpars)
                             for tval in tvals]

            SELOGGER.info("generating corresponding DT ends")
            newDTseqs = []
            for i, echoose in enumerate(endchoosersDT):
                newDTseqs.append(
                    sd.easyends(
                        'DT',
                        5,
                        number=len(newDT),
                        oldends=oldDTseqs,
                        energetics=energetics,
                        interaction=targetint,
                        echoose=echoose,
                        **sdopts))
                pl.update(i+trials)

            TDarr = [[oldTDarray.concat(x),
                      oldDTarray.concat(y)]
                     for x, y in zip(newTDseqs, newDTseqs)]

            scores = [sd.multimodel.deviation_score(list(e), all_energetics)
                      for e in TDarr]

            sort = np.argsort(scores)

            newTDseqs = newTDseqs[sort[0]].tolist()[len(oldTDseqs):]
            newDTseqs = newDTseqs[sort[0]].tolist()[len(oldDTseqs):]
            info['score'] = scores[sort[0]]
            info['maxscore'] = scores[sort[-1]]
            info['meanscore'] = np.mean(scores)

        # FIXME: move to stickydesign
        assert len(newTDseqs) == len(newTD)
        assert len(newDTseqs) == len(newDT)

        # Shuffle the lists of end sequences, to ensure that they're
        # random order, and that ends used earlier in the set are not
        # always better than those used later. But only shuffle if
        # there were no templates:
        if not TDtemplates:
            shuffle(newTDseqs)
        if not DTtemplates:
            shuffle(newDTseqs)

        # Make sure things are consistent if there are templates:
        if TDtemplates:
            for t, s in zip(TDtemplates, newTDseqs):
                seq.merge(t, s)
        if DTtemplates:
            for t, s in zip(DTtemplates, newDTseqs):
                seq.merge(t, s)
            
        for end, s in zip(newDT, newDTseqs):
            ends[end.name].fseq = s
        for end, s in zip(newTD, newTDseqs):
            ends[end.name].fseq = s

        ends.check_consistent()

        # Ensure that the old and new sets have consistent end definitions,
        # and that the tile definitions still fit.
        tileset.ends.merge(ends)
        newtileset.tiles.endlist().merge(ends)

        # Apply new sequences to tile system.
        newtileset.ends = ends
        if 'info' not in newtileset.keys():
            newtileset['info'] = {}
        if 'end_design' not in newtileset['info'].keys():
            newtileset['info']['end_design'] = []
        if isinstance('end_design', dict):  # convert old
            newtileset['info']['end_design'] = [newtileset['info']['end_design']]
        newtileset['info']['end_design'].append(info)

        return (newtileset, [e.name for e in newTD] + [e.name for e in newDT])

    def reorder_ends(tileset,
                     newends=[],
                     hightemp=0.1,
                     lowtemp=1e-7,
                     steps=45000,
                     update=1000,
                     energetics=None):
        """Given a tileset dictionary that includes sticky end sequences, reorder these
        to try to optimize error rates.
        """
        from . import endreorder
        from . import anneal

        if energetics is None:
            energetics = DEFAULT_ENERGETICS

        tset = tileset.copy()

        if 'info' not in tset.keys():
            tset['info'] = {}

        reordersys = endreorder.EndSystemFseq(
            tset, newends, energetics=energetics)

        # FIXME: better parameter control here.
        annealer = anneal.Annealer(reordersys.score, reordersys.mutate)

        newstate = annealer.anneal(reordersys.initstate, hightemp, lowtemp, steps,
                                   update)

        # Now take that new state, and apply it to the new tileset.
        seqs = reordersys.slowseqs(newstate[0])
        for end in tset.ends:
            if end.etype in ['DT', 'TD']:
                eloc = reordersys.enlocs[end['name']]
                end.fseq = seqs[eloc[1]].tolist()[eloc[0]]

        ri = {}

        ri['score'] = reordersys.score(newstate[0])

        tset['info']['reorder'] = ri

        # Ensure that only ends in newends moved: that all others remain mergeable:
        if newends:
            old_ends_from_new_set = EndList(end for end in tset.ends
                                            if end['name'] not in newends)
            tileset.ends.merge(old_ends_from_new_set)

        # Ensure system consistency
        tset.check_consistent()
        return tset


    def create_strand_sequences(tileset,
                                basename,
                                includes=None,
                                spurious_pars="verboten_weak=1.5",
                                *options):
        """Given a tileset dictionary with sticky ends sequences, create core sequences
    for tiles.
        """

        newtileset = copy.deepcopy(tileset)

        newtileset.create_pepper_input_files(basename)

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

        tileset_with_strands = newtileset.load_pepper_output_files(basename)

        # Ensure:
        tileset.ends.merge(
            tileset_with_strands.tiles.endlist())  # Ends still fit
        for tile in tileset_with_strands.tiles:
            oldtile = tileset.tiles[tile.name]
            if 'fullseqs' in oldtile.keys():
                for old, new in zip(oldtile['fullseqs'], tile['fullseqs']):
                    seq.merge(old, new)  # old tile sequences remain
            assert oldtile.ends == tile.ends

        # Check that old end sequences remain
        tileset['ends'].merge(tileset_with_strands['ends'])

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

        for tile in tileset.tiles:
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
            tiletype = tile['structure']
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

        for tile in tset.tiles:
            pepperstrands = re.compile('strand ' + tile['name'] +
                                       '-([^ ]+) = ([^\n]+)').findall(seqsstring)
            tile['fullseqs'] = tilestructures.order_pepper_strands(pepperstrands)

        for adapter in tset['seed']['adapters']:
            pepperstrands = re.compile('strand ' + adapter['name'] +
                                       '-([^ ]+) = ([^\n]+)').findall(seqsstring)
            adapter['fullseqs'] = tilestructures.order_pepper_strands(pepperstrands)

        return tset

    def create_guard_strand_sequences(tileset):
        tset = copy.deepcopy(tileset)

        for guard in tset['guards']:
            tile = tilestructures.gettile(tset, guard[0])
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

    def generate_xgrow_dict(ts,
                            perfect=False,
                            rotate=False,
                            energetics=None):

        # Combine ends and tile-specified adjacents
        newtiles = []
        newends = []
        doubleends = []
        doubles = []
        vdoubleends = []
        vdoubles = []
        ts = copy.deepcopy(ts)
        newtiles.append({'name': 'origami',
                         'edges': ['origami', 'origami', 'origami', 'origami'],
                         'stoic': 0,
                         'color': 'white'})

        atiles = [None]*16
        to_use = []
        # If we have use_adapters, use that, otherwise use every adapter:
        if 'use_adapters' in ts['seed']:
            for tilename in ts['seed']['use_adapters']:
                try:
                    tile = [x for x in ts['seed']['adapters']
                            if x.get('name') == tilename][0]
                    to_use.append(tile)
                except IndexError as e:
                    raise Exception("Can't find {}".format(tilename)) from e
        else:
            to_use = ts['seed']['adapters']

        for tile in to_use:
            newtile = {}
            newtile['edges'] = [ 'origami' ] +  [ re.sub('/','_c',x) for x in tile['ends'] ] + [ 'origami' ]
            newtile['name'] = tile.get('name','')
            newtile['stoic'] = 0
            newtile['color'] = 'white'
            atiles[tile['loc']-1] = newtile
        for tile in atiles:
            if tile:
                newtiles.append(tile)
            else:
                newtiles.append( { 'name': 'emptyadapt', 'edges': ['origami',0,0,'origami'], 'stoic': 0, 'color': 'white'} )


        if rotate:
            rotatedtiles = []
            for tile in ts.tiles:
                if tile['structure'] == 'tile_daoe_3up' or tile['structure'] == 'tile_daoe_5up':
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_lrf'
                    newtile['ends']=[tile['ends'][x] for x in (1,0,3,2)]
                    rotatedtiles.append(newtile)
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_udf'
                    newtile['structure']='tile_daoe_'+{'5up':'3up','3up':'5up'}[tile['structure'][-3:]]
                    newtile['ends']=[tile['ends'][x] for x in (3,2,1,0)]
                    rotatedtiles.append(newtile)
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_bf'
                    newtile['structure']='tile_daoe_'+{'5up':'3up','3up':'5up'}[tile['structure'][-3:]]
                    newtile['ends']=[tile['ends'][x] for x in (2,3,0,1)]
                    rotatedtiles.append(newtile)
                elif tile['structure'] == 'tile_daoe_doublehoriz_35up':
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_lrf'
                    newtile['structure']='tile_daoe_doublevert_53up'
                    newtile['ends']=[tile['ends'][x] for x in (2,1,0,5,4,3)]
                    rotatedtiles.append(newtile)
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_udf'
                    newtile['structure']='tile_daoe_doublevert_53up'
                    newtile['ends']=[tile['ends'][x] for x in (5,4,3,2,1,0)]
                    rotatedtiles.append(newtile)
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_bf'
                    newtile['ends']=[tile['ends'][x] for x in (3,4,5,0,1,2)]
                    rotatedtiles.append(newtile)
                elif tile['structure'] == 'tile_daoe_doublevert_35up':
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_lrf'
                    newtile['structure']='tile_daoe_doublehoriz_53up'
                    newtile['ends']=[tile['ends'][x] for x in (2,1,0,5,4,3)]
                    rotatedtiles.append(newtile)
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_udf'
                    newtile['structure']='tile_daoe_doublehoriz_53up'
                    newtile['ends']=[tile['ends'][x] for x in (5,4,3,2,1,0)]
                    rotatedtiles.append(newtile)
                    newtile = copy.deepcopy(tile)
                    newtile['name']+='_bf'
                    newtile['ends']=[tile['ends'][x] for x in (3,4,5,0,1,2)]
                    rotatedtiles.append(newtile)

            ts.tiles += rotatedtiles

        for tile in ts.tiles:
            if tile['structure'] == 'tile_daoe_3up' or tile['structure'] == 'tile_daoe_5up':
                newtile = {}
                newtile['edges'] = [ re.sub('/','_c',x) for x in tile['ends'] ]
                if 'name' in tile: newtile['name'] = tile['name']
                if 'conc' in tile: newtile['stoic'] = tile['conc']
                if 'color' in tile: newtile['color'] = tile['color']
                newtiles.append(newtile)

            if tile['structure'] == 'tile_daoe_doublehoriz_35up' or tile['structure'] == 'tile_daoe_doublehoriz_53up':
                newtile1 = {}
                newtile2 = {}
                newtile1['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][0:1] ] \
                    + [ tile['name']+'_db' ] \
                    + [ re.sub('/','_c',x) for x in tile['ends'][4:] ]
                newtile2['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][1:4] ] \
                    + [ tile['name']+'_db' ]            
                newtile1['name'] = tile['name']+'_left'
                newtile2['name'] = tile['name']+'_right'

                doubleends.append( tile['name']+'_db' )
                doubles.append( (newtile1['name'], newtile2['name']) )

                if 'conc' in tile: 
                    newtile1['stoic'] = tile['conc']
                    newtile2['stoic'] = tile['conc']

                if 'color' in tile: 
                    newtile1['color'] = tile['color']
                    newtile2['color'] = tile['color']

                newtiles.append(newtile1)
                newtiles.append(newtile2)
            if tile['structure'] == 'tile_daoe_doublevert_35up' or tile['structure'] == 'tile_daoe_doublevert_53up':
                newtile1 = {}
                newtile2 = {}
                newtile1['edges'] = [ re.sub('/','_c',x) for x in tile['ends'][0:2] ] \
                    + [ tile['name']+'_db' ] \
                    + [ re.sub('/','_c',x) for x in tile['ends'][5:] ]
                newtile2['edges'] = [ tile['name']+'_db' ] + [ re.sub('/','_c',x) for x in tile['ends'][2:5] ] 
                newtile1['name'] = tile['name']+'_top'
                newtile2['name'] = tile['name']+'_bottom'

                vdoubleends.append( tile['name']+'_db' )
                vdoubles.append( (newtile1['name'], newtile2['name']) )

                if 'conc' in tile: 
                    newtile1['stoic'] = tile['conc']
                    newtile2['stoic'] = tile['conc']

                if 'color' in tile: 
                    newtile1['color'] = tile['color']
                    newtile2['color'] = tile['color']

                newtiles.append(newtile1)
                newtiles.append(newtile2)

        newends.append( { 'name': 'origami', 'strength': 100 } )

        for end in doubleends:
            newends.append( { 'name': end, 'strength': 10 } )
        for end in vdoubleends:
            newends.append( { 'name': end, 'strength': 10 } )

        gluelist = []
        if not perfect: 
            glueends = {'DT': [], 'TD': []}
            for end in ts['ends']:
                newends.append( { 'name': end['name'], 'strength': 0 } )
                newends.append( { 'name': end['name']+'_c', 'strength': 0 } )
                if (end['type'] == 'TD') or (end['type'] == 'DT'):
                    glueends[end['type']].append((end['name'],end['fseq']))

            if energetics:
                ef = energetics
            else:
                ef = DEFAULT_ENERGETICS

            eavg = {}
            for t in ['DT', 'TD']:
                names, fseqs = zip(*glueends[t])
                ea = sd.endarray(fseqs, t)
                eavg[t] = np.average( ef.matching_uniform( ea ) )
            eavg_combined = ( eavg['DT'] + eavg['TD'] ) / 2.0

            for t in ['DT','TD']:
                names, fseqs = zip(*glueends[t])
                allnames = names + tuple( x+'_c' for x in names )
                ea = sd.endarray(fseqs, t)
                ar = sd.energy_array_uniform(ea,ef) / eavg_combined
                for i1,n1 in enumerate(names):
                    for i2,n2 in enumerate(allnames):
                        gluelist.append([n1,n2,max(float(ar[i1,i2]),0.0)])


        else:
            if 'ends' not in ts.keys():
                ts['ends']=[]
            endsinlist = set( e['name'] for e in ts['ends'] )
            endsintiles = set()
            for tile in ts.tiles:
                endsintiles.update( re.sub('/','',e) for e in tile['ends'] if e != 'hp')
            for end in ts['ends'] + list({'name': e} for e in endsintiles):
                newends.append( { 'name': end['name'], 'strength': 0 } )
                newends.append( { 'name': end['name']+'_c', 'strength': 0 } )
                gluelist.append([end['name'],end['name']+'_c',1.0]) 



        newends.append( {'name': 'hp', 'strength': 0} )

        xga = {}
        xga['doubletiles'] = [ list(x) for x in doubles ]
        xga['vdoubletiles'] = [ list(x) for x in vdoubles ]
        xga.update( ts['xgrow_options'] )
        xga.update( ts['xgrow_options'] )
        #if not perfect:
        #    xga['gse_calc_avg'] = eavg_combined


        sts = { 'tiles': newtiles, 'bonds': newends, 'xgrowargs': xga, 'glues': gluelist }

        return sts

    
    def plot_se_hists(tileset,
                      all_energetics=None,
                      energetics_names=None,
                      title=None, **kwargs):

        if all_energetics is None:
            all_energetics = DEFAULT_MULTIMODEL_ENERGETICS

        if energetics_names is None:
            energetics_names = DEFAULT_MM_ENERGETICS_NAMES

        if 'ends' in tileset.keys():
            ends = tileset['ends']
        else:
            ends = tileset

        if title is None:
            # FIXME
            title = 'Title'

        td = sd.endarray([x['fseq'] for x in ends
                          if x['type'] == 'TD'], 'TD')

        dt = sd.endarray([x['fseq'] for x in ends
                          if x['type'] == 'DT'], 'DT')

        sdplots.hist_multi([td, dt],
                           all_energetics,
                           energetics_names,
                           title, **kwargs)


    def plot_adjacent_regions(tileset,
                              energetics=None):

        if energetics is None:
            energetics = DEFAULT_REGION_ENERGETICS

        regions = [t.structure._side_bound_regions
                   for t in tileset.tiles]
        regions = [[x.lower() for x in y] for y in regions]
        allregions = sum(regions, [])
        count = [[Counter(x) for x in y] for y in regions]
        gc_count = [[x['g']+x['c'] for x in c] for c in count]
        gc_counts = sum(gc_count, [])

        ens = energetics.matching_uniform(sd.endarray(allregions, 'DT'))

        pylab.figure(figsize=(10, 4))
        pylab.subplot(121)
        pylab.hist(gc_counts,
                   bins=np.arange(min(gc_counts)-0.5, max(gc_counts)+0.5))
        pylab.title('G/C pairs in arms')
        pylab.ylabel('# of 15 nt arms')
        pylab.xlabel('# of G/C pairs')
        pylab.subplot(122)
        pylab.hist(ens)
        pylab.title('ΔG, T=33, no coaxparams/danglecorr')
        pylab.ylabel('# of 8 nt regions')
        pylab.xlabel('stickydesign ΔG')
        pylab.suptitle('8 nt end-adjacent region strengths')

    def plot_side_strands(tileset,
                          energetics=None):

        if energetics is None:
            energetics = DEFAULT_REGION_ENERGETICS

        regions = [t.structure._short_bound_full
                   for t in tileset.tiles]
        regions = [[x.lower() for x in y] for y in regions]
        allregions = sum(regions, [])
        count = [[Counter(x) for x in y] for y in regions]
        gc_count = [[x['g']+x['c'] for x in c] for c in count]
        gc_counts = sum(gc_count, [])

        ens = energetics.matching_uniform(sd.endarray(allregions, 'DT'))

        pylab.figure(figsize=(10, 4))
        pylab.subplot(121)
        pylab.hist(gc_counts,
                   bins=np.arange(min(gc_counts)-0.5, max(gc_counts)+0.5))
        pylab.title('G/C pairs in arms')
        pylab.ylabel('# of 15 nt arms')
        pylab.xlabel('# of G/C pairs')
        pylab.subplot(122)
        pylab.hist(ens)
        pylab.title('ΔG, T=33, no coaxparams/danglecorr')
        pylab.ylabel('# of 16 nt regions')
        pylab.xlabel('stickydesign ΔG')
        pylab.suptitle('16 nt arm region strengths')


RoundTripRepresenter.add_representer(TileSet,
                                     RoundTripRepresenter.represent_dict)

