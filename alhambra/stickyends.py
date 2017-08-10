from .tilesets import TileSet
from . import tiletypes
from . import util

import stickydesign as sd
import numpy as np
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


def create_sequences(tileset, method='default', energetics=None,
                     trials=100, sdopts={}, ecpars={}):
    """Create sticky end sequences for a tileset, using stickydesign.  This new
version should be more flexible, and should be able to keep old sticky ends,
accounting for them in creating new ones.

Parameters
----------

tileset: the tileset to create sticky ends sequences for.  This will be copied
     and returned, not modified.

method: [default 'default'] if 'default', use the default, single-model
sequence design.  If 'multimodel', use multimodel end choice.

energetics: the energetics instance to use for the design.  If None (default),
will use alhambra.designer.DEFAULT_ENERGETICS.  Note that if
method='multimodel', energetics *must* be specified, and should be a list of
energy models.  The first model in the list will be used as primary.

Outputs (tileset, new_ends) where new_ends is a list of new end names that
were designed.
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
    newtileset = TileSet(tileset).copy()

    # Build a list of ends from the endlist in the tileset.  Do this
    # by creating a named_list, then merging them into it.
    ends = util.named_list()

    if 'ends' in newtileset.keys():
        ends = util.merge_endlists(
            ends, newtileset['ends'], fail_immediate=False, in_place=True)

    # This is the endlist from the tiles themselves.
    if 'tiles' in newtileset.keys():  # maybe you just want ends?
        # this checks for end/complement usage, and whether any
        # previously-describedends are unused
        # FIXME: implement
        # tiletypes.check_end_usage(newtileset['tiles'], ends)

        endlist_from_tiles = tiletypes.endlist_from_tilelist(
            newtileset['tiles'])

    ends = util.merge_endlists(ends, endlist_from_tiles, in_place=True)

    # Ensure that if there are any resulting completely-undefined ends, they
    # have their sequences removed.
    for end in ends:
        if 'fseq' in end.keys() and end['fseq'] == 'nnnnnnn':
            del (end['fseq'])

    # Build inputs suitable for stickydesign: lists of old sequences for TD/DT,
    # and numbers of new sequences needed.
    oldDTseqs = [
        end['fseq'] for end in ends
        if 'fseq' in end.keys() and end['type'] == 'DT'
    ]
    oldTDseqs = [
        end['fseq'] for end in ends
        if 'fseq' in end.keys() and end['type'] == 'TD'
    ]

    newTDnames = [
        end['name'] for end in ends
        if 'fseq' not in end.keys() and end['type'] == 'TD'
    ]
    newDTnames = [
        end['name'] for end in ends
        if 'fseq' not in end.keys() and end['type'] == 'DT'
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
            energetics.matching_uniform(sd.endarray(oldDTseqs, 'DT')))
    if len(oldTDseqs) > 0:
        targets.append(
            energetics.matching_uniform(sd.endarray(oldTDseqs, 'TD')))
    targetint = np.average(targets)

    if method == 'multimodel' and (len(oldDTseqs) > 0 or len(oldTDseqs) > 0):
        raise NotImplementedError()

    if method == 'default':
        # Create new sequences.
        newTDseqs = sd.easyends(
            'TD',
            5,
            number=len(newTDnames),
            energetics=energetics,
            interaction=targetint,
            **sdopts).tolist()

        newDTseqs = sd.easyends(
            'DT',
            5,
            number=len(newDTnames),
            energetics=energetics,
            interaction=targetint,
            **sdopts).tolist()

    elif method == 'multimodel':
        SELOGGER.info(
            "starting multimodel sticky end generation " +
            "of TD ends for {} DT and {} TD ends, {} trials.".format(
                len(newDTnames), len(newTDnames), trials))

        endchooser = sd.multimodel.endchooser(all_energetics, **ecpars)

        newTDseqs = []
        pl = util.ProgressLogger(SELOGGER, trials*2)
        for i in range(0, trials):
            newTDseqs.append(
                sd.easyends(
                    'TD',
                    5,
                    number=len(newTDnames),
                    energetics=energetics,
                    interaction=targetint,
                    echoose=endchooser,
                    **sdopts))
            pl.update(i)

        tvals = [[e.matching_uniform(x[0:1])
                 for e in all_energetics] for x in newTDseqs]
        endchoosers = [sd.multimodel.endchooser(all_energetics,
                                                target_vals=tval,
                                                **ecpars)
                       for tval in tvals]

        SELOGGER.info("generating corresponding DT ends")
        newDTseqs = []
        for i, echoose in enumerate(endchoosers):
            newDTseqs.append(
                sd.easyends(
                    'DT',
                    5,
                    number=len(newDTnames),
                    energetics=energetics,
                    interaction=targetint,
                    echoose=echoose,
                    **sdopts))
            pl.update(i+trials)

        scores = [sd.multimodel.deviation_score(list(e), all_energetics)
                  for e in zip(newTDseqs, newDTseqs)]

        sort = np.argsort(scores)

        newTDseqs = newTDseqs[sort[0]].tolist()
        newDTseqs = newDTseqs[sort[0]].tolist()
        info['score'] = scores[sort[0]]
        info['maxscore'] = scores[sort[-1]]
        info['meanscore'] = np.mean(scores)

    # FIXME: move to stickydesign
    assert len(newTDseqs) == len(newTDnames)
    assert len(newDTseqs) == len(newDTnames)

    # Shuffle the lists of end sequences, to ensure that they're random order,
    # and that ends used earlier in the set are not always better than those
    # used later.
    shuffle(newTDseqs)
    shuffle(newDTseqs)

    for name, s in zip(newDTnames, newDTseqs):
        ends[name]['fseq'] = s
    for name, s in zip(newTDnames, newTDseqs):
        ends[name]['fseq'] = s

    ends.check_consistent()

    # Ensure that the old and new sets have consistent end definitions,
    # and that the tile definitions still fit.
    tiletypes.merge_endlists(tileset['ends'], ends)
    tiletypes.merge_endlists(
        tiletypes.endlist_from_tilelist(newtileset['tiles']), ends)

    # Apply new sequences to tile system.
    newtileset['ends'] = ends
    if 'info' not in newtileset.keys():
        newtileset['info'] = {}
    newtileset['info']['end_design'] = info

    return (newtileset, newTDnames + newDTnames)


def reorder(tileset,
            newends=[],
            hightemp=0.1,
            lowtemp=1e-7,
            steps=45000,
            update=1000,
            energetics=None):
    """Given a tileset dictionary that includes sticky end sequences, reorder these
    to try to optimize error rates.
    """
    # from . import endreorder
    from . import endreorder_fast
    from . import anneal

    if energetics is None:
        energetics = DEFAULT_ENERGETICS

    tset = TileSet(tileset).copy()

    if 'info' not in tset.keys():
        tset['info'] = {}

    reordersys = endreorder_fast.EndSystemFseq(
        tset, newends, energetics=energetics)

    # reordersys_old = endreorder.EndSystemFseq( tset, energetics=energetics )

    # FIXME: better parameter control here.
    annealer = anneal.Annealer(reordersys.score, reordersys.mutate)

    newstate = annealer.anneal(reordersys.initstate, hightemp, lowtemp, steps,
                               update)

    # Now take that new state, and apply it to the new tileset.
    seqs = reordersys.slowseqs(newstate[0])
    for end in tset['ends']:
        if end['type'] in ['DT', 'TD']:
            eloc = reordersys.enlocs[end['name']]
            end['fseq'] = seqs[eloc[1]].tolist()[eloc[0]]

    ri = {}

    ri['score'] = reordersys.score(newstate[0])

    tset['info']['reorder'] = ri

    # Ensure that only ends in newends moved: that all others remain mergeable:
    if newends:
        old_ends_from_new_set = util.named_list(end for end in tset['ends']
                                                if end['name'] not in newends)
        util.merge_endlists(tileset['ends'], old_ends_from_new_set)

    # Ensure system consistency
    tset.check_consistent()
    return tset
