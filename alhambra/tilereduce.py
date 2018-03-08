from .util import base, comp
from itertools import product, combinations
from .endreduce import equate_pair, latticedefects
import re
import logging
from .tilesets import TileSet
from .tiles import TileList
from random import shuffle

log = logging.getLogger(__name__)

RS = [2, 3, 0, 1]


def inputtilepairset(tileset, tile, rti=None):
    tpg = []
    for s, e, ip in zip(tile.structure._dirs, tile['ends'], tile['input']):
        if not ip:
            continue
        tt = []
        for tile2 in tileset.tiles:
            if comp(e) not in tile2['ends']:
                continue
            for e2, s2, ip2 in zip(tile2['ends'], tile2.structure._dirs,
                                   tile2['input']):
                if (not ip2) and (RS[s] == s2) and (comp(e) == e2):
                    if (rti is None) or (tile2.name not in rti):
                        tt.append(tile2.name)
                    else:
                        tt.append(rti[0])
                        tt.append(rti[1])
        tpg.append(tt)
    return set(product(*tpg))


def equate_multiends(ts, pairs, eqret=False):
    rpairs = []
    pairs = list(pairs)
    while len(pairs) > 0:
        e1, e2 = pairs.pop()
        swap = not ((e1[-1] == '/' and e2[-1] == '/') or
                    (e1[-1] != '/' and e2[-1] != '/'))

        def rfunc(match):
            if swap:
                return comp(base(e1) + match.group(1))
            else:
                return base(e1) + match.group(1)

        r = re.compile(base(e2) + r'(/?)$')
        rpairs.append((e1, e2))
        for i, pair in enumerate(pairs):
            pair = [r.sub(rfunc, x) for x in pair]
            # we would set a pair = its complement
            if pair[0] == comp(pair[1]):
                if not eqret:
                    return False
                else:
                    return (False, False)
            # we need to update the pair in the list: FIXME is this true?
            pairs[i] = pair
            
        pairs = [p for p in pairs if p[0] != p[1]]
    for p in rpairs:
        ts = equate_pair(ts, p, unsafe=True)
    if not eqret:
        return ts
    else:
        return ts, rpairs


def equate_tiles(tileset, tilepair, eqret=False):
    if eqret:
        fret = (False, False)
    else:
        fret = False
    T1, T2 = tilepair
    if (len(T1.ends) != 4) or (len(T2.ends) != 4):
        return fret
    for tp1, tp2 in zip(T1.structure._endtypes, T2.structure._endtypes):
        if tp1 != tp2:
            return fret
    endcombs, endreps = equate_multiends(
        tileset, zip(T1.ends, T2.ends), eqret=True)
    if endcombs is False:
        return fret
    del (endcombs.tiles[T2.name])
    if not eqret:
        return endcombs
    else:
        return endcombs, endreps


def tryreducerot(ts, rti, rot=0, checkld=False):
    if rot == 0:
        rts, pr = equate_tiles(
            ts, (ts.tiles[rti[0]], ts.tiles[rti[1]]), eqret=True)
    else:
        rts, pr = equate_tiles(
            ts, (ts.tiles[rti[0]], ts.tiles[rti[1]].rotations[rot - 1]),
            eqret=True)
        if rts is False:
            return False
        # add fake tile: this is a really stupid method.
        fakets = TileSet({'tiles': TileList([ts.tiles[rti[1]]].copy())})
        for p in pr:
            fakets = equate_pair(fakets, p, unsafe=True)

        rts.tiles.append(fakets.tiles[0])
    if rts is False:
        return False
    allpairset = set.union(*(inputtilepairset(ts, x) for x in ts.tiles))
    for x in ts.tiles:
        if (rot != 0):
            oldps = inputtilepairset(ts, x)
            newps = set.intersection(
                inputtilepairset(rts, rts.tiles[x.name]), allpairset)
        elif (x.name not in rti):
            oldps = inputtilepairset(ts, x)
            newps = set.intersection(
                inputtilepairset(rts, rts.tiles[x.name], rti=rti), allpairset)
        else:
            oldps = set.union(
                inputtilepairset(ts, ts.tiles[rti[0]]),
                inputtilepairset(ts, ts.tiles[rti[1]]))
            newps = set.intersection(
                inputtilepairset(rts, rts.tiles[rti[0]], rti=rti), allpairset)
        if oldps != newps:
            return False
    # Now note the fake tile
    if rot > 0:
        rts.tiles[rti[1]]['fake'] = 1
    if not check_changes_multi(ts, rts, pr, checkld=checkld):
        return False
    log.debug("Changes: {}".format(pr))
    if rts and rts.seed and (rot == 0):
        for t in rts.seed['adapters']:
            if rti[1] == t['tilebase']:
                t['tilebase'] = rti[0]
    return rts


def check_changes_multi(oldts,
                        newts,
                        equatedpairs,
                        oldsc=None,
                        newsc=None,
                        checkld=False):
    if oldsc is None:
        oldsc = oldts.sensitivity_classes()
    if newsc is None:
        newsc = newts.sensitivity_classes()

    reps = equatedpairs

    for pair in newsc['2GO']:
        if pair not in oldsc['2GO']:
            for eremain, ereplace in reps:
                if eremain in pair:
                    pair = frozenset.union(pair - {eremain}, {ereplace})
                elif comp(eremain) in pair:
                    pair = frozenset.union(pair - {comp(eremain)},
                                           {comp(ereplace)})
            if pair not in oldsc['2GO']:
                return False

    # FIXME: add lattice defect check
    if checkld:
        oldld = latticedefects(oldts, depth=checkld)
        newld = latticedefects(newts, depth=checkld)
        if (len(newld) > len(oldld)) or (len(
            [x for x in newld if x not in oldld]) > 0):
            log.debug("lattice defect: {}".format(equatedpairs))
            return False
    return True


def reduce_tiles(tileset,
                 colors=None,
                 rotation=True,
                 checkld=True,
                 _unsafe=True,
                 update=1000):
    """Attempt maximal end reduction on a tileset.
    
    Parameters
    ----------

    tileset : alhambra.TileSet
        The TileSet to reduce.  Note that reduce_tiles does not currently preprocess
        to avoid clashes with manual cleverness, so you should avoid using tilesets
        that have had their ends reduced, or that reuse the same ends for multiple
        purposes.

    colors : None, dict, or other mapping
    """

    # Start by making a copy of the set, to make absolutely sure we don't mangle it:
    ts = tileset.copy()
    oldts = ts  # FIXME remove when confident
    pairstodo = list(combinations(ts.tiles, 2))
    shuffle(pairstodo)

    pairspassed = []
    removedtiles = []
    reminfo = []
    idone = 0
    if rotation:
        dirs = [0, 1, 2, 3]
    else:
        dirs = [0]

    log.info(
        "Started tile reduction: {} tile combinations".format(len(pairstodo)))

    while len(pairstodo) > 0:
        t1, t2 = pairstodo.pop()

        if ('fake' in t1.keys()) or ('fake' in t2.keys()):
            continue
        if (t1.name in removedtiles) or (t2.name in removedtiles):
            log.warning(
                "{} or {} in removed tiles, but should have been removed from todo".
                format(t1.name, t2.name))

        if (idone > 0) and (idone % update == 0):
            log.info(
                "Finished {} pairs, {} left.".format(idone, len(pairstodo)))

        shuffle(dirs)
        for r in dirs:
            trialts = tryreducerot(
                ts, (t1.name, t2.name), rot=r, checkld=checkld)
            if trialts is not False:
                ts = trialts
                removedtiles.append(t2.name)
                reminfo.append((t1.name, t2.name, r))

                ts.tiles[t1.name]['functionsas'] = t1.get(
                    'functionsas', []) + t2.get('functionsas', []) + [t2.name]

                # Remove pairs that have the removed tile.
                pairstodo = [
                    pair for pair in pairstodo
                    if (t2.name != pair[0].name) and (t2.name != pair[1].name)
                ]
                pairspassed = [
                    pair for pair in pairspassed
                    if (t2.name != pair[0].name) and (t2.name != pair[1].name)
                ]

                # We've removed a tile, so now we'll need to add the passed pairs back into the todo list.  We'll shuffle, too:
                shuffle(pairspassed)
                pairstodo += pairspassed
                pairspassed = []

                log.info("{} -> {} (rot {}), {} removed, {} pairs left".format(
                    t2.name, t1.name, r, len(removedtiles), len(pairstodo)))
                break
        idone += 1
    log.info(
        "Finished tile reduction: {} -> {} ({} reduction), checked {} pairs.".
        format(
            len(tileset.tiles),
            len(tileset.tiles) - len(removedtiles), len(removedtiles), idone))

    # End by ensuring that unsafe didn't do anything bad:
    if oldts != tileset:
        log.warning("ts != tileset in reduce_tiles")
        return ts, oldts

    ts.add_info('tilereduce',
               {'rotation': rotation,
                'checkld': checkld,
                'removed': reminfo})

    return ts
