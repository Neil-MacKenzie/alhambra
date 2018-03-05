from .util import base, comp
from itertools import product, combinations
from .endreduce import equate_pair, latticedefects
import re
import logging
from .tilesets import TileSet
from .tiles import TileList

log = logging.getLogger(__name__)

RS = [2, 3, 0, 1]


def inputtilepairset(tileset, tile, rti=None):
    if len(tile.ends) != 4:
        return set()
    tpg = []
    for s, (e, ip) in enumerate(zip(tile.ends, tile['input'])):
        if not ip:
            continue
        tt = []
        for tile2 in tileset.tiles:
            if len(tile.ends) != 4:
                continue
            if (tile2.ends[RS[s]] == comp(e)) and (not tile2['input'][RS[s]]):
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
        for pair in pairs:
            pair = [r.sub(rfunc, x) for x in pair]
            if pair[0] == comp(pair[1]):
                if not eqret:
                    return False
                else:
                    return (False, False)
        pairs = [p for p in pairs if p[0] != p[1]]
    for p in rpairs:
        ts = equate_pair(ts, p)
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


def tryreducenpr(ts, rti):
    rts, pr = equate_tiles(
        ts, (ts.tiles[rti[0]], ts.tiles[rti[1]]), eqret=True)
    if rts is False:
        return False
    allpairset = set.union(*(inputtilepairset(ts, x) for x in ts.tiles))
    for x in ts.tiles:
        if x.name not in rti:
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
    if not check_changes_multi(ts, rts, pr):
        return False
    return rts


def tryreducerot(ts, rti, rot=0):
    if rot == 0:
        rts, pr = equate_tiles(
            ts, (ts.tiles[rti[0]], ts.tiles[rti[1]]), eqret=True)
        if rts and rts.seed:
            for t in rts.seed['adapters']:
                if rti[1] == t['tilebase']:
                    t['tilebase'] = rti[0]
    else:
        rts, pr = equate_tiles(
            ts, (ts.tiles[rti[0]], ts.tiles[rti[1]].rotations[rot - 1]),
            eqret=True)
        if rts is False:
            return False
        # add fake tile: this is a really stupid method.
        fakets = TileSet({
            'tiles': TileList([ts.tiles[rti[1]]].copy())
        })
        for p in pr:
            fakets = equate_pair(fakets, p)

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
    if not check_changes_multi(ts, rts, pr):
        return False
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
