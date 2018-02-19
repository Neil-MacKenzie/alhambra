# New sensitivity code.
from .tilestructures import tile_daoe_single
from .tiles import TileList, Tile
from collections import Counter
from .util import comp


class CounterSet(Counter):
    def add(self, val):
        self[val] += 1


def _fakesingle(tile):
    if not tile.structure.double:
        return TileList([tile])
    ret = TileList()
    for es, c in zip(tile.structure.singleends, ['', '/']):
        ft = Tile()
        ft.structure = tile_daoe_single()
        ft.structure._endlocs = [None, None, None, None]
        ft['ends'] = []
        if 'input' in tile.keys():
            ft['input'] = []
        ft.structure._endtypes = []
        for i in es:
            if i is not None:
                ft['ends'].append(tile.ends[i])
                ft.structure._endtypes.append(tile.structure._endtypes[i])
                if 'input' in tile.keys():
                    ft['input'].append(tile['input'][i])
            else:
                ft['ends'].append(tile.name + '_db' + c)
                ft.structure._endtypes.append('fakedouble')
                if 'input' in tile.keys():
                    ft['input'].append(0)
        ret.append(ft)
    return ret


def _fakesingles(tiles):
    return sum((_fakesingle(x) for x in tiles), TileList())



_rev = [2, 3, 0, 1]


def sensitivity_classes(tileset, count=False):

    singles = _fakesingles(tileset.tiles)
    rotatedsingles = singles + _fakesingles(
        sum([x.rotations for x in tileset.tiles], TileList()))

    if count:
        sclass = CounterSet
    else:
        sclass = set
    
    spairs = {'1GO': sclass(), '1NGO': sclass(), '2NGO': sclass(), '2GO': sclass()}

    for t1 in singles:
        for i in range(0, 4):
            if t1.structure._endtypes[i] in {'fakedouble', 'hairpin'}:
                continue
            for t2 in rotatedsingles:
                if t2.ends[i] != t1.ends[i]:
                    continue
                for j in range(0, 4):
                    if i == j:
                        continue
                    if t1.ends[j] == t2.ends[j]:
                        continue
                    if t1.structure._endtypes[j] != t2.structure._endtypes[j]:
                        continue
                    if t1.structure._endtypes[j] in {'fakedouble', 'hairpin'}:
                        continue
                    spairs['1NGO'].add(frozenset((t2.ends[j], t1.ends[j])))
                    spairs['1NGO'].add(frozenset((comp(t2.ends[j]), comp(t1.ends[j]))))
                    go1 = False
                    if t1['input'][i] and t1['input'][j]:
                        go1 = True
                        spairs['1GO'].add(frozenset((t2.ends[j], t1.ends[j])))
                        spairs['1GO'].add(frozenset((comp(t2.ends[j]), comp(t1.ends[j]))))
                    for k in (set(range(0, 4)) - {i, j}):
                        if (t1.structure._endtypes[k] == 'hairpin') or (
                                t2.structure._endtypes[k] == 'hairpin'):
                            continue
                        ec1 = comp(t1.ends[k])
                        ec2 = comp(t2.ends[k])
                        kc = _rev[k]
                        for t12 in singles:
                            if t12.ends[kc] != ec1:
                                continue
                            for t22 in rotatedsingles:
                                if t22.ends[kc] != ec2:
                                    continue
                                for m in range(0, 4):
                                    if m == kc:
                                        continue
                                    if t22.ends[m] != t12.ends[m]:
                                        continue
                                    if t12.structure._endtypes[m] in {
                                            'fakedouble', 'hairpin'
                                    }:
                                        continue
                                    spairs['2NGO'].add(
                                        frozenset((t2.ends[j], t1.ends[j])))
                                    spairs['2NGO'].add(
                                        frozenset((comp(t2.ends[j]), comp(t1.ends[j]))))
                                    if t12['input'][m] and t12['input'][kc] and go1:
                                        spairs['2GO'].add(
                                            frozenset((t2.ends[j], t1.ends[j]
                                                       )))
                                        spairs['2GO'].add(
                                            frozenset((comp(t2.ends[j]), comp(t1.ends[j])
                                                       )))

    return spairs
