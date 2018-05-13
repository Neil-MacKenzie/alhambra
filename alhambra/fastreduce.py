from collections import namedtuple
from copy import copy
import numpy as np
from . import tilestructures as ts
from .tiles import TileList
from .ends import End
from random import shuffle

FTile = namedtuple("FTile", ("color", "use", "glues", "name", "used",
                             "structure", "dfake"))

FTilesArray = namedtuple("FTilesArray", ("color", "use", "glues", "name",
                                         "used", "structure", "dfake"))

TAU = 2


class FGlueList():
    def __init__(self, glues):
        self.name = []
        self.strength = []
        self.structure = []
        self.complement = []
        self.tonum = {}
        # self.use = []
        for i, g in enumerate(glues):
            self.name.append(g.name)
            self.name.append(g.name + '/')
            self.complement.append(2 * i + 1)
            self.complement.append(2 * i)
            self.structure.append(g.etype)
            self.structure.append(g.etype)
            self.tonum.update({g.name: 2 * i, g.name + '/': 2 * i + 1})
            self.strength.append(g.strength)
            self.strength.append(g.strength)
        self.name = np.array(self.name)
        self.strength = np.array(self.strength)
        self.structure = np.array(self.structure)
        self.complement = np.array(self.complement)

    def blankequiv(self):
        return np.arange(0, len(self.name))

    def iseq(self, equiv, a, b):
        return equiv[a] == equiv[b]

    def domerge(self, equiv, a, b):
        if self.structure[a] != self.structure[b]:
            raise ValueError("structure")
        elif self.strength[a] != self.strength[b]:
            raise ValueError("strength")
        elif equiv[a] == equiv[self.complement[b]]:
            raise ValueError("self-comp")
        else:
            equiv = copy(equiv)
            newg, oldg = sorted((equiv[a], equiv[b]))
            newc, oldc = sorted((equiv[self.complement[a]],
                                 equiv[self.complement[b]]))
            equiv[equiv == oldg] = newg
            equiv[equiv == oldc] = newc
        return equiv


U_UNUSED = -1
U_INPUT = 1
U_OUTPUT = 0


class FTileList():
    def __init__(self, tiles, gluelist):
        self.tiles = []
        self.totile = {}
        for t in tiles:
            glues = [gluelist.tonum[x] for x in t.ends]
            if 'input' not in t.keys():
                used = False
                use = [U_UNUSED for _ in t.ends]
            else:
                used = True
                use = t['input']
            color = 'label' in t.keys()
            self.tiles.append(
                FTile(
                    name=t.name,
                    color=color,
                    use=[int(x) for x in use],
                    glues=glues,
                    used=used,
                    structure=t.structure,
                    dfake=False))
            if used:
                self.totile[t.name] = self.tiles[-1]
        stiles = []
        htiles = []
        vtiles = []

        for t in self.tiles:
            if isinstance(t.structure, ts.tile_daoe_single):
                stiles.append(t)
            elif isinstance(t.structure, ts.tile_daoe_doublehoriz):
                stiles += _ffakesingle(t, gluelist)
                htiles.append(t)
            elif isinstance(t.structure, ts.tile_daoe_doublevert):
                stiles += _ffakesingle(t, gluelist)
                vtiles.append(t)
            else:
                raise NotImplementedError
        self.stiles = _ft_to_fta(stiles)

        for tile in stiles:
            x, y = _ffakedouble(tile, self.stiles, gluelist)
            htiles += x
            vtiles += y

        self.htiles = _ft_to_fta(htiles)
        self.vtiles = _ft_to_fta(vtiles)


RSEL = (2, 3, 0, 1)
FTI = (1, 0, 1, 0)
FTS = (ts.tile_daoe_doublevert(), ts.tile_daoe_doublehoriz(),
       ts.tile_daoe_doublevert(), ts.tile_daoe_doublehoriz())


def _fdg(dir, gs1, gs2):
    if dir == 0:
        return [gs2[0], gs2[1], gs1[1], gs1[2], gs1[3], gs2[3]]
    elif dir == 1:
        return [gs1[0], gs2[0], gs2[1], gs2[2], gs1[2], gs1[3]]
    elif dir == 2:
        return [gs1[0], gs1[1], gs2[1], gs2[2], gs2[3], gs1[3]]
    elif dir == 3:
        return [gs2[0], gs1[0], gs1[1], gs1[2], gs2[2], gs2[3]]
    else:
        raise ValueError(dir)


def _ffakedouble(tile, sta, gluelist):
    faketiles = ([], [])
    for dir in range(0, 4):
        if tile.use[dir] == 1:
            continue
        oti = np.nonzero(
            (gluelist.complement[tile.glues[dir]] == sta.glues[:, RSEL[dir]]) &
            (sta.use[:, RSEL[dir]] == 1) & (sta.used))
        for i in oti[0]:
            faketiles[FTI[dir]].append(
                FTile(
                    color=False,
                    used=True,
                    name=tile.name + '_{}_'.format(dir) + sta.name[i],
                    structure=FTS[dir],
                    glues=_fdg(dir, tile.glues, sta.glues[i]),
                    use=_fdg(dir, tile.use, sta.use[i]),
                    dfake=True))
    return faketiles


# THESE SELECTORS ARE 1-INDEXED!!! 0 corresponds to fake double tile bond.
HSEL = ((1, 0, 5, 6), (2, 3, 4, 0))
VSEL = ((1, 2, 0, 6), (0, 3, 4, 5))


def _ffakesingle(ftile, gluelist):
    # FIXME: should be more generalized.  Currently only tau=2
    if isinstance(ftile.structure, ts.tile_daoe_doublehoriz):
        sel = HSEL
    elif isinstance(ftile.structure, ts.tile_daoe_doublevert):
        sel = VSEL
    else:
        raise NotImplementedError

    # Start by making the tiles, then change around the inputs
    fdb = gluelist.tonum['fakedouble']
    glues = [[([fdb] + ftile.glues)[x] for x in y] for y in sel]
    fuse = [[([-1] + ftile.use)[x] for x in y] for y in sel]
    use = []
    used = []
    names = [ftile.name + '_fakedouble_a', ftile.name + 'fakedouble_b']
    for gu in zip(glues, fuse):
        if sum(gluelist.strength[g] for g, u in zip(*gu) if u == 1) >= TAU:
            use.append(gu[1])
            used.append(True)
        else:
            use.append([-1, -1, -1, -1])
            used.append(False)
    return [
        FTile(
            color=ftile.color,
            use=u,
            glues=g,
            name=n,
            used=ud,
            structure=ts.tile_daoe_single(),
            dfake=False) for u, g, n, ud in zip(use, glues, names, used)
    ]


def _ft_to_fta(ftiles):
    return FTilesArray(
        color=np.array([x.color for x in ftiles]),
        name=np.array([x.name for x in ftiles]),
        glues=np.array([x.glues for x in ftiles]),
        use=np.array([x.use for x in ftiles]),
        used=np.array([x.used for x in ftiles]),
        structure=np.array([x.structure.name for x in ftiles]),
        dfake=np.array([x.dfake for x in ftiles]))


class FTileSystem():
    def __init__(self, tilesystem):
        self.gluelist = FGlueList(tilesystem.allends + [
            End({
                'name': 'hp',
                'type': 'hairpin',
                'strength': 0
            }),
            End({
                'name': 'fakedouble',
                'type': 'fakedouble',
                'strength': 0
            })
        ])
        self.tilelist = FTileList(tilesystem.tiles + sum(
            [x.named_rotations()
             for x in tilesystem.tiles], TileList()), self.gluelist)

    def applyequiv(self, ts, equiv):
        ts = ts.copy()
        alreadythere = []
        for tile in ts.tiles:
            tile.ends = [
                self.gluelist.name[equiv[self.gluelist.tonum[e]]]
                for e in tile.ends
            ]
            if (tile.ends, 'label' in tile.keys()) in alreadythere:
                tile['fake'] = True
                continue
            rs = [tile] + tile.rotations
            alreadythere += [(t.ends, 'label' in t.keys()) for t in rs]
        if 'seed' in ts.keys():
            for t in ts.seed['adapters']:
                if 'ends' in t.keys():
                    t['ends'] = [
                        self.gluelist.name[equiv[self.gluelist.tonum[e]]]
                        for e in t['ends']
                    ]
        return ts


def ptins(fts, equiv, tau=2):
    """Calculate potential tile attachments to input neighborhoods"""
    ptins = []
    for ta in [fts.tilelist.stiles, fts.tilelist.htiles, fts.tilelist.vtiles]:
        ptin = []
        for ti in np.arange(0, len(ta.used))[ta.used]:
            # Only iterate through used tiles
            isel = ta.use[ti] == 1  # used edges
            gsel = equiv[ta.glues[ti, isel]] == equiv[ta.glues[:, isel]]
            # matching glues
            strs = np.sum(
                fts.gluelist.strength[ta.glues[ti, isel]] * gsel,
                axis=1)  # strengths of matching
            matches = np.nonzero((strs >= tau) & (~ta.dfake))
            # indices of matching
            # (excludes fake doubles, which can't actually attach
            # (because they are
            # actually two singles) to their own local neighborhoods)
            ptin.append(matches)
        ptins.append(ptin)
    return ptins


def isatamequiv(fts, equiv, initptins=None):
    if initptins is None:
        initptins = ptins(fts, fts.gluelist.blankequiv())
    npt = ptins(fts, equiv)
    for x, y, tl in zip(initptins, npt, [
            fts.tilelist.stiles, fts.tilelist.htiles, fts.tilelist.vtiles
    ]):
        for xx, yy in zip(x, y):
            if len(xx[0]) == len(yy[0]) and (np.all(
                    xx[0] == yy[0])):  # No change
                continue
            elif len(xx[0]) == 1:  # Deterministic start
                mm = ~np.all(
                    equiv[tl.glues[xx]] == equiv[tl.glues[yy]], axis=1)
                if np.any(mm):
                    return False, (tl.name[xx[0][0]], tl.name[yy[0][mm][0]])
            elif len(xx[0]) == 0:
                return False, None
            else:
                raise NotImplementedError
    return True, None


def tilemerge(fts, equiv, t1, t2):
    if t1.structure.name != t2.structure.name:
        raise ValueError
    if t1.color != t2.color:
        raise ValueError
    for g1, g2 in zip(t1.glues, t2.glues):
        equiv = fts.gluelist.domerge(equiv, g1, g2)
    return equiv


def _findpotentialtilemerges(fts, equiv):
    ppairs = []
    for ti in range(0, len(fts.tilelist.tiles)):
        t1 = fts.tilelist.tiles[ti]
        if not t1.used:
            continue
        ppairs += [(t1, t) for t in fts.tilelist.tiles[ti:]
                   if (t1.color == t.color) and (
                       t1.structure.name == t.structure.name)]
    shuffle(ppairs)
    return ppairs


def _findpotentialgluemerges(fts, equiv):
    ppairs = []
    for g1 in np.arange(0, len(fts.gluelist.strength)):
        g2s = g1 + 1 + np.nonzero(
            (fts.gluelist.strength[g1 + 1:] == fts.gluelist.strength[g1]) &
            (fts.gluelist.structure[g1 + 1:] == fts.gluelist.structure[g1]))[0]
        ppairs += [(g1, g2) for g2 in g2s]
    shuffle(ppairs)
    return ppairs


def _recfix(fts, equiv, tp, initptins):
    equiv = tilemerge(fts, equiv, fts.tilelist.totile[tp[0]],
                      fts.tilelist.totile[tp[1]])
    ae, badpair = isatamequiv(fts, equiv, initptins=initptins)
    if ae:
        return equiv
    elif badpair is None:
        raise ValueError
    else:
        return _recfix(fts, equiv, badpair, initptins)


def _tilereduce(fts, equiv, initptins=None):
    todo = _findpotentialtilemerges(fts, equiv)
    if initptins is None:
        initptins = ptins(fts, equiv)
    for t1, t2 in todo:
        try:
            nequiv = tilemerge(fts, equiv, t1, t2)
        except ValueError:
            continue
        ae, badpair = isatamequiv(fts, nequiv, initptins=initptins)
        if ae:
            equiv = nequiv
        elif badpair is None:
            continue
        else:
            try:
                equiv = _recfix(fts, equiv, badpair, initptins)
            except ValueError:
                continue
            except KeyError:
                continue
    return equiv


def _gluereduce(fts, equiv, initptins=None):
    todo = _findpotentialgluemerges(fts, equiv)
    if initptins is None:
        initptins = ptins(fts, equiv)
    for g1, g2 in todo:
        try:
            nequiv = fts.gluelist.domerge(equiv, g1, g2)
        except ValueError:
            continue
        ae, badpair = isatamequiv(fts, nequiv, initptins=initptins)
        if ae:
            equiv = nequiv
        elif badpair is None:
            continue
        else:
            try:
                equiv = _recfix(fts, equiv, badpair, initptins)
            except ValueError:
                continue
            except KeyError:
                continue
    return equiv
