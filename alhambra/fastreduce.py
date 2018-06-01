from collections import namedtuple
from copy import copy
import numpy as np
from . import tilestructures as ts
from .tiles import TileList
from .ends import End
from random import shuffle
from . import util
from . import fastlatticedefect as fld
import logging

FTile = namedtuple("FTile", ("color", "use", "glues", "name", "used",
                             "structure", "dfake", "sfake"))

FTilesArray = namedtuple("FTilesArray",
                         ("color", "use", "glues", "name", "used", "structure",
                          "dfake", "sfake"))

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
            glues = np.array([gluelist.tonum[x] for x in t.ends])
            if 'fake' in t.keys():
                continue
            if 'input' not in t.keys():
                used = False
                use = np.array([U_UNUSED for _ in t.ends])
            else:
                used = True
                use = np.array([int(x) for x in t['input']])
            color = 'label' in t.keys()
            self.tiles.append(
                FTile(
                    name=t.name,
                    color=color,
                    use=use,
                    glues=glues,
                    used=used,
                    structure=t.structure,
                    dfake=0,
                    sfake=False))
            assert t.name not in self.totile.keys()
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


def _ffakedouble(tile, sta, gluelist, outputonly=True, dir4=False, equiv=None):
    if equiv is None:
        equiv = gluelist.blankequiv()
    if not dir4:
        faketiles = ([], [])
        fti = FTI
    else:
        faketiles = ([], [], [], [])
        fti = (0, 1, 2, 3)
    for dir in range(0, 4):
        if (tile.use[dir] == 1) and outputonly:
            continue
        if outputonly:
            oti = np.nonzero((equiv[gluelist.complement[tile.glues[dir]]] ==
                              equiv[sta.glues[:, RSEL[dir]]]) &
                             (sta.use[:, RSEL[dir]] == 1) & (sta.used))
        else:
            oti = np.nonzero((equiv[gluelist.complement[tile.glues[dir]]] ==
                              equiv[sta.glues[:, RSEL[dir]]]))
        for i in oti[0]:
            faketiles[fti[dir]].append(
                FTile(
                    color=False,
                    used=True,
                    name=tile.name + '_{}_'.format(dir) + sta.name[i],
                    structure=FTS[dir],
                    glues=np.array(_fdg(dir, tile.glues, sta.glues[i])),
                    use=np.array(_fdg(dir, tile.use, sta.use[i])),
                    dfake=dir + 1,
                    sfake=(tile.sfake or sta.sfake[i])))
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
    glues = [[([fdb] + list(ftile.glues))[x] for x in y] for y in sel]
    fuse = [[([-1] + list(ftile.use))[x] for x in y] for y in sel]
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
            dfake=0,
            sfake=True) for u, g, n, ud in zip(use, glues, names, used)
    ]


def _ft_to_fta(ftiles):
    return FTilesArray(
        color=np.array([x.color for x in ftiles]),
        name=np.array([x.name for x in ftiles]),
        glues=np.array([x.glues for x in ftiles]),
        use=np.array([x.use for x in ftiles]),
        used=np.array([x.used for x in ftiles]),
        structure=np.array([x.structure.name for x in ftiles]),
        dfake=np.array([x.dfake for x in ftiles]),
        sfake=np.array([x.sfake for x in ftiles]))


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
        ts['info'] = ts.get('info', dict())
        ts['info']['fgluemerge'] = ts['info'].get('fgluemerge', list())
        ts['info']['fgluemerge'].append([int(x) for x in equiv])
        return ts

    def togluemergespec(self, ts, equiv):
        gms = util.GlueMergeSpec()
        for i in range(0, len(equiv)):
            if i != equiv[i]:
                gms.add(self.gluelist.name[equiv[i]], self.gluelist.name[i])
        return gms


def ptins(fts, equiv=None, tau=2):
    """Calculate potential tile attachments to input neighborhoods"""
    ptins = []
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
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
            matches = np.nonzero((strs >= tau) & (ta.dfake == 0))
            # indices of matching
            # (excludes fake doubles, which can't actually attach
            # (because they are
            # actually two singles) to their own local neighborhoods)
            ptin.append(matches)
        ptins.append(ptin)
    return ptins


def is_2go_single_nn(fts,
                     tn,
                     un,
                     equiv,
                     tau=2,
                     in2go=None,
                     retall=False,
                     also22go=False):
    # Before starting, check to see if t and u are actually the same:
    if np.all(equiv[fts.tilelist.stiles.glues[fts.tilelist.stiles.used][tn]] ==
              equiv[fts.tilelist.stiles.glues[un]]):
        if not also22go:
            if not retall:
                return False
            else:
                return False, None
        else:
            if not retall:
                return False, False
            else:
                return False, False, None, None

    if in2go is None:
        in2go = _ffakedouble(
            fta_to_ft(fts.tilelist.stiles, tn, fts.tilelist.stiles.used),
            fts.tilelist.stiles,
            fts.gluelist,
            outputonly=True,
            dir4=True)
    incorrect = _ffakedouble(
        fta_to_ft(fts.tilelist.stiles, un),
        fts.tilelist.stiles,
        fts.gluelist,
        outputonly=False,
        dir4=True,
        equiv=equiv)
    is2go = None
    for tvs, uws in zip(in2go, incorrect):
        for tv in tvs:
            for uw in uws:
                if np.sum(fts.gluelist.strength[tv.glues[(
                        equiv[tv.glues] == equiv[uw.glues]) & (
                            tv.use == 1)]]) >= tau:
                    if also22go and uw.sfake:
                        if not retall:
                            return True, True
                        else:
                            return True, True, (ts, uw), None
                    elif also22go:
                        is2go = (tv, uw)
                        for dirr, (tvs2,
                                   uws2) in enumerate(zip(in2go, incorrect)):
                            if dirr + 1 == tv.dfake:
                                continue
                            for tv2 in tvs2:
                                for uw2 in uws2:
                                    if np.sum(fts.gluelist.strength[tv2.glues[(
                                            equiv[tv2.glues] == equiv[uw2.glues]
                                    ) & (tv2.use == 1)]]) >= tau:
                                        if not retall:
                                            return True, True
                                        else:
                                            return True, True, (tv, uw), (tv2,
                                                                          uw2)
                    elif not retall:
                        return True
                    else:
                        return True, (tv, uw)
    if also22go:
        if retall:
            return (not (is2go is None)), False, is2go, None
        else:
            return (not (is2go is None)), False
    elif not retall:
        return False
    else:
        return False, None


def gen_2go_single_ins(fts, tau=2):
    ssel = np.nonzero(fts.tilelist.stiles.used)
    ins = []
    for tn in ssel[0]:
        ins.append(
            _ffakedouble(
                fta_to_ft(fts.tilelist.stiles, tn),
                fts.tilelist.stiles,
                fts.gluelist,
                outputonly=True,
                dir4=True))
    return ins


def gen_2go_profile(fts, equiv=None, ins2go=None, also22go=False):
    if ins2go is None:
        ins2go = gen_2go_single_ins(fts, tau=2)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    sens1s = ptins(fts, tau=1, equiv=equiv)[0]
    sens2s = []
    sens22s = []
    for tn, (uns, in2go) in enumerate(zip(sens1s, ins2go)):
        x2 = []
        x22 = []
        for un in uns[0]:
            if also22go:
                s2, s22 = is_2go_single_nn(
                    fts, tn, un, equiv, tau=2, in2go=in2go, also22go=True)
            else:
                s2 = is_2go_single_nn(
                    fts, tn, un, equiv, tau=2, in2go=in2go, also22go=False)
            if s2:
                x2.append(un)
            if also22go and s22:
                x22.append(un)
        sens2s.append(x2)
        if also22go:
            sens22s.append(x22)

    if also22go:
        return sens2s, sens22s
    else:
        return sens2s


def is_2go_equiv(fts, equiv=None, ins2go=None, origsens=None):
    if ins2go is None:
        ins2go = gen_2go_single_ins(fts, tau=2)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    if origsens is None:
        origsens = gen_2go_profile(fts)
    sens1s = ptins(fts, tau=1, equiv=equiv)[0]

    for tn, (uns, in2go, os) in enumerate(zip(sens1s, ins2go, origsens)):
        for un in uns[0]:
            r, p = is_2go_single_nn(
                fts, tn, un, equiv, tau=2, in2go=in2go, retall=True)
            if (r and (un not in os)):
                # Next if checks to see if sensitive tile is actually identical
                # to an already sensitive tile:
                # FIXME: this would not be necessary if we stored whether a
                # tile was a duplicate.  However, that would break the
                # glue-equiv-only method...
                if not np.any(
                        np.all(
                            equiv[fts.tilelist.stiles.glues[os]] ==
                            equiv[fts.tilelist.stiles.glues[un]],
                            axis=1)):
                    return False, (
                        fts.tilelist.stiles.name[fts.tilelist.stiles.used][tn],
                        fts.tilelist.stiles.name[un]), p
    return True, None, None


def is_22go_equiv(fts, equiv=None, ins2go=None, orig22go=None):
    if ins2go is None:
        ins2go = gen_2go_single_ins(fts, tau=2)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    if orig22go is None:
        _, orig22go = gen_2go_profile(fts, also22go=True)
    sens1s = ptins(fts, tau=1, equiv=equiv)[0]

    for tn, (uns, in2go, os22) in enumerate(zip(sens1s, ins2go, orig22go)):
        for un in uns[0]:
            r2, r22, p2, p22 = is_2go_single_nn(
                fts,
                tn,
                un,
                equiv,
                tau=2,
                in2go=in2go,
                also22go=True,
                retall=True)
            if (r22 and (un not in os22)):
                # Next if checks to see if sensitive tile is actually identical
                # to an already sensitive tile:
                # FIXME: this would not be necessary if we stored whether a
                # tile was a duplicate.  However, that would break the
                # glue-equiv-only method...
                if not np.any(
                        np.all(
                            equiv[fts.tilelist.stiles.glues[os22]] ==
                            equiv[fts.tilelist.stiles.glues[un]],
                            axis=1)):
                    return False, (
                        fts.tilelist.stiles.name[fts.tilelist.stiles.used][tn],
                        fts.tilelist.stiles.name[un]), p22
    return True, None, None


def fta_to_ft(stiles, un, sel=None):
    if sel is None:
        sel = np.ones_like(stiles.used, dtype=bool)
    return FTile(
        color=stiles.color[sel][un],
        use=stiles.use[sel, :][un, :],
        glues=stiles.glues[sel, :][un, :],
        name=stiles.name[sel][un],
        used=stiles.used[sel][un],
        structure=stiles.structure[sel][un],
        dfake=stiles.dfake[sel][un],
        sfake=stiles.sfake[sel][un])


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
                    (equiv[tl.glues[xx]] == equiv[tl.glues[yy]]),
                    axis=1) | ~(tl.color[xx] == tl.color[yy])
                if np.any(mm):
                    return False, (tl.name[xx[0][0]], tl.name[yy[0][mm][0]])
            elif len(xx[0]) == 0:
                return False, None
            else:
                raise NotImplementedError(xx)
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
        ppairs += [
            (t1, t) for t in fts.tilelist.tiles[ti:]
            if (t1.color == t.color) and (
                t1.structure.name == t.structure.name) and (t1.name != t.name)
        ]
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


def _recfix(fts,
            equiv,
            tp,
            initptins,
            check2go=False,
            ins2go=None,
            orig2go=None,
            orig22go=None,
            check22go=False,
            checkld=False,
            chain=None):
    log = logging.getLogger(__name__)
    if chain is None:
        chain = []
    equiv = tilemerge(fts, equiv, fts.tilelist.totile[tp[0]],
                      fts.tilelist.totile[tp[1]])
    ae, badpair = isatamequiv(fts, equiv, initptins=initptins)
    if check2go and ae:
        ae, badpair, _ = is_2go_equiv(
            fts, equiv, ins2go=ins2go, origsens=orig2go)
    if ae and check22go:
        ae, badpair, _ = is_22go_equiv(
            fts, equiv, ins2go=ins2go, orig22go=orig22go)
    if ae and checkld:
        ld_e = fld.latticedefects(fts, 'e', equiv=equiv)
        ld_w = fld.latticedefects(fts, 'w', equiv=equiv)
        if (len(ld_e) > 0) or (len(ld_w) > 0):
            raise ValueError
    if ae:
        log.debug(
            "Recfix succeeds with {}, {}, {}".format(tp[0], tp[1], chain))
        return equiv
    elif badpair is None:
        raise ValueError
    else:
        chain.append(tp)
        return _recfix(
            fts,
            equiv,
            badpair,
            initptins,
            check2go,
            ins2go,
            orig2go,
            orig22go=orig22go,
            check22go=check22go,
            checkld=checkld,
            chain=chain)


def _tilereduce(fts,
                equiv=None,
                check2go=False,
                initptins=None,
                check22go=False,
                checkld=False):
    log = logging.getLogger(__name__)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    todo = _findpotentialtilemerges(fts, equiv)
    if initptins is None:
        initptins = ptins(fts, equiv)
    if check2go or check22go:
        ins2go = gen_2go_single_ins(fts)
        origsens, orig22go = gen_2go_profile(fts, ins2go=ins2go, also22go=True)
    else:
        ins2go = None
        orig22go = None
        origsens = None
    for todoi, (t1, t2) in enumerate(todo):
        try:
            nequiv = tilemerge(fts, equiv, t1, t2)
        except ValueError:
            continue
        ae, badpair = isatamequiv(fts, nequiv, initptins=initptins)
        if ae and check2go:
            ae, badpair, _ = is_2go_equiv(
                fts, nequiv, ins2go=ins2go, origsens=origsens)
        if ae and check22go:
            ae, badpair, _ = is_22go_equiv(
                fts, nequiv, ins2go=ins2go, orig22go=orig22go)
        if ae and checkld:
            ld_e = fld.latticedefects(fts, 'e', equiv=nequiv)
            ld_w = fld.latticedefects(fts, 'w', equiv=nequiv)
            if (len(ld_e) > 0) or (len(ld_w) > 0):
                continue
        if ae:
            equiv = nequiv
            log.debug("Reduced {}, {} ({} of {} pairs done)".format(
                t1.name, t2.name, todoi, len(todo)))
        elif badpair is None:
            continue
        else:
            try:
                equiv = _recfix(
                    fts,
                    equiv,
                    badpair,
                    initptins,
                    check2go,
                    ins2go,
                    origsens,
                    check22go=check22go,
                    orig22go=orig22go,
                    checkld=checkld,
                    chain=[(t1.name, t2.name)])
            except ValueError:
                continue
            except KeyError:
                continue
    return equiv


def _gluereduce(fts,
                equiv=None,
                check2go=False,
                check22go=False,
                checkld=False,
                initptins=None):
    log = logging.getLogger(__name__)
    if equiv is None:
        equiv = fts.gluelist.blankequiv()
    todo = _findpotentialgluemerges(fts, equiv)
    if initptins is None:
        initptins = ptins(fts)
    if check2go or check22go:
        ins2go = gen_2go_single_ins(fts)
        origsens, orig22go = gen_2go_profile(fts, ins2go=ins2go, also22go=True)
    else:
        ins2go = None
        origsens = None
        orig22go = None
    for todoi, (g1, g2) in enumerate(todo):
        try:
            nequiv = fts.gluelist.domerge(equiv, g1, g2)
        except ValueError:
            continue
        ae, badpair = isatamequiv(fts, nequiv, initptins=initptins)
        if ae and check2go:
            ae, badpair, _ = is_2go_equiv(
                fts, nequiv, ins2go=ins2go, origsens=origsens)
        if ae and check22go:
            ae, badpair, _ = is_22go_equiv(
                fts, nequiv, ins2go=ins2go, orig22go=orig22go)
        if ae and checkld:
            ld_e = fld.latticedefects(fts, 'e', equiv=nequiv)
            ld_w = fld.latticedefects(fts, 'w', equiv=nequiv)
            if (len(ld_e) > 0) or (len(ld_w) > 0):
                continue
        if ae:
            equiv = nequiv
            log.debug("Glue reduction: {}, {} ({}/{} done)".format(
                fts.gluelist.name[g1], fts.gluelist.name[g2], todoi, len(
                    todo)))
        elif badpair is None:
            continue
        else:
            try:
                equiv = _recfix(
                    fts,
                    equiv,
                    badpair,
                    initptins,
                    check2go,
                    ins2go,
                    origsens,
                    check22go=check22go,
                    checkld=checkld,
                    chain=[(fts.gluelist.name[g1], fts.gluelist.name[g2])])
            except ValueError:
                continue
            except KeyError:
                continue
    return equiv
