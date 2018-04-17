from .util import GlueMergeSpec, TileMergeSpec, comp
from .tiles import TileList
import random
from . import endreduce as er
from itertools import combinations
from . import latticedefect as ld
from . import sensitivityprofiles as sp
import logging

log = logging.getLogger(__name__)


def eq_noio(t1, t2, m=GlueMergeSpec([])):
    return (t1.structure.name == t2.structure.name) and all(
        m.eq(x, y) for x, y in zip(t1.ends, t2.ends))


def input_eq_siminput(a, b, m=GlueMergeSpec([])):
    # FIXME: COULD BE BETTER
    for d, g, io in zip(a.structure._dirs, a.ends, a['input']):
        if not io:
            continue
        else:
            if not any([
                    m.eq(g, g2) for d2, g2 in zip(b.structure._dirs, b.ends)
                    if d == d2
            ]):
                return False
    return True


def mergeproblem(ts, m):
    for a in ts.tiles:
        for b in ts.tiles + sum([x.rotations for x in ts.tiles], TileList()):
            if input_eq_siminput(a, b):
                continue
            if input_eq_siminput(a, b, m):
                if eq_noio(a, b, m) and (('label' in a) == ('label' in b)):
                    continue
                else:
                    return (a, b)
    return None


def applymerge(ts, m, tm=TileMergeSpec([])):
    for i in range(0, len(m._ecs) // 2):
        for x in m._ecs[2 * i]:
            ts = er.equate_pair(ts, ('merged{}'.format(i), x), doseed=True)
        for x in m._ecs[2 * i + 1]:
            ts = er.equate_pair(ts, ('merged{}/'.format(i), x), doseed=True)
    for tc in tm._ecs:
        z = list(tc)
        for tn in z[1:]:
            try:
                ts.tiles[tn]['fake'] = z[0]
            except:
                continue
    return ts


def mergetiles(ts,
               t,
               m2,
               tm2,
               checkprofiles=None,
               oldclasses=None,
               checkld=False):
    t1, t2 = t
    if t1.structure.name != t2.structure.name:
        return None, None
    if t1.structure._endtypes != t2.structure._endtypes:
        return None, None
    if ('label' in t1.keys()) != ('label' in t2.keys()):
        return None, None
    try:
        for e1, e2, gt1, gt2 in zip(t1['ends'], t2['ends'],
                                    t1.structure._endtypes,
                                    t2.structure._endtypes):
            assert gt1 == gt2
            m2 = m2.copyadd(e1, e2)
    except ValueError:
        return None, None
    tm2 = tm2.copyadd(t1.name, t2.name)
    mp = mergeproblem(ts, m2)
    if (not mp) and checkprofiles:
        mp = sp.sensitivity_profiles_fakesingles(ts, 3, m2, tm2, oldclasses,
                                                 checkprofiles, True)
    if not mp:
        if checkld:
            oldld = ld.latticedefects(ts, direction='e', depth=2)
            newld = ld.latticedefects(ts, direction='e', depth=2, gms=m2)
            if (len(newld) > len(oldld)) or (len(
                [x for x in newld if x not in oldld]) > 0):
                log.debug("lattice defect (e): {}".format((t1.name, t2.name)))
                return False, False

            oldld = ld.latticedefects(ts, direction='w', depth=2)
            newld = ld.latticedefects(ts, direction='w', depth=2, gms=m2)
            if (len(newld) > len(oldld)) or (len(
                [x for x in newld if x not in oldld]) > 0):
                log.debug("lattice defect (w): {}".format((t1.name, t2.name)))
                return False, False
        return m2, tm2
    else:
        return mergetiles(ts, mp, m2, tm2, checkprofiles, oldclasses, checkld)


def newtilereduce(ts, checkprofiles=None, oldclasses=None, checkld=False):
    rti = ts.tiles + sum([x.rotations for x in ts.tiles], TileList())
    tpairs = [(x, y) for x in ts.tiles for y in rti]
    random.shuffle(tpairs)
    m = GlueMergeSpec([])
    tm = TileMergeSpec([])
    for t1, t2 in tpairs:
        if tm.eq(t1.name, t2.name):
            continue
        m2, tm2 = mergetiles(
            ts, (t1, t2), m, tm, checkprofiles, oldclasses, checkld=checkld)
        if m2:
            m = m2
            tm = tm2
            log.debug("tile reduction: {} {}".format(t1.name, t2.name))
            tm.add(t1.name, t2.name)
    log.info("tilereduced {} - {} tiles, {} - {} ends".format(
        len(ts.tiles),
        sum(len(x) - 1 for x in tm._ecs),
        len(ts.allends), sum(len(x) - 1 for x in m._ecs) // 2))
    return m, tm


def newgluereduce(ts,
                  m=GlueMergeSpec([]),
                  tm=TileMergeSpec([]),
                  checkprofiles=None,
                  oldclasses=None,
                  checkld=False):
    pts = - sum(len(x)-1 for x in tm._ecs)
    pgs = - sum(len(x)-1 for x in m._ecs)//2
    potentials = list(er.find_nonsens_pairs(ts))
    random.shuffle(potentials)
    while len(potentials) > 0:
        pp = tuple(potentials.pop())
        if m.eq(pp[0], pp[1]):
            continue
        try:
            m2 = m.copyadd(*pp)
        except ValueError:
            continue
        mp = mergeproblem(ts, m2)
        if (not mp) and checkprofiles:
            mp = sp.sensitivity_profiles_fakesingles(ts, 3, m2, tm, oldclasses,
                                                     checkprofiles, True)
        if (not mp):
            if checkld:
                oldld = ld.latticedefects(ts, direction='e', depth=2)
                newld = ld.latticedefects(ts, direction='e', depth=2, gms=m2)
                if (len(newld) > len(oldld)) or (len(
                    [x for x in newld if x not in oldld]) > 0):
                    log.debug("lattice defect (e): {}".format(pp))
                    continue

                oldld = ld.latticedefects(ts, direction='w', depth=2)
                newld = ld.latticedefects(ts, direction='w', depth=2, gms=m2)
                if (len(newld) > len(oldld)) or (len(
                    [x for x in newld if x not in oldld]) > 0):
                    log.debug("lattice defect (w): {}".format(pp))
                    continue
            log.debug(str(pp))
            m = m2
        else:
            m2, tm2 = mergetiles(ts, mp, m2, tm, checkprofiles, oldclasses,
                                 checkld)
            if m2:
                log.debug(str(m2))
                m = m2
                tm = tm2
    log.info("tilereduced {} - {} (-{} p) tiles, {}  - {} (-{} p) ends".format(
        len(ts.tiles),
        sum(len(x) - 1 for x in tm._ecs),
        pts,
        len(ts.allends),
        sum(len(x) - 1 for x in m._ecs) // 2,
        pgs))
    return m, tm
