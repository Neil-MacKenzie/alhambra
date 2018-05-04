from . import sensitivitynew as sens
from .util import comp, base
from random import shuffle
import logging
from .seq import _WC
from .latticedefect import latticedefects
import re
import copy
from .tiles import TileList

log = logging.getLogger(__name__)


def getsimadjs(ts, adjlike):
    endtomimic = ts.allends[adjlike]
    assert endtomimic.fseq
    adjs = (endtomimic.fseq[0], endtomimic.fseq[-1])
    allowed = []
    for end in ts.allends:
        if not end.fseq:
            continue
        if end.etype != endtomimic.etype:
            continue
        if (end.fseq[0], end.fseq[-1]) == adjs:
            allowed.append(end.name)
        if (_WC[end.fseq[-1]], _WC[end.fseq[0]]) == adjs:
            allowed.append(end.name + '/')
    return allowed


def find_potential_end_removal(ts, end, sc=None, adjlike=None):
    allpotential = find_nonsens_pairs(ts, sc)

    # Note: this is a bit of a hack, returning a set of the partner
    # in each potential pair including end.  It uses min to return
    # the sole element of a frozenset.
    partners = {
        min(pair.difference({end}))
        for pair in allpotential if end in pair
    }
    if adjlike is not None:
        allowed = getsimadjs(ts, adjlike)
        partners = {partner for partner in partners if partner in allowed}

    pairlist = [(partner, end) for partner in partners]
    shuffle(pairlist)
    return pairlist



def find_nonsens_pairs(ts, sc=None, _1go_only=False):
    if sc is None:
        sc = ts.sensitivity_classes()

    if not _1go_only:
        allsenspairs = set.union(*sc.values())
    else:
        allsenspairs = set(sc['1GO'])

    ae = ts.allends

    potentialpairs = set.union({
        frozenset((x.name, y.name))
        for x in ae for y in ae if x != y and x.etype == y.etype
    }, {
        frozenset((x.name, comp(y.name)))
        for x in ae for y in ae if x != y and x.etype == y.etype
    })

    potentialpairs.difference_update(allsenspairs)

    return potentialpairs


def check_changes(oldts,
                  newts,
                  equatedpair,
                  oldsc=None,
                  newsc=None,
                  checkld=False,
                  _classes=('2GO',)):
    if oldsc is None:
        oldsc = oldts.sensitivity_classes()
    if newsc is None:
        newsc = newts.sensitivity_classes()
    equatedpair = tuple(equatedpair)
    eremain, ereplace = equatedpair

    for sc in _classes:
        for pair in newsc[sc]:
            if pair not in oldsc[sc]:
                if eremain in pair:
                    pair = frozenset.union(pair - {eremain}, {ereplace})
                    if pair not in oldsc[sc]:
                        return False
                elif comp(eremain) in pair:
                    pair = frozenset.union(pair - {comp(eremain)},
                                           {comp(ereplace)})
                    if pair not in oldsc[sc]:
                        return False
                else:
                    return False

    # FIXME: add lattice defect check
    if checkld:
        oldld = latticedefects(oldts, direction='e', depth=checkld)
        newld = latticedefects(newts, direction='e', depth=checkld)
        if (len(newld) > len(oldld)) or (len(
            [x for x in newld if x not in oldld]) > 0):
            log.debug("lattice defect (e): {}".format(equatedpair))
            return False

        oldld = latticedefects(oldts, direction='w', depth=checkld)
        newld = latticedefects(newts, direction='w', depth=checkld)
        if (len(newld) > len(oldld)) or (len(
            [x for x in newld if x not in oldld]) > 0):
            log.debug("lattice defect (w): {}".format(equatedpair))
            return False
    return True


def equate_pair(ts, pair, unsafe=False, doseed=False):
    pair = tuple(pair)  # ensure that the pair is actually ordered.
    e1, e2 = pair

    swap = not ((e1[-1] == '/' and e2[-1] == '/') or
                (e1[-1] != '/' and e2[-1] != '/'))
    if not unsafe:
        newts = ts.copy()
    else:
        newts = copy.copy(ts)
        newts['tiles'] = TileList([copy.copy(t) for t in ts['tiles']])
        if doseed and ('seed' in ts.keys()):
            newts['seed'] = ts['seed'].copy()
            newts['seed']['adapters'] = newts['seed']['adapters'].copy()
    r = re.compile(r'^' + base(e2) + r'(/?)$')

    def rfunc(match):
        if swap:
            return comp(base(e1) + match.group(1))
        else:
            return base(e1) + match.group(1)

    for t in newts.tiles:
        t['ends'] = [r.sub(rfunc, x) for x in t['ends']]
    if doseed and newts.seed:
        for t in newts.seed['adapters']:
            if 'ends' in t.keys():
                t['ends'] = [r.sub(rfunc, x) for x in t['ends']]
    return newts


def reduce_ends(ts, checkld=False, _wraparound=False, _classes=('2GO',), _smo=2, _unsafe=False, _1go_only=False):
    
    if ('22GO' in _classes) or ('22NGO' in _classes):
        _smo = 3
    oldts = ts.copy()
    sc = ts.sensitivity_classes(_maxorder=_smo)
    potentials = list(find_nonsens_pairs(oldts, sc, _1go_only))
    removedpairs = []
    shuffle(potentials)

    log.info("start removal: {} ends, {} potential pairs".format(
        len(ts.allends), len(potentials)))

    if checkld and (checkld < 2):
        log.info("Checkld should be depth (int>=2), \
not bool or 1. Setting to 2")
        checkld = 2
    
    changed = False
    while len(potentials) > 0:
        pair = potentials.pop()
        trialts = equate_pair(oldts, pair, unsafe=_unsafe)
        trialsc = trialts.sensitivity_classes(_maxorder=_smo)
        if check_changes(oldts, trialts, pair, sc, trialsc, checkld=checkld, _classes=_classes):
            oldts = trialts
            sc = trialsc
            # Reset potentials to be only things we haven't visited yet
            # that are still available:
            newpairs = find_nonsens_pairs(oldts, sc)
            potentials = [
                x for x in potentials if x in newpairs
            ]
            if _wraparound:
                changed = True
            removedpairs.append(tuple(pair))
            if log.isEnabledFor(logging.INFO):
                log.info("removed {}. {} ends remain, {} potential pairs".
                         format(set(pair), len(oldts.allends), len(newpairs)))
        # But if this isn't anything, then reset:
        if (len(potentials) == 0) and changed:
            log.info("reset potentials")
            potentials = newpairs

    oldts.add_info('endreduce', {'checkld': checkld, 'wrap': _wraparound, 'removed': removedpairs})
            
    return oldts.copy()


def attempt_end_removal(ts, end, adjlike=None, checkld=False):
    oldts = ts
    sc = ts.sensitivity_classes()
    potentials = find_potential_end_removal(ts, end, sc, adjlike)
    newts = None
    while len(potentials) > 0:
        pair = potentials.pop()
        trialts = equate_pair(ts, pair)
        if check_changes(oldts, trialts, pair, sc, checkld=checkld):
            newts = trialts
            break

    if newts is None:
        pair = None

    return newts, pair
