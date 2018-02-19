from . import sensitivitynew as sens
from .util import comp, base
from random import shuffle
import logging
from .seq import _WC
import re

log = logging.getLogger()


def getsimadjs(ts, adjlike):
    endtomimic = ts.allends[adjlike]
    assert endtomimic.fseq
    adjs = (endtomimic.fseq[0], endtomimic.fseq[-1])
    logging.debug(adjs)
    allowed = []
    for end in ts.allends:
        if not end.fseq:
            continue
        if end.etype != endtomimic.etype:
            continue
        if (end.fseq[0], end.fseq[-1]) == adjs:
            allowed.append(end.name)
        if (_WC[end.fseq[-1]], _WC[end.fseq[0]]) == adjs:
            allowed.append(end.name+'/')
    logging.debug(allowed)
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


#def check_pair_sensitivity(ts, sc=None):
#    pass


def find_nonsens_pairs(ts, sc=None):
    if sc is None:
        sc = ts.sensitivity_classes()

    allsenspairs = set.union(*sc.values())

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


def check_changes(oldts, newts, equatedpair, oldsc=None, newsc=None):
    if oldsc is None:
        oldsc = oldts.sensitivity_classes()
    if newsc is None:
        newsc = newts.sensitivity_classes()
    equatedpair = tuple(equatedpair)
    eremain, ereplace = equatedpair

    for pair in newsc['2GO']:
        if pair not in oldsc['2GO']:
            if eremain in pair:
                log.debug("in replace: {}".format(pair))
                pair = frozenset.union(pair-{eremain},{ereplace})
                log.debug("now: {}".format(pair))
                if pair not in oldsc['2GO']:
                    return False
                else:
                    log.debug("caught by replace")
            elif comp(eremain) in pair:
                log.debug("in creplace: {}".format(pair))
                pair = frozenset.union(pair-{comp(eremain)},{comp(ereplace)})
                log.debug("now: {}".format(pair))
                if pair not in oldsc['2GO']:
                    return False
                else:
                    log.debug("caught by compreplace")
            else:
                log.debug((pair,ereplace,eremain))
                return False

    # FIXME: add lattice defect check
    return True


def equate_pair(ts, pair):
    pair = tuple(pair)  # ensure that the pair is actually ordered.
    e1, e2 = pair

    swap = not ((e1[-1] == '/' and e2[-1] == '/') or
                (e1[-1] != '/' and e2[-1] != '/'))
    newts = ts.copy()
    r = re.compile(base(e2) + r'(/?)$')

    def rfunc(match):
        if swap:
            return comp(base(e1) + match.group(1))
        else:
            return base(e1) + match.group(1)

    for t in newts.tiles:
        t.ends = [r.sub(rfunc, x) for x in t.ends]
    for t in newts.seed['adapters']:
        t['ends'] = [r.sub(rfunc, x) for x in t['ends']]
    return newts

def _update_nonsens_pairs(potentials, pair, trialsc):
    def rename(x):
        if pair[1] in x:
            return frozenset.union(x-{pair[1]},{pair[0]})
        elif comp(pair[1]) in x:
            return frozenset.union(x-{comp(pair[1])},{comp(pair[0])})
        else:
            return x
    p1 = [rename(x) for x in potentials]
    return [x for x in p1 if x not in set.union(*trialsc.values()) and
            len({base(y) for y in x}) == 2]

def fast_reduce_ends(ts):
    oldts = ts
    sc = ts.sensitivity_classes()
    potentials = find_nonsens_pairs(oldts, sc)
    removedpairs = []

    log.info("start removal: {} ends, {} potential pairs".format(
        len(ts.allends), len(potentials)))

    while len(potentials) > 0:
        pair = tuple(potentials.pop())
        trialts = equate_pair(oldts, pair)
        trialsc = trialts.sensitivity_classes()
        if check_changes(oldts, trialts, pair, sc, trialsc):
            oldts = trialts
            sc = trialsc
            potentials = _update_nonsens_pairs(potentials, pair, trialsc)
            removedpairs.append(pair)
            if log.isEnabledFor(logging.INFO):
                log.info("removed {}. {} ends remain, {} potential pairs".format(
                    pair, len(oldts.allends), len(potentials)))

    return oldts, removedpairs


def reduce_ends(ts):
    oldts = ts
    sc = ts.sensitivity_classes()
    potentials = find_nonsens_pairs(oldts, sc)
    removedpairs = []

    log.info("start removal: {} ends, {} potential pairs".format(
        len(ts.allends), len(potentials)))

    while len(potentials) > 0:
        pair = potentials.pop()
        trialts = equate_pair(oldts, pair)
        trialsc = trialts.sensitivity_classes()
        if check_changes(oldts, trialts, pair, sc, trialsc):
            oldts = trialts
            sc = trialsc
            potentials = find_nonsens_pairs(oldts, sc)
            removedpairs.append(pair)
            if log.isEnabledFor(logging.INFO):
                log.info("removed {}. {} ends remain, {} potential pairs".format(
                    pair, len(oldts.allends), len(potentials)))

    return oldts, removedpairs

def attempt_end_removal(ts, end, adjlike=None):
    oldts = ts
    sc = ts.sensitivity_classes()
    potentials = find_potential_end_removal(ts, end, sc, adjlike)
    newts = None
    log.debug(potentials)
    while len(potentials) > 0:
        pair = potentials.pop()
        logging.debug(pair)
        trialts = equate_pair(ts, pair)
        if check_changes(oldts, trialts, pair, sc):
            newts = trialts
            break

    if newts is None:
        pair = None

    return newts, pair
