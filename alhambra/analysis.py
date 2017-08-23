from . import stickyends
from . import tiletypes
import stickydesign.plots as sdplots
import stickydesign as sd
from collections import Counter
from stickydesign import EnergeticsDAOE
from matplotlib import pylab
import numpy as np

DEFAULT_REGION_ENERGETICS = EnergeticsDAOE(temperature=33,
                                           coaxparams=False,
                                           danglecorr=False)


def plot_se_hists(tileset,
                  all_energetics=None,
                  energetics_names=None,
                  title=None, **kwargs):

    if all_energetics is None:
        all_energetics = stickyends.DEFAULT_MULTIMODEL_ENERGETICS

    if energetics_names is None:
        energetics_names = stickyends.DEFAULT_MM_ENERGETICS_NAMES

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

    regions = [tiletypes.tfactory.parse(t)._side_bound_regions
               for t in tileset['tiles']]
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

    regions = [tiletypes.tfactory.parse(t)._short_bound_full
               for t in tileset['tiles']]
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
