from . import stickyends
from . import util
from .tilesets import TileSet
import stickydesign.plots as sdplots
import stickydesign as sd
from collections import Counter


def plot_se_hists(tileset,
                  all_energetics=None,
                  energetics_names=None,
                  title=None):

    if all_energetics is None:
        all_energetics = stickyends.DEFAULT_MULTIMODEL_ENERGETICS

    if energetics_names is None:
        energetics_names = stickyends.DEFAULT_MM_ENERGETICS_NAMES

    if title is None:
        # FIXME
        title = 'Title'

    td = sd.endarray([x['fseq'] for x in tileset['ends']
                      if x['type'] == 'TD'], 'TD')

    dt = sd.endarray([x['fseq'] for x in tileset['ends']
                      if x['type'] == 'DT'], 'DT')

    sdplots.hist_multi([td, dt],
                       all_energetics,
                       energetics_names,
                       title)

    
def plot_side_arms(tileset,
                   energetics=None):
    
