# End reordering code to optimize placement.

import anneal
import sensitivity as sens
import stickydesign as sd
import numpy.random as random
import numpy as np
import random as pyrand
from copy import deepcopy

def flatten(seq):
    for item in seq:
        if is_string_like(item) or not iterable(item):
            yield item
        else:
            for subitem in flatten(item):
                yield subitem

def ecomp(x):
    if x[-1]=='/':
        return x[:-1]
    else:
        return x+'/'

class FseqState:
    def __init__(self, seqs=None):
        if not seqs:
            self.seqs = {}
        else:
            self.seqs = seqs
    def copy(self):
        return FseqState({'DT': self.seqs['DT'].copy(), 'TD': self.seqs['TD'].copy()})

class EndSystemFseq:
    def __init__(self, tilesys, pairs=None):
        # Set up variables, etc.
        self.ef = sd.energyfuncs_santalucia(mismatchtype='max')
        tilesys = deepcopy(tilesys)
        self.ends = tilesys['ends']
        self.tiles = tilesys['tiles']
        self.tilesystem = tilesys
        
        if not pairs:
            pairs = sens.consolidate_pairs( sens.senspairs(tilesys), comcomp=1, onlytop=True )

        self.initstate = FseqState()
        self.names = {}
        fseqsTD, self.names['TD'] = (list(x) for x in zip(*[ [end['fseq'].lower(),end['name']] for end in self.ends if end['type'] == 'TD']))
        fseqsDT, self.names['DT'] = (list(x) for x in zip(*[ [end['fseq'].lower(),end['name']] for end in self.ends if end['type'] == 'DT']))
        self.initstate.seqs['TD'] = sd.endarray(fseqsTD,'TD')
        self.initstate.seqs['DT'] = sd.endarray(fseqsDT,'DT')
        self.enlocs = {}
        for i, endn in enumerate(self.names['TD']):
            self.enlocs[endn] = (i,'TD')
        for i, endn in enumerate(self.names['DT']):
            self.enlocs[endn] = (i,'DT')        
        
        # Get the mean non-spurious interaction
        self.meangse = 0.5*( np.mean(self.ef.matching_uniform( self.initstate.seqs['TD'] ))+np.mean(self.ef.matching_uniform( self.initstate.seqs['DT'] )) )
        self.mult = {'1NGO': np.exp(-2.0*self.meangse), '2NGO': np.exp(-1.65*self.meangse), '1GO': np.exp(-1.5*self.meangse), '2GO': np.exp(-1.1*self.meangse)}
        
        self.pairdict = {}
        for pairclass,memberset in pairs.items():
            for x,y in memberset:
                self.pairdict[(x,ecomp(y))] = pairclass
                self.pairdict[(y,ecomp(x))] = pairclass
        
        
    def mutate(self, state):
        # Start by deciding to swap TD or DT ends.
        if random.rand() > 1.0*len(state.seqs['TD'])/(len(state.seqs['DT'])+len(state.seqs['TD'])):
            t = 'DT'
        else:
            t = 'TD'
        
        en = len(state.seqs[t])
        
        a = random.randint( 0, en )
        b = random.randint( 0, en )
        state.seqs[t][[a,b],:] = state.seqs[t][[b,a],:]

    def score(self, state):
        
        sc = 0.0
        
        for (xn,yn),pairclass in self.pairdict.items():
            
            # set comp flags
            xc = False
            yc = False
            if xn[-1] == '/':
                xc = True
                xn = xn[:-1]
            if yn[-1] == '/':
                yc = True
                yn = yn[:-1]
                
            # get end indexes and types
            xi,xt = self.enlocs[xn]
            yi,yt = self.enlocs[yn]
            #print "%s, %s: (%s %s)" % (xn,yn,xt,yt)
            # skip if not same type
            if xt != yt: continue
            
            if yc and xc:
                val = self.ef.uniform(state.seqs[xt][xi:xi+1].comps,state.seqs[yt][yi:yi+1].comps)[0]
            elif xc:
                val = self.ef.uniform(state.seqs[xt][xi:xi+1].comps,state.seqs[yt][yi:yi+1].ends)[0]
            elif yc:
                val = self.ef.uniform(state.seqs[xt][xi:xi+1].ends,state.seqs[yt][yi:yi+1].comps)[0]
            else:
                val = self.ef.uniform(state.seqs[xt][xi:xi+1].ends,state.seqs[yt][yi:yi+1].ends)[0]
            
            sc += self.mult[pairclass]*np.exp( val )
    
        return sc   

            
wcd = {   'a': 't',
         'b': 'v',
         'c': 'g',
         'd': 'h',
         'g': 'c',
         'h': 'd',
         'k': 'm',
         'm': 'k',
         'n': 'n',
         's': 's',
         't': 'a',
         'v': 'b',
         'w': 'w' }
         
def wc(seqstr):
    return ''.join(wcd[x] for x in reversed(seqstr))
