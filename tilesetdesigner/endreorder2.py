# End reordering code to optimize placement.

import anneal
import sensitivity as sens
import stickydesign as sd
import numpy.random as random
import numpy as np
import random as pyrand
from pylab import flatten
from copy import deepcopy

def ecomp(x):
    if x[-1]=='/':
        return x[:-1]
    else:
        return x+'/'

class State:
    def __init__(self, seqs={}):
        self.seqs = seqs
    def copy(self):
        return State({'DT': list(self.seqs['DT']), 'TD': list(self.seqs['TD'])})

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
        
                
    def create_avails(self):
        self.avail = {}
        self.avail['TD'] = sd.get_accept_set( 'TD', 5, 8.35, 0.1, 0.4, energyfuncs=self.ef )
        self.avail['DT'] = sd.get_accept_set( 'DT', 5, 8.35, 0.1, 0.4, energyfuncs=self.ef )
        
    def nsmutate(self,state):
        # Start by deciding to swap TD or DT ends.
        if random.rand() > 1.0*len(state.seqs['TD'])/(len(state.seqs['DT'])+len(state.seqs['TD'])):
            t = 'DT'
        else:
            t = 'TD'
        
        en = len(state.seqs[t])
        a = random.randint( 0, en )
        b = random.randint( 0, len(self.avail[t]) )
        
        state.seqs[t][a,:] = self.avail[t][b]

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
            
    def nsscore(self, state):
        
        sc = 0.0
        
        for t in ['TD','DT']:
            #print "0",sc
            match = self.ef.matching_uniform(state.seqs[t])
            sc += sum(np.exp( abs( match-np.mean(match) ) )-1)
            #print "A",sc
            ea = sd.energy_array_uniform( state.seqs[t], self.ef )
            ea1 = ea[len(ea)/2:,:]/8.35
            ea2 = ea[0:len(ea)/2,0:len(ea)/2]/8.35
            np.fill_diagonal(ea1,0)
            sc += 10*sum(sum(ea2>0.5))
            sc += 10*sum(sum(ea1>0.5))
            #print "B",sc
        
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

class EndSystem:
    mult = {'1NGO': 1.0, '2NGO': 2.0, '1GO': 5.0, '2GO': 10.0}
    def __init__(self, tilesys, tileadj=False, pairs=None):
        
        # Set up variables, etc.
        self.ef = sd.energyfuncs_santalucia(mismatchtype='max')
        self.tileadj = tileadj
        self.ends = tilesys['ends']
        self.tiles = tilesys['tiles']
        self.tilesystem = tilesys
        
        # Create a dict to give us pair classes from an input of two interacting ends 
        # (not the usual two ends that can go in each other's place!)
        if not pairs:
            pairs = sens.consolidate_pairs( sens.senspairs(tilesys), comcomp=1, onlytop=True )
                
        state = State()
        self.names = {}
        # Now pull sequences from the ends:
        state.seqs['TD'], self.names['TD'] = (list(x) for x in zip(*[ [end['seq'].lower(),end['name']] for end in self.ends if end['type'] == 'TD']))
        state.seqs['DT'], self.names['DT'] = (list(x) for x in zip(*[ [end['seq'].lower(),end['name']] for end in self.ends if end['type'] == 'DT']))
        self.adjs = {'DT':[], 'TD':[] }
        self.cadjs = {'DT':[], 'TD':[] }
        self.enlocs = {}
        for i, endn in enumerate(self.names['TD']):
            self.enlocs[endn] = (i,'TD')
        for i, endn in enumerate(self.names['DT']):
            self.enlocs[endn] = (i,'DT')
        te = list(flatten([ x['ends'] for x in self.tiles+tilesys['seed']['adapters'] ]))
        ta = list(flatten([ x['adjs'] for x in self.tiles+tilesys['seed']['adapters'] ]))
        
        self.pairdict = {}
        for pairclass,memberset in pairs.items():
            for x,y in memberset:
                self.pairdict[(x,ecomp(y))] = pairclass
                self.pairdict[(y,ecomp(x))] = pairclass
        
        for t in [ 'TD', 'DT' ]:
            for end in self.names[t]:
                self.adjs[t].append( set([ adj.lower() for adj,name in zip(ta,te) if name == end ]) )
                self.cadjs[t].append( set([ adj.lower() for adj,name in zip(ta,te) if name == ecomp(end) ]) )
        self.initstate = state
            
        
    def mutate(self, state):
        # Start by deciding to swap TD or DT ends.
        if random.rand() > 1.0*len(state.seqs['TD'])/(len(state.seqs['DT'])+len(state.seqs['TD'])):
            t = 'DT'
        else:
            t = 'TD'
        
        en = len(state.seqs[t])
        
        a = random.randint( 0, en-1 )
        b = random.randint( 0, en-1 )
        state.seqs[t][a],state.seqs[t][b] = state.seqs[t][b],state.seqs[t][a]
    
    def randomize(self, state):
        pyrand.shuffle(state.seqs['DT'])
        pyrand.shuffle(state.seqs['TD'])
        
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
            
            # create lists
            if xt == 'TD': 
                if not xc:
                    xsq = sd.endarray([ state.seqs['TD'][xi]+z for z in self.adjs['TD'][xi] ], 'TD')
                else:
                    xsq = sd.endarray([ wc(state.seqs['TD'][xi])+z for z in self.cadjs['TD'][xi] ], 'TD')
            else:
                if not xc:
                    xsq = sd.endarray([ z+state.seqs['DT'][xi] for z in self.adjs['DT'][xi] ], 'DT')
                else:
                    xsq = sd.endarray([ z+wc(state.seqs['DT'][xi]) for z in self.cadjs['DT'][xi] ], 'DT')
            if yt == 'TD': 
                if not yc:
                    ysq = sd.endarray([ state.seqs['TD'][yi]+z for z in self.adjs['TD'][yi] ], 'TD')
                else:
                    ysq = sd.endarray([ wc(state.seqs['TD'][yi])+z for z in self.cadjs['TD'][yi] ], 'TD')
            else:
                if not yc:
                    ysq = sd.endarray([ z+state.seqs['DT'][yi] for z in self.adjs['DT'][yi] ], 'DT')
                else:
                    ysq = sd.endarray([ z+wc(state.seqs['DT'][yi]) for z in self.cadjs['DT'][yi] ], 'DT')
            
            # now create an array and average it.
            enarray = self.ef.uniform( np.repeat(xsq,ysq.shape[0],0), np.tile(ysq,(xsq.shape[0],1)) )
            #print "%s %s, %s: %g (%g)" % (xn, yn, pairclass, np.mean(enarray), np.mean(enarray)*self.mult[pairclass])
            sc += np.mean(enarray)*self.mult[pairclass]
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