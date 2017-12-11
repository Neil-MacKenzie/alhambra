import copy
from . import seq
from ruamel.yaml.comments import CommentedSeq
from ruamel.yaml.representer import RoundTripRepresenter

class NamedList(CommentedSeq):
    """\
A class for a list of dicts, where some dicts have a 'name' item which should be unique, but others might not.
Indexing works with either number or name.
Note that updating dicts may make the list inconsistent."""
    def __init__(self,x=[]):
        CommentedSeq.__init__(self,x)
        
    def __getitem__(self, i):
        if isinstance(i,str):
            r = [ x for x in self if x.get('name',None)==i ]
            if len(r)>1:
                raise KeyError("There are {} elements named {}.".format(len(r),i))
            elif len(r)==0:
                raise KeyError("No element named {} found.".format(i))
            else:
                return r[0]
        else:
            return CommentedSeq.__getitem__(self,i)

    def check_consistent(self):
        """\
Checks that each name appears only once.  On failure, returns a ValueError
with (message, {failed_name: count}).  Otherwise, return with no output.
"""
        names = [ v['name'] for v in self if 'name' in v.keys() ]

        # if we have *no* names, the tests will fail, but we are obviously
        # consistent, so:
        if not names:
            return
        
        from collections import Counter
        namecounts = Counter(names)
       
        if max(namecounts.values()) > 1:
            badcounts = { n:v for n,v in namecounts.items() if v>1 }
            raise ValueError("Inconsistent NamedList.", badcounts)
        
        
    def __setitem__(self, i, v):
        if isinstance(i,str):
            r = [ (ii,x) for ii,x in enumerate(self) if x.get('name',None)==i ]
            if len(r)>1:
                raise KeyError("There are {} elements named {}.".format(len(r),i))
            elif len(r)==0:
                self.append(v)
            else:
                CommentedSeq.__setitem__(self, r[0][0], v)
        else:
            CommentedSeq.__setitem__(self,i,v)

    def keys(self):
        return [ x['name'] for x in self if 'name' in x.keys() ]

    def __deepcopy__(self, memo):
        # FIXME: this is here waiting on ruamel.yaml bugfix.
        res = self.__class__()
        memo[id(self)] = res
        for k in self:
            res.append(copy.deepcopy(k))
            self.copy_attributes(res, deep=True)
        return res

RoundTripRepresenter.add_representer(NamedList,
                                     RoundTripRepresenter.represent_list)

lton = { 'a': (0,),
         'b': (1, 2, 3),
         'c': (1,),
         'd': (0, 2, 3),
         'g': (2,),
         'h': (0, 1, 3),
         'k': (2, 3),
         'm': (0, 1),
         'n': (0, 1, 2, 3),
         's': (1, 2),
         't': (3,),
         'v': (0, 1, 2),
         'w': (0, 3) }
wc = {   'a': 't',
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
ntol={v:i for i,v in lton.items()}

def merge_ends(end1, end2):
    """\
Given ends end1 and end2, assuming they describe the same sticky end, merge them into a single end, combining
information from each and enforcing that the two input ends consistently make a single output end.
    """
    # Of things in the ends: fseq might need special care, everything else just needs to match.
    out = copy.deepcopy(end1)
    for i,v in end2.items():
        if i in out.keys() and i=='fseq':
            # we merge sequences    
            out[i] = seq.merge(out[i],v)
        elif i in out.keys() and out[i]!=v:
            # we insisted all others must be equal
            raise ValueError(\
                             "end1 has {}={}, end2 has {}={}".format(i,out[i],i,v))
        elif i not in out.keys():
            out[i] = copy.deepcopy(v)
    return out


import time


class ProgressLogger(object):
    def __init__(
            self,
            logger,
            N,
            seconds_interval=60):
        self.logger = logger
        stime = time.perf_counter()
        self.stime = stime
        self.ltime = stime
        self.li = 0
        self
        self.seconds_interval = seconds_interval
        self.N = N
        self.logger.info("starting {} tasks".format(self.N))
        
    def update(self, i):
        ctime = time.perf_counter()
        if ctime - self.ltime > self.seconds_interval:
            self.logger.info(
                "finished {}/{}, {} s elapsed, {} s est remaining".format(
                    i, self.N, int(ctime - self.stime),
                    int((self.N-i)*(ctime-self.stime)/i)))
            self.ltime = ctime
            self.li = i
