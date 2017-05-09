import copy
from . import tiletypes
from . import seq
from ruamel.yaml.comments import CommentedSeq
from ruamel.yaml.representer import RoundTripRepresenter

class named_list(CommentedSeq):
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
            raise ValueError("Inconsistent named_list.", badcounts)
        
        
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

RoundTripRepresenter.add_representer(named_list,
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


def merge_endlists(endlist1, endlist2, fail_immediate=False, 
                   in_place=False ):
    """\
Given end lists `endlist1` and `endlist2`, merge the two lists, using
`merge_ends` to merge any named ends that are present in both.

Parameters
----------

endlist1: named_list of sticky ends.

endlist2: named_list OR just a list of sticky ends.  If it just a list,
    which may have multiple copies of the same named sticky end, each
    sticky end in order is merged.

fail_immediate: (default False) if True, fail immediately on a failed
merge, passing through the ValueError from merge_ends.  If False, finish
merging the two lists, then raise a ValueError listing *all* ends that
failed to merge.

in_place: (default False) if True, do merging in place in endlist1.  Note
   the merged and added ends from endlist2 will be copies regardless.

output: a merged named_list of sticky ends.


Exceptions
----------

ValueError: In the event of a failed merge, when named ends cannot be
    merged.  If fail_immediate is True, then this is passed through from
    merge_ends.  If fail_immediate is False, then the ValueError has the
    following args: 
    ("message", [exceptions], [failed_name,failed_end1,failed_end2) ...],out)
"""
    endlist1 = named_list(endlist1)

    # Check consistency of each named_list
    endlist1.check_consistent()
    try:
        endlist2.check_consistent()
    except AttributeError:
        pass
    
    if not in_place:
        out = copy.deepcopy(endlist1)
    else:
        out = endlist1
    
    exceptions = []
    errors = []
    
    for end in endlist2:
        if end['name'] in out.keys():
            try:
                out[end['name']] = merge_ends( out[end['name']], end)
            except ValueError as e:
                if fail_immediate:
                    raise e
                else:
                    exceptions.append(e)
                    errors.append((end['name'],out[end['name']],end))
        else:
            out.append(copy.deepcopy(end))
    
    if errors:
        errorstring = " ".join( [e[0] for e in errors] )
        raise ValueError(\
              "Errors merging {}".format(errorstring), exceptions, errors, out)
    
    return out
