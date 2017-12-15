# Utility functions for sequences
import string

_L_TO_N = { 'a': frozenset((0,)),
            'b': frozenset((1, 2, 3)),
            'c': frozenset((1,)),
            'd': frozenset((0, 2, 3)),
            'g': frozenset((2,)),
            'h': frozenset((0, 1, 3)),
            'k': frozenset((2, 3)),
            'm': frozenset((0, 1)),
            'n': frozenset((0, 1, 2, 3)),
            'r': frozenset((0,2)),
            's': frozenset((1, 2)),
            't': frozenset((3,)),
            'v': frozenset((0, 1, 2)),
            'w': frozenset((0, 3)),
            'y': frozenset((1,3))}
# FIXME: X (synonym for N) is not allowed.

_N_TO_L = { v: i for i,v in _L_TO_N.items() }

_WC = {k: _N_TO_L[frozenset(3-v for v in s)] for k, s in _L_TO_N.items()}
_WC_WITH_PUNC = {**_WC, **{'-': '-', '+': '+', ' ': ' '}}

_AMBBASES = frozenset( _L_TO_N.keys() - {'a','c','g','t'} )

validnts = _L_TO_N.keys()

def revcomp(s):
    return "".join(_WC_WITH_PUNC[l] for l in reversed(s))

def is_null(seq):
    """Return True if a sequence consists only of Ns, or is empty. 
    Return False otherwise."""
    if not seq:
        return True
    check_bases(seq)
    return set(seq.lower()).issubset(set('n').union(set(string.whitespace)))

def is_definite(seq):
    """Return True if a sequence consists only of defined bases.  Return
    False otherwise.  If blank, return False.
    """
    if not seq:
        return False
    check_bases(seq)
    if set(seq.lower()).issubset(set(string.whitespace)):
        return False
    return set(seq.lower()).issubset({'a','g','c','t'}.union(set(string.whitespace))) 

def check_bases(seq):
    if not set(seq.lower()).issubset(set(_L_TO_N.keys()).union(set(string.whitespace))):
        raise ValueError("Sequence has unknown bases.")

def count_ambiguous(seq):
    """Return the number of ambiguous bases in a sequence."""
    check_bases(seq)
    return sum( 1 for x in seq.lower() if x in _AMBBASES )

def length(seq):
    """Return the length of a sequence, stripping whitespace.  This does not handle
    extended labels, etc."""
    check_bases(seq)
    return len( seq.translate( str.maketrans('','',string.whitespace) ) )


def merge(seq1, seq2):
    """Merge two sequences together, returning a single sequence that
    represents the constraint of both sequences.  If the sequences
    can't be merged, raise a MergeConflictError.

    FIXME: this needs to intelligently handle case and whitespace.
    """
    
    check_bases(seq1)
    check_bases(seq2)
    if len(seq1) != len(seq2):
        raise MergeConflictError(seq1,seq2,'length',len(seq1),len(seq2))
    
    out = []
    for i,(n1,n2) in enumerate(zip(seq1.lower(),seq2.lower())):
        try:
            out.append(_N_TO_L[frozenset.intersection(_L_TO_N[n1], _L_TO_N[n2])])
        except KeyError as e:
            if e.args[0] == frozenset():
                raise MergeConflictError(seq1,seq2,i,n1,n2) from None
            else:
                raise e
    return ''.join(out)


class MergeConflictError(ValueError):
    """
    Merge of items failed because of conflicting information.
    Arguments are (item1, item2, location or property, value1, value2)
    """
    
class MergeConflictsError(ValueError):
    """
    Merge of multiple items failed because individual merges
    raised MergeConflictErrors.
    Arguments are ([list of MergeConflictErrors])
    """
